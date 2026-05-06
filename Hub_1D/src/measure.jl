import Pkg

required_packages = ["ITensors", "ITensorMPS", "JLD2"]

installed_packages = [pkg.name for pkg in values(Pkg.dependencies())]

for pkg in required_packages
    if !(pkg in installed_packages)
        println("installing... [$pkg]")
        Pkg.add(pkg)
    end
end
using ITensors
using ITensorMPS
using LinearAlgebra
using Statistics
using JLD2 

# === Folded chain ordering for better PBC convergence ===
# Maps physical site index → MPS site index
# Standard: 1,2,3,...,N  →  MPS 1,2,3,...,N
# Folded:   1,N,2,N-1,3,...  →  phys_to_mps reorders so PBC bond is short-range
function folded_permutation(N::Int)
    mps_to_phys = zeros(Int, N)
    left, right = 1, N
    for m in 1:N
        if isodd(m)
            mps_to_phys[m] = left
            left += 1
        else
            mps_to_phys[m] = right
            right -= 1
        end
    end
    # Invert: phys_to_mps[p] = m
    phys_to_mps = zeros(Int, N)
    for m in 1:N
        phys_to_mps[mps_to_phys[m]] = m
    end
    return phys_to_mps
end

identity_permutation(N::Int) = collect(1:N)

function nearest_neighbor_std(mat, N::Int)
    vals = [mat[i, mod1(i + 1, N)] for i in 1:N]
    return std(vals)
end

function pbc_convergence_stats(cdagc_up, cdagc_dn, N::Int=size(cdagc_up, 1))
    return (
        nn_std_up=nearest_neighbor_std(cdagc_up, N),
        nn_std_dn=nearest_neighbor_std(cdagc_dn, N),
        boundary_diff_up=abs(cdagc_up[1, 2] - cdagc_up[N, 1]),
        boundary_diff_dn=abs(cdagc_dn[1, 2] - cdagc_dn[N, 1]),
    )
end

# === Custom Observer for PBC convergence check per sweep ===
mutable struct PBCObserver <: AbstractObserver
    N::Int
    separate_chains::Bool
    perm::Vector{Int}
    pbc_threshold::Float64
    energies::Vector{Float64}
    max_stds_up::Vector{Float64}
    max_stds_dn::Vector{Float64}
    boundary_diffs_up::Vector{Float64}
    boundary_diffs_dn::Vector{Float64}
    converged::Bool
end

function PBCObserver(N::Int, separate_chains::Bool, perm::Vector{Int}; pbc_threshold::Float64=1e-4)
    return PBCObserver(N, separate_chains, perm, pbc_threshold,
        Float64[], Float64[], Float64[], Float64[], Float64[], false)
end

function ITensorMPS.measure!(obs::PBCObserver; kwargs...)
    half_sweep = kwargs[:half_sweep]
    bond = kwargs[:bond]
    # Only measure at the last bond of the second half-sweep (end of full sweep)
    if half_sweep != 2 || bond != 1
        return nothing
    end

    psi = kwargs[:psi]
    sweep = kwargs[:sweep]
    energy = kwargs[:energy]
    sites = siteinds(psi)

    push!(obs.energies, energy)

    # Measure cdagc (in physical site order)
    cdagc_up, cdagc_dn = measure_cdagc(psi, sites, obs.N, obs.separate_chains, obs.perm)

    stats = pbc_convergence_stats(cdagc_up, cdagc_dn, obs.N)
    push!(obs.max_stds_up, stats.nn_std_up)
    push!(obs.max_stds_dn, stats.nn_std_dn)
    push!(obs.boundary_diffs_up, stats.boundary_diff_up)
    push!(obs.boundary_diffs_dn, stats.boundary_diff_dn)

    println("  [PBC Sweep $sweep] E=$energy | NN std(up)=$(round(stats.nn_std_up; sigdigits=4)) NN std(dn)=$(round(stats.nn_std_dn; sigdigits=4)) | Δbnd(up)=$(round(stats.boundary_diff_up; sigdigits=4)) Δbnd(dn)=$(round(stats.boundary_diff_dn; sigdigits=4))")

    obs.converged = (stats.nn_std_up < obs.pbc_threshold && stats.nn_std_dn < obs.pbc_threshold)
end

function ITensorMPS.checkdone!(obs::PBCObserver; kwargs...)
    return obs.converged
end

function build_sites(N::Int, separate_chains::Bool, Nup, Ndn, perm::Vector{Int})
    nup = isnothing(Nup) ? N ÷ 2 : Nup
    ndn = isnothing(Ndn) ? N ÷ 2 : Ndn

    if separate_chains
        # Fermion mode: 2N sites
        sites = siteinds("Fermion", 2N; conserve_qns=true)
        init_state = fill("0", 2N)
        for i in 1:nup
            init_state[2perm[i] - 1] = "1"  # up site at MPS position
        end
        for i in 1:ndn
            init_state[2perm[i]] = "1"       # dn site at MPS position
        end
    else
        # Electron mode: N sites, dim=4
        sites = siteinds("Electron", N; conserve_qns=true)
        init_state = fill("Emp", N)
        ndouble = min(nup, ndn)
        nup_remaining = nup - ndouble
        ndn_remaining = ndn - ndouble
        idx = 1
        for _ in 1:ndouble
            init_state[perm[idx]] = "UpDn"
            idx += 1
        end
        for _ in 1:nup_remaining
            init_state[perm[idx]] = "Up"
            idx += 1
        end
        for _ in 1:ndn_remaining
            init_state[perm[idx]] = "Dn"
            idx += 1
        end
    end

    return sites, init_state, nup, ndn
end

function build_hamiltonian(sites, N::Int, t::Float64, tp::Float64, U::Float64, separate_chains::Bool, perm::Vector{Int})
    os = OpSum()

    if !separate_chains
        # === Electron mode ===
        # NN hopping (PBC)
        for i in 1:N
            j = mod1(i + 1, N)
            mi, mj = perm[i], perm[j]
            os += -t, "Cdagup", mi, "Cup", mj
            os += -t, "Cdagup", mj, "Cup", mi
            os += -t, "Cdagdn", mi, "Cdn", mj
            os += -t, "Cdagdn", mj, "Cdn", mi
        end
        # NNN hopping (PBC)
        for i in 1:N
            j = mod1(i + 2, N)
            mi, mj = perm[i], perm[j]
            os += -tp, "Cdagup", mi, "Cup", mj
            os += -tp, "Cdagup", mj, "Cup", mi
            os += -tp, "Cdagdn", mi, "Cdn", mj
            os += -tp, "Cdagdn", mj, "Cdn", mi
        end
        # Hubbard U
        for i in 1:N
            os += U, "Nupdn", perm[i]
        end
    else
        # === Fermion mode: 2N sites ===
        # NN hopping (t)
        for i in 1:N
            inext = mod1(i + 1, N)
            # up-up hopping
            a = 2perm[i] - 1
            b = 2perm[inext] - 1
            os += -t, "Cdag", a, "C", b
            os += -t, "Cdag", b, "C", a
            # dn-dn hopping
            a = 2perm[i]
            b = 2perm[inext]
            os += -t, "Cdag", a, "C", b
            os += -t, "Cdag", b, "C", a
        end
        # NNN hopping (t')
        for i in 1:N
            inext2 = mod1(i + 2, N)
            # up-up
            a = 2perm[i] - 1
            b = 2perm[inext2] - 1
            os += -tp, "Cdag", a, "C", b
            os += -tp, "Cdag", b, "C", a
            # dn-dn
            a = 2perm[i]
            b = 2perm[inext2]
            os += -tp, "Cdag", a, "C", b
            os += -tp, "Cdag", b, "C", a
        end
        # Hubbard U: up and dn at same physical site
        for i in 1:N
            os += U, "N", 2perm[i] - 1, "N", 2perm[i]
        end
    end

    H = MPO(os, sites)
    return H
end

function run_dmrg(H, sites, init_state; χ::Int=10, nsweeps::Int=20, cutoff::Float64=1e-10, noise=nothing, observer=nothing, psi_init::Union{MPS,Nothing}=nothing)
    psi0 = isnothing(psi_init) ? randomMPS(sites, init_state; linkdims=min(χ, 10)) : psi_init

    # Build maxdim schedule: ramp up to χ (skip ramp if warm-starting)
    if !isnothing(psi_init)
        maxdim_schedule = fill(χ, nsweeps)
    else
        ramp = filter(x -> x < χ, [10, 20, 50, 100, 200, 500])
        maxdim_schedule = vcat(ramp, fill(χ, max(nsweeps - length(ramp), 0)))
        if length(maxdim_schedule) > nsweeps
            maxdim_schedule = maxdim_schedule[1:nsweeps]
        end
    end

    kwargs = Dict{Symbol,Any}(
        :maxdim => maxdim_schedule,
        :cutoff => cutoff,
        :nsweeps => nsweeps,
    )

    if !isnothing(noise)
        noise_schedule = length(noise) >= nsweeps ? noise[1:nsweeps] : vcat(noise, fill(0.0, nsweeps - length(noise)))
        kwargs[:noise] = noise_schedule
    end

    if !isnothing(observer)
        kwargs[:observer] = observer
    end

    energy, psi = dmrg(H, psi0; kwargs...)

    println("DMRG energy = $energy")
    return energy, psi
end

function measure_cdagc(psi, sites, N::Int, separate_chains::Bool, perm::Vector{Int})
    if !separate_chains
        # Electron mode: correlation_matrix returns in MPS order
        cdagc_mps_up = correlation_matrix(psi, "Cdagup", "Cup")
        cdagc_mps_dn = correlation_matrix(psi, "Cdagdn", "Cdn")
        # Permute back to physical order: cdagc_phys[i,j] = cdagc_mps[perm[i], perm[j]]
        cdagc_up = zeros(Float64, N, N)
        cdagc_dn = zeros(Float64, N, N)
        for i in 1:N, j in 1:N
            cdagc_up[i, j] = cdagc_mps_up[perm[i], perm[j]]
            cdagc_dn[i, j] = cdagc_mps_dn[perm[i], perm[j]]
        end
    else
        # Fermion mode: full correlation matrix in MPS order, then extract
        cdagc_full = correlation_matrix(psi, "Cdag", "C")
        cdagc_up = zeros(Float64, N, N)
        cdagc_dn = zeros(Float64, N, N)
        for i in 1:N, j in 1:N
            # Physical site i → MPS up index 2*perm[i]-1, dn index 2*perm[i]
            cdagc_up[i, j] = cdagc_full[2perm[i] - 1, 2perm[j] - 1]
            cdagc_dn[i, j] = cdagc_full[2perm[i], 2perm[j]]
        end
    end

    return cdagc_up, cdagc_dn
end


function check_pbc_convergence(io::IO, cdagc_up, cdagc_dn, N::Int)
    stats = pbc_convergence_stats(cdagc_up, cdagc_dn, N)

    println(io, "\n=== PBC Convergence Check ===")
    println(io, "[up] NN std = $(round(stats.nn_std_up; sigdigits=4))")
    println(io, "[up] |⟨c†_1 c_2⟩ - ⟨c†_N c_1⟩| = $(round(stats.boundary_diff_up; sigdigits=4))")
    println(io, "[dn] NN std = $(round(stats.nn_std_dn; sigdigits=4))")
    println(io, "[dn] |⟨c†_1 c_2⟩ - ⟨c†_N c_1⟩| = $(round(stats.boundary_diff_dn; sigdigits=4))")
    println(io, "=============================\n")
end

function check_pbc_convergence(cdagc_up, cdagc_dn, N::Int)
    check_pbc_convergence(stdout, cdagc_up, cdagc_dn, N)
end

function get_Hubbard1D_chain(N, t, tp, U;
    χ::Int=20,
    separate_chains::Bool=false,
    Nup=nothing,
    Ndn=nothing,
    nsweeps::Int=20,
    cutoff::Float64=1e-10,
    noise=nothing,
    save_path=nothing,
    pbc_threshold::Float64=1e-3,
    auto_chi::Bool=true,
    chi_max::Int=500,
    chi_step::Float64=1.25,
    folded::Bool=false,
)
    perm = folded ? folded_permutation(N) : identity_permutation(N)

    println("=== 1D Hubbard Model (N=$N, t=$t, t'=$tp, U=$U) ===")
    println("Mode: $(separate_chains ? "Fermion (2N=$(2N) sites)" : "Electron (N=$N sites)")$(folded ? " [folded]" : "")")

    # 1. Build sites and initial state
    sites, init_state, nup, ndn = build_sites(N, separate_chains, Nup, Ndn, perm)
    println("Half-filling: Nup=$nup, Ndn=$ndn")

    # 2. Build Hamiltonian with PBC
    H = build_hamiltonian(sites, N, Float64(t), Float64(tp), Float64(U), separate_chains, perm)

    # 3. Run DMRG (with optional auto-χ)
    if !auto_chi
        println("Bond dimension χ=$χ, nsweeps=$nsweeps")
        obs = PBCObserver(N, separate_chains, perm; pbc_threshold=pbc_threshold)
        energy, psi = run_dmrg(H, sites, init_state; χ=χ, nsweeps=nsweeps, cutoff=cutoff, noise=noise, observer=obs)
    else
        # Auto-χ: increase χ until PBC converges
        step = chi_step > 0 ? chi_step : 1.5  # default step = 1.5
        current_chi = χ
        psi = nothing
        obs = nothing
        energy = 0.0
        all_sweep_data = Dict{String,Any}()
		stack = 0.0
        while current_chi <= chi_max
            println("\n--- Auto-χ: running with χ=$current_chi ---")
            obs = PBCObserver(N, separate_chains, perm; pbc_threshold=pbc_threshold)
            energy, psi = run_dmrg(H, sites, init_state; χ=current_chi, nsweeps=nsweeps, cutoff=cutoff, noise=noise, observer=obs, psi_init=psi)

            cdagc_up_chi, cdagc_dn_chi = measure_cdagc(psi, sites, N, separate_chains, perm)
            all_sweep_data["chi_$current_chi"] = Dict(
                "energies" => copy(obs.energies),
                "max_stds_up" => copy(obs.max_stds_up),
                "max_stds_dn" => copy(obs.max_stds_dn),
                "boundary_diffs_up" => copy(obs.boundary_diffs_up),
                "boundary_diffs_dn" => copy(obs.boundary_diffs_dn),
                "pbc_converged" => obs.converged,
                "cdagc_up" => copy(cdagc_up_chi),
                "cdagc_dn" => copy(cdagc_dn_chi),
            )

            if obs.converged
                println("Auto-χ: PBC converged at χ=$current_chi (NN std threshold)")
                break
            end

            current_chi += Int(round(step*current_chi))
			if stack > 0 
				println("Auto-χ: reached chi_max=500 without PBC convergence")
				break
			end

            if current_chi > chi_max
				stack +=1
				current_chi = chi_max
            end
        end
        χ = min(current_chi, chi_max)
    end

    nsweeps_done = length(obs.energies)
    if obs.converged
        println("PBC converged (NN std threshold=$pbc_threshold) at sweep $nsweeps_done")
    else
        println("WARNING: PBC did NOT converge after $nsweeps_done sweeps (NN std threshold=$pbc_threshold)")
    end

    # 4. Measure final cdagc (in physical site order)
    cdagc_up, cdagc_dn = measure_cdagc(psi, sites, N, separate_chains, perm)

    # 5. Final PBC convergence summary
    check_pbc_convergence(cdagc_up, cdagc_dn, N)

    # 6. Save if requested
    if !isnothing(save_path)
        params = Dict(
            "N" => N, "t" => t, "tp" => tp, "U" => U,
            "chi" => χ, "separate_chains" => separate_chains,
            "Nup" => nup, "Ndn" => ndn,
            "folded" => folded,
        )
        sweep_data = if auto_chi
            all_sweep_data
        else
            Dict(
                "energies" => obs.energies,
                "max_stds_up" => obs.max_stds_up,
                "max_stds_dn" => obs.max_stds_dn,
                "boundary_diffs_up" => obs.boundary_diffs_up,
                "boundary_diffs_dn" => obs.boundary_diffs_dn,
                "pbc_converged" => obs.converged,
            )
        end
        @save save_path cdagc_up cdagc_dn energy params sweep_data
        println("Results saved to $save_path")
    end

    return cdagc_up, cdagc_dn
end
