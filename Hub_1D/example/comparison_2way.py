import numpy as np
import matplotlib.pyplot as plt
from juliacall import Main as jl

jl.include("../src/measure.jl")

N_sites = 30     
t_hop = 1.0     
tp_hop = 1/6    
U_int = 4.0     

print("\n--- Running DMRG with Folded Chain ---")
cdagc_up_folded, _ = jl.get_Hubbard1D_chain(
    N_sites, t_hop, tp_hop, U_int, 
    χ=20,                  
    nsweeps=N_sites,            
    auto_chi=True,       
    chi_max=200,        
    folded=True,          
    save_path="data_folded.jld2",
    pbc_threshold=1e-3
)

print("\n--- Running DMRG with Identity Chain ---")
cdagc_up_identity, _ = jl.get_Hubbard1D_chain(
    N_sites, t_hop, tp_hop, U_int, 
    χ=20,                  
    nsweeps=N_sites,            
    auto_chi=True,       
    chi_max=200,        
    folded=False,          
    save_path="data_identity.jld2",
    pbc_threshold=1e-3
)

corr_folded = [cdagc_up_folded[0, i] for i in range(N_sites)]
corr_identity = [cdagc_up_identity[0, i] for i in range(N_sites)]

plt.figure(figsize=(10, 6))

plt.plot(range(N_sites), corr_folded, marker='o', linestyle='-', color='blue', label="Folded Chain")
plt.plot(range(N_sites), corr_identity, marker='s', linestyle='--', color='red', label="Identity Chain")

plt.xlabel("Site index (j)", fontsize=12)
plt.ylabel("$\\langle c^\\dagger_0 c_j \\rangle$", fontsize=12)
plt.title(f"Real-space Correlation Function: Folded vs Identity\n(N={N_sites}, t={t_hop}, t'={tp_hop:.3f}, U={U_int})", fontsize=14)
plt.legend(fontsize=12)
plt.grid(True, linestyle=':', alpha=0.7)
plt.tight_layout()

plt.show()
