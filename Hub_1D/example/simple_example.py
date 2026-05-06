import numpy as np
import matplotlib.pyplot as plt
from juliacall import Main as jl

jl.include("../src/measure.jl")

N_sites = 20     
t_hop = -0.136     #(NN)
tp_hop = 0.020    #(NNN) 
U_int = 5.42



cdagc_up, cdagc_dn = jl.get_Hubbard1D_chain(
    N_sites, t_hop, tp_hop, U_int, # necessary
    χ=20,                  
    nsweeps=N_sites,           
    auto_chi=True,       
    chi_max=500,        
    folded=True,          
	save_path=f"../data/data_N{N_sites:d}_t{t_hop:.3f}_tp{tp_hop:.3f}_U{U_int:.3f}.jld2",
	pbc_threshold=1e-3
)

print("\n=== result ===")
print(f"type of cdagc_up: {type(cdagc_up)}")
print(f"size : {np.shape(cdagc_up)}")

corr_with_first_site = [cdagc_up[0, i] for i in range(N_sites)]


#simple plot
plt.plot(range(N_sites), corr_with_first_site, marker='o', label="Up Spin Correlation")
plt.xlabel("Site index (j)")
plt.ylabel("$\\langle c^\\dagger_0 c_j \\rangle$")
plt.title("Real-space Correlation Function")
plt.legend()
plt.grid(True)
plt.show()
