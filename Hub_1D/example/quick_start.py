import matplotlib.pyplot as plt
from juliacall import Main as jl

jl.include("../src/measure.jl")

N = 20
t = 1.0
tp = 1/6
U = 4.0
cdagc_up, cdagc_dn = jl.get_Hubbard1D_chain(N, t, tp, U)

plt.plot(cdagc_up[0, :], 'bo-', label="Up Spin")
plt.plot(cdagc_dn[0, :], 'rs--', label="Down Spin")

plt.title("1D Hubbard DMRG Correlation (N=20, U=4.0)")
plt.xlabel("Site index (j)")
plt.ylabel("$\\langle c^\\dagger_0 c_j \\rangle$")
plt.legend()
plt.grid(True)
plt.show()
