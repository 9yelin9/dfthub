# dfthub
First-Principles Hubbard Parameter Extraction.

## Systems
You can choose one of the systems in `data/` directory.
| Directory | Lattice constant | Spin polarization | Description |
| :--- | :---: | :---: | :---: |
| Li | $a_0$ | O | - |
| Li-a1 | $a_0$ | X | - |
| Li-a2 | $2a_0$ | X | - |

## Prerequisite
- [ABINIT v10.6+](https://github.com/abinit/abinit/releases/tag/10.6.5)
- [Wannier90 v3.1+](https://wannier.org/download/)
- Python v3.12+
- Python dependencies: **abipy**, numpy, scipy, matplotlib

## Usage
### 1. Hopping parameter ($t$) extraction
1. Navigate to the target system's directory:
```bash
cd data/Li-a2_init
```
2. Run ABINIT-Wannier90 interface (⚠️Use single processor):
```bash
mpirun -np 1 abinit wann_.abi
```
3. Check `w90_hr.dat` (Columns: R1, R2, R3, i, j, Re(t), Im(t)):
```bash
...
 -4    0    0    1    1    0.003014    0.000000
 -3    0    0    1    1   -0.005475    0.000000
 -2    0    0    1    1    0.020240    0.000000
 -1    0    0    1    1   -0.136003   -0.000000
  0    0    0    1    1   -2.898835    0.000000
  1    0    0    1    1   -0.136003    0.000000
  2    0    0    1    1    0.020240   -0.000000
  3    0    0    1    1   -0.005475   -0.000000
  4    0    0    1    1    0.003014   -0.000000
...
```

### 2. Interaction parameters ($U$ and $J$)
1. Navigate to the target system's directory:
```bash
cd data/Li-a2_init
```
2. Run ABINIT-cRPA calculation (⚠️Use single processor):
```bash
mpirun -np 1 abinit crpa_.abi
```
3. Check `crpa_.abo`:
```bash
...
== Calculation of the screened interaction on the correlated orbital U m == 

= Start loop over frequency

--- For frequency w =   1  -------------


Diagonal cRPA interaction
1  5.425

U'=U(m1,m2,m1,m2) for the cRPA interaction
-      1    
1  5.425

Hubbard cRPA interaction for w =  1, U=1/(2l+1)**2 \sum U(m1,m2,m1,m2)=    5.4251    0.0000

(Hubbard cRPA interaction for w =   1, U=1/(2l+1) \sum U(m1,m1,m1,m1)=    5.4251    0.0000)

Hund coupling J=U(m1,m1,m2,m2) for the cRPA interaction
-      1    
1  5.425

cRPA interaction value of J=U-1/((2l+1)(2l)) \sum_{m1,m2} (U(m1,m2,m1,m2)-U(m1,m2,m2,m1))=    0.0000    0.0000


Hund coupling J2=U(m1,m2,m2,m1) for the cRPA interaction
-      1    
1  5.425
...
```

### 3. Corrected density calculation
1. Run the Hubbard model simulation code:
```bash
cd Hub_1D/example
python3 simple_example.py
```

2. Run the density correction script with the simulation result: 
```bash
./dfthub.py data/Li-a2_init -d Hub_1D/data/data_N20_t-0.136_tp0.020_U5.420.jld2
```

3. Check the output directory and the corrected density `dfthub_DEN`:
```bash
ls data/Li-a2_t-0.136_tp0.020_U5.420
```
*Note: The output directory name is automatically generated based on the simulation result name.*

Expected files:
```bash
crpa_.abi  dfthub_DEN  w90.win  wann_.abi
```

4. Now repeat the workflow!
