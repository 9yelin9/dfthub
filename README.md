# dfthub
First-Principles Hubbard Parameter Extraction

## Systems
| Directory | Lattice constant | Spin polarization | Description |
| :--- | :---: | :---: | :--- |
| Li | $a_0$ | O | |
| Li-a1 | $a_0$ | X | |
| Li-a2 | $2a_0$ | X | |

## Usage
### 1. Hopping parameter ($t$) extraction:
1. Navigate to target system's directory:
```bash
cd Li-a2
```
2. Execute ABINIT-Wannier90 interface (⚠️ Use single processor):
```bash
mpirun -np 1 abinit Li_wann_.abi
```
3. Check `w90_hr.dat` (Columns: $R_1, R_2, R_3, i, j, \text{Re}(t), \text{Im}(t)):
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

### 2. Interaction parameters ($U$ and $J$):
1. Navigate to target system's directory:
```bash
cd Li-a2
```
2. Execute ABINIT-cRPA calculation (⚠️ Use single processor):
```bash
mpirun -np 1 abinit Li_crpa_.abi
```
3. Check `Li_crpa_.abo`:
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
1. Run density correction script: 
```bash
./dfthub.py Li-a2 -d
```
