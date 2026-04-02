# dfthub
First-Principles Hubbard Parameter Extraction

## Usage
### 1. Hopping parameter ($t$) extraction:
1. Navigate to target material's directory:
```bash
cd [mat_dir]
```
2. Execute ABINIT-Wannier90 interface (⚠️ Use single processor):
```bash
mpirun -np 1 abinit [mat]_wann_.abi
```
3. Check `w90_hr.dat`

### 2. Interaction parameters ($U$ and $J$):
1. Navigate to target material's directory:
```bash
cd [mat_dir]
```
2. Execute ABINIT-cRPA calculation (⚠️ Use single processor):
```bash
mpirun -np 1 abinit [mat]_crpa_.abi
```
3. Check `[mat]_crpa_.abo`

### 3. Corrected density calculation
1. Run density correction script: 
```bash
./dfthub.py [mat_dir] -d
```
