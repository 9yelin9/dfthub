# 1D Hubbard Python Interface

This project lets you run an ITensor-based DMRG calculation for the 1D Hubbard model from Python through Julia. The numerical backend is implemented in Julia and loaded through `juliacall`, so you can work from Python while using ITensor for the tensor-network computation.

If you want the shortest path to a first run, start with [QUICKSTART.md](./QUICKSTART.md).

## Overview

The main entry point is:

```python
get_Hubbard1D_chain(N, t, tp, U; ...)
```

It solves a 1D Hubbard chain with periodic boundary conditions and returns real-space single-particle correlation matrices.

For many user workflows, you can treat it as a black-box solver: choose the model parameters, set the filling you want through `Nup` and `Ndn` when needed, and give a reasonable `chi_max` while keeping `auto_chi=True`.

## What The Code Returns

The calculation returns two matrices:

- `cdagc_up[i, j] = <c^\dagger_{i,up} c_{j,up}>`
- `cdagc_dn[i, j] = <c^\dagger_{i,dn} c_{j,dn}>`

Both outputs are `N x N` correlation matrices in physical site order. In typical Python use through `juliacall`, you can treat them like array-like matrix objects and convert them to NumPy arrays if needed.

## Requirements

You need:

- a local copy of this repository
- Python 3
- the ability to install the required Python packages: `juliacall`, `numpy`, and `matplotlib`

The Julia side uses `ITensors`, `ITensorMPS`, and `JLD2`. The source file installs missing Julia packages on the first run, so the first execution can take noticeably longer than later runs.

## Setup

You can run this project in any Python environment that already has `juliacall`, `numpy`, and `matplotlib` installed.

Using a virtual environment is recommended, but not required.

If you want an isolated environment, from the project directory:

```bash
cd /path/to/1DHubbard_python
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

`venv` is not macOS-specific. It is part of the Python standard library and works on Linux, macOS, and Windows. Only the activation command changes by platform:

```bash
# macOS / Linux
source venv/bin/activate
```

```powershell
# Windows PowerShell
venv\Scripts\Activate.ps1
```

```cmd
# Windows Command Prompt
venv\Scripts\activate.bat
```

If you already have a working Python environment with the required packages, you can skip virtual-environment creation and use that environment directly.

The repository includes a `requirements.txt` file for the Python-side dependencies.

## Quick Example

Run this from the repository root:

```python
from juliacall import Main as jl

jl.include("src/measure.jl")

cdagc_up, cdagc_dn = jl.get_Hubbard1D_chain(
    20, 1.0, 1 / 6, 4.0,
    χ=20,
    nsweeps=20,
    auto_chi=True,
    chi_max=200,
    folded=True,
    save_path="data.jld2",
    pbc_threshold=1e-3,
)
```

This computes the spin-up and spin-down correlation matrices for a 20-site 1D Hubbard chain with nearest-neighbor hopping `t`, next-nearest-neighbor hopping `tp`, and on-site interaction `U`.

## Main Function

The public user-facing function is:

```python
cdagc_up, cdagc_dn = jl.get_Hubbard1D_chain(
    N, t, tp, U,
    χ=20,
    separate_chains=False,
    Nup=None,
    Ndn=None,
    nsweeps=20,
    cutoff=1e-10,
    noise=None,
    save_path=None,
    pbc_threshold=1e-3,
    auto_chi=True,
    chi_max=500,
    chi_step=1.25,
    folded=False,
)
```

If `Nup` and `Ndn` are not provided, the calculation uses half filling by default.

### Default Behavior

If you call `get_Hubbard1D_chain(N, t, tp, U)` with only the four required model parameters, the code uses these defaults:

- `χ=20`
- `separate_chains=False`
- `Nup=None`, `Ndn=None` (half filling)
- `nsweeps=20`
- `cutoff=1e-10`
- `noise=None`
- `save_path=None`
- `pbc_threshold=1e-3`
- `auto_chi=True`
- `chi_max=500`
- `chi_step=1.25`
- `folded=False`

### Recommended Usage

For most users, the defaults are a good starting point.

In particular, it is usually fine to keep most options unchanged and let the code handle bond-dimension growth automatically with `auto_chi=True`.

The main option to revisit is `chi_max`. If you increase `N`, you should usually increase `chi_max` as well.

The code prints periodic-boundary-condition convergence messages during the run. If you see output such as `WARNING: PBC did NOT converge ...` or `Auto-χ: reached chi_max=500 without PBC convergence`, increase `chi_max` and run again.

PBC convergence is measured from the nearest-neighbor correlation values in the correlation matrices. During each sweep, the code computes the standard deviation of `<c^\dag c_{i+1}>` around the ring for both spin sectors and declares convergence when both `NN std(up)` and `NN std(dn)` are smaller than `pbc_threshold`. The printed `Δbnd` values are an additional boundary diagnostic comparing `<c^\dagger_1 c_2>` and `<c^\dag_N c_1>`.

### Black-Box Usage

If you do not want to tune many DMRG details, the simplest workflow is:

- choose `N`, `t`, `tp`, and `U`
- keep the default half filling, or set `Nup` and `Ndn` for the filling you want
- leave `auto_chi=True`
- set a reasonable `chi_max`

In that mode, the code automatically increases the bond dimension up to `chi_max` and stops early if the periodic-boundary convergence check succeeds.

## Key Arguments

| Argument | Meaning |
| --- | --- |
| `N` | Number of physical sites in the 1D chain |
| `t` | Nearest-neighbor hopping amplitude |
| `tp` | Next-nearest-neighbor hopping amplitude |
| `U` | On-site Hubbard interaction |
| `χ` | Starting bond dimension used by DMRG |
| `nsweeps` | Number of DMRG sweeps for each run |
| `auto_chi` | If `True` (recommended), the code increases bond dimension automatically until the periodic-boundary convergence criterion is met or `chi_max` is reached |
| `chi_max` | Maximum bond dimension allowed when `auto_chi=True` |
| `chi_step` | Multiplicative growth factor for the next bond dimension in `auto_chi` mode |
| `folded` | If `True`, uses a reordered site layout that helps periodic-boundary-condition convergence |
| `separate_chains` | Chooses the site representation: `False` for Electron mode, `True` for Fermion mode |
| `Nup`, `Ndn` | Number of spin-up and spin-down particles; omitted means half filling |
| `cutoff` | SVD truncation cutoff used in DMRG |
| `noise` | Optional sweep-dependent noise schedule |
| `pbc_threshold` | Convergence threshold used in the periodic-boundary check |
| `save_path` | Optional `.jld2` output file path |

### Electron Mode vs Fermion Mode

- `separate_chains=False`: Electron mode with one site per physical lattice site
- `separate_chains=True`: Fermion mode with separate up and down chains in the underlying representation

For most users, Electron mode is the simpler default unless you specifically want to compare the two representations.

## Saved Output

If you pass `save_path="something.jld2"`, the code writes a JLD2 file containing:

- `cdagc_up`
- `cdagc_dn`
- `energy`
- `params`
- `sweep_data`

This is useful when you want to keep the final correlation matrices together with the run parameters and convergence history.

## Example Scripts

The repository already includes runnable Python examples:

- `example/quick_start.py`: the shortest example
- `example/simple_example.py`: a configurable run with explicit keyword arguments
- `example/comparison_2way.py`: compares folded and non-folded runs

Important: these example scripts load `../src/measure.jl`, so run them from inside the `example/` directory:

```bash
cd /path/to/1DHubbard_python/example
python quick_start.py
```

## Notes

- The first run can be slow because Julia packages may need to be installed or initialized.
- The code uses periodic boundary conditions and includes a convergence check based on nearest-neighbor statistics.
- In `auto_chi` mode, the code can save intermediate sweep information for different bond dimensions in `sweep_data`.
- The returned matrices are already reordered back to physical site order, even when `folded=True`.
