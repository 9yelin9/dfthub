# Quick Start

This guide is the fastest way to run the project and get the `cdagc_up` and `cdagc_dn` correlation matrices from Python.

## Prerequisites

You need:

- a local copy of this repository
- Python 3
- the ability to install the required Python packages: `juliacall`, `numpy`, and `matplotlib`

## Minimal Setup

You can run the project in any Python environment that already has `juliacall`, `numpy`, and `matplotlib` installed.

Using a virtual environment is recommended, but not required.

If you want an isolated environment, from the project root:

```bash
cd /path/to/1DHubbard_python
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

`venv` is not macOS-specific. It works on Linux, macOS, and Windows. If you are not on macOS, use the activation command that matches your shell:

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

If you already have a working Python environment with the required packages, you can skip virtual-environment creation.

If the virtual environment already exists, just activate it:

```bash
cd /path/to/1DHubbard_python
source venv/bin/activate
```

If you prefer not to use a virtual environment, install the same packages into your current Python environment:

```bash
pip install -r requirements.txt
```

## Run The Example

The provided example scripts use `../src/measure.jl`, so run them from inside the `example/` directory:

```bash
cd /path/to/1DHubbard_python/example
python quick_start.py
```

## What You Get

`quick_start.py` does the following:

- loads `../src/measure.jl` through `juliacall`
- runs `get_Hubbard1D_chain(N, t, tp, U)` with `N=20`, `t=1.0`, `tp=1/6`, and `U=4.0`
- returns `cdagc_up` and `cdagc_dn`
- plots the first-row correlations for the up-spin and down-spin sectors

The first run can be slower than usual because Julia packages may need to initialize or be installed automatically.

In many cases, you can treat the solver as a black box: choose the model parameters, keep half filling or set `Nup` and `Ndn` for a different filling, and set a reasonable `chi_max` if you want to control the maximum bond dimension used by `auto_chi`.

By default, calling `get_Hubbard1D_chain(N, t, tp, U)` uses half filling, `auto_chi=True`, and `chi_max=500` together with the other built-in defaults from [README.md](./README.md).

For most users, the defaults are a good starting point. The main option to revisit is `chi_max`: if you increase `N`, you should usually increase `chi_max` as well. If the run prints messages such as `WARNING: PBC did NOT converge ...` or `Auto-χ: reached chi_max=500 without PBC convergence`, increase `chi_max` and try again.

PBC convergence is checked from the spread of nearest-neighbor correlations around the ring. In the printed output, `NN std(up)` and `NN std(dn)` are the quantities used for the convergence test against `pbc_threshold`, while `Δbnd` is an additional boundary diagnostic.

## Next Steps

- Read [README.md](./README.md) for the full argument list and output description.
- Run `python simple_example.py` from the same `example/` directory for a more configurable workflow.
- Run `python comparison_2way.py` from the same `example/` directory to compare folded and non-folded runs.
