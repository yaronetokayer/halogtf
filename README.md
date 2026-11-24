# HaloGTF

**HaloGTF** is a modern Python implementation of a 1D Lagrangian gravothermal fluid code. It simulates the dynamical evolution of self-interacting dark matter halos using the fluid approximation, based on a Fortran code originally developed by Prof. Frank van den Bosch (Yale University).

See [the runnable example notebook](examples/demo.ipynb), the
[scattering conventions](docs/scattering.md), and the
[output and units reference](docs/output-reference.md).

Contact yarone.tokayer-at-yale.edu with any questions/comments.

---

## Scattering choices

HaloGTF supports elastic hard-sphere scattering (either model with
`w=float("inf")`), screened Rutherford scattering (`model="rutherford"`),
and Møller scattering (`model="moller"`). Both velocity-dependent models
are normalized so that their low-speed limit corresponds to the same
`sigma_m_0` in cm²/g as the hard-sphere model, with no additional scaling
factor for the user to apply.

The [scattering reference](docs/scattering.md) explains the prescriptions,
the meaning of `sigma_m_0` and `w`, the viscosity-cross-section convention,
and the SMFP normalization.

## Overview

The code is organized around **two user-facing classes**:

### 1. `Config`

Stores static input parameters, grouped into modules:

* `io`: Output paths and snapshot cadence
* `grid`: Grid resolution and domain
* `init`: Initial profile (NFW, truncated NFW, or $\alpha$-$\beta$-$\gamma$)
* `scatter`: Scattering model, cross section, velocity scale, and conduction calibration
* `num`: Precision tolerances and iteration limits
* `stop`: Absolute time and central-density stopping limits
* `cosmo`: Cosmological parameters

### 2. `State`

Holds the dynamically evolving quantities:

* Radial grid and shell midpoints: `r`, `rmid`
* Physical variables: `m`, `rho`, `v2` (pressure is `rho*v2`; specific internal energy is `1.5*v2`)
* Diagnostics: `kn`, `mfp`, `kn_cond`, `mfp_cond`, `kn_c`, `kn_cond_c`, `ltemp`, `t_cool`, `t_dyn`, `du_max`
* Time tracking: `t`, `dt`, `step_count`, `snapshot_index`
* Characteristic scales (derived from `Config`)
* A `run()` method to evolve the system until collapse or a stopping condition is reached

Construct an initialized state with `State.from_config(config)`. Optionally pass
`ic_filepath` to start a new run from a saved profile. `State.from_dir()` is not
implemented and raises an error; it cannot resume a saved simulation.

`r` and `m` are edge arrays of length N+1. `rho`, `v2`, and the cell diagnostics
have length N and are associated with `rmid`.

In addition to these two classes, there are three plotting functions that are automatically imported, each of which also exists as a method of `State`:

### 1. `plot_time_evolution()`

Plot the time evolution of any one of a number of system-wide quantities.  Multiple simulations can be passed to compare their evolutions on the same plot.

### 2. `plot_snapshots()`

Plot up to three profiles of one simulation at specified points in time.

### 3. `make_movie()`

Animate the full evolution of up to three profiles.  This requires ffmpeg to be installed and callable with `ffmpeg` from the working directory.

Use `help()` for documentation on any of these classes, methods, and functions.

---

## Getting started

### Installation

```bash
git clone https://github.com/yaronetokayer/halogtf.git
cd halogtf
pip install -e .
```

Dependencies: Python 3.12+, `numpy`, `scipy`, `numba`, `matplotlib`, `tqdm`

### Example usage

```python
import halogtf as gtf

config = gtf.Config()
state = gtf.State.from_config(config)
state.run(steps=10000)  # Bounded demonstration; omit steps for a full run
```

To read an existing simulation without constructing a state:

```python
from halogtf.io.read import load_snapshot_bundle
snapshot = load_snapshot_bundle('./sims/Model00002')  # Latest saved snapshot
snapshot = load_snapshot_bundle('./sims/Model00002', snapshot=59)
```

Alternatively, you can call `from halogtf import Config, State`.  In that case, plotting functions will not be automatically imported.

We can also run for a specified duration:
```python
state.run(steps=100) # Run for 100 simulation steps
state.run(stoptime=55.0) # Run for a duration of 55.0 simulation time units
state.run(rho_c=5000.0) # Run until the central density exceeds 5000.0
```

The first satisfied user or configuration stopping condition ends the run. The configured density limit is checked only after time > 50; `run(rho0=...)` applies immediately after a step. Multiple `run()` commands can be executed in succession, and each will continue from the current state.  Call `state.reset()` to reset the state to its initial condition and reset the step counter.

To customize defaults:

```python
# Customize initial profile
config.init = "abg"                                 # Use ABG with default params
config.init = ("abg", {"alpha": 3.5, "beta": 4.5})  # Custom ABG

# Customize other configuration parameters
config.grid.drfrac_init = 0.05                    # Initial dr/r sets the cell count
config.io.model_no = 42
config.io.base_dir = "./sims"                    # Default is the current working directory
config.scatter.sigma_m_0 = 1.0                         # Low-velocity amplitude [cm^2/g]
config.scatter.w = 50.0                                # Velocity scale [km/s]; default inf is constant scattering
config.scatter.smfp_order = 2                          # Normalized second order (default); 1 is also supported

# Switch to a truncated NFW
config.init = ("truncated_nfw", {"Zt": 0.05, "deltaP": 1e-4})

# Turn off chatter
config.io.chatter = False
```

If you don’t explicitly assign `config.io.model_no`, it is automatically set to the first available model number in `config.io.base_dir` (starting at zero) when the model number is first accessed. You can explicitly assign a `model_no` with `config = gtf.Config(io={"model_no": 5})` or with `config.io.model_no = 5` once config is instantiated.  Note that in that case, existing outputs may be overwritten. Set `config.io.overwrite = False` to reject initialization in an existing model directory.

### Scattering, calibration, and units

Set these parameters **before** constructing the state:

```python
config = gtf.Config(
    scatter={"model": "rutherford", "sigma_m_0": 10.0, "w": 50.0,
             "smfp_order": 2, "alph": 0.551, "beta": 0.842},
    num={"eps_du": 1e-4},
    stop={"t_halt": 1000.0, "rho0_halt": 1e5},
    io={"base_dir": "./sims", "overwrite": False},
)
```

`sigma_m_0` is in cm²/g and `w` is in km/s. See the
[scattering reference](docs/scattering.md) for model choices, low-speed
limits, and normalization conventions. The default `smfp_order=2` uses the
normalized second-order transport factor; `smfp_order=1` is also supported.

The current calibration defaults are `alph=0.552` and `c=0.842` (`c` is the
coefficient called beta in our calibration discussion). These are the adopted
Rutherford, second-order calibration, based on mean-core-density arrival
times at 10^2.6. They are not a separately established calibration for Møller
or first-order transport. The historical `(0.5, 0.9)` baseline can be selected
explicitly with `alph=0.5, c=0.9`; it is not automatically imposed at large `w`.
The production implementation evaluates transport moments locally as the
halo evolves.

The common reference dispersion is fixed at `0.30*v_s` for both transport
moments. With `F = sqrt(K_L_ref*K_S_ref)`, the dimensionless cross section is
`hat_sigma_eff = sigma_m_0*rho_s*r_s*F` (using consistent physical units), and
`t_s = 1/(a*sigma_m_0*rho_s*v_s*K_L_ref)`. Consequently saved `time` is
physical time divided by this new `t_s`, not the old amplitude-only time unit.
`time_Gyr` remains physical time. In the hard-sphere limit both reference
moments and `F` equal one. See the [units reference](docs/output-reference.md#units-and-normalization)
for conversions and the corresponding `state.char` attributes.

`stop.t_halt` is an **absolute** time in the new units; `run(stoptime=...)`
is an additional duration in those units.

### Reusing initial conditions

```python
from pathlib import Path

config.init = ("truncated_nfw", {"Mvir": 1e12, "cvir": 30.0, "Zt": 0.00807})
ic_path = Path("./ics/nfw_c30_Zt000807.dat")
ic_path.parent.mkdir(parents=True, exist_ok=True)
if not ic_path.exists():
    gtf.State.make_ic_file(config, str(ic_path))
state = gtf.State.from_config(config, ic_filepath=str(ic_path))
```

`Mvir` is in Msun/h. Reuse an IC file only with matching halo/profile and grid
assumptions; keep its generating configuration with it. Loading the file
starts a new run and relaxes the saved initial profile; this is not a restart.
Set scattering and output parameters independently for each new run.

### Plotting

There are three plotting functions that are imported with the `halogtf` package:

The `plot_time_evolution()` function plots the evolution of system-wide parameters over time.  It can plot any of the columns in the `time_evolution.txt` output.

```python
import halogtf as gtf

# Compare the central density evolution of two different simulations
# The default plot is rho0 (the innermost-cell density)
gtf.plot_time_evolution([state1, state2], show=True)

# Alternatively, the simulations can be called by their Config objects or by their model numbers
# v_max is in characteristic velocity units v_s
gtf.plot_time_evolution([config1, config2], quantity="v_max", show=True)

# base_dir needs to be specified if simulations are called by model number:
gtf.plot_time_evolution([5, 6], quantity="kn_cond_c", base_dir='./', show=True) # This is useful for simulations run in a different session

# The plot can be saved to a file
# Use 'show' to show the figure in standard output as well
gtf.plot_time_evolution(state1, filepath='./rho0_vs_time.png', show=True)
```

The `plot_snapshots()` function plots up to three profiles in separate panels for one or multiple snapshots of the simulation.  Snapshots are specified by the index of the `profile_x.dat` file.  Like `plot_time_evolution()`, the State object, Config object, or model number can be used to specify the simulation you wish to plot.

```python
import halogtf as gtf

# Plot the initial density profile
gtf.plot_snapshots(state, show=True)

# Plot the mass profile at a specified snapshot
gtf.plot_snapshots(config, snapshots=50, profiles='m', show=True)

# Plot the density, v^2, and Knudsen number profiles, comparing several snapshots
gtf.plot_snapshots(4, snapshots=[0, 50, 100], profiles=['rho', 'v2', 'kn_cond'], base_dir='./', show=True)

# The plot can be saved to a file
# Use 'show=True' to show the figure in standard output as well
gtf.plot_snapshots(state, filepath='./initial_rho.png')
```

The `make_movie()` function generates animations of up to three profiles in separate panels for a simulation.  Like `plot_time_evolution()`, the State object, Config object, or model number can be used to specify the simulation you wish to plot.  Only the snapshots of the most recent run for the simulation model_no will be included, even if profiles with higher indices are in the directory from previous simulation runs.  You can check the current version of the `snapshot_conversion.txt` file for all snapshots that will be included in the animation

```python
import halogtf as gtf

# Plot the density profile
gtf.make_movie(state)

# All keyword arguments available in plot_snapshots function, other than 'snapshots', can be used here
gtf.make_movie(2, base_dir='./', profiles=['v2', 'p'], grid=True)
```

The default movie profiles are `['rho', 'v2']`, with a `rho0` inset in the
first panel. Set `parallel=False` for serial rendering; `make_movie_deluxe()`
is an alias for `make_movie()`.

```python
state.make_movie(profiles=['rho', 'v2'], insets=False, parallel=False)
state.make_movie(profiles=['rho', 'v2'], insets=['rho0', None])
```

`insets=False` disables all insets, `insets=None` retains the default, and a
list selects an inset (or no inset) for each panel. Without insets or radius
annotations, a movie needs only the profile files and `snapshot_conversion.txt`.

The plotting functions also exist as methods to the `State` object:

```python
import halogtf as gtf

config = gtf.Config()
state = gtf.State.from_config(config)
state.run()

state.plot_time_evolution(show=True)         # Accepts all keyword arguments, but cannot compare between simulations when used this way

state.plot_snapshots(show=True)              # Defaults to the latest state
state.plot_snapshots(snapshots=0, filepath="./initial_profs.png")   # Plot and save initial profiles
# Note that while the standalone function defaults to the initial profile, the `State` method defaults to the current state.

state.make_movie(profiles=['rho', 'kn', 'v2'])
```

---

## Output files

Outputs live in `config.io.base_dir/ModelXXXXX/`. The complete
[output reference](docs/output-reference.md) documents every column of
`time_evolution.txt`, `profile_*.dat`, and `logfile.txt`, including units,
core definitions, transport diagnostics, and undefined boundary values.
It also describes `snapshot_conversion.txt`, `char_params.txt`, and
`model_metadata.txt`.

Profiles and histories are normally triggered by fractional changes in
**innermost-cell density** (`io.drho_prof` and `io.drho_tevol`). Initialization
and run endpoints write outputs too. `io.nlog` controls step-based logging;
`io.chatter` also sends log entries to the console. Saved `rho0` is the
innermost-cell density; saved `rho_c` is the mean density enclosed by the
half-central-density radius.

---

## Package Layout

```
halogtf/
├── config.py               # Defines 'Config' class
├── state.py                # Defines 'State' class
│
├── parameters/             # Parameter subclasses for 'Config' attributes
│   ├── char_params.py
│   ├── cosmo_params.py
│   ├── constants.py
│   ├── grid_params.py
│   ├── init_params.py
│   ├── io_params.py
│   ├── num_params.py
│   ├── scatter_params.py
│   └── stop_params.py
│
├── profiles/               # Profile specific tools to set initial conditions of 'State'
│   ├── abg.py
│   ├── nfw.py
│   ├── truncated_nfw.py
│   └── profile_routines.py
│
├── evolve/                 # Integration and solver, used by 'State' methods
│   ├── integrator.py
│   ├── transport.py
│   └── hydrostatic.py
│
├── io/                     # I/O routines
│   ├── write.py
│   └── read.py
│
├── plot/                   # Plotting routines
│   ├── time_evolution.py
│   └── snapshot.py
```

---

## License

MIT License. See [LICENSE](./LICENSE) for details.

---

## References

Listed in chronological order by publication year.

1. **Lynden-Bell, D. & Eggleton, P. P. (1980).** *On the consequences of the gravothermal catastrophe.* MNRAS, **191**, 483–498. [DOI: 10.1093/mnras/191.3.483](https://doi.org/10.1093/mnras/191.3.483).

   Develops the conducting-gas description of gravothermal collapse in star clusters that underlies the fluid approach used by HaloGTF.

2. **Balberg, S., Shapiro, S. L. & Inagaki, S. (2002).** *Self-Interacting Dark Matter Halos and the Gravothermal Catastrophe.* ApJ, **568**, 475–487. [DOI: 10.1086/339038](https://doi.org/10.1086/339038).

   Adapts the gravothermal fluid model to SIDM halos, including evolution between long- and short-mean-free-path regimes, providing the foundation for the code's conduction model.

3. **Koda, J. & Shapiro, P. R. (2011).** *Gravothermal collapse of isolated self-interacting dark matter haloes: N-body simulation versus the fluid model.* MNRAS, **415**, 1125–1137. [DOI: 10.1111/j.1365-2966.2011.18684.x](https://doi.org/10.1111/j.1365-2966.2011.18684.x).

   Compares conducting-fluid calculations with Monte Carlo N-body simulations, establishing benchmarks for calibrating the fluid conductivity and testing halo evolution.

4. **Nishikawa, H., Boddy, K. K. & Kaplinghat, M. (2020).** *Accelerated core collapse in tidally stripped self-interacting dark matter halos.* Physical Review D, **101**, 063009. [DOI: 10.1103/PhysRevD.101.063009](https://doi.org/10.1103/PhysRevD.101.063009).

   Uses the gravothermal fluid model to show how outer-halo stripping accelerates collapse.

5. **Outmezguine, N. J., Boddy, K. K., Gad-Nasr, S., Kaplinghat, M. & Sagunski, L. (2023).** *Universal gravothermal evolution of isolated self-interacting dark matter halos for velocity-dependent cross-sections.* MNRAS, **523**, 4786–4800. [DOI: 10.1093/mnras/stad1705](https://doi.org/10.1093/mnras/stad1705).

   Develops velocity-dependent transport averages and approximate evolutionary scaling relations that inform the code's thermal moments and reference-unit rescaling.

6. **Gad-Nasr, S., Boddy, K. K., Kaplinghat, M., Outmezguine, N. J. & Sagunski, L. (2024).** *On the late-time evolution of velocity-dependent self-interacting dark matter halos.* JCAP, **2024**(05), 131. [DOI: 10.1088/1475-7516/2024/05/131](https://doi.org/10.1088/1475-7516/2024/05/131).

   Studies deep core collapse with velocity-dependent scattering and supplies the thermal-moment fitting prescription and second-order SMFP treatment used by the code, with HaloGTF's normalization described above.

7. **van den Bosch, F. C. & Tokayer, Y. M. (2026).** *Not so Fast: Revised Calibration of SIDM Core Collapse.* MNRAS, stag2015. [DOI: 10.1093/mnras/stag2015](https://doi.org/10.1093/mnras/stag2015).

   Calibrates gravothermal fluid evolution against SphDSMC simulations and motivates the gentler LMFP–SMFP interpolation and constant-cross-section baseline implemented in HaloGTF.
