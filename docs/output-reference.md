# Output reference

This reference describes the current velocity-dependent HaloGTF output schema.
Column numbers below are **one-based**. Read columns by header name when
possible: older runs can have different schemas. The first line is an
uncommented header; subsequent lines contain whitespace-separated values.
Log headers have spaces inside angle brackets, so they require special care.

## Units and normalization

Unless explicitly indicated, output values are dimensionless:

| Quantity | Multiply by to obtain physical units |
| --- | --- |
| Radius, mean free path, temperature scale length | `state.char.r_s` [Mpc] |
| Enclosed mass | `state.char.m_s` [Msun] |
| Density | `state.char.rho_s` [Msun/Mpc³] |
| 1D velocity dispersion | `state.char.v_s` [km/s] |
| Squared dispersion (`v2`) | `state.char.v_s**2` [(km/s)²] |
| Time, timestep, core crossing time | `state.char.t_s` [seconds] |

The definitions are `rho_s = m_s/(4*pi*r_s**3)` and
`v_s**2 = G*m_s/r_s`. Therefore dimensionless mean density inside radius `r`
is `3*m/r**3`, with no extra 4*pi. Pressure is `rho*v2` in units of
`rho_s*v_s**2`; specific internal energy is `1.5*v2` in units of `v_s**2`.

Let `hat_sigma_0 = sigma_m_0*rho_s*r_s`, with physical factors converted to
consistent units. The common reference 1D dispersion is `v_ref = 0.33*v_s`.
`K_L_ref` and `K_S_ref` are the selected model/order's moments there, and
`F = sqrt(K_L_ref*K_S_ref)`. The code uses:

- `sigma_m_s_hardsphere = 1/(rho_s*r_s)` in cm²/g.
- `sigma_m_s = sigma_m_s_hardsphere/F` in cm²/g.
- `sigma_m_0_char_hardsphere = hat_sigma_0`.
- `sigma_m_0_char = hat_sigma_eff = hat_sigma_0*F`.
- `t_s = 1/(a*sigma_m_0*rho_s*v_s*K_L_ref)` in seconds.
- `dynamical_time_factor = (r_s/v_s)/t_s = a*hat_sigma_0*K_L_ref`.

The unit-defining moments stay fixed during a run. Local moments still evolve
with temperature; transport uses `K_L/K_L_ref` and `K_S/K_S_ref`. The saved
moment columns contain the local factors themselves, **not** these ratios.
`K_L = K5`. First-order `K_S = K5`; second-order `K_S` is the combination of
K3, K5, K7 and K9 returned by `factors`, normalized by 45/44 to equal one at
constant scattering. Both scattering models tend to hard-sphere scattering
with amplitude `sigma_m_0` as `w` tends to infinity. No additional 44/45 or
sqrt(2) conversion is needed.

In the old amplitude-only time unit `t0 = 1/(a*sigma_m_0*rho_s*v_s)`, a given
physical time had value `time_old`. Its new value is
`time = K_L_ref*time_old`. Use `time_Gyr` for comparisons across runs without
manually converting time units.

## Core definitions

`rho0` is the density of the innermost cell. `r_c` is the first radius where
interpolated density falls to `rho0/2` (local power-law interpolation between
cell midpoints). **`rho_c` is the mean enclosed density `3*m_c/r_c**3`, not
`rho0/2`.** The latter is the density at the core boundary.

`r_m2` is the first downward crossing of logarithmic density slope -2,
interpolating the slopes in log radius. Both radius finders return the
outermost edge if no crossing is found. Enclosed masses use power-law
interpolation of the edge mass array. `v2_c` and `v2_m2` are arithmetic,
mass-weighted means of the squared 1D dispersion, including partial shells.
`kn_c` and `kn_cond_c` are mass-weighted **geometric** means within `r_c`.

## `time_evolution.txt`

Density-triggered rows use `io.drho_tevol`; profile snapshots and run endpoints
also write rows. `io.t_evol=False` disables only density-triggered history
output. Successive `run()` calls retain earlier rows and regenerate the final
slope columns.

| Column | Name | Meaning |
| --- | --- | --- |
| 1 | `step` | Cumulative accepted integration-step count. |
| 2 | `dt` | Most recent timestep, in `t_s` units (initial value is a seed). |
| 3 | `n` | Current number of radial cells. |
| 4 | `time` | Time since initialization, in `t_s` units. |
| 5 | `time_Gyr` | Same time in physical Gyr. |
| 6 | `rho0` | Innermost-cell density, in `rho_s` units. |
| 7 | `v_max` | Maximum **1D dispersion** over cells, in `v_s` units; not circular velocity. |
| 8 | `kn_c` | Core geometric mean of amplitude-reference Knudsen number. |
| 9 | `kn_cond_c` | Core geometric mean of conductivity-effective Knudsen number. |
| 10 | `r_c` | Half-central-density radius, in `r_s` units. |
| 11 | `rho_c` | Mean density inside `r_c`, in `rho_s` units. |
| 12 | `m_c` | Enclosed mass at `r_c`, in `m_s` units. |
| 13 | `v2_c` | Mass-weighted squared 1D dispersion inside `r_c`, in `v_s²` units. |
| 14 | `x_c` | `sqrt(v2_c)/w_char`. |
| 15 | `K_L_c` | LMFP factor evaluated at `v2_c`. |
| 16 | `K_S_c` | Normalized SMFP factor evaluated at `v2_c`. |
| 17 | `r_m2` | Density-slope -2 radius, in `r_s` units. |
| 18 | `rho_m2` | Mean enclosed density `3*m_m2/r_m2**3`. |
| 19 | `m_m2` | Mass enclosed within `r_m2`. |
| 20 | `v2_m2` | Mass-weighted squared 1D dispersion inside `r_m2`. |
| 21 | `x_m2` | `sqrt(v2_m2)/w_char`. |
| 22 | `K_L_m2` | LMFP factor evaluated at `v2_m2`. |
| 23 | `K_S_m2` | Normalized SMFP factor evaluated at `v2_m2`. |
| 24 | `drfrac_max` | Maximum shell width/geometric-midpoint radius, excluding the central cell. |
| 25 | `tsc_c` | `r_c/sqrt(v2_c)*dynamical_time_factor`, in `t_s` units. |
| 26 | `dlnmc_dlnvc` | Local fitted `d ln(m_c)/d ln(v_c)`, where `v_c=sqrt(v2_c)`. |
| 27 | `dlnrhocdlnvc` | Local fitted `d ln(rho_c)/d ln(v_c)`. |
| 28 | `zeta_balb` | `1 + d ln(v2_c)/d ln(m_c)`, fitted independently. |

The core moment values are evaluated at the mean dispersion; they are not
averages of local K values. Columns 26–28 are appended on run finalization,
using local log-log regressions over up to 31 saved rows, shifted at the
endpoints. A degenerate fit returns 0 for the first two slopes and 1 for
`zeta_balb`; these fallbacks are not physical measurements. Sparse histories
and turning points can give unreliable slopes.

## `profile_<index>.dat`

Each row describes one cell, indexed from zero. There are `n` rows; the
origin edge is omitted. `log_r` and `m` refer to that cell's **outer edge**;
cell quantities use `log_rmid`. The `_e` conductivity columns also refer to
that outer edge, and the outermost boundary entries are NaN.

| Column | Name | Meaning |
| --- | --- | --- |
| 1 | `i` | Zero-based cell index. |
| 2 | `log_r` | Base-10 log of outer-edge radius in `r_s` units. |
| 3 | `log_rmid` | Base-10 log of cell-midpoint radius in `r_s` units. |
| 4 | `m` | Enclosed mass at outer edge, in `m_s` units. |
| 5 | `rho` | Cell density, in `rho_s` units. |
| 6 | `v2` | Cell squared 1D dispersion, in `v_s²` units. |
| 7 | `kn` | `1/(hat_sigma_0*sqrt(rho*v2))`. |
| 8 | `ltemp` | Temperature-gradient length `v2/abs(dv2/dr)`, in `r_s` units. |
| 9 | `mfp` | Amplitude-reference mean free path `1/(hat_sigma_0*rho)`, in `r_s` units. |
| 10 | `drfrac` | Shell width/geometric-midpoint radius; central cell NaN. |
| 11 | `drltemp` | Shell width divided by `ltemp`; central cell NaN. |
| 12 | `mfpltemp` | `mfp/ltemp`. |
| 13 | `k_sc` | SMFP radial conduction coefficient at cell midpoint. |
| 14 | `k_lc` | LMFP radial conduction coefficient at cell midpoint. |
| 15 | `k_totc` | Interpolated total radial conduction coefficient at cell midpoint. |
| 16 | `k_se` | SMFP radial conduction coefficient at outer cell face. |
| 17 | `k_le` | LMFP radial conduction coefficient at outer cell face. |
| 18 | `k_tote` | Interpolated total radial conduction coefficient at outer cell face. |
| 19 | `krat_c` | `k_sc/k_lc`. |
| 20 | `krat_e` | `k_se/k_le`. |
| 21 | `dttcool` | Most recent `dt/t_cool`. |
| 22 | `tdyntcool` | Local dynamical time divided by `t_cool`. |
| 23 | `s` | Entropy proxy `1.5*ln(v2)-ln(rho)` using dimensionless quantities. |
| 24 | `dsdr` | Radial derivative of `s` with respect to `r/r_s`. |
| 25 | `dlnrhodlnp` | Radial logarithmic density/pressure slope ratio, with `p=rho*v2`. |
| 26 | `kn_cond` | `kn/sqrt(K_L*K_S)`. |
| 27 | `mfp_cond` | `mfp/sqrt(K_L*K_S)`, in `r_s` units. |
| 28 | `x` | Local `sqrt(v2)/w_char`; zero for infinite `w`. |
| 29 | `K_L` | Local LMFP transport factor. |
| 30 | `K_S` | Local normalized SMFP transport factor for the selected order. |
| 31 | `mfp_cond_ltemp` | `mfp_cond/ltemp`. |

The `kn`/`mfp` columns use the low-speed amplitude, without velocity
averaging. The conductivity-effective columns combine the two transport
moments; they are not literal collision mean free paths. In particular,
`k_S/k_L = b/(a*c)*kn_cond**2` at the same location.

The `k_*` columns are dimensionless **radial flux coefficients**, including
the factor `1.5*r**2`, rather than dimensional material conductivities.
Writing `T=v2`, `L=K_L/K_L_ref` and `S=K_S/K_S_ref`, they are
`k_L = 1.5*r**2*sqrt(T)*c*rho*T*L` and
`k_S = 1.5*r**2*sqrt(T)*b/(a*hat_sigma_eff**2*S)`.
The total is `(k_L**(-alph)+k_S**(-alph))**(-1/alph)`.
Faces use arithmetic-mean neighboring temperatures and the solver's
interpolated density. Boundary conditions set the outer flux separately.

`t_cool` measures the absolute fractional conduction change in squared
dispersion from the most recent conduction solve, not a signed radiative
cooling time. `t_dyn = dynamical_time_factor/sqrt(rho)`.
Initial timescale diagnostics can be placeholders. On ordinary grids, the first
two cells have undefined `ltemp`, `dsdr`, and `dlnrhodlnp`. Undefined derivatives,
central geometric ratios and outer face coefficients can be NaN; flat
temperature gradients can give infinite `ltemp`. The timescale writer returns
zero when its denominator is zero. Do not interpret these special values as
resolved physical measurements.

## `logfile.txt` and console log

The same columns are printed to the console when `io.chatter=True`. Logging
occurs every `io.nlog` accepted steps and at run endpoints. Angle-bracket
quantities average the accepted steps since the previous log entry; they are
`N/A` when no steps have accumulated. Counters reset after each entry.

| Column | Name | Meaning |
| --- | --- | --- |
| 1 | `step` | Cumulative accepted step count. |
| 2 | `time` | Current time in `t_s` units. |
| 3 | `<dt>` | Mean accepted timestep in `t_s` units. |
| 4 | `n` | Current number of cells. |
| 5 | `rho0` | Current innermost-cell density in `rho_s` units. |
| 6 | `v_max` | Current maximum 1D dispersion in `v_s` units. |
| 7 | `kn_c` | Current core geometric mean of `kn`. |
| 8 | `kn_cond_c` | Current core geometric mean of `kn_cond`. |
| 9 | `eps_du_eff` | Current energy-change tolerance, including the low-Knudsen boost. |
| 10 | `<du lim>` | Mean of each accepted step's maximum fractional conduction change divided by its effective tolerance. |
| 11 | `<dr lim>` | Mean final hydrostatic fractional-radius correction divided by `num.eps_dr`; not total shell displacement. |
| 12 | `<n_retry_du>` | Mean conduction retries beyond the first trial per accepted step. |
| 13 | `<n_iter_dr>` | Mean additional hydrostatic solves beyond the first per accepted step. |
| 14 | `n_split` | Number of split operations since the previous log entry. |
| 15 | `n_merge` | Number of merge operations since the previous log entry. |

`eps_du_eff` is evaluated at the current state, whereas `<du lim>` averages
the actual tolerance ratios used on preceding steps. It is not the ratio of
the average change to the currently printed tolerance. These columns diagnose
solver behavior, not physical observables. Log column widths are configurable
in `write_log_entry()`; names and values are never truncated.

To change the layout, edit the `columns` list in `halogtf/io/write.py`:
each entry is `(header, value, minimum_width)`, for example
`('time', state.t, 13)`. The third value sets the minimum width, which expands
when required by the header or value; floating-point values retain six decimal
places in scientific notation. This affects both file and console output.

## Other files and reading data

`snapshot_conversion.txt` has three columns: `index` (snapshot suffix),
`time` (in `t_s` units), and `step` (accepted-step count). Use this mapping to
select snapshots, rather than assuming equally spaced times. Initial and
run-endpoint profiles are written even when `io.profiles=False`.

`model_metadata.txt` records the configuration groups, including `scatter`,
`num`, and `stop`. It is configuration metadata, not a restart checkpoint.
`char_params.txt` records characteristic scales and reference factors by name:
`r_s`, `m_s`, `rho_s`, `v_s`, `t_s`, the two cross-section scales and amplitudes
listed above, `w_char=w/v_s`, `k_l_ref`, `k_s_ref`,
`v_ref_over_vs=0.33`, and `dynamical_time_factor`. It also contains `fc`
(the NFW mass normalization), `chi` (the ABG mass normalization, when used),
and `mtot_m200` (truncated-NFW total/virial mass ratio, when available).
Unavailable values are NaN. Long names in this file's fixed-width header may
abut; do not assume a simple whitespace header parser can recover all names.
Use `state.char` directly while the state is available.

```python
from pathlib import Path
import numpy as np
from halogtf.io.read import load_snapshot_bundle

model_dir = Path(state.config.io.base_dir) / state.config.io.model_dir
history = np.atleast_1d(np.genfromtxt(model_dir / 'time_evolution.txt', names=True))
core_density = history['rho_c']
physical_time_Gyr = history['time_Gyr']
snapshot = load_snapshot_bundle(str(model_dir))  # latest saved snapshot
```

For `logfile.txt`, use `np.genfromtxt(path, skip_header=1)` for the numerical
body (`N/A` becomes NaN), then assign the 15 names above manually. Whitespace
splitting the header would incorrectly split `<du lim>` and `<dr lim>`.
`State.from_dir()` is not implemented: reading outputs does not resume a run.
