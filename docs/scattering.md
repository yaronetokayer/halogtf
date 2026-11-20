# Scattering models and normalization

HaloGTF supports three elastic-scattering choices:

| Scattering prescription | Configuration |
| --- | --- |
| Constant, isotropic ("hard-sphere") scattering | Either model with `w=float("inf")` |
| Screened Rutherford scattering | `model="rutherford"` with finite `w` |
| Møller scattering | `model="moller"` with finite `w` |

## Meaning of the input parameters

The parameter `sigma_m_0` specifies $\sigma_0/m$ in **cm²/g**, and `w`
specifies the characteristic velocity scale in **km/s**. "Rutherford" here
means the screened, Yukawa-like prescription with a finite low-speed limit,
rather than unscreened Coulomb scattering.

**Both velocity-dependent prescriptions use a normalization chosen to recover
hard-sphere scattering with amplitude $\sigma_0/m$ at low velocities.** Any
model-dependent constant prefactor in the underlying scattering prescription
is absorbed into this definition. Consequently, the same `sigma_m_0` has the
same low-speed meaning for Rutherford, Møller, and hard-sphere scattering;
users should not multiply it by an additional normalization factor when
switching models.

## Low-speed and hard-sphere limits

The fluid solver incorporates velocity dependence through normalized thermal
transport moments. These satisfy

$$
K_p(v/w)\longrightarrow 1
\qquad\text{as}\qquad v/w\longrightarrow 0,
$$

where $v$ is the **single-particle, one-dimensional velocity dispersion**.
Thus, setting `w=float("inf")` recovers the constant-scattering model exactly,
while finite `w` approaches that model when the local dispersion is much
smaller than `w`. At finite $v/w$, Rutherford and Møller generally predict
different transport rates. No extra factor of $\sqrt{2}$ should be inserted
into the moment argument: the thermal averaging already accounts for
relative particle velocities.

Set the scattering parameters before constructing a `State`:

```python
import halogtf as gtf

config = gtf.Config()

# Hard-sphere elastic scattering with sigma_0/m = 10 cm^2/g
config.scatter.model = "rutherford"  # "moller" gives the same result at w=inf
config.scatter.sigma_m_0 = 10.0
config.scatter.w = float("inf")

# Velocity-dependent screened Rutherford scattering
config.scatter.w = 50.0  # km/s

# Velocity-dependent Møller scattering, with the same low-speed amplitude
config.scatter.model = "moller"
```

The assignments above illustrate switching prescriptions; constructing a
state after the final line selects Møller scattering with `w=50.0` km/s.

## Viscosity cross section versus input amplitude

The normalization convention should be distinguished from the definition of
the **viscosity cross section**,

$$
\sigma_V=\int \sin^2\theta\,d\sigma.
$$

For isotropic scattering with total cross section $\sigma_0$, this angular
integral gives $\sigma_V=2\sigma_0/3$; its isotropic-equivalent value is
therefore $3\sigma_V/2=\sigma_0$. HaloGTF's input amplitude and normalized
transport moments account for this distinction. The thermal moments follow
the viscosity-weighted convention described by
[Gad-Nasr et al. (2024)](https://doi.org/10.1088/1475-7516/2024/05/131).

## SMFP normalization and simulation units

The second-order SMFP factor is also normalized to approach one in the
hard-sphere limit, preserving the meaning of the conductivity coefficient
`b`. Its $45/44$ normalization is already included internally and requires
**no adjustment to `sigma_m_0`**. First-order SMFP transport is also available
through `config.scatter.smfp_order=1`.

This physical amplitude convention is separate from the reference-moment
rescaling of dimensionless cross section and time. `sigma_m_0` always remains
a physical input in cm²/g; the code performs the unit conversion internally.
See the [output and units reference](output-reference.md#units-and-normalization)
for the conversion formulas, and the
[README](../README.md#scattering-calibration-and-units) for the adopted
calibration and its scope.
