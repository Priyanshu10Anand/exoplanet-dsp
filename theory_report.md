# Theoretical Framework: The Physics and Mathematics of Exoplanet Detection

This report details the mathematical and physical principles utilized in the extraction and characterization of exoplanetary signals from raw photometric telemetry. The process is divided into two domains: **Signal Mathematics** (extracting the signal) and **Orbital Physics** (characterizing the planet).

---

## Part 1: The Mathematics of Signal Extraction

Spacecraft photometry is dominated by stochastic noise, instrumental drift, and stellar variability. To extract a weak exoplanet transit, we apply rigorous Digital Signal Processing (DSP) techniques.

### 1.1 Baseline Detrending (Savitzky-Golay Filtering)

Stellar pulsation and thermal drifts manifest as low-frequency (1/f) colored noise. We isolate this using a Savitzky-Golay filter, which performs a local polynomial regression to smooth the data without destroying high-frequency features (like the sharp edges of a transit).

For a data point `y_i`, the filter considers a symmetric window of `2m+1` points (from `i-m` to `i+m`). It finds a polynomial of degree `k` (where `k < 2m+1`) that minimizes the least-squares error:
`E = Σ [p(j) - y_{i+j}]^2`
Where `p(j) = a_0 + a_1*j + a_2*j^2 + ... + a_k*j^k`.
By solving this over the sliding window, we extract the low-frequency stellar continuum, `F_trend(t)`. The pre-whitened, normalized flux is then:
`F_norm(t) = F_raw(t) / F_trend(t)`

### 1.2 Matched Filtering (Box Least Squares - BLS)

Because planetary transits have a very small duty cycle (the transit lasts only a tiny fraction of the total orbit), standard Fourier transforms fail to capture the energy efficiently. Instead, we use the Box Least Squares (BLS) algorithm, which acts as a matched filter for a periodic step-function.

The signal is modeled as a boxcar with two discrete states:

- **Out of transit (High state):** H
- **In transit (Low state):** L

For a trial period `P`, trial duration `τ`, and reference epoch `t0`, the algorithm computes the fractional transit duration `q = τ / P`. The BLS algorithm minimizes the squared residuals between the data and this boxcar model. The detection statistic (BLS Power) is equivalent to the **Signal Residue (SR)** reduction:
`SR = max[ s^2 / (q * (1 - q)) ]`
Where `s` is the sum of the weights of the data points falling inside the trial transit window. The highest peak in the resulting BLS periodogram defines the true orbital period `P`.

### 1.3 Coherent Integration (Phase Folding)

To boost the Signal-to-Noise Ratio (SNR) of the minuscule dip, we exploit the periodicity of the orbit. We map the linear time domain `t` into a circular orbital phase `φ ∈ [-0.5, 0.5]` using modulo arithmetic:
`φ_i = [ (t_i - t0) / P + 0.5 ] mod 1 - 0.5`
By stacking `N` transits, independent Gaussian white noise `σ` cancels out, amplifying the signal SNR by a factor of `√N`:
`SNR_folded ≈ SNR_single * √N`

---

## Part 2: The Physics of Transit Photometry

Once the periodic signal is cleanly extracted, we transition from mathematics to geometry to determine the physical size of the planet.

### 2.1 The Geometric Transit Model (Transit Depth)

When an exoplanet crosses in front of its host star, it eclipses a tiny fraction of the stellar disk.
Assuming the star is a uniform sphere of radius `R*` and the planet is an opaque sphere of radius `Rp`, the projected two-dimensional areas they present to the telescope are `A* = π * R*^2` and `Ap = π * Rp^2`.

The fractional drop in measured light, known as the transit depth (`δ`), is simply the ratio of these areas:
`δ = ΔF / F = Ap / A* = (π * Rp^2) / (π * R*^2) = (Rp / R*)^2`
By algebraically rearranging this, we can solve for the planetary radius relative to the star:
`Rp = R* * √δ`
_(Note: In reality, stars exhibit "limb darkening"—they are brighter in the center than at the edges. This causes the bottom of the transit dip to be slightly U-shaped rather than a perfectly flat box, but the geometric area approximation remains highly accurate)._

---

## Part 3: Orbital Mechanics and Astrophysics

Finally, we use the detected period (`P`) and the known mass of the host star (`M*`) to deduce the planet's orbital dynamics using classical Newtonian mechanics.

### 3.1 Kepler's Third Law of Planetary Motion

To find the distance between the planet and the star (the semi-major axis, `a`), we equate Newton's Law of Universal Gravitation with the centripetal force required to maintain a circular orbit.

Gravitational Force:
`Fg = (G * M* * mp) / a^2`
Centripetal Force:
`Fc = (mp * v^2) / a`
Where `G` is the gravitational constant, `mp` is the planet's mass, and `v` is the orbital velocity. Setting `Fg = Fc`:
`(G * M* * mp) / a^2 = (mp * v^2) / a => v^2 = (G * M*) / a`
We know that the velocity of an object in a circular orbit is the circumference divided by the period:
`v = (2 * π * a) / P`
Substituting this into the velocity-squared equation gives:
`((2 * π * a) / P)^2 = (G * M*) / a`
`(4 * π^2 * a^2) / P^2 = (G * M*) / a`
Solving for `a^3` yields **Newton's derivation of Kepler's Third Law**:
`a^3 = (G * M* * P^2) / (4 * π^2)`
By taking the cube root, we successfully calculate the distance from the planet to the star purely from the timing of the light dips.

### 3.2 Mean Orbital Velocity

With the semi-major axis (`a`) and the orbital period (`P`) now mathematically proven, calculating the physical speed of the planet through space is a straightforward kinematic equation (assuming an eccentricity `e ≈ 0` for a circular orbit):
`v_orbit = (2 * π * a) / P`
For extreme systems like Kepler-10b, this results in velocities exceeding 200 km/s.

---

## Summary of the Theoretical Flow

1. **t -> F(t)**: The telescope measures photons over time.
2. **F(t) -> P, t0, δ**: DSP mathematics (Savitzky-Golay, BLS, Phase-Folding) extracts the orbital period, epoch, and transit depth.
3. **δ -> Rp**: Photometric geometry translates the depth into the planet's physical radius.
4. **P -> a, v**: Newtonian physics translates the period into the orbital distance and velocity.
