# Theoretical Framework: The Physics and Mathematics of Exoplanet Detection

This report details the mathematical and physical principles utilized in the extraction and characterization of exoplanetary signals from raw photometric telemetry (the process of collecting an unprocessed, continuous digital stream of light-intensity measurements from a telescope). The process is divided into two domains: **Signal Mathematics** (extracting the signal) and **Orbital Physics** (characterizing the planet).

---

## Part 1: The Mathematics of Signal Extraction

Spacecraft photometry is dominated by stochastic (unpredictable, random background fluctuation) noise, instrumental drift, and stellar variability. To extract a weak exoplanet transit, we apply rigorous Digital Signal Processing (DSP) techniques.

### 1.1 Baseline Detrending (Savitzky-Golay Filtering)

Stellar pulsation and thermal drifts manifest as low-frequency ($1/f$) colored noise. We isolate this using a Savitzky-Golay filter (a digital mathematical tool used to smooth noisy data without flattening important peaks or valleys), which performs a local polynomial regression to smooth the data without destroying high-frequency features (like the sharp edges of a transit).

For a data point $y_i$, the filter considers a symmetric window of $2m + 1$ points (from $i - m$ to $i + m$). It finds a polynomial of degree $k$ (where $k < 2m + 1$) that minimizes the least-squares error:

$$E = \sum [p(j) - y_{i+j}]^2$$

where:
$$p(j) = a_0 + a_1 j + a_2 j^2 + \dots + a_k j^k$$

By solving this over the sliding window, we extract the low-frequency stellar continuum, $F_{\text{trend}}(t)$. The pre-whitened, normalized flux is then:

$$F_{\text{norm}}(t) = \frac{F_{\text{raw}}(t)}{F_{\text{trend}}(t)}$$

### 1.2 Matched Filtering (Box Least Squares - BLS)

Because planetary transits have a very small duty cycle (the transit lasts only a tiny fraction of the total orbit), standard Fourier transforms fail to capture the energy efficiently. Instead, we use the Box Least Squares (BLS) algorithm, which acts as a matched filter for a periodic step-function.

The signal is modeled as a boxcar with two discrete states:

- **Out of transit (High state):** $H$
- **In transit (Low state):** $L$

For a trial period $P$, trial duration $\tau$, and reference epoch $t_0$, the algorithm computes the fractional transit duration $q = \tau / P$. The BLS algorithm minimizes the squared residuals between the data and this boxcar model. The detection statistic (BLS Power) is equivalent to the **Signal Residue (SR)** reduction:

$$SR = \max \left[ \frac{s^2}{q(1 - q)} \right]$$

Where $s$ is the sum of the weights of the data points falling inside the trial transit window. The highest peak in the resulting BLS periodogram defines the true orbital period $P$.

### 1.3 Coherent Integration (Phase Folding)

To boost the Signal-to-Noise Ratio (SNR) of the minuscule dip, we exploit the periodicity of the orbit. We map the linear time domain $t$ into a circular orbital phase $\phi \in [-0.5, 0.5]$ using modulo arithmetic:

$$\phi_i = \left[ \frac{t_i - t_0}{P} + 0.5 \right] \bmod 1 - 0.5$$

By stacking $N$ transits, independent Gaussian white noise $\sigma$ cancels out, amplifying the signal SNR by a factor of $\sqrt{N}$:

$$\text{SNR}_{\text{folded}} \approx \text{SNR}_{\text{single}} \times \sqrt{N}$$

---

## Part 2: The Physics of Transit Photometry

Once the periodic signal is cleanly extracted, we transition from mathematics to geometry to determine the physical size of the planet.

### 2.1 The Geometric Transit Model (Transit Depth)

When an exoplanet crosses in front of its host star, it eclipses a tiny fraction of the stellar disk. Assuming the star is a uniform sphere of radius $R_*$ and the planet is an opaque sphere of radius $R_p$, the projected two-dimensional areas they present to the telescope are $A_* = \pi R_*^2$ and $A_p = \pi R_p^2$.

The fractional drop in measured light, known as the transit depth ($\delta$), is simply the ratio of these areas:

$$\delta = \frac{\Delta F}{F} = \frac{A_p}{A_*} = \frac{\pi R_p^2}{\pi R_*^2} = \left(\frac{R_p}{R_*}\right)^2$$

By algebraically rearranging this, we can solve for the planetary radius relative to the star:

$$R_p = R_* \sqrt{\delta}$$

> **Note:** In reality, stars exhibit "limb darkening"—they are brighter in the center than at the edges. This causes the bottom of the transit dip to be slightly U-shaped rather than a perfectly flat box, but the geometric area approximation remains highly accurate.

---

## Part 3: Orbital Mechanics and Astrophysics

Finally, we use the detected period ($P$) and the known mass of the host star ($M_*$) to deduce the planet's orbital dynamics using classical Newtonian mechanics.

### 3.1 Orbital Dynamics and Kepler's Third Law

To find the distance between the planet and the star (the semi-major axis, $a$), we equate Newton's Law of Universal Gravitation with the centripetal force required to maintain a circular orbit.

**Gravitational Force:**
$$F_g = \frac{G M_* m_p}{a^2}$$

**Centripetal Force:**
$$F_c = \frac{m_p v^2}{a}$$

Where $G$ is the gravitational constant, $m_p$ is the planet's mass, and $v$ is the orbital velocity. Setting $F_g = F_c$ allows the planet's mass to cancel out:

$$\frac{G M_* m_p}{a^2} = \frac{m_p v^2}{a} \implies v^2 = \frac{G M_*}{a}$$

The velocity of an object in a stable circular orbit is its orbital circumference divided by its period $P$:

$$v = \frac{2 \pi a}{P}$$

Substituting this expression into the velocity-squared equation yields:

$$\left( \frac{2 \pi a}{P} \right)^2 = \frac{G M_*}{a}$$

$$\frac{4 \pi^2 a^2}{P^2} = \frac{G M_*}{a}$$

Solving for $a^3$ yields **Newton's derivation of Kepler's Third Law**:

$$a^3 = \frac{G M_* P^2}{4 \pi^2}$$

By taking the cube root, we successfully calculate the semi-major axis $a$ purely from the timing of the light dips:

$$a = \left( \frac{G M_* P^2}{4 \pi^2} \right)^{1/3}$$

### 3.2 Mean Orbital Velocity

With the semi-major axis ($a$) and the orbital period ($P$) now mathematically proven, calculating the physical speed of the planet through space is a straightforward kinematic equation (assuming an eccentricity $e \approx 0$ for a circular orbit):

$$v_{\text{orbit}} = \frac{2 \pi a}{P}$$

For extreme systems like **Kepler-10b**, this results in velocities exceeding 200 km/s.

---

## Summary of the Theoretical Flow

1. **$t \rightarrow F(t)$**: The telescope measures photons over time.
2. **$F(t) \rightarrow P, t_0, \delta$**: DSP mathematics (Savitzky-Golay, BLS, Phase-Folding) extracts the orbital period, epoch, and transit depth.
3. **$\delta \rightarrow R_p$**: Photometric geometry translates the depth into the planet's physical radius.
4. **$P \rightarrow a, v$**: Newtonian physics translates the period into the orbital distance and velocity.
