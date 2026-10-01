# Exoplanet Transit Detection via Digital Signal Processing

## 1. Project Abstract

This project implements an end-to-end digital signal processing (DSP) pipeline to detect and physically characterize exoplanets from noisy spacecraft telemetry. By treating astronomical photometry as a classic signal-to-noise extraction problem, the pipeline identifies periodic transit dips buried in 1/f colored noise.

Using NASA's Kepler Space Telescope data for the target system **Kepler-10**, the pipeline successfully performs pre-whitening, matched filtering, and coherent stacking to extract the 0.016% signal drop caused by the rocky exoplanet **Kepler-10b**. Finally, the pipeline derives the planet's physical metrics using Kepler's laws of planetary motion.

## 2. Technical Stack

- **Language:** Python 3.9+
- **Astronomy Libraries:** `lightkurve` (NASA MAST API interface), `astropy` (time and physical constants)
- **DSP & Math Libraries:** `numpy`, `scipy` (Savitzky-Golay filtering, robust statistics)
- **Visualization:** `matplotlib` (Publication-quality dashboards)

---

## 3. The DSP Pipeline Architecture

Detecting an Earth-sized exoplanet is fundamentally a problem of detecting a weak, low-duty-cycle periodic pulse train in non-white noise. The pipeline is divided into 6 distinct stages.

### Step 1: Data Acquisition & Telemetry Parsing

The pipeline queries the **Mikulski Archive for Space Telescopes (MAST)** to download Quarter 3 observations of Kepler-10.

- **Cadence:** We use Long Cadence (30-minute exposure) data, yielding a continuous ~90-day baseline.
- **Signal Type:** We utilize **PDCSAP Flux** (Pre-search Data Conditioning Simple Aperture Photometry). This is the raw photon count that has been partially corrected for spacecraft pointing jitter and thermal focus drift using cotrending basis vectors.

### Step 2: Data Cleaning & Robust Masking

Spacecraft data is contaminated by cosmic ray hits, thruster firings, and data link dropouts.

- **NaN & Quality Filtering:** We drop missing cadences and filter by `quality == 0` to discard known hardware anomalies (e.g., reaction wheel desaturations).
- **Asymmetric Sigma-Clipping:** Standard standard deviation (σ) is heavily skewed by outliers. We instead compute the **Median Absolute Deviation (MAD)** to find the robust statistical spread.
  - _The EE Concept:_ Cosmic rays cause artificial _positive_ charge spikes, while exoplanets cause physical _negative_ dips. We apply an asymmetric threshold (+4σ upper, -8σ lower). This rigorously scrubs cosmic rays without accidentally amputating the true exoplanet signal.

### Step 3: Baseline Detrending (Pre-Whitening)

Stars are not static; they pulsate, have magnetic starspot cycles, and the spacecraft undergoes slow thermal expansion. This manifests as high-amplitude 1/f (colored) baseline drift.

- **The Filter:** We apply a **Savitzky-Golay polynomial filter** (window size ~3 days).
- **The DSP Concept:** This acts as a specialized low-pass filter. The window is deliberately chosen to be much wider than the planetary transit duration (~1.6 hours). Therefore, the polynomial perfectly tracks the slow stellar drift but glides right over the rapid transit dips. Dividing the raw flux by this low-pass trend yields a flattened, pre-whitened residual signal centered strictly at 1.0.

### Step 4: Matched Filtering via Box Least Squares (BLS)

Standard Fourier Transforms (FFT) fail here because a transit has a very small duty cycle (~1%). Spreading that square pulse over hundreds of sinusoidal harmonics destroys the Signal-to-Noise Ratio (SNR).

- **The BLS Algorithm:** We use the Box Least Squares periodogram, which acts as a **Matched Filter**.
- **Execution:** We define a mathematical "boxcar" template (a U-shaped dip). The algorithm sweeps across a massive grid of trial orbital periods (P), transit durations (τ), and reference epochs (t0). It calculates the cross-correlation power for each combination. The highest peak in the BLS spectrum reveals the true orbital period of the planet.

### Step 5: Coherent Phase Folding & Stacking

A single transit of Kepler-10b drops the light by only ~160 parts per million (ppm), which is almost invisible to the naked eye against the instrument noise floor.

- **Phase Wrapping:** Using modulo arithmetic (φ = [(t - t0) / P] mod 1), we fold the linear time series into a circular orbital phase.
- **Synchronous Detection:** By stacking all ~100 individual transits on top of each other and averaging them into 10-minute phase bins, uncorrelated Gaussian white noise cancels itself out. The SNR is boosted by √N (where N is the number of transits), revealing a pristine, high-fidelity transit profile.

---

## 4. Astrophysical Parameter Extraction (Step 6)

From the DSP pipeline, we extract three critical empirical observables:

1. **Period (P):** 0.8375 days
2. **Transit Depth (δ):** 162.0 ppm (0.0162%)
3. **Transit Duration (τ):** 1.66 hours

Using these observables, the pipeline mathematically derives the physical reality of the alien solar system:

#### A. Planet-to-Star Radius Ratio (Rp / R\*)

The fraction of light blocked (δ) is directly proportional to the area of the planet's disk relative to the star's disk.
`δ = (π * Rp^2) / (π * R*^2)` => `Rp = R* * √δ`
_Result:_ Kepler-10b has a radius of **1.48 Earth Radii**. (NASA literature confirms 1.47 ± 0.03 R⊕).

#### B. Semi-Major Axis (a)

Using Newton's derivation of **Kepler's Third Law**, we calculate the distance between the planet and its star using the orbital period and the star's mass (M*):
`a = ((G * M\* _ P^2) / (4 _ π^2))^(1/3)`
_Result:_ The planet orbits at **0.0169 AU** (~2.5 million km). This is more than 20 times closer to its star than Mercury is to our Sun, categorizing it as an ultra-short-period "Lava World".

#### C. Mean Orbital Velocity (v)

Assuming a circular orbit, we calculate the planet's velocity:
`v = (2 * π * a) / P`
_Result:_ The planet travels at a blistering **218.9 km/s** (roughly 490,000 mph).

---

## 5. Conclusion & Project Value

This project demonstrates the powerful intersection of digital signal processing and modern astrophysics. By systematically addressing non-Gaussian outliers, 1/f baseline drift, and low-SNR pulse detection, the pipeline successfully uncovers a rocky world hundreds of light-years away using nothing but open-source Python toolkits and fundamental math.
