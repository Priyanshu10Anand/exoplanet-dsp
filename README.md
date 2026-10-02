# Exoplanet Transit Detection via Digital Signal Processing

## 1. Project Abstract

This project implements an end-to-end **digital signal processing (DSP) pipeline** to detect and physically characterize exoplanets from noisy spacecraft telemetry. By treating astronomical photometry as a classic signal-to-noise extraction problem, the pipeline identifies periodic transit dips buried in $1/f$ colored noise.

Using NASA's **Kepler Space Telescope** data for the target system **Kepler-10**, the pipeline successfully performs **pre-whitening, matched filtering, and coherent stacking** to extract the $0.016\%$ signal drop caused by the rocky exoplanet **Kepler-10b**. Finally, the pipeline derives the planet's physical metrics using Kepler's laws of planetary motion.

---

## 2. Technical Stack

| Category                | Technology                                                                      |
| :---------------------- | :------------------------------------------------------------------------------ |
| **Language**            | Python 3.9+                                                                     |
| **Astronomy Libraries** | `lightkurve` (NASA MAST API interface), `astropy` (time and physical constants) |
| **DSP & Mathematics**   | `numpy`, `scipy` (Savitzky-Golay filtering, statistics)                         |
| **Visualization**       | `matplotlib` (dashboards)                                                       |

---

## 3. DSP Pipeline Architecture

Detecting an Earth-sized exoplanet is fundamentally a problem of detecting a **weak, low-duty-cycle periodic pulse train in non-white noise**.

The pipeline is divided into **six distinct stages**:

```text
            Kepler Telemetry
                   │
                   ▼
      ┌──────────────────────────┐
      │ 1. Data Acquisition      │
      │    & Telemetry Parsing   │
      └────────────┬─────────────┘
                   ▼
      ┌──────────────────────────┐
      │ 2. Data Cleaning         │
      │    & Robust Masking      │
      └────────────┬─────────────┘
                   ▼
      ┌──────────────────────────┐
      │ 3. Baseline Detrending   │
      │    / Pre-Whitening       │
      └────────────┬─────────────┘
                   ▼
      ┌──────────────────────────┐
      │ 4. Matched Filtering     │
      │    via BLS               │
      └────────────┬─────────────┘
                   ▼
      ┌──────────────────────────┐
      │ 5. Phase Folding         │
      │    & Coherent Stacking   │
      └────────────┬─────────────┘
                   ▼
      ┌──────────────────────────┐
      │ 6. Astrophysical         │
      │    Parameter Extraction  │
      └──────────────────────────┘
```

---

## 4. Step 1: Data Acquisition & Telemetry Parsing

The pipeline queries the **Mikulski Archive for Space Telescopes (MAST)** to download observations of Kepler-10.

### Dataset Characteristics

- **Target:** Kepler-10
- **Mission:** NASA Kepler Space Telescope
- **Quarter:** Quarter 3
- **Cadence:** Long Cadence ($\sim 30$-minute exposure)
- **Observation baseline:** Approximately 90 days
- **Signal:** PDCSAP Flux

**PDCSAP Flux (Pre-search Data Conditioning Simple Aperture Photometry)** is used as the input photometric signal. It is already partially corrected for instrumental systematics such as spacecraft pointing jitter and thermal focus drift using cotrending basis vectors.

---

## 5. Step 2: Data Cleaning & Robust Masking

Spacecraft photometric data can contain several types of artifacts, including:

- Cosmic-ray hits
- Thruster firings
- Data dropouts
- Instrumental anomalies
- Reaction-wheel desaturation events

### 5.1 NaN & Quality Filtering

Missing cadences are removed and only measurements satisfying:

```python
quality == 0
```

are retained. This removes cadences associated with known spacecraft or instrumental anomalies.

### 5.2 Asymmetric Sigma Clipping

Standard deviation-based outlier rejection can be strongly affected by extreme outliers. Instead, the pipeline uses the **Median Absolute Deviation (MAD)** as a robust estimator of statistical spread.

The key observation is that different artifacts have different signatures:

- **Cosmic rays:** typically produce artificial positive flux spikes.
- **Exoplanet transits:** produce physical negative flux dips.

Therefore, an asymmetric clipping strategy is applied:

| Direction | Threshold  | Purpose                               |
| :-------- | :--------- | :------------------------------------ |
| Positive  | $+4\sigma$ | Aggressively remove cosmic-ray spikes |
| Negative  | $-8\sigma$ | Preserve genuine transit dips         |

This reduces the risk of accidentally removing the planetary signal during outlier rejection.

---

## 6. Step 3: Baseline Detrending / Pre-Whitening

Stars can exhibit intrinsic variability due to pulsations, magnetic activity and starspots. In addition, spacecraft thermal effects and other instrumental systematics can introduce slow baseline variations.

These effects manifest as high-amplitude, low-frequency drift, commonly modeled as **colored $1/f$ noise**.

### 6.1 Savitzky-Golay Filtering

The pipeline uses a **Savitzky-Golay polynomial filter** with a window of approximately **3 days**.

The window is deliberately much wider than the planetary transit duration of approximately **1.6 hours**. This allows the filter to track slow stellar and instrumental variations while largely preserving the short-duration transit signal.

The filter operates over a symmetric window of $2m+1$ points around each data point $i$. A polynomial of degree $k$ is fitted by minimizing the least-squares error:

$$E = \sum_j [p(j) - y_{i+j}]^2$$

where:

$$p(j) = a_0 + a_1 j + a_2 j^2 + \dots + a_k j^k$$

The resulting low-frequency stellar/instrumental continuum is $F_{\text{trend}}(t)$.

The normalized, pre-whitened flux is then obtained by dividing the raw flux by this trend:

$$F_{\text{norm}}(t) = \frac{F_{\text{raw}}(t)}{F_{\text{trend}}(t)}$$

The resulting signal is centered approximately around:

$$F_{\text{norm}} \approx 1.0$$

with the planetary transits appearing as small downward excursions.

---

## 7. Step 4: Matched Filtering via Box Least Squares (BLS)

A conventional Fourier Transform is not ideal for detecting planetary transits because a transit is a **short-duration, low-duty-cycle signal**.

A transit typically occupies only a small fraction of the orbital period. Representing such a square-like pulse using sinusoidal harmonics spreads the signal's energy across many frequencies, reducing detection efficiency.

### 7.1 Box Least Squares

The **Box Least Squares (BLS)** algorithm acts as a specialized matched filter for periodic transit-like signals.

The algorithm sweeps over a grid of trial:

- Orbital periods $P$
- Transit durations $\tau$
- Reference epochs $t_0$

For each combination, it evaluates how well a box-shaped transit model matches the observed data.

The transit is represented by two discrete flux states:

- **Out of transit:** High state $H$
- **In transit:** Low state $L$

For a trial period $P$ and duration $\tau$, the fractional transit duration is:

$$q = \frac{\tau}{P}$$

The BLS algorithm evaluates the reduction in squared residuals produced by placing a box-shaped dip at a particular phase.

One form of the BLS detection statistic is the **Signal Residue (SR)**:

$$SR = \max \left[ \frac{s^2}{q(1-q)} \right]$$

where $s$ represents the weighted sum associated with data points falling within the trial transit window.

The strongest peak in the BLS periodogram identifies the candidate orbital period.

### 7.2 Detection Result

For Kepler-10b, the recovered orbital period is approximately:

$$P \approx 0.8375\ \text{days}$$

---

## 8. Step 5: Coherent Phase Folding & Stacking

A single transit of Kepler-10b produces a flux decrease of only approximately **160 parts per million (ppm)**, making the individual event difficult to distinguish visually from the instrument noise.

### 8.1 Phase Folding

Once the orbital period $P$ and reference epoch $t_0$ are known, the linear time series is folded into a single orbital phase.

The phase is calculated as:

$$\phi_i = \left[ \frac{t_i - t_0}{P} + 0.5 \right] \bmod 1 - 0.5$$

This maps all observations from multiple orbital cycles onto a common phase interval.

### 8.2 Coherent Stacking

All detected transits are then aligned in phase and averaged into approximately **10-minute phase bins**.

Because the planetary signal is coherent while much of the measurement noise is uncorrelated, stacking improves the signal-to-noise ratio approximately as:

$$\text{SNR}_{\text{folded}} \approx \text{SNR}_{\text{single}} \sqrt{N}$$

where $N$ is the number of observed transits.

For approximately 100 individual transits, the coherent stacking process dramatically improves the visibility of the transit profile.

---

## 9. Step 6: Astrophysical Parameter Extraction

The DSP pipeline extracts three primary empirical observables:

| Observable                    | Measured Value |
| :---------------------------- | :------------- |
| **Orbital Period ($P$)**      | $0.8375$ days  |
| **Transit Depth ($\delta$)**  | $162.0$ ppm    |
| **Transit Duration ($\tau$)** | $1.66$ hours   |

These observables can then be converted into physical properties of the planetary system.

### 9.1 Planet-to-Star Radius Ratio

For a simplified transit model, the fractional loss of stellar flux is approximately equal to the ratio of the projected areas of the planet and star:

$$\delta = \frac{A_p}{A_*} = \frac{\pi R_p^2}{\pi R_*^2}$$

Therefore:

$$\delta = \left( \frac{R_p}{R_*} \right)^2$$

and:

$$R_p = R_* \sqrt{\delta}$$

Using the measured transit depth and the stellar parameters for Kepler-10 gives:

$$R_p \approx 1.48\,R_\oplus$$

This is consistent with published values of approximately:

$$R_p \approx 1.47 \pm 0.03\,R_\oplus$$

### 9.2 Semi-Major Axis

Assuming a circular orbit, Newton's formulation of Kepler's Third Law can be derived by equating gravitational and centripetal forces.

The gravitational force is:

$$F_g = \frac{G M_* m_p}{a^2}$$

while the centripetal force is:

$$F_c = \frac{m_p v^2}{a}$$

Equating the two:

$$\frac{G M_* m_p}{a^2} = \frac{m_p v^2}{a}$$

which gives:

$$v^2 = \frac{G M_*}{a}$$

For a circular orbit:

$$v = \frac{2\pi a}{P}$$

Substituting this expression and solving for $a$:

$$a = \left( \frac{G M_* P^2}{4\pi^2} \right)^{1/3}$$

Using the stellar mass and measured orbital period produces:

$$a \approx 0.0169\ \text{AU}$$

or approximately:

$$a \approx 2.5 \times 10^6\ \text{km}$$

This extremely small orbital separation places Kepler-10b much closer to its host star than Mercury is to the Sun.

### 9.3 Mean Orbital Velocity

Assuming an approximately circular orbit, the orbital velocity is:

$$v_{\text{orbit}} = \frac{2\pi a}{P}$$

Using the derived semi-major axis and measured period:

$$v_{\text{orbit}} \approx 218.9\ \text{km/s}$$

Thus, Kepler-10b travels around its host star at a mean orbital speed of approximately **219 km/s**.

---

## 10. Results Summary

The complete DSP pipeline recovers the following key characteristics of Kepler-10b:

| Parameter                 | Pipeline Result              |
| :------------------------ | :--------------------------- |
| **Orbital Period**        | $0.8375$ days                |
| **Transit Depth**         | $162.0$ ppm                  |
| **Transit Depth**         | $0.0162\%$                   |
| **Transit Duration**      | $1.66$ hours                 |
| **Planet Radius**         | $\approx 1.48\,R_\oplus$     |
| **Semi-Major Axis**       | $\approx 0.0169$ AU          |
| **Orbital Distance**      | $\approx 2.5 \times 10^6$ km |
| **Mean Orbital Velocity** | $\approx 218.9$ km/s         |

---

## 11. Signal Processing Interpretation

The project demonstrates how concepts from classical DSP can be applied directly to modern astronomical observations.

| Astronomical Challenge          | DSP Solution                         |
| :------------------------------ | :----------------------------------- |
| Cosmic-ray spikes               | Robust MAD-based asymmetric clipping |
| Instrumental anomalies          | Quality-flag filtering               |
| Slow stellar/instrumental drift | Savitzky-Golay detrending            |
| Colored $1/f$ noise             | Baseline removal / pre-whitening     |
| Low-duty-cycle transit          | Box Least Squares matched filtering  |
| Weak individual events          | Coherent phase folding               |
| Uncorrelated noise              | Multi-transit stacking               |
| Physical interpretation         | Keplerian orbital mechanics          |

The key idea is that **the transit signal is weak but highly structured**. Its periodicity, duration and shape provide a strong prior that can be exploited using matched-filter techniques.

---

## 12. Conclusion & Project Value

This project demonstrates the powerful intersection of **digital signal processing, computational astronomy, and fundamental physics**.

By systematically addressing:

1. Non-Gaussian outliers,
2. Instrumental and stellar baseline drift,
3. Colored $1/f$ noise,
4. Low-SNR periodic pulse detection,
5. Coherent signal accumulation,

the pipeline is able to extract a planetary transit signal from spacecraft photometry and convert it into meaningful physical parameters.

The project therefore illustrates a complete scientific workflow:

```text
         Raw Spacecraft Telemetry
                    ↓
              Data Cleaning
                    ↓
            Noise Suppression
                    ↓
              Pre-Whitening
                    ↓
            Matched Filtering
                  (BLS)
                    ↓
             Period Detection
                    ↓
              Phase Folding
                    ↓
            Coherent Stacking
                    ↓
         Transit Characterization
                    ↓
           Keplerian Dynamics
                    ↓
      Physical Planetary Parameters
```

Ultimately, the pipeline shows that the detection of distant worlds can be framed as a rigorous **weak-signal detection problem**, where DSP techniques provide the tools necessary to recover signals that are otherwise buried beneath instrumental and astrophysical noise.

---

## 13. Core Equations

For quick reference, the principal equations used in the pipeline are:

### Savitzky-Golay Detrending:

$$F_{\text{norm}}(t) = \frac{F_{\text{raw}}(t)}{F_{\text{trend}}(t)}$$

### Transit Depth:

$$\delta = \left( \frac{R_p}{R_*} \right)^2$$

### Planet Radius:

$$R_p = R_* \sqrt{\delta}$$

### Kepler's Third Law:

$$a = \left( \frac{G M_* P^2}{4\pi^2} \right)^{1/3}$$

### Orbital Velocity:

$$v_{\text{orbit}} = \frac{2\pi a}{P}$$

### Phase Folding:

$$\phi_i = \left[ \frac{t_i - t_0}{P} + 0.5 \right] \bmod 1 - 0.5$$

### Coherent Stacking SNR:

$$\text{SNR}_{\text{folded}} \approx \text{SNR}_{\text{single}} \sqrt{N}$$

---

## 14. Glossary

## Signal Processing & Telemetry

### Digital Signal Processing (DSP)

The conversion of raw analog signals from celestial objects into numerical data streams that computers can clean, filter and analyze.

### Telemetry

The automated collection of measurements and data from remote spacecraft and its wireless transmission back to Earth for monitoring and analysis.

### White Noise

A random signal that has equal intensity across all frequencies, resulting in a flat and constant power spectral density.

### $1/f$ Colored Noise (Pink Noise)

A random, time-correlated signal fluctuation where the power spectral density $S(f)$ is inversely proportional to the frequency ($f$) of the signal:
$$S(f) \propto \frac{1}{f}$$
It is referred to as "colored" because it exhibits greater power at lower frequencies.

### Pre-whitening

An iterative data-processing technique used to extract individual periodic oscillation frequencies (such as stellar pulsation modes) from noisy time-series data.

### Matched Filtering

A signal processing technique used to detect faint, known signal patterns hidden inside loud observational noise by correlating the noisy data with a theoretical template of the signal.

### Coherent Stacking

A data processing technique where multiple individual signals or observations are aligned by both their amplitude and phase before being added together constructively.

---

## Photometry & Astronomical Observations

### Photometry

The technique of measuring the flux or intensity of light radiated by celestial objects.

### Cadence

The time interval defining how frequently an instrument observes a specific target.

### Simple Aperture Photometry (SAP)

The measurement of a target star's brightness calculated by summing the raw light values of all the pixels within a designated aperture on a digital sensor.

### Pre-search Data Conditioning Simple Aperture Photometry (PDCSAP)

A corrected version of SAP data that has been algorithmically cleaned to remove instrumental noise, spacecraft artifacts, and long-term systematic errors while preserving short-term astrophysical variations like exoplanet transits.

### Epoch

A specific moment in time used as a reference point to measure and calculate the changing positions, coordinates or orbital paths of celestial objects.

### Starspots (Star sports)

Dark or bright temporary patches on the surface of a star caused by magnetic field concentrations.

### Transit Depth

The fractional decrease in a star's brightness observed when an exoplanet passes directly in front of it.

---

## Instrumentation & Systematics

### Pointing Jitter

Rapid, uncommanded micro-vibrations of a telescope's line-of-sight over short time intervals (typically in s or ms).

### Thermal Focus Drift

The gradual loss of sharp focus in a telescope's optics caused by changing thermal conditions during an observation session.

### Cotrending Basis Vectors (CBVs)

A set of mathematical reference profiles used in astronomy to model and remove systematic instrumental noise from space telescope light curves.

### Reaction-Wheel Desaturation Event

_(Also known as momentum dumping or momentum unloading)_  
A scheduled operational procedure where a spacecraft fires thrusters to slow down its rapidly spinning reaction wheels, shedding accumulated angular momentum.

---

## Time-Series Analysis & Filtering

### Savitzky-Golay Polynomial Filter

A digital data-smoothing method used to reduce high-frequency noise in spectra, light curves and time-series measurements without distorting critical signal features like peak heights and line widths.

### Periodogram

A graph and statistical tool used to identify hidden repeating cycles or periodicities within time-series data.

### Phase Folding

A data analysis technique that stacks multiple repeating cycles of a periodic signal on top of each other by plotting brightness against orbital phase rather than absolute time.

---

## Scientific Python Libraries

- **`lightkurve`**: An open-source Python package designed to analyze astronomical flux time-series data. It provides a user-friendly API to download, inspect and analyze data collected by NASA’s Kepler, K2 and TESS exoplanet missions.
- **`astropy`**: An open-source Python library used for core astronomical data analysis, coordinate transformations, physical calculations and research.
- **`numpy`**: An open-source foundational library for multi-dimensional array operations, mathematical calculations and numeric data analysis in Python.
- **`scipy`**: An open-source Python library built on NumPy used for advanced scientific calculations, technical computing, optimization and signal processing algorithms.
- **`matplotlib`**: An open-source Python visualization library used for creating static, animated and interactive publication-quality data plots.
