"""
Step 4: Matched Filtering via Box Least Squares (BLS)
- Sweeps trial orbital periods P and transit durations tau
- Computes likelihood power spectrum (equivalent to matched filter response)
- Detects the orbital period and mid-transit epoch (t0)
"""

import astropy.units as u
from astropy.timeseries import BoxLeastSquares
import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np

def clean_and_detrend(lc):
    """Steps 2 & 3: Clean telemetry, normalize, and detrend via SavGol."""
    clean = lc.remove_nans()
    clean = clean[clean.quality == 0]

    flux = clean.flux.value
    med = np.nanmedian(flux)
    mad = np.nanmedian(np.abs(flux - med))
    sigma = 1.4826 * mad
    mask = (flux <= med + 4.0 * sigma) & (flux >= med - 8.0 * sigma)
    clean = clean[mask].normalize()

    flat_lc = clean.flatten(window_length=151, polyorder=2)
    return flat_lc

def main():
    print("[+] Loading and preprocessing Kepler-10 (Q3)...")
    raw_lc = lk.search_lightcurve(
        "Kepler-10",
        mission="Kepler",
        quarter=3,
        author="Kepler",
        exptime=1800
    ).download()

    flat_lc = clean_and_detrend(raw_lc)

    time = flat_lc.time.value
    flux = flat_lc.flux.value

    # Initialize BLS Matched Filter model
    bls = BoxLeastSquares(time, flux)

    # Define search parameter grid
    # We test periods from 0.4 days (ultra-short) to 5.0 days
    # Durations from 0.03 days (~45 mins) to 0.12 days (~2.9 hours)
    durations = np.linspace(0.03, 0.12, 15)
    periods = np.linspace(0.4, 5.0, 10000)

    print(f"[+] Running BLS matched filter across {len(periods)} trial periods...")
    periodogram = bls.power(periods, durations)

    # Locate peak in power spectrum
    best_idx = np.argmax(periodogram.power)
    best_period = periodogram.period[best_idx]
    best_t0 = periodogram.transit_time[best_idx]
    best_duration = periodogram.duration[best_idx]
    best_depth = periodogram.depth[best_idx]
    max_power = periodogram.power[best_idx]

    print("\n--- Detection Results ---")
    print(f"Detected Orbital Period: {best_period:.5f} days ({best_period * 24:.2f} hours)")
    print(f"Mid-transit Epoch (t0):  {best_t0:.4f} BJD")
    print(f"Transit Duration (tau):  {best_duration * 24:.2f} hours")
    print(f"Transit Depth (delta):   {best_depth * 1e6:.1f} ppm ({best_depth * 100:.4f}%)")
    print(f"Peak BLS Power:          {max_power:.2f}")

    # Plot BLS Periodogram (Matched Filter Frequency Response)
    fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
    ax.plot(periodogram.period, periodogram.power, "navy", lw=0.9, label="BLS Spectrum")
    ax.axvline(best_period, color="crimson", linestyle="--", lw=1.5,
               label=f"Detected Planet: P = {best_period:.4f} d")

    # Mark harmonics (half period / 2x period) to demonstrate filter harmonic structure
    ax.axvline(best_period * 2.0, color="orange", linestyle=":", lw=1.2, label=f"2x Harmonic ({best_period * 2:.3f} d)")
    ax.axvline(best_period * 0.5, color="green", linestyle=":", lw=1.2, label=f"0.5x Sub-Harmonic ({best_period * 0.5:.3f} d)")

    ax.set_title("Box Least Squares (BLS) Matched Filter Spectrum")
    ax.set_xlabel("Trial Period (Days)")
    ax.set_ylabel("Filter Power (Signal Residue Reduction)")
    ax.set_xlim(0.4, 5.0)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="upper right")

    plt.show()

if __name__ == "__main__":
    main()