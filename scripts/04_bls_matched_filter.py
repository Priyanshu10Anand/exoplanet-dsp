"""
Step 4: Matched Filtering via Box Least Squares (BLS)
- Sweeps trial orbital periods P and transit durations tau
- Computes likelihood power spectrum (equivalent to matched filter response)
- Detects the orbital period and mid-transit epoch (t0)
"""

import warnings
warnings.filterwarnings(
    "ignore", 
    category = UserWarning, 
    module = "lightkurve"
)
import astropy.units as u
from astropy.timeseries import BoxLeastSquares
import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def clean_and_detrend(lc):
    """Steps 2 & 3: Clean telemetry, normalize, and detrend via SavGol."""
    # 1. Strip NaNs and non-zero hardware anomaly flags
    clean = lc.remove_nans()
    clean = clean[clean.quality == 0]

    # 2. Outlier-immune noise estimation via MAD
    flux = clean.flux.value
    med = np.nanmedian(flux)
    mad = np.nanmedian(np.abs(flux - med))
    sigma = 1.4826 * mad

    # 3. Asymmetric gate: reject cosmic rays (+4s), preserve transit dips (-8s)
    mask = (flux <= med + 4.0 * sigma) & (flux >= med - 8.0 * sigma)
    clean = clean[mask].normalize()

    # 4. Remove low-frequency stellar baseline using SavGol filter
    flat_lc = clean.flatten(
        window_length = 151, 
        polyorder = 2
    )
    return flat_lc

def main():
    # 1. Fetch Quarter 3 data and condition the signal
    print("[+] Loading and preprocessing Kepler-10 (Q3)...")
    raw_lc = lk.search_lightcurve(
        "Kepler-10",
        mission = "Kepler",
        quarter = 3,
        author = "Kepler",
        exptime = 1800
    ).download()

    flat_lc = clean_and_detrend(raw_lc)

    time = flat_lc.time.value
    flux = flat_lc.flux.value

    # 2. Configure Box Least Squares matched filter
    bls = BoxLeastSquares(time, flux)

    # 3. Grid trial transit durations (~45m to ~2.9h) and periods (0.4d to 5.0d)
    durations = np.linspace(0.03, 0.12, 15)
    periods = np.linspace(0.4, 5.0, 10000)

    # 4. Compute power spectrum across parameter space
    print(f"[+] Running BLS matched filter across {len(periods)} trial periods...")
    periodogram = bls.power(periods, durations)

    # 5. Extract optimal transit parameters from maximum spectral peak
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

    # 6. Plot BLS periodogram with fundamental peak and harmonic markers
    fig, ax = plt.subplots(
        figsize = (11, 5), 
        constrained_layout = True
    )
    ax.plot(
        periodogram.period, 
        periodogram.power, 
        "navy", 
        lw = 0.9, 
        label = "BLS Spectrum"
    )
    ax.axvline(
        best_period, 
        color = "crimson", 
        linestyle = "--", 
        lw = 1.5,
        label = f"Detected Planet: P = {best_period:.4f} d"
    )

    # Mark filter response harmonics (integer multiples and sub-multiples)
    ax.axvline(
        best_period * 2.0, 
        color = "orange", 
        linestyle = ":", 
        lw = 1.2, 
        label = f"2x Harmonic ({best_period * 2:.3f} d)"
    )
    ax.axvline(
        best_period * 0.5, 
        color = "green", 
        linestyle = ":", 
        lw = 1.2, 
        label = f"0.5x Sub-Harmonic ({best_period * 0.5:.3f} d)"
    )

    ax.set_title("Box Least Squares (BLS) Matched Filter Spectrum")
    ax.set_xlabel("Trial Period (Days)")
    ax.set_ylabel("Filter Power (Signal Residue Reduction)")
    ax.set_xlim(0.4, 5.0)
    ax.grid(
        True, 
        linestyle = "--", 
        alpha = 0.5
    )
    ax.legend(loc = "upper right")

    # 7. Export high-resolution figure prior to canvas render
    out_file = "04_bls_periodogram.png"
    print(f"\n[+] Saving plot to {out_file}...")
    plt.savefig(out_file, dpi = 300, bbox_inches = "tight")
    print(f"[+] Saved {out_file} successfully!")

    plt.show()

if __name__ == "__main__":
    main()