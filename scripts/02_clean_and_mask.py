"""
Step 2: Data Cleaning & Quality Masking
- Drops missing values and hardware-flagged cadences
- Estimates robust noise spread via MAD
- Rejects cosmic-ray spikes while preserving physical transit dips
- Rescales continuum to unity
"""

import warnings
warnings.filterwarnings(
    "ignore", 
    category = UserWarning, 
    module = "lightkurve"
)
import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def clean_lightcurve(lc):
    """Clean telemetry, reject non-Gaussian outliers, and normalize baseline."""
    # 1. Strip missing/unrecorded cadences
    clean = lc.remove_nans()

    # 2. Retain only optimal cadences (bitmask 0 = nominal spacecraft operations)
    clean = clean[clean.quality == 0]

    # 3. Compute robust noise scale via Median Absolute Deviation (outlier-immune)
    flux = clean.flux.value
    time = clean.time.value

    median_flux = np.nanmedian(flux)
    # Scale factor 1.4826 aligns MAD with standard deviation for Gaussian noise
    mad = np.nanmedian(np.abs(flux - median_flux))
    sigma_mad = 1.4826 * mad

    # 4. Asymmetric thresholding: aggressively clip positive cosmic rays, preserve negative transit dips
    upper_threshold = median_flux + 4.0 * sigma_mad   # Tight upper gate for cosmic-ray hits
    lower_threshold = median_flux - 8.0 * sigma_mad   # Wide lower gate to protect transit depths

    mask = (flux <= upper_threshold) & (flux >= lower_threshold)
    clipped_lc = clean[mask]

    # 5. Rescale baseline to 1.0 for fractional transit depth analysis
    normalized_lc = clipped_lc.normalize()

    stats = {
        "raw_count": len(lc),
        "nan_dropped": len(lc) - len(clean),
        "outliers_clipped": len(clean) - len(clipped_lc),
        "final_count": len(normalized_lc),
        "median_flux": median_flux,
        "sigma_mad": sigma_mad,
    }
    return normalized_lc, stats

def main():
    # 1. Fetch raw Quarter 3 target pixel telemetry
    print("[+] Fetching Kepler-10 Quarter 3 data...")
    lc = lk.search_lightcurve(
        "Kepler-10",
        mission="Kepler",
        quarter=3,
        author="Kepler",
        exptime=1800
    ).download()

    # 2. Execute cleaning pipeline and generate data accounting metrics
    print("[+] Cleaning and performing asymmetric sigma clipping...")
    clean_lc, stats = clean_lightcurve(lc)

    print("\n--- Cleaning Summary ---")
    print(f"Original points:      {stats['raw_count']}")
    print(f"NaN / quality flags:  -{stats['nan_dropped']}")
    print(f"Outliers clipped:     -{stats['outliers_clipped']}")
    print(f"Final usable samples:  {stats['final_count']}")
    print(f"Photometric scatter:   {stats['sigma_mad'] / stats['median_flux'] * 1e6:.1f} ppm")

    # 3. Comparative diagnostic plots: raw telemetry vs conditioned baseline
    fig, axes = plt.subplots(
        2, 1, figsize = (11, 6), 
        sharex = True, 
        constrained_layout = True
    )

    # Panel 1: Raw flux showing sensor systematics and cosmic-ray spikes
    axes[0].scatter(
        lc.time.value, lc.flux.value, s = 1.5, 
        color = "gray", 
        alpha = 0.7, 
        label = "Raw PDCSAP Flux"
    )
    axes[0].set_ylabel("Flux (e⁻ / s)")
    axes[0].set_title("Kepler-10 (Q3): Raw Telemetry")
    axes[0].legend(loc = "upper right")

    # Panel 2: Normalized flux isolated around unity (1.0)
    axes[1].scatter(
        clean_lc.time.value, 
        clean_lc.flux.value, 
        s = 1.5, 
        color = "midnightblue", 
        label = "Cleaned & Normalized Flux"
    )
    axes[1].axhline(
        1.0, color = "crimson", 
        linestyle = "--", 
        linewidth = 0.8, 
        alpha = 0.8, 
        label = "Baseline (1.0)"
    )
    axes[1].set_ylabel("Normalized Flux")
    axes[1].set_xlabel("Time (BJD - 2454833)")
    axes[1].set_title("Cleaned, Masked & Normalized Time Series")
    axes[1].legend(loc = "upper right")

    plt.show()

if __name__ == "__main__":
    main()