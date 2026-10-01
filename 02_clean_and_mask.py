"""
Step 2: Data Cleaning & Quality Masking
- NaN drop
- Bitmask telemetry filtering
- Robust Median Absolute Deviation (MAD) calculation
- Asymmetric sigma-clipping (protects negative transit dips)
- Baseline normalization
"""

import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np

def clean_lightcurve(lc):
    """Clean and normalize a raw Kepler light curve."""
    # 1. Drop NaN cadences
    clean = lc.remove_nans()

    # 2. Kepler telemetry quality flag filtering
    # quality == 0 isolates cadences unaffected by known spacecraft anomalies
    clean = clean[clean.quality == 0]

    # 3. Robust statistics using Median Absolute Deviation (MAD)
    flux = clean.flux.value
    time = clean.time.value

    median_flux = np.nanmedian(flux)
    # 1.4826 converts MAD to equivalent standard deviation for a normal distribution
    mad = np.nanmedian(np.abs(flux - median_flux))
    sigma_mad = 1.4826 * mad

    # Asymmetric clipping thresholds
    upper_threshold = median_flux + 4.0 * sigma_mad   # Strict on cosmic ray spikes
    lower_threshold = median_flux - 8.0 * sigma_mad   # Lenient on transit dips

    mask = (flux <= upper_threshold) & (flux >= lower_threshold)
    clipped_lc = clean[mask]

    # 4. Normalize median flux to 1.0
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
    print("[+] Fetching Kepler-10 Quarter 3 data...")
    lc = lk.search_lightcurve(
        "Kepler-10",
        mission="Kepler",
        quarter=3,
        author="Kepler",
        exptime=1800
    ).download()

    print("[+] Cleaning and performing asymmetric sigma clipping...")
    clean_lc, stats = clean_lightcurve(lc)

    print("\n--- Cleaning Summary ---")
    print(f"Original points:      {stats['raw_count']}")
    print(f"NaN / quality flags:  -{stats['nan_dropped']}")
    print(f"Outliers clipped:     -{stats['outliers_clipped']}")
    print(f"Final usable samples:  {stats['final_count']}")
    print(f"Photometric scatter:   {stats['sigma_mad'] / stats['median_flux'] * 1e6:.1f} ppm")

    # Plot raw vs cleaned comparison
    fig, axes = plt.subplots(2, 1, figsize=(11, 6), sharex=True, constrained_layout=True)

    # Raw plot
    axes[0].scatter(lc.time.value, lc.flux.value, s=1.5, color="gray", alpha=0.7, label="Raw PDCSAP Flux")
    axes[0].set_ylabel("Flux (e⁻ / s)")
    axes[0].set_title("Kepler-10 (Q3): Raw Telemetry")
    axes[0].legend(loc="upper right")

    # Cleaned & Normalized plot
    axes[1].scatter(clean_lc.time.value, clean_lc.flux.value, s=1.5, color="midnightblue", label="Cleaned & Normalized Flux")
    axes[1].axhline(1.0, color="crimson", linestyle="--", linewidth=0.8, alpha=0.8, label="Baseline (1.0)")
    axes[1].set_ylabel("Normalized Flux")
    axes[1].set_xlabel("Time (BJD - 2454833)")
    axes[1].set_title("Cleaned, Masked & Normalized Time Series")
    axes[1].legend(loc="upper right")

    plt.show()

if __name__ == "__main__":
    main()