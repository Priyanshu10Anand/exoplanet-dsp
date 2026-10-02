"""
Step 3: Baseline Detrending (Pre-Whitening)
- Removes 1/f stellar variability and spacecraft thermal drift
- Employs Savitzky-Golay polynomial filtering (preserves transit edges)
- Pre-whitens time series to unity baseline for matched filtering
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

def clean_and_normalize(lc):
    """Step 2 pre-filter: reject NaNs, quality flags, and asymmetric outliers."""
    # 1. Strip missing records and non-zero hardware anomaly flags
    clean = lc.remove_nans()
    clean = clean[clean.quality == 0]

    # 2. Compute outlier-resistant noise spread via MAD
    flux = clean.flux.value
    med = np.nanmedian(flux)
    mad = np.nanmedian(np.abs(flux - med))
    sigma = 1.4826 * mad

    # 3. Gate data: clip positive cosmic rays (+4s), preserve negative transits (-8s)
    mask = (flux <= med + 4.0 * sigma) & (flux >= med - 8.0 * sigma)
    return clean[mask].normalize()

def main():
    # 1. Ingest Quarter 3 calibrated telemetry
    print("[+] Loading Kepler-10 (Q3) data...")
    raw_lc = lk.search_lightcurve(
        "Kepler-10",
        mission = "Kepler",
        quarter = 3,
        author = "Kepler",
        exptime = 1800
    ).download()

    clean_lc = clean_and_normalize(raw_lc)

    # 2. Savitzky-Golay filter geometry:
    # 151 cadences ~ 3.1 days (Window >> Transit Duration ~ 1.8h)
    window_length = 151   # Must be an odd integer for symmetric window centering
    polyorder = 2         # 2nd-order local polynomial preserves transit morphology

    print(f"[+] Applying Savitzky-Golay detrending (window = {window_length}, order = {polyorder})...")
    # Low-pass filter baseline and divide it out to isolate flat residuals
    flat_lc, trend_lc = clean_lc.flatten(
        window_length = window_length,
        polyorder = polyorder,
        return_trend = True
    )

    # 3. Compute baseline noise-floor reduction (ppm)
    raw_std = np.std(clean_lc.flux.value) * 1e6
    flat_std = np.std(flat_lc.flux.value) * 1e6

    print("\n--- Detrending Performance ---")
    print(f"Filter window:          {window_length} cadences (~{window_length * 29.4 / 60 / 24:.2f} days)")
    print(f"Pre-filter scatter:     {raw_std:.1f} ppm")
    print(f"Post-filter scatter:    {flat_std:.1f} ppm")
    print(f"Noise floor reduction:  {((raw_std - flat_std) / raw_std) * 100:.1f}%")

    # 4. Multi-tier diagnostic visualization
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), constrained_layout=True)

    # Panel 1: Conditioned flux overlaid with low-pass SavGol trend
    axes[0].scatter(
        clean_lc.time.value, clean_lc.flux.value, s = 1.2, 
        color = "gray", 
        alpha = 0.6, 
        label = "Cleaned Flux"
    )
    axes[0].plot(
        trend_lc.time.value, trend_lc.flux.value, 
        color = "crimson", 
        lw = 1.8, 
        label = "SavGol Baseline Trend"
    )
    axes[0].set_ylabel("Normalized Flux")
    axes[0].set_title("1. Low-Frequency Baseline Extraction")
    axes[0].legend(loc = "upper right")

    # Panel 2: Pre-whitened residuals (continuum leveled to 1.0)
    axes[1].scatter(
        flat_lc.time.value, 
        flat_lc.flux.value, 
        s = 1.2, 
        color = "teal", 
        alpha = 0.7, 
        label = "Detrended Flux (Residuals)"
    )
    axes[1].axhline(
        1.0, color = "black", 
        linestyle = "--", 
        lw = 0.8, 
         alpha = 0.7
    )
    axes[1].set_ylabel("Normalized Flux")
    axes[1].set_title("2. Pre-Whitened Signal (Stellar Drift Removed)")
    axes[1].legend(loc = "upper right")

    # Panel 3: Zoomed 5-day window revealing periodic transit dips to the eye
    t_start, t_end = 280.0, 285.0
    zoom_mask = (flat_lc.time.value >= t_start) & (flat_lc.time.value <= t_end)

    axes[2].scatter(
        flat_lc.time.value[zoom_mask], 
        flat_lc.flux.value[zoom_mask], 
        s = 12, 
        color = "navy", 
        label = "Zoomed 5-Day Segment"
    )
    axes[2].plot(
        flat_lc.time.value[zoom_mask], 
        flat_lc.flux.value[zoom_mask], 
        color = "navy", 
        lw = 0.7, 
        alpha = 0.3
    )
    axes[2].axhline(
        1.0, color = "crimson", 
        linestyle = "--", 
        lw = 0.8, 
        label = "Baseline (1.0)"
    )
    axes[2].set_xlim(t_start, t_end)
    axes[2].set_ylim(0.9990, 1.0008)
    axes[2].set_ylabel("Normalized Flux")
    axes[2].set_xlabel("Time (BJD - 2454833)")
    axes[2].set_title("3. Micro-View: Periodic Transit Dips Visible to the Eye")
    axes[2].legend(loc = "lower right")

    # 5. Export high-resolution graphic before rendering interactive canvas
    out_file = "03_savgol_detrending.png"
    print(f"\n[+] Saving plot to {out_file}...")
    plt.savefig(out_file, dpi = 300, bbox_inches = "tight")
    print(f"[+] Saved {out_file} successfully!")

    plt.show()

if __name__ == "__main__":
    main()