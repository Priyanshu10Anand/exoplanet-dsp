"""
Step 5: Coherent Phase Folding & Stacking
- Folds time-series into phase space phi in [-0.5, +0.5]
- Coherently stacks ~100 transits across the 90-day baseline
- Suppresses uncorrelated white noise by sqrt(N) via phase binning
- Overlays best-fit boxcar matched-filter model
"""

from astropy.timeseries import BoxLeastSquares
import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np

def clean_and_detrend(lc):
    """Steps 2 & 3: Telemetry cleaning, asymmetric MAD clipping, and SavGol detrending """
    # 1. Filter missing cadences and spacecraft telemetry anomalies
    clean = lc.remove_nans()
    clean = clean[clean.quality == 0]

    # 2. Asymmetric noise gating: drop cosmic rays (+4s), retain transit dips (-8s)
    flux = clean.flux.value
    med = np.nanmedian(flux)
    mad = np.nanmedian(np.abs(flux - med))
    sigma = 1.4826 * mad
    mask = (flux <= med + 4.0 * sigma) & (flux >= med - 8.0 * sigma)
    clean = clean[mask].normalize()

    # 3. Pre-whiten baseline to 1.0 using Savitzky-Golay polynomial filter
    return clean.flatten(window_length=151, polyorder=2)

def main():
    # 1. Ingest Quarter 3 telemetry and condition signal
    print("[+] Loading Kepler-10 (Q3) and preprocessing...")
    raw_lc = lk.search_lightcurve(
        "Kepler-10",
        mission = "Kepler",
        quarter = 3,
        author = "Kepler",
        exptime = 1800
    ).download()

    flat_lc = clean_and_detrend(raw_lc)

    # 2. Refined BLS search around the known peak to lock period and epoch
    time = flat_lc.time.value
    flux = flat_lc.flux.value
    bls = BoxLeastSquares(time, flux)

    durations = np.linspace(0.03, 0.10, 10)
    periods = np.linspace(0.80, 0.86, 3000)  # Narrow search around detected peak
    periodogram = bls.power(periods, durations)

    # Extract best-fit transit parameters
    best_idx = np.argmax(periodogram.power)
    period = periodogram.period[best_idx]
    t0 = periodogram.transit_time[best_idx]
    duration = periodogram.duration[best_idx]
    depth = periodogram.depth[best_idx]

    print(f"\n--- Matched Filter Solution ---")
    print(f"Period (P):    {period:.5f} days")
    print(f"Epoch (t0):    {t0:.4f} BJD")
    print(f"Duration (tau): {duration * 24:.2f} hours")
    print(f"Depth (delta): {depth * 1e6:.1f} ppm")

    # 3. Phase wrapping: map linear time into circular orbital phase centered at t0
    folded_lc = flat_lc.fold( period = period, epoch_time = t0)

    # 4. Phase binning: 10-minute uniform bins boost SNR by averaging uncorrelated noise
    bin_size_days = 10.0 / (60.0 * 24.0)
    binned_lc = folded_lc.bin(time_bin_size = bin_size_days)

    # 5. Synthesize analytic matched-filter boxcar template
    phase_dense = np.linspace(-0.25, 0.25, 1000)
    box_model = np.ones_like(phase_dense)
    in_transit = np.abs(phase_dense) <= (duration / 2.0)
    box_model[in_transit] = 1.0 - depth

    # 6. Diagnostic plotting: folded scatter vs. high-SNR binned profile
    fig, axes = plt.subplots(
        2, 
        1, 
        figsize = (11, 7), 
        sharex = True, 
        constrained_layout = True
    )

    # Panel 1: Complete coherently stacked cadences (~100 transits superimposed)
    axes[0].scatter(
        folded_lc.time.value, 
        folded_lc.flux.value, 
        s = 1.5, 
        color = "teal", 
        alpha = 0.35, 
        label = "Phase-Folded Cadences (All ~100 Transits)"
    )
    axes[0].plot(
        phase_dense, 
        box_model, 
        color = "crimson", 
        lw = 2.0, 
        label = "BLS Boxcar Model"
    )
    axes[0].set_ylabel("Normalized Flux")
    axes[0].set_title(f"Kepler-10b: Phase-Folded Light Curve (P = {period:.5f} d)")
    axes[0].set_ylim(0.9988, 1.0012)
    axes[0].legend(loc = "lower right")

    # Panel 2: Binned transit profile resolving steep ingress/egress walls
    axes[1].scatter(
        folded_lc.time.value, 
        folded_lc.flux.value, 
        s = 1.0, 
        color = "gray", 
        alpha = 0.2, 
        label = "Folded Scatter"
    ) 
    axes[1].errorbar(
        binned_lc.time.value,
        binned_lc.flux.value,
        yerr = binned_lc.flux_err.value if binned_lc.flux_err is not None else None,
        fmt = "o",
        color = "navy",
        markersize = 3.5,
        elinewidth = 0.8,
        capsize = 1.5,
        label = "Binned (10-min bins)"
    )
    axes[1].plot(
        phase_dense, 
        box_model, 
        color = "crimson", 
        lw = 2.0, 
        label = f"Box Model (Depth = {depth*1e6:.0f} ppm)"
    )
    axes[1].set_ylabel("Normalized Flux")
    axes[1].set_xlabel("Time from Mid-Transit (Days)")
    axes[1].set_xlim(-0.2, 0.2)
    axes[1].set_ylim(0.9992, 1.0006)
    axes[1].axvline(
        0.0, 
        color = "black", 
        linestyle = ":", 
        lw = 0.8, 
        alpha = 0.6
    )
    axes[1].axhline(
        1.0, 
        color = "black", 
        linestyle = "--", 
        lw = 0.8, 
        alpha = 0.6
    )
    axes[1].set_title("Coherent Stack & Binned Profile (Noise Averaged Out)")
    axes[1].legend(loc = "lower right")

    # 7. Export publication-quality graphic before displaying GUI canvas
    out_file = "05_phase_folded_transit.png"
    print(f"\n[+] Saving plot to {out_file}...")
    plt.savefig(out_file, dpi = 300, bbox_inches = "tight")
    print(f"[+] Saved {out_file} successfully!")

    plt.show()

if __name__ == "__main__":
    main()