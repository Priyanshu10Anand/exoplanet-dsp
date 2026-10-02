"""
Step 6: Physical Metric Calculation & Final Validation Dashboard
- Extracts planet-to-star radius ratio (Rp / R*)
- Computes planetary radius in Earth radii (R_Earth)
- Calculates semi-major axis (AU) and orbital speed (km/s) using Kepler's Third Law
- Generates a 4-panel DSP validation report
"""

import warnings
warnings.filterwarnings(
    "ignore", 
    category = UserWarning, 
    module = "lightkurve"
)
import astropy.constants as const
import astropy.units as u
from astropy.timeseries import BoxLeastSquares
import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Stellar priors for Kepler-10 host star (Batalha et al., NASA Kepler Team)
R_STAR = 1.065 * u.R_sun       # Stellar radius
M_STAR = 0.910 * u.M_sun       # Stellar mass

def clean_and_detrend(lc):
    """Steps 2 & 3: Drop flags, apply asymmetric MAD mask, and flatten with SavGol."""
    # 1. Reject NaN telemetry and nonzero hardware fault flags
    clean = lc.remove_nans()
    clean = clean[clean.quality == 0]

    # 2. Outlier rejection: gate cosmic rays (+4s) while protecting transit dips (-8s)
    flux = clean.flux.value
    med = np.nanmedian(flux)
    mad = np.nanmedian(np.abs(flux - med))
    sigma = 1.4826 * mad
    mask = (flux <= med + 4.0 * sigma) & (flux >= med - 8.0 * sigma)

    # 3. Pre-whiten baseline using a 3.1-day Savitzky-Golay polynomial filter
    return clean[mask].normalize().flatten(window_length = 151, polyorder = 2)

def main():
    # 1. Ingest Quarter 3 calibrated telemetry and isolate residuals
    print("[+] Loading Kepler-10 (Q3) and preparing validation report...")
    raw_lc = lk.search_lightcurve(
        "Kepler-10",
        mission = "Kepler",
        quarter = 3,
        author = "Kepler",
        exptime = 1800
    ).download()

    flat_lc = clean_and_detrend(raw_lc)

    # 2. Matched filter optimization: narrow sweep to isolate fundamental transit peak
    time_val = flat_lc.time.value
    flux_val = flat_lc.flux.value
    bls = BoxLeastSquares(time_val, flux_val)

    periods = np.linspace(0.80, 0.86, 3000)
    durations = np.linspace(0.03, 0.10, 10)
    periodogram = bls.power(periods, durations)

    # 3. Extract optimal transit observables
    best_idx = np.argmax(periodogram.power)
    raw_period = periodogram.period[best_idx]
    period_days = (raw_period.value if hasattr(raw_period, "value") else raw_period) * u.day
    
    raw_t0 = periodogram.transit_time[best_idx]
    t0_days = raw_t0.value if hasattr(raw_t0, "value") else raw_t0

    raw_dur = periodogram.duration[best_idx]
    duration_days = (raw_dur.value if hasattr(raw_dur, "value") else raw_dur) * u.day

    raw_depth = periodogram.depth[best_idx]
    depth = raw_depth.value if hasattr(raw_depth, "value") else raw_depth

    # 4. Geometric translation: transit depth to physical planet radius
    rp_over_rs = np.sqrt(depth)
    r_planet = (rp_over_rs * R_STAR).to(u.R_earth)

    # 5. Orbital mechanics: solve Kepler's Third Law for orbital radius (a)
    period_seconds = period_days.to(u.s)
    m_star_kg = M_STAR.to(u.kg)
    a_cubed = (const.G * m_star_kg * period_seconds**2) / (4 * np.pi**2)
    a_semi_major = (a_cubed**(1/3)).to(u.AU)

    # 6. Kinematics: compute mean circular orbital speed
    v_orbit = ((2 * np.pi * a_semi_major) / period_seconds).to(u.km / u.s)

    # 7. Print physical solution report
    print("\n" + "=" * 50)
    print("       EXOPLANET DETECTION & PHYSICAL REPORT")
    print("=" * 50)
    print(f"Target System:           Kepler-10")
    print(f"Detected Orbital Period: {period_days.value:.5f} days ({period_days.to(u.hour).value:.2f} hours)")
    print(f"Mid-Transit Epoch (t0):  {t0_days:.4f} BJD")
    print(f"Transit Duration (tau):  {duration_days.to(u.hour).value:.2f} hours")
    print(f"Observed Transit Depth:  {depth * 1e6:.1f} ppm ({depth * 100:.4f}%)")
    print("-" * 50)
    print(f"Planet-to-Star Ratio:    Rp / R* = {rp_over_rs:.4f}")
    print(f"Planetary Radius:        {r_planet.value:.2f} Earth Radii (R_earth)")
    print(f"Semi-Major Axis (a):     {a_semi_major.value:.4f} AU (~{a_semi_major.to(u.km).value:.1e} km)")
    print(f"Mean Orbital Velocity:   {v_orbit.value:.1f} km/s")
    print("=" * 50)

    # 8. Render publication-style 4-panel DSP validation dashboard
    fig = plt.figure(
        figsize = (12, 10), 
        constrained_layout = True
    )
    gs = fig.add_gridspec(3, 2)

    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    ax4 = fig.add_subplot(gs[2, :])

    # Panel A: Complete pre-whitened ~90-day time series
    ax1.scatter(
        time_val, 
        flux_val, 
        s = 1, 
        color = "slategray", 
        alpha = 0.6
    )
    ax1.axhline(
        1.0, 
        color = "crimson", 
        linestyle = "--", 
        lw = 0.8
    )
    ax1.set_ylabel("Normalized Flux")
    ax1.set_xlabel("Time (BJD - 2454833)")
    ax1.set_title("A. Quarter 3 Full Pre-Whitened Photometric Baseline (~90 Days)")

    # Panel B: Matched-filter periodogram with detection peak
    period_plot = periodogram.period.value if hasattr(periodogram.period, "value") else periodogram.period
    ax2.plot(
        period_plot, 
        periodogram.power, 
        color = "navy", 
        lw = 1.0
    )  
    ax2.axvline(
        period_days.value, 
        color = "crimson", 
        linestyle = "--", 
        label = f"P = {period_days.value:.4f} d")
    ax2.set_xlabel("Trial Period (Days)")
    ax2.set_ylabel("Filter Power")
    ax2.set_title("B. BLS Matched Filter Periodogram")
    ax2.legend(loc = "upper right")
    ax2.grid(
        True, 
        linestyle = ":", 
        alpha = 0.5
    )

    # Panel C: Synthesized square-wave matched filter pulse template
    phase_dense = np.linspace(-0.15, 0.15, 500)
    box = np.ones_like(phase_dense)
    box[np.abs(phase_dense) <= (duration_days.value / 2.0)] = 1.0 - depth
    ax3.plot(
        phase_dense, 
        box, 
        "crimson", 
        lw = 2, 
        label = "Fitted Matched Box"
    )
    ax3.axvline(
        0.0, 
        color = "black", 
        linestyle = ":", 
        lw = 0.7
    )
    ax3.set_xlabel("Phase (Days from Mid-Transit)")
    ax3.set_ylabel("Relative Flux")
    ax3.set_title("C. Analytical Transit Pulse Template")
    ax3.set_ylim(0.9996, 1.0002)
    ax3.legend(loc = "lower right")

    # Panel D: Coherent phase fold and 10-minute noise-suppression binning
    folded_lc = flat_lc.fold(
        period = period_days.value, 
        epoch_time = t0_days
    )
    binned_lc = folded_lc.bin(time_bin_size = 10.0 / (60.0 * 24.0))

    ax4.scatter(
        folded_lc.time.value, 
        folded_lc.flux.value, 
        s = 1.2, 
        color = "lightsteelblue", 
        alpha = 0.4, 
        label = "Folded Samples"
        )
    ax4.errorbar(
        binned_lc.time.value,
        binned_lc.flux.value,
        yerr = binned_lc.flux_err.value if binned_lc.flux_err is not None else None,
        fmt = "o",
        color = "midnightblue",
        markersize = 3.5,
        elinewidth = 0.8,
        label = "Binned (10 min)"
    )
    ax4.plot(
        phase_dense, 
        box, 
        "crimson", 
        lw = 2.0, 
        label = f"Transit Depth = {depth*1e6:.0f} ppm"
    )
    ax4.set_xlim(-0.25, 0.25)
    ax4.set_ylim(0.9993, 1.0006)
    ax4.set_xlabel("Orbital Phase (Days from Mid-Transit)")
    ax4.set_ylabel("Normalized Flux")
    ax4.set_title(f"D. Coherent Detection: Kepler-10b (Rp = {r_planet.value:.2f} R_Earth, a = {a_semi_major.value:.4f} AU)")
    ax4.legend(loc = "lower right")

    out_file = "06_physical_metrics.png"
    print(f"\n[+] Saving plot to {out_file}...")
    plt.savefig(out_file, dpi = 300, bbox_inches = "tight")
    print(f"[+] Saved {out_file} successfully!")

    plt.show()

if __name__ == "__main__":
    main()