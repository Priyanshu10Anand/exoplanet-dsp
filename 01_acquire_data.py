"""
Step 1: Environment & Data Acquisition
Target: Kepler-10 (Quarter 3)
Mission: NASA Kepler Space Telescope
"""

import astropy.units as u
import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np

def main():
    print("[+] Querying NASA MAST archive for Kepler-10 (Quarter 3)...")
    # Filter explicitly by exptime=1800 (Long Cadence, 30 min) to get the full quarter
    search_q3 = lk.search_lightcurve(
        "Kepler-10", 
        mission="Kepler", 
        quarter=3, 
        author="Kepler",
        exptime=1800
    )
    print(search_q3)

    print("\n[+] Downloading single Long Cadence Q3 file...")
    lc = search_q3.download()

    # Calculate sampling cadence using numpy difference on time values (in days)
    dt_days = np.nanmedian(np.diff(lc.time.value))
    dt_minutes = (dt_days * u.day).to(u.minute)

    # Telemetry and metadata inspection
    print("\n--- LightCurve Telemetry Summary ---")
    print(f"Time format:       {lc.time.format} ({lc.time.scale})")
    print(f"Total data points: {len(lc)}")
    print(f"Sampling cadence:  ~{dt_minutes:.1f}")
    print(f"Baseline span:     {lc.time.value[-1] - lc.time.value[0]:.2f} days")
    print(f"Available columns: {lc.colnames}")

    # Plot raw instrument-corrected flux (PDCSAP)
    fig, ax = plt.subplots(figsize=(11, 4))
    lc.plot(
        ax=ax,
        column="pdcsap_flux",
        label="PDCSAP Flux (Pre-search Data Conditioned)",
        color="black",
        lw=0.6
    )
    ax.set_title("Kepler-10 Quarter 3 Raw Photometric Time Series")
    ax.set_ylabel("Flux (e⁻ / s)")
    ax.set_xlabel("Time (BJD - 2454833)")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()