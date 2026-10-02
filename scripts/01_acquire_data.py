"""
Step 1: Environment & Data Acquisition
Target: Kepler-10 (Quarter 3)
Mission: NASA Kepler Space Telescope
"""

import warnings
warnings.filterwarnings(
    "ignore", 
    category = UserWarning, 
    module = "lightkurve"
)
import astropy.units as u
import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def main():
    # 1. Query MAST archive for 30-min long-cadence observations
    print("[+] Querying NASA MAST archive for Kepler-10 (Quarter 3)...")
    search_q3 = lk.search_lightcurve(
        "Kepler-10", 
        mission = "Kepler", 
        quarter = 3, 
        author = "Kepler",
        exptime = 1800
    )
    print(search_q3)

    # 2. Download calibrated FITS data product
    print("\n[+] Downloading single Long Cadence Q3 file...")
    lc = search_q3.download()

    # 3. Compute empirical sampling interval across cadences
    dt_days = np.nanmedian(np.diff(lc.time.value))
    dt_minutes = (dt_days * u.day).to(u.minute)

    # 4. Display telemetry parameters and coordinate baseline
    print("\n--- LightCurve Telemetry Summary ---")
    print(f"Time format:       {lc.time.format} ({lc.time.scale})")
    print(f"Total data points: {len(lc)}")
    print(f"Sampling cadence:  ~{dt_minutes:.1f}")
    print(f"Baseline span:     {lc.time.value[-1] - lc.time.value[0]:.2f} days")
    print(f"Available columns: {lc.colnames}")

    # 5. Plot spacecraft-corrected flux (PDCSAP removes pointing jitter & thermal drift)
    fig, ax = plt.subplots(figsize=(11, 4))
    lc.plot(
        ax = ax,
        column = "pdcsap_flux",
        label = "PDCSAP Flux (Pre-search Data Conditioned)",
        color = "black",
        lw = 0.6
    )
    ax.set_title("Kepler-10 Quarter 3 Raw Photometric Time Series")
    ax.set_ylabel("Flux (e⁻ / s)")
    ax.set_xlabel("Time (BJD - 2454833)")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()