"""
tle_calibration_mvp.py

MVP prototype:
- Creates 15 "public" (slightly erroneous) TLEs for LEO satellites (Starlink-like)
- Creates 15 "truth" orbits (the simulated real LEO constellation)
- Simulates SDR observations (position vectors) for satellite passes over a ground station
- Matches observations to public TLEs (RMS distance)
- Performs a simple calibration (adjust mean motion and mean anomaly) to produce "calibrated" TLEs
- Reports RMS before/after and writes calibrated TLEs to output

Author: ChatGPT (prototype for user)
Date: 2025-10-15
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from sgp4.api import Satrec, jday
import math
import sys

# -------------------------
# Configuration
# -------------------------
np.random.seed(42)  # reproducible
N_SATS = 15
GROUND_LAT = 28.6139
GROUND_LON = 77.2090  # New Delhi
START_TIME = datetime(2025, 10, 15, 12, 0, 0)  # epoch for simulation
NUM_OBS_PER_SAT = 80  # number of discrete observation timestamps to sample (across window)
OBS_INTERVAL_SEC = 30  # seconds between observation timestamps
NOISE_SIGMA_KM = 2.0  # SDR position noise (km) (typical LEO observation noise for demo)

# Physical constants
EARTH_RADIUS_KM = 6378.137

# -------------------------
# Helpers
# -------------------------
def tle_epoch_string(dt: datetime):
    """
    Format TLE epoch as YYDDD.FFFFF where DDD is day-of-year and fraction.
    Returns string like '25288.50000000' where '25' -> 2025 and 288th day.
    """
    year_short = dt.year % 100
    day_of_year = dt.timetuple().tm_yday
    # fraction of day
    seconds_in_day = dt.hour*3600 + dt.minute*60 + dt.second + dt.microsecond/1e6
    frac = seconds_in_day / 86400.0
    epoch = f"{year_short:02d}{day_of_year:03d}.{frac:8.8f}"
    return epoch

def format_tle_lines(norad, name, i_deg, raan_deg, ecc, argp_deg, m_deg, mm_rev_per_day, epoch_dt):
    """
    Build simple two-line TLE strings (synthetic) with given classical elements.
    Note: This is a formatted TLE string intended for use with sgp4.twoline2rv.
    It attempts to respect column widths but is still synthetic (OK for prototyping).
    """
    # TLE Line 1 components
    epoch_str = tle_epoch_string(epoch_dt)
    line1 = f"1 {norad:05d}U 23001A   {epoch_str}  .00000000  00000-0  00000-0 0  999{norad%10}"

    # TLE Line 2: inclination (8.4), raan (8.4), ecc (7 digits no decimal), arg perigee (8.4), mean anomaly (8.4), mean motion (11.8)
    ecc_str = f"{int(round(ecc * 1e7)):07d}"  # ecc as 7-digit integer in TLE format
    line2 = (f"2 {norad:05d} "
             f"{i_deg:8.4f} "
             f"{raan_deg:8.4f} "
             f"{ecc_str} "
             f"{argp_deg:8.4f} "
             f"{m_deg:8.4f} "
             f"{mm_rev_per_day:11.8f}    {norad%100:02d}")
    return name, line1, line2

def propagate_positions_from_tle(line1, line2, times):
    """Propagate a Satrec built from TLE lines to produce ECI position vectors (km) for given datetimes."""
    sat = Satrec.twoline2rv(line1, line2)
    positions = []
    errors = []
    for t in times:
        jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute, t.second + t.microsecond/1e-6)
        e, r, v = sat.sgp4(jd, fr)
        if e != 0:
            # propagate failure: fill NaNs
            positions.append([np.nan, np.nan, np.nan])
            errors.append(e)
        else:
            positions.append(r)
            errors.append(0)
    return np.array(positions), np.array(errors)

def rms_of_differences(a, b):
    """Compute RMS distance between arrays of shape (N,3). Ignores NaNs pairwise."""
    diff = a - b
    # mask rows with any nan
    valid = ~np.isnan(diff).any(axis=1)
    if valid.sum() == 0:
        return np.nan
    norms = np.linalg.norm(diff[valid], axis=1)
    return np.sqrt(np.mean(norms**2))

# -------------------------
# 1) Build "truth" constellation (N_SATS) -- slightly ideal orbits
# -------------------------
base_incl = 53.0  # degrees
base_raan = 150.0  # degrees
base_ecc = 0.00012
base_argp = 90.0
base_mean_anomaly = 0.0
base_mean_motion = 15.05  # revs per day ~ typical LEO

truth_tles = []
for i in range(N_SATS):
    norad = 43000 + i
    # introduce small systematic spacing in RAAN and mean anomaly to mimic plane / phasing
    i_deg = base_incl + np.random.normal(0, 0.05)  # tiny inclination jitter
    raan_deg = (base_raan + i * 2.0 + np.random.normal(0, 0.5)) % 360
    ecc = base_ecc + np.random.normal(0, 5e-6)
    argp = base_argp + np.random.normal(0, 0.5)
    m_deg = (base_mean_anomaly + i * (360.0/N_SATS) + np.random.normal(0, 1.0)) % 360
    mm = base_mean_motion + np.random.normal(0, 0.002)  # small mean motion jitter
    name = f"TRUTH-{norad}"
    name, l1, l2 = format_tle_lines(norad, name, i_deg, raan_deg, ecc, argp, m_deg, mm, START_TIME)
    truth_tles.append((norad, name, l1, l2))

# -------------------------
# 2) Build "public" TLEs (slightly erroneous) that the engine will receive
#    (public TLEs are perturbed versions of truth)
# -------------------------
public_tles = []
for (norad, name, l1, l2) in truth_tles:
    # perturb the TLE elements more strongly to simulate "public" inaccuracies
    # parse line2 components we wrote earlier (we kept them in consistent format)
    # We'll extract numeric fields by splitting and re-build with perturbation.
    parts = l2.split()
    i_deg = float(parts[1]) + np.random.normal(0, 0.05) * 5.0  # larger jitter
    raan_deg = (float(parts[2]) + np.random.normal(0, 0.2) * 5.0) % 360
    # ecc was encoded; convert back:
    ecc = int(parts[3]) / 1e7 + np.random.normal(0, 2e-5)
    argp = float(parts[4]) + np.random.normal(0, 0.5) * 2.0
    m_deg = (float(parts[5]) + np.random.normal(0, 2.0)) % 360
    mm = float(parts[6]) + np.random.normal(0, 0.01)  # public TLE has up to 0.01 rev/day error in this synthetic demo
    pub_name = name.replace("TRUTH", "PUBLIC")
    pub_name, pub_l1, pub_l2 = format_tle_lines(norad, pub_name, i_deg, raan_deg, ecc, argp, m_deg, mm, START_TIME)
    public_tles.append((norad, pub_name, pub_l1, pub_l2))

# -------------------------
# 3) Create observation timestamps and "SDR" observations
#    We'll sample a wide window so many passes are present; the "truth" orbits produce the 'observations'
# -------------------------
# Build a global timestamp list for the observation window
total_obs_points = NUM_OBS_PER_SAT
times = [START_TIME + timedelta(seconds=OBS_INTERVAL_SEC * t) for t in range(total_obs_points)]

# For each truth satellite, propagate to times and add noise (SDR observations)
truth_positions_all = {}
obs_positions_all = {}

for idx, (norad, name, l1, l2) in enumerate(truth_tles):
    positions, errs = propagate_positions_from_tle(l1, l2, times)
    # Add SDR noise to simulated observed position (km)
    noise = np.random.normal(0, NOISE_SIGMA_KM, size=positions.shape)
    obs_positions = positions + noise
    truth_positions_all[norad] = positions
    obs_positions_all[norad] = obs_positions

# -------------------------
# 4) Matching engine:
#    For each observed "sat pass" (here we treat each truth sat obs as an observation to be matched),
#    we will compare predicted tracks from all public TLEs and pick the best-match TLE by RMS.
# -------------------------
match_results = []  # list of dicts with norad_truth, best_public_norad, rms_before

for true_norad, obs_positions in obs_positions_all.items():
    best_rms = None
    best_pub_norad = None
    best_pub_entry = None
    # For matching, use the same set of times to predict public TLE positions and compute RMS
    for pub_norad, pub_name, pub_l1, pub_l2 in public_tles:
        pred_pub_pos, errp = propagate_positions_from_tle(pub_l1, pub_l2, times)
        rms = rms_of_differences(obs_positions, pred_pub_pos)
        if np.isnan(rms):
            continue
        if best_rms is None or rms < best_rms:
            best_rms = rms
            best_pub_norad = pub_norad
            best_pub_entry = (pub_norad, pub_name, pub_l1, pub_l2)
    match_results.append({
        "truth_norad": true_norad,
        "matched_public_norad": best_pub_norad,
        "rms_before_km": best_rms,
        "public_entry": best_pub_entry
    })

# -------------------------
# 5) Calibration engine (very simple heuristic)
#    For each matched public TLE:
#      - compute RMS before,
#      - compute simple correction:
#         * adjust mean motion by a small fraction proportional to RMS error
#         * adjust mean anomaly by projecting the mean difference vector onto the velocity direction (approx)
#      - produce new TLE lines (calibrated)
# -------------------------
calibrated_tles = []
summary_rows = []

for res in match_results:
    true_norad = res["truth_norad"]
    pub_entry = res["public_entry"]
    if pub_entry is None:
        # matching failed: skip
        continue
    pub_norad, pub_name, pub_l1, pub_l2 = pub_entry
    # propagate public predicted and truth
    pred_pub_pos, _ = propagate_positions_from_tle(pub_l1, pub_l2, times)
    truth_pos = truth_positions_all[true_norad]
    rms_before = rms_of_differences(truth_pos, pred_pub_pos)

    # Estimate correction factors
    # simple rule: fraction = (rms_before / (EARTH_RADIUS_KM + 500 km)) clipped
    orbit_radius_mean = EARTH_RADIUS_KM + 500.0  # typical LEO radius used for scaling
    fraction = float(rms_before) / float(orbit_radius_mean) if not np.isnan(rms_before) else 0.0
    fraction = np.clip(fraction, 0.0, 0.05)  # don't change mean motion by >5% in this toy demo

    # parse public TLE line2 to extract mean motion and mean anomaly
    parts = pub_l2.split()
    i_deg = float(parts[1])
    raan_deg = float(parts[2])
    ecc = int(parts[3]) / 1e7
    argp_deg = float(parts[4])
    mean_anom_deg = float(parts[5])
    mean_motion = float(parts[6])

    # Correct mean motion (simple proportional correction: increase slightly if errors positive)
    # we pick a small signed adjustment based on average radial projection sign:
    # compute average difference vector over times
    diffs = truth_pos - pred_pub_pos
    valid_mask = ~np.isnan(diffs).any(axis=1)
    delta_sign = 0.0
    if valid_mask.sum() > 0:
        avg_diff = np.mean(diffs[valid_mask], axis=0)
        avg_norm = np.linalg.norm(avg_diff)
        if avg_norm > 0:
            # sign with projection along position vector (approx along-track if dot with velocity > 0)
            # approximate by checking dot product with position (if truth is ahead of pred, mean anomaly should increase)
            avg_pos = np.mean(pred_pub_pos[valid_mask], axis=0)
            dot = np.dot(avg_pos, avg_diff)
            delta_sign = np.sign(dot)
    # correction magnitude
    mm_correction = mean_motion * (fraction * 0.2 * delta_sign)  # scaled small
    new_mean_motion = mean_motion + mm_correction

    # adjust mean anomaly slightly opposite to along-track offset
    # approximate ΔM (deg) = -k * avg_alongtrack_km / orbit_radius_km * (180/pi)
    k_anom = 0.7
    avg_alongtrack_km = 0.0
    if valid_mask.sum() > 0:
        # project avg_diff onto avg_pos unit to get approx radial; for along-track we do cross heuristic
        avg_diff = np.mean(diffs[valid_mask], axis=0)
        avg_pos = np.mean(pred_pub_pos[valid_mask], axis=0)
        avg_pos_unit = avg_pos / np.linalg.norm(avg_pos)
        # alongtrack approx via cross: assume velocity ~ perpendicular to position in circular orbit
        # so alongtrack component approx = component perpendicular to radial -> use difference minus radial projection
        radial_comp = np.dot(avg_diff, avg_pos_unit) * avg_pos_unit
        alongtrack_vec = avg_diff - radial_comp
        avg_alongtrack_km = np.linalg.norm(alongtrack_vec)
        # sign for anomaly change: if alongtrack_vec dot (some approx velocity direction) positive, increase M
        # approximate velocity direction as cross(k, r) where k = [0,0,1]; so v approx [-y, x, 0]
        approx_vel = np.array([-avg_pos[1], avg_pos[0], 0.0])
        sign_anom = np.sign(np.dot(alongtrack_vec, approx_vel))
        deltaM_deg = -k_anom * (avg_alongtrack_km / orbit_radius_mean) * (180.0 / math.pi) * sign_anom
    else:
        deltaM_deg = 0.0

    new_mean_anom = (mean_anom_deg + deltaM_deg) % 360.0

    # Build new calibrated TLE lines
    cal_name = pub_name.replace("PUBLIC", "CAL")
    cal_name, cal_l1, cal_l2 = format_tle_lines(pub_norad, cal_name, i_deg, raan_deg, ecc, argp_deg, new_mean_anom, new_mean_motion, START_TIME)

    # compute RMS after applying calibrated TLE
    cal_pos, _ = propagate_positions_from_tle(cal_l1, cal_l2, times)
    rms_after = rms_of_differences(truth_pos, cal_pos)

    calibrated_tles.append((pub_norad, cal_name, cal_l1, cal_l2, rms_before, rms_after))
    summary_rows.append({
        "truth_norad": true_norad,
        "public_norad": pub_norad,
        "rms_before_km": rms_before,
        "rms_after_km": rms_after,
        "delta_rms_km": (rms_before - rms_after) if not (np.isnan(rms_before) or np.isnan(rms_after)) else np.nan
    })

# -------------------------
# 6) Reporting and outputs
# -------------------------
df_summary = pd.DataFrame(summary_rows)
df_summary = df_summary.sort_values("truth_norad").reset_index(drop=True)

print("\nSummary (first rows):")
print(df_summary.head(20).to_string(index=False))

# Save calibrated TLEs to file
with open("calibrated_tles_output.txt", "w") as fh:
    fh.write("# Calibrated TLEs generated by MVP prototype\n\n")
    for norad, name, l1, l2, rms_b, rms_a in calibrated_tles:
        fh.write(f"{name}\n")
        fh.write(f"{l1}\n")
        fh.write(f"{l2}\n\n")
print("\nCalibrated TLEs written to 'calibrated_tles_output.txt'")

# Show improvement statistics
valid = df_summary[~df_summary["rms_before_km"].isna()]
if len(valid) > 0:
    mean_before = valid["rms_before_km"].mean()
    mean_after = valid["rms_after_km"].mean()
    print(f"\nMean RMS before calibration: {mean_before:.3f} km")
    print(f"Mean RMS after calibration:  {mean_after:.3f} km")
    print(f"Mean improvement (km):       {mean_before - mean_after:.3f} km")

# -------------------------
# 7) Visualize a few sample satellites (predicted vs observed vs calibrated)
# -------------------------
plot_examples = min(3, len(calibrated_tles))
fig = plt.figure(figsize=(12, 4 * plot_examples))
for p in range(plot_examples):
    norad, cal_name, cal_l1, cal_l2, rms_b, rms_a = calibrated_tles[p]
    # find truth norad matching in summary
    rec = df_summary[df_summary["public_norad"] == norad]
    if rec.empty:
        continue
    truth_norad = int(rec["truth_norad"].iloc[0])
    truth_pos = truth_positions_all[truth_norad]
    obs_pos = obs_positions_all[truth_norad]
    pub_entry = next((x for x in public_tles if x[0] == norad), None)
    if pub_entry is None:
        continue
    _, _, pub_l1, pub_l2 = pub_entry
    pred_pub_pos, _ = propagate_positions_from_tle(pub_l1, pub_l2, times)
    cal_pos, _ = propagate_positions_from_tle(cal_l1, cal_l2, times)

    ax = fig.add_subplot(plot_examples, 1, p+1, projection='3d')
    ax.plot(truth_pos[:,0], truth_pos[:,1], truth_pos[:,2], label='Truth Orbit (noiseless)', linewidth=1)
    ax.plot(obs_pos[:,0], obs_pos[:,1], obs_pos[:,2], label='SDR Observations (noisy)', linestyle=':', linewidth=1)
    ax.plot(pred_pub_pos[:,0], pred_pub_pos[:,1], pred_pub_pos[:,2], label='Public TLE pred', linestyle='--', linewidth=1)
    ax.plot(cal_pos[:,0], cal_pos[:,1], cal_pos[:,2], label='Calibrated TLE pred', linestyle='-.', linewidth=1)
    ax.set_title(f"Satellite NORAD {norad} — RMS before {rms_b:.3f} km, after {rms_a:.3f} km")
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.legend()

plt.tight_layout()
plt.show()
