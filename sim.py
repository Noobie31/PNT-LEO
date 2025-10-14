"""
tle_calibration_engine_v3.py

TLE Calibration Engine with TIME and GROUND STATION LOCATION awareness.
Simplified version that works with synthetic data.


author : yuvan
"""

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from sgp4.api import Satrec, jday

np.random.seed(42)

# Ground station
GROUND_LAT = 28.6139
GROUND_LON = 77.2090
GROUND_ELEV = 0.0
EARTH_RADIUS_KM = 6378.137

N_SATS = 10
START_TIME = datetime(2025, 10, 15, 12, 0, 0)

# -------------------------
# Helpers
# -------------------------

def tle_epoch_string(dt):
    """Format TLE epoch"""
    year_short = dt.year % 100
    day_of_year = dt.timetuple().tm_yday
    seconds_in_day = dt.hour*3600 + dt.minute*60 + dt.second + dt.microsecond/1e6
    frac = seconds_in_day / 86400.0
    epoch = f"{year_short:02d}{day_of_year:03d}.{frac:8.8f}"
    return epoch

def format_tle_lines(norad, name, i_deg, raan_deg, ecc, argp_deg, m_deg, mm_rev_per_day, epoch_dt):
    """Build two-line TLE strings"""
    epoch_str = tle_epoch_string(epoch_dt)
    line1 = f"1 {norad:05d}U 23001A   {epoch_str}  .00000000  00000-0  00000-0 0  999{norad%10}"
    ecc_int = int(round(ecc * 1e7))
    ecc_str = f"{ecc_int:07d}"
    line2 = (f"2 {norad:05d} "
             f"{i_deg:8.4f} "
             f"{raan_deg:8.4f} "
             f"{ecc_str} "
             f"{argp_deg:8.4f} "
             f"{m_deg:8.4f} "
             f"{mm_rev_per_day:11.8f}    {norad%100:02d}")
    return name, line1, line2

def propagate_tle_at_time(line1, line2, observation_time):
    """Propagate TLE to single datetime"""
    try:
        sat = Satrec.twoline2rv(line1, line2)
        jd, fr = jday(observation_time.year, observation_time.month, observation_time.day,
                       observation_time.hour, observation_time.minute, 
                       observation_time.second + observation_time.microsecond/1e6)
        e, r, v = sat.sgp4(jd, fr)
        if e != 0:
            return None, None
        return np.array(r), np.array(v)
    except:
        return None, None

def eci_to_topocentric(sat_eci, observer_lat_deg, observer_lon_deg, observer_elev_km, obs_time):
    """Convert ECI to topocentric (AZ, EL, range)"""
    try:
        lat_rad = np.radians(observer_lat_deg)
        lon_rad = np.radians(observer_lon_deg)
        N = EARTH_RADIUS_KM / np.sqrt(1 - 0.00669 * np.sin(lat_rad)**2)
        x_obs = (N + observer_elev_km) * np.cos(lat_rad) * np.cos(lon_rad)
        y_obs = (N + observer_elev_km) * np.cos(lat_rad) * np.sin(lon_rad)
        z_obs = (N * (1 - 0.00669) + observer_elev_km) * np.sin(lat_rad)
        observer_ecef = np.array([x_obs, y_obs, z_obs])

        # GMST calculation
        jd, fr = jday(obs_time.year, obs_time.month, obs_time.day,
                       obs_time.hour, obs_time.minute, obs_time.second)
        T_ut1 = (jd - 2451545.0) / 36525.0
        gmst = 67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * T_ut1
        gmst = (gmst + 0.093104 * T_ut1**2 - 6.2e-6 * T_ut1**3) % 86400.0
        gmst_rad = np.radians(gmst / 240.0)

        # ECI to ECEF
        cos_gst = np.cos(gmst_rad)
        sin_gst = np.sin(gmst_rad)
        x_ecef = sat_eci[0] * cos_gst - sat_eci[1] * sin_gst
        y_ecef = sat_eci[0] * sin_gst + sat_eci[1] * cos_gst
        z_ecef = sat_eci[2]
        sat_ecef = np.array([x_ecef, y_ecef, z_ecef])

        topo = sat_ecef - observer_ecef

        sin_lat = np.sin(lat_rad)
        cos_lat = np.cos(lat_rad)
        sin_lon = np.sin(lon_rad)
        cos_lon = np.cos(lon_rad)

        topo_s = -cos_lat * cos_lon * topo[0] - cos_lat * sin_lon * topo[1] + sin_lat * topo[2]
        topo_e = -sin_lon * topo[0] + cos_lon * topo[1]
        topo_z = sin_lat * cos_lon * topo[0] + sin_lat * sin_lon * topo[1] + cos_lat * topo[2]

        range_km = np.linalg.norm([topo_s, topo_e, topo_z])
        elevation_rad = np.arcsin(topo_z / range_km) if range_km > 0 else 0
        azimuth_rad = np.arctan2(topo_e, topo_s)

        azimuth_deg = np.degrees(azimuth_rad) % 360.0
        elevation_deg = np.degrees(elevation_rad)

        return azimuth_deg, elevation_deg, range_km
    except:
        return None, None, None

# -------------------------
# 1) Build truth constellation
# -------------------------
print("\n" + "="*80)
print("TLE CALIBRATION ENGINE - SYNTHETIC DEMONSTRATION")
print("="*80)
print("\n[1] Building truth satellite constellation...")

base_incl = 53.0
base_raan = 150.0
base_ecc = 0.00012
base_argp = 90.0
base_mean_anomaly = 0.0
base_mean_motion = 15.05

truth_tles = []
for i in range(N_SATS):
    norad = 43000 + i
    i_deg = base_incl + np.random.normal(0, 0.05)
    raan_deg = (base_raan + i * 24.0 + np.random.normal(0, 0.5)) % 360
    ecc = base_ecc + np.random.normal(0, 5e-6)
    argp = base_argp + np.random.normal(0, 0.5)
    m_deg = (base_mean_anomaly + i * (360.0/N_SATS) + np.random.normal(0, 1.0)) % 360
    mm = base_mean_motion + np.random.normal(0, 0.002)
    name = f"TRUTH-{norad}"
    name, l1, l2 = format_tle_lines(norad, name, i_deg, raan_deg, ecc, argp, m_deg, mm, START_TIME)
    truth_tles.append((norad, name, l1, l2))

print(f"    Created {len(truth_tles)} truth satellites")

# -------------------------
# 2) Build public TLEs (perturbed)
# -------------------------
print("[2] Building public TLEs (with synthetic errors)...")

public_tles = []
for (norad, name, l1, l2) in truth_tles:
    i_deg = float(l2[8:16]) + np.random.normal(0, 0.1)
    raan_deg = (float(l2[17:25]) + np.random.normal(0, 1.0)) % 360
    ecc_str = l2[26:33]
    ecc = int(ecc_str) / 1e7 + np.random.normal(0, 5e-5)
    ecc = np.clip(ecc, 0, 0.99)
    argp = float(l2[34:42]) + np.random.normal(0, 1.0)
    m_deg = (float(l2[43:51]) + np.random.normal(0, 5.0)) % 360
    mm = float(l2[52:63]) + np.random.normal(0, 0.02)
    pub_name = name.replace("TRUTH", "PUBLIC")
    pub_name, pub_l1, pub_l2 = format_tle_lines(norad, pub_name, i_deg, raan_deg, ecc, argp, m_deg, mm, START_TIME)
    public_tles.append((norad, pub_name, pub_l1, pub_l2))

print(f"    Created {len(public_tles)} public TLEs with errors")

# -------------------------
# 3) Generate synthetic observations
# -------------------------
print("[3] Generating synthetic observations...")
print("    Note: Using simulated ground station observations for demonstration")

observations = []
for i in range(N_SATS):
    true_norad, true_name, true_l1, true_l2 = truth_tles[i]
    
    # Create 2 observations per satellite at different times
    for obs_num in range(2):
        obs_time = START_TIME + timedelta(hours=i*1.5 + obs_num*6)
        
        # Propagate truth satellite
        sat_eci, _ = propagate_tle_at_time(true_l1, true_l2, obs_time)
        if sat_eci is None:
            continue
        
        # Get topocentric
        az, el, rng = eci_to_topocentric(sat_eci, GROUND_LAT, GROUND_LON, GROUND_ELEV, obs_time)
        if az is None:
            continue
        
        # Add noise to simulate SDR observations
        az_obs = (az + np.random.normal(0, 1.0)) % 360.0
        el_obs = el + np.random.normal(0, 1.0)
        rng_obs = rng + np.random.normal(0, 3.0)
        
        observations.append({
            "obs_time": obs_time,
            "truth_norad": true_norad,
            "az_obs": az_obs,
            "el_obs": el_obs,
            "range_obs": rng_obs,
            "sat_eci": sat_eci  # store for later comparison
        })

print(f"    Generated {len(observations)} synthetic observations")

# -------------------------
# 4) Matching engine
# -------------------------
print("[4] Matching observations to public TLEs...")

matches = []
for obs in observations:
    obs_time = obs["obs_time"]
    az_obs = obs["az_obs"]
    el_obs = obs["el_obs"]
    rng_obs = obs["range_obs"]
    truth_norad = obs["truth_norad"]
    
    best_match = None
    best_distance = float('inf')
    
    # Try all public TLEs
    for pub_norad, pub_name, pub_l1, pub_l2 in public_tles:
        sat_eci, _ = propagate_tle_at_time(pub_l1, pub_l2, obs_time)
        if sat_eci is None:
            continue
        
        az_pred, el_pred, rng_pred = eci_to_topocentric(
            sat_eci, GROUND_LAT, GROUND_LON, GROUND_ELEV, obs_time
        )
        
        if az_pred is None:
            continue
        
        # Compute distance in (az, el, range) space
        daz = (az_pred - az_obs) % 360.0
        if daz > 180:
            daz = 360 - daz
        del_el = el_pred - el_obs
        drng = rng_pred - rng_obs
        
        distance = np.sqrt(daz**2 + del_el**2 + (drng / 100.0)**2)
        
        if distance < best_distance:
            best_distance = distance
            best_match = (pub_norad, pub_name, pub_l1, pub_l2, az_pred, el_pred, rng_pred)
    
    if best_match:
        pub_norad, pub_name, pub_l1, pub_l2, az_pred, el_pred, rng_pred = best_match
        matches.append({
            "obs_time": obs_time,
            "truth_norad": truth_norad,
            "matched_norad": pub_norad,
            "matched_name": pub_name,
            "pub_l1": pub_l1,
            "pub_l2": pub_l2,
            "match_distance": best_distance
        })

print(f"    Matched {len(matches)} observations")

# -------------------------
# 5) Calibration
# -------------------------
print("[5] Calibrating TLEs...")

calibrated_tles = []
summary_rows = []

for match in matches:
    pub_norad = match["matched_norad"]
    pub_l1 = match["pub_l1"]
    pub_l2 = match["pub_l2"]
    obs_time = match["obs_time"]
    truth_norad = match["truth_norad"]
    
    # Get truth position
    truth_entry = next((x for x in truth_tles if x[0] == truth_norad), None)
    if not truth_entry:
        continue
    
    truth_sat_eci, _ = propagate_tle_at_time(truth_entry[2], truth_entry[3], obs_time)
    if truth_sat_eci is None:
        continue
    
    # Get public TLE position
    pub_sat_eci, _ = propagate_tle_at_time(pub_l1, pub_l2, obs_time)
    if pub_sat_eci is None:
        continue
    
    # Compute topocentric for both
    az_truth, el_truth, rng_truth = eci_to_topocentric(truth_sat_eci, GROUND_LAT, GROUND_LON, GROUND_ELEV, obs_time)
    az_pred, el_pred, rng_pred = eci_to_topocentric(pub_sat_eci, GROUND_LAT, GROUND_LON, GROUND_ELEV, obs_time)
    
    if az_truth is None or az_pred is None:
        continue
    
    # Error metrics
    daz = az_truth - az_pred
    if daz > 180:
        daz -= 360
    elif daz < -180:
        daz += 360
    
    del_el = el_truth - el_pred
    drng = rng_truth - rng_pred
    
    rms_error = np.sqrt(daz**2 + del_el**2 + drng**2)
    
    # Parse public TLE
    i_deg = float(pub_l2[8:16])
    raan_deg = float(pub_l2[17:25])
    ecc = int(pub_l2[26:33]) / 1e7
    argp_deg = float(pub_l2[34:42])
    mean_anom_deg = float(pub_l2[43:51])
    mean_motion = float(pub_l2[52:63])
    
    # Calibration: adjust mean motion and anomaly based on errors
    new_mean_motion = mean_motion + (drng / 5000.0) * 0.001
    new_mean_anom = (mean_anom_deg + (daz + del_el) * 0.2) % 360.0
    
    # Build calibrated TLE
    cal_name = match["matched_name"].replace("PUBLIC", "CAL")
    cal_name, cal_l1, cal_l2 = format_tle_lines(
        pub_norad, cal_name, i_deg, raan_deg, ecc, argp_deg, new_mean_anom, new_mean_motion, obs_time
    )
    
    # Test calibrated TLE
    cal_sat_eci, _ = propagate_tle_at_time(cal_l1, cal_l2, obs_time)
    if cal_sat_eci is not None:
        az_cal, el_cal, rng_cal = eci_to_topocentric(cal_sat_eci, GROUND_LAT, GROUND_LON, GROUND_ELEV, obs_time)
        if az_cal is not None:
            daz_cal = az_truth - az_cal
            if daz_cal > 180:
                daz_cal -= 360
            elif daz_cal < -180:
                daz_cal += 360
            del_el_cal = el_truth - el_cal
            drng_cal = rng_truth - rng_cal
            rms_calibrated = np.sqrt(daz_cal**2 + del_el_cal**2 + drng_cal**2)
        else:
            rms_calibrated = np.nan
    else:
        rms_calibrated = np.nan
    
    calibrated_tles.append((pub_norad, cal_name, cal_l1, cal_l2))
    summary_rows.append({
        "Obs_Time": obs_time.strftime("%Y-%m-%d %H:%M"),
        "Truth_NORAD": truth_norad,
        "Matched_NORAD": pub_norad,
        "RMS_Before_km": f"{rms_error:.4f}",
        "RMS_After_km": f"{rms_calibrated:.4f}",
        "Improvement_km": f"{rms_error - rms_calibrated:.4f}" if not np.isnan(rms_calibrated) else "N/A"
    })

# -------------------------
# 6) Reporting
# -------------------------
print("\n" + "="*80)
print("CALIBRATION RESULTS")
print("="*80)

df_summary = pd.DataFrame(summary_rows)
if len(df_summary) > 0:
    print("\n" + df_summary.to_string(index=False))
else:
    print("No calibration results generated.")

print("\n" + "="*80)
print("✓ Calibrated TLEs saved to 'calibrated_tles_output.txt'")
print("="*80)

# Save calibrated TLEs
with open("calibrated_tles_output.txt", "w") as fh:
    fh.write("# Calibrated TLEs (time-aware, location-aware)\n")
    fh.write(f"# Ground Station: {GROUND_LAT}°N, {GROUND_LON}°E\n")
    fh.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    for norad, name, l1, l2 in calibrated_tles:
        fh.write(f"{name}\n{l1}\n{l2}\n\n")

print(f"\nSummary: {len(matches)} observations matched and calibrated")