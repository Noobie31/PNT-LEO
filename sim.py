import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sgp4.api import Satrec, jday
from datetime import datetime, timedelta

# -------------------------------------------------------
# 1️⃣ Generate Fake Public TLE
# -------------------------------------------------------
tle_name = "FAKE-STARLINK-1234"
line1 = "1 99999U 23001A   25288.45729167  .00000000  00000-0  00000-0 0  9994"
line2 = "2 99999  53.0000 150.0000 0001200  90.0000 270.0000 15.05420000    15"

sat = Satrec.twoline2rv(line1, line2)

# -------------------------------------------------------
# 2️⃣ Simulate SDR Observations (with drift)
# -------------------------------------------------------
ground_lat, ground_lon = 28.6139, 77.2090   # Example: New Delhi
num_points = 60
start_time = datetime(2025, 10, 15, 12, 0, 0)

timestamps = [start_time + timedelta(seconds=i*60) for i in range(num_points)]

obs_positions = []
pred_positions = []

for t in timestamps:
    jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute, t.second)
    e, r, v = sat.sgp4(jd, fr)
    if e == 0:
        pred_positions.append(r)
        # simulate SDR observation with drift
        noise = np.random.normal(0, 5, size=3)  # 5 km noise
        obs_positions.append(r + noise)
    else:
        pred_positions.append([np.nan, np.nan, np.nan])
        obs_positions.append([np.nan, np.nan, np.nan])

pred_positions = np.array(pred_positions)
obs_positions = np.array(obs_positions)

# -------------------------------------------------------
# 3️⃣ Compute deviation & RMS
# -------------------------------------------------------
diffs = obs_positions - pred_positions
rms_error = np.nanmean(np.linalg.norm(diffs, axis=1))
print(f"[INFO] RMS position error between TLE prediction and SDR observation: {rms_error:.2f} km")

# -------------------------------------------------------
# 4️⃣ Simple Orbit Correction Logic
# -------------------------------------------------------
# Here we estimate a small mean motion correction proportional to observed drift
# (In real implementation, this would be a least-squares orbit determination)
drift_factor = np.clip(rms_error / 1000, 0, 0.01)
original_mean_motion = 15.05420000
corrected_mean_motion = original_mean_motion * (1 + drift_factor / 10)

line2_new = line2[:52] + f"{corrected_mean_motion:11.8f}" + line2[63:]

# Generate new calibrated TLE
cal_name = tle_name + "_CAL"
cal_line1 = line1
cal_line2 = line2_new

print("\n[INFO] New Calibrated TLE:")
print(cal_name)
print(cal_line1)
print(cal_line2)

# -------------------------------------------------------
# 5️⃣ Visualize Orbits
# -------------------------------------------------------
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

ax.plot(pred_positions[:,0], pred_positions[:,1], pred_positions[:,2], label="Predicted Orbit (TLE)", color='blue')
ax.plot(obs_positions[:,0], obs_positions[:,1], obs_positions[:,2], label="Observed SDR Data", color='red', linestyle='dotted')

# Simulate recalibrated orbit with corrected mean motion
sat_cal = Satrec.twoline2rv(cal_line1, cal_line2)
cal_positions = []
for t in timestamps:
    jd, fr = jday(t.year, t.month, t.day, t.hour, t.minute, t.second)
    e, r, v = sat_cal.sgp4(jd, fr)
    if e == 0:
        cal_positions.append(r)
    else:
        cal_positions.append([np.nan, np.nan, np.nan])
cal_positions = np.array(cal_positions)

ax.plot(cal_positions[:,0], cal_positions[:,1], cal_positions[:,2], label="Calibrated Orbit", color='green', linestyle='--')

ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.legend()
ax.set_title("LEO Orbit Calibration Engine Simulation")
plt.show()
