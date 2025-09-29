import math

# --- Constants ---
mu = 398600.4418  # km^3/s^2
Re = 6378.137  # km

# --- Orbit parameters (public TLE made-up) ---
alt_km = 550.0
a = Re + alt_km
n = math.sqrt(mu / a**3)  # rad/s

incl_deg = 53.0
raan_deg = 120.0
argp_deg = 0.0

# 5 satellites along-track
num_sats = 5
mean_anomaly_offsets_deg = [i*2 for i in range(num_sats)]  # small spacing

# Epochs (seconds)
epochs = [0, 60, 120, 180, 240, 300]

def rot_x(vec, ang_rad):
    x, y, z = vec
    ca = math.cos(ang_rad); sa = math.sin(ang_rad)
    return (x, ca*y - sa*z, sa*y + ca*z)

def rot_z(vec, ang_rad):
    x, y, z = vec
    ca = math.cos(ang_rad); sa = math.sin(ang_rad)
    return (ca*x - sa*y, sa*x + ca*y, z)

def orbital_to_eci(a, incl_rad, raan_rad, argp_rad, M_rad):
    nu = M_rad
    r_pf = (a*math.cos(nu), a*math.sin(nu), 0.0)
    r1 = rot_z(r_pf, argp_rad)
    r2 = rot_x(r1, incl_rad)
    r3 = rot_z(r2, raan_rad)
    return r3

# --- Generate public TLE positions ---
public_tle = []
public_positions = []
for i in range(num_sats):
    M0 = math.radians(mean_anomaly_offsets_deg[i])
    mean_motion_rpd = n*86400/(2*math.pi)
    public_tle.append({
        "sat_id": f"SAT{i+1}",
        "a_km": a,
        "incl_deg": incl_deg,
        "raan_deg": raan_deg,
        "argp_deg": argp_deg,
        "mean_anomaly_deg": mean_anomaly_offsets_deg[i],
        "mean_motion_rev_per_day": mean_motion_rpd
    })
    positions = []
    for t in epochs:
        M_t = M0 + n*t
        r = orbital_to_eci(a, math.radians(incl_deg), math.radians(raan_deg), math.radians(argp_deg), M_t)
        positions.append(r)
    public_positions.append(positions)

# --- Simulate drifted satellite positions (SAT1) ---
drift_rpd = 0.0005  # rev/day drift
delta_n = drift_rpd*2*math.pi/86400
drift_positions = []
for t in epochs:
    M0 = math.radians(mean_anomaly_offsets_deg[0])
    M_t = M0 + (n + delta_n)*t
    r = orbital_to_eci(a, math.radians(incl_deg), math.radians(raan_deg), math.radians(argp_deg), M_t)
    drift_positions.append(r)

# --- Engine: estimate corrected mean motion from drift positions ---
# Use first and last position to estimate mean motion (simplified)
def eci_to_angle(r):
    x, y, z = r
    return math.atan2(y, x)

angle0 = eci_to_angle(drift_positions[0])
angleN = eci_to_angle(drift_positions[-1])
# unwrap
while angleN - angle0 < 0:
    angleN += 2*math.pi
n_est = (angleN - angle0) / (epochs[-1] - epochs[0])
mean_motion_rpd_est = n_est*86400/(2*math.pi)

corrected_tle = public_tle[0].copy()
corrected_tle["mean_motion_rev_per_day"] = round(mean_motion_rpd_est, 9)

# --- Terminal Output ---
print("=== Public TLE (simulated) ===")
for sat in public_tle:
    print(sat)

print("\n=== Drifted SAT1 Positions (ECI km) ===")
for t,r in zip(epochs, drift_positions):
    print(f"t={t}s: {r}")

print("\n=== Corrected TLE (SAT1) ===")
print(corrected_tle)
