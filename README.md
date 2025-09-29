# LEO Satellite Drift Simulation & TLE Correction Engine

This README explains the **mathematical approach** behind simulating a LEO constellation, introducing a drift in a satellite’s orbit, and estimating a corrected TLE from observed signals.

---

## 1. Circular Orbit Basics

For a satellite in a **circular orbit**, the radius from the Earth's center is:

\[
r = R_e + h
\]

Where:  
- \(R_e\) = Earth radius (~6378.137 km)  
- \(h\) = satellite altitude above Earth (~550 km for LEO)  

The **mean motion** \(n\) (rad/s) is:

\[
n = \sqrt{\frac{\mu}{r^3}}
\]

Where:  
- \(\mu = 398600.4418\ \text{km}^3/\text{s}^2\) (Earth’s gravitational parameter)

The **orbital period** \(T\) is:

\[
T = \frac{2 \pi}{n}
\]

---

## 2. Generating “Public TLE” (Fake)

Each satellite is initialized with:  
- Semi-major axis \(a = r\)  
- Inclination \(i\)  
- RAAN \(\Omega\)  
- Argument of perigee \(\omega\)  
- Mean anomaly \(M_0\)  

For multiple satellites along the same orbital plane, we offset the **mean anomaly** by a small angle:

\[
M_0^{(i)} = M_0 + \Delta M_i
\]

Where \(\Delta M_i\) is the along-track spacing in degrees.

The **position in the perifocal frame** is:

\[
\mathbf{r}_{pf} = 
\begin{bmatrix} 
a \cos \nu \\ 
a \sin \nu \\ 
0 
\end{bmatrix}
\]

Where \(\nu\) is the **true anomaly**, equal to the mean anomaly \(M\) for circular orbits:

\[
\nu(t) = M_0 + n \cdot t
\]

Transform to **ECI frame** via rotation matrices:

\[
\mathbf{r}_{eci} = R_z(\Omega) \cdot R_x(i) \cdot R_z(\omega) \cdot \mathbf{r}_{pf}
\]

Where:  
- \(R_x\) and \(R_z\) are standard rotation matrices about X and Z axes.  

---

## 3. Simulating Satellite Drift

To simulate a real-world scenario where the SDR observes a slight orbital drift, we modify the **mean motion** of a satellite:

\[
n_\text{drift} = n + \Delta n
\]

Where \(\Delta n\) represents a **small drift**, e.g., from orbital decay or perturbations.  

Then the **drifted mean anomaly** becomes:

\[
M_\text{drift}(t) = M_0 + n_\text{drift} \cdot t
\]

The satellite’s ECI position is then recalculated as:

\[
\mathbf{r}_{eci, drift}(t) = R_z(\Omega) \cdot R_x(i) \cdot R_z(\omega) \cdot
\begin{bmatrix} a \cos M_\text{drift}(t) \\ a \sin M_\text{drift}(t) \\ 0 \end{bmatrix}
\]

---

## 4. Estimating Corrected Mean Motion

Given a set of observed positions from SDR (\(\mathbf{r}_{eci, drift}\)) and **predicted public TLE**, we can estimate the corrected mean motion (\(n_\text{est}\)).

1. Compute the **orbital angle** in the perifocal plane for first and last observation:

\[
\theta_0 = \arctan2(y_0, x_0), \quad
\theta_N = \arctan2(y_N, x_N)
\]

2. Unwrap the angle difference:

\[
\Delta \theta = \theta_N - \theta_0
\]

3. Compute **corrected mean motion**:

\[
n_\text{est} = \frac{\Delta \theta}{t_N - t_0}
\]

Convert to **revolutions per day** (TLE convention):

\[
\text{mean motion (rev/day)} = n_\text{est} \cdot \frac{86400}{2\pi}
\]

---

## 5. Corrected TLE

The updated TLE parameters for the drifted satellite are:

| Field | Value |
|-------|-------|
| Semi-major axis \(a\) | Same as public TLE |
| Inclination \(i\) | Same as public TLE |
| RAAN \(\Omega\) | Same as public TLE |
| Argument of perigee \(\omega\) | Same as public TLE |
| Mean anomaly \(M_0\) | Same initial value |
| Mean motion | \(n_\text{est}\) (from drift fitting) |

---

## 6. Summary of Steps

1. Generate **fake public TLE** with initial orbital elements.  
2. Sample positions at discrete times (ECI).  
3. Simulate a **drift** in mean motion for a satellite.  
4. Observe satellite positions via “SDR” (simulated drifted positions).  
5. Estimate corrected mean motion using first/last position angular difference.  
6. Update TLE with **corrected mean motion**, forming **new TLE**.  

This forms the **proof-of-concept pipeline**: feeding public TLE + observed drift → identify satellite → correct TLE → propagate positions.

---

## 7. Notes

- The simulation assumes **circular orbits**; eccentricity can be added for higher fidelity.  
- Drift estimation is **simplified** (linear fit of angle vs time). More accurate methods can use **least-squares on all sampled positions**.  
- This method is exactly what the **Python engine provided** in this project implements.
