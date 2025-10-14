# Satellite Calibration and TLE Refinement Simulation

This project is a **single-file Python simulation** demonstrating how a ground-based calibration engine can identify and refine satellite orbital data (TLEs) using simulated observations and nominal TLEs.

---

## Overview

The simulation has **three main components**:

1. **Public TLE Generator**  
   - Simulates **15 nominal satellite TLEs** representing predicted orbits.  
   - Acts as reference data for the calibration engine.

2. **Satellite Pass Simulation (SDR Input)**  
   - Simulates **15 satellite passes** as if observed by a ground SDR receiver.  
   - Each pass incorporates **drift and observational errors** to mimic real-world deviations from the public TLEs.

3. **Predictive and Calibration Engine**  
   - Ground-based engine initialized with **accurate time and location**.  
   - Observes the simulated satellite passes.  
   - Searches through the **public TLEs** using a **position and time tolerance** to find the best match.  
   - Once a satellite is correctly mapped, the engine **refines the TLE** using the observed trajectory.

---

## Working Principle

1. **Public TLEs** are generated first — providing nominal satellite orbit predictions.  
2. **Satellite passes** are then simulated, including drift and observational errors.  
3. The **engine** observes each simulated pass.  
4. It **searches the public TLEs** to find the one matching the observed satellite within a specified tolerance.  
5. Once a match is found, the engine **maps the observed pass to the TLE**.  
6. The TLE is then **updated/refined** based on the observed trajectory, correcting for drift and errors.

---

## How to Run

### 1. Install Dependencies
```bash
pip install numpy pandas matplotlib sgp4
