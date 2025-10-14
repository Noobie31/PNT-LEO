Satellite Calibration and TLE Refinement Simulation

This project is a single-file Python simulation demonstrating how a ground-based calibration engine can identify and refine satellite orbital data (TLEs) using simulated observations and nominal TLEs.

Overview

The simulation has three main components:

Public TLE Generator

Simulates 15 nominal satellite TLEs representing predicted orbits.

Acts as reference data for the calibration engine.

Satellite Pass Simulation (SDR Input)

Simulates 15 satellite passes as if observed by a ground SDR receiver.

Each pass incorporates drift and observational errors to mimic real-world deviations from the public TLEs.

Predictive and Calibration Engine

Ground-based engine initialized with accurate time and location.

Observes the simulated satellite passes.

Searches through the public TLEs using a position and time tolerance to find the best match.

Once a satellite is correctly mapped, the engine refines the TLE using the observed trajectory.

Working Principle

Public TLEs are generated first — providing nominal satellite orbit predictions.

Satellite passes are then simulated, including drift and observational errors.

The engine observes each simulated pass.

It searches the public TLEs to find the one matching the observed satellite within a specified tolerance.

Once a match is found, the engine maps the observed pass to the TLE.

The TLE is then updated/refined based on the observed trajectory, correcting for drift and errors.

How to Run
1. Install Dependencies
pip install numpy pandas matplotlib sgp4


If you encounter permission issues:

pip install --user numpy pandas matplotlib sgp4

2. (Recommended) Use a Virtual Environment
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install numpy pandas matplotlib sgp4

3. Run the Simulation
python .\sim.py

Notes

Entire simulation is implemented in a single Python file (sim.py).

Public TLEs and satellite passes are fully simulated; no real-world data is required.

The engine acts like a cognitive calibration node, capable of identifying, mapping, and refining TLEs based on observed trajectories.

Demonstrates how a network of ground stations could cooperatively improve public TLE accuracy in real time.
