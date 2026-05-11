# AlGaAs barriers: electron Monte Carlo simulation
This Julia source code was developed for the ensemble Monte Carlo (MC) Boltzmann transport simulation of electrons in an AlGaAs solid-state thermionic device. Unique features include a three-valley model with ellipsoidal, anisotropic electron valleys and detailed tracking of phonon features. Results generated from this simulation have been included in the following article:

Franceschetti, L., Shin, S., and Kaviany, M., “Phonon valleytronics: Enhanced phonon absoprtion by electron valley-energy filtering,” under review.

## Structure
The simulation is initialized from MonteCarlo.jl, which calls the primary MC temporal loop (executes one timestep). Within that timestep, every particle is drifted and potentially scattered. At the end of the timestep, particles are injected from the contacts when nescessary. The Poisson equation is then evaluated, and data is collected. This continues until the end of the simulation with a specified duration. The detailed algorthim, which follows a rather standard MC for ensemble electron transport, is illustrated in the flowchart below.
<img src="https://github.com/lorfrances/algaas-barriers-MC/blob/main/images/flowchart.png" width="600">

## Usage
A sample input file (INPUT.txt) has been included with a few customizable parameters. The program is run by executing _MonteCarlo.jl_ followed by a nescessary input file from the command line (not the Julia REPL).
