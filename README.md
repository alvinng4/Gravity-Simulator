# Gravity Simulator
This is a 3D N-body gravity simulator projected on the xy plane.
* With the plotting module, you may setup your own system and plot the simulation results.
* With the interactive module, you may interact with the real-time gravity simulator.

Available integrators:
* Euler
* Euler_cromer
* Fourth Order Runge-Kutta
* Leapfrog
* Runge–Kutta–Fehlberg 4(5)
* Dormand–Prince method (DOPRI)
* Verner's method (DVERK)
* Runge–Kutta–Fehlberg 7(8)

![Screenshot](gravity_plot/examples/solar_system.png)

## Quick Start
### Python version
This program requires Python version 3.11. 

### Installation
Download the source file, or clone this repository by running the following command in terminal:
```
git clone https://github.com/alvinng4/Gravity-Simulator
```
Install the required packages by
```
pip install -r requirements.txt
```
## Plotting module
Once you have downloaded the source files, navigate to the directory of the source files in terminal and run
```
python gravity_plot
```
If you want to setup your own system, choose the "custom" option.
Note that the default unit is in one solar mass, AU and day.

The system data will be saved once all the required information has been entered.
If you wish to make any changes, you can access the file at gravity_sim/gravity_plot/customized_systems.csv and edit it accordingly.
## Interactive module
Once you have downloaded the source files, navigate to the directory of the source files in terminal and run
```
python gravity_sim
```
### Control
Move camera: `W` `A` `S` `D`\
Menu: `Esc`\
Pause: `P`\
Toggle full-screen mode: `F`\
Hide user interface: `H`\
Reset parameters: `R`\
Create new star: 
Hold the right mouse button to create a star + drag the mouse to give it an initial boost.\
Adjust parameter values: Left-click on the parameters panel to select a parameter + scroll to change its value.\
Switch integrators: Left-click the integrator on the integrators panel.

Warning: switching integrators in the middle of simulation may produce numerical error.
### Changing the resolution
The default resolution is set to the user's screen size. However, you can set your own resolution by the following command:
```
python3 gravity_sim -r <width> <height>
```
## Data References
1. Park, R.S., et al., 2021, “The JPL Planetary and Lunar Ephemerides DE440 and DE441”, https://ssd.jpl.nasa.gov/doc/Park.2021.AJ.DE440.pdf, *Astronomical Journal*, 161:105.
2. Horizons System, Jet Propulsion Laboratory, https://ssd.jpl.nasa.gov/horizons/

## Bibliography
1. Roa, Javier, et al. *Moving Planets Around: An Introduction to N-Body Simulations Applied to Exoplanetary Systems*, MIT Press, 2020
