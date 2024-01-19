# N-body-gravity-simulator
Welcome to my Python project: N-body gravity simulator. This project is currently in development. More functions are going to be implemented in the future.

## Progress
* Fixed step-size Euler, Euler-Cromer, RK2, RK4 and leapfrog methods for N-body are implemented.
* 4 modes are available: Void, Solar System, Figure-8 Orbit and Pythagorean Three-Body Orbit.
* RKF4(5) in development

## Quick Start
### Install required packages
Before running the program, make sure that you have the following packages installed:
1. pygame 
2. numpy
3. scipy
4. numba

You can install them by running the following commands in terminal:
```
pip3 install pygame
pip3 install numpy
pip3 install scipy
pip3 install numba
```
If you do not want to install scipy and numba, you can still run the program by removing the following lines in `gravity_sim/simulator.py`.
Warning: This may leads to a significant drop in performance.
```
import numba as nb
```
```
@nb.njit
```

### Running the program
Once you have downloaded the source files, open a terminal window, navigate to the directory of the source files and run:
```
python3 gravity_sim
```
Note: the buttons may need some time to react.

## Control
Move camera: `W` `A` `S` `D`\
Menu: `Esc`\
Pause: `P`\
Toggle full-screen mode: `F`\
Create new star: 
Hold the right mouse button to create a star + drag the mouse to give it an initial boost.\
Adjust parameter values: Left-click on the parameters panel to select a parameter + scroll to change its value.\
Switch integrators: Left-click the integrator on the integrators panel.

## Changing the resolution
The default resolution is set to 1920 x 1080. However, you can set your own resolution by the following command:
```
python3 gravity_sim -r <width> <height>
```

## Functions in development:
1. More advanced simulation: RKF4(5) method.

## Bibliography
1. Roa, Javier, et al. *Moving Planets Around: An Introduction to N-Body Simulations Applied to Exoplanetary Systems*, MIT Press, 2020
2. Matthes, Eric. *Python Crash Course, 3rd Edition: A Hands-On, Project-Based Introduction to Programming*, No Starch Press, 2023