# 2D-N-body-gravity-simulator
Welcome to my CS50 Python project: 2D N-body gravity simulator. This project is currently in the initial stage of development. More functions are going to be implemented in the future.

## Quick Start
### Install required packages

Before running the program, make sure that you have the pygame and numpy packages installed. You can install them by running the following commands in terminal:
```
pip3 install pygame
pip3 install numpy
```
### Run the program
Once you have downloaded the source files, open a terminal window, navigate to the directory of the source files and run:
```
python3 gravity_sim
```
Please note that due to unknown reasons, the menu buttons may not work perfectly on the first execution of the program.

## Control
`W` `A` `S` `D`: move camera\
`Esc`: Open the menu 

## Changing the settings
When starting the program, you can customize the settings using command-line arguments.

### Resolution and full screen
The default resolution is 1366 x 768. However, you can set your own resolution by the following commands:
```
python3 gravity_sim -r <width> <height>
```
```
python3 gravity_sim --resolution <width> <height>
```
For example, if you want to set the resolution to 1920 x 1080, you may enter
```
python3 gravity_sim -r 1920 1080
```
If you want to set it to full screen, simply enter 0 for both width and height.
```
python3 gravity_sim -r 0 0
```

### Image scale
All the images of the gravitational objects are scaled by some ratio for better visibility. You can change the scale by the following commands:
```
python3 gravity_sim -i <sun_img_scale> <obj_img_scale>
```
```
python3 gravity_sim --img_scale <sun_img_scale> <obj_img_scale>
```
For example, if you want the sun and the objects to be scaled by 10x and 30x, you may enter
```
python3 gravity_sim -i 10 30
```
You may also change both the resolution and image scale together:
```
python3 gravity_sim -r 1920 1080 -i 10 30
```

## Functions in development:
1. Simulation function (Main focus)
2. `Left click`: Create new gravitational object
3. `Arrow up and down`: Change the rate of mass change when creating new object
4. `P`: Pause


## Bibliography
1. Roa, Javier, et al. Moving Planets Around: An Introduction to N-Body Simulations Applied to Exoplanetary Systems, MIT Press, 2020
2. Matthes, Eric. Python Crash Course, 3rd Edition: A Hands-On, Project-Based Introduction to Programming, No Starch Press, 2023