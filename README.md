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
`Esc`: Open the menu\
`F`: Full screen mode

## Changing the settings
When starting the program, you can customize the settings using command-line arguments.
```
python3 gravity_sim [-r | --resolution] [-i | --img_scale] [-d | --distance_scale] <args> 
```

### Resolution
The default resolution is set to 1920 x 1080. However, you can set your own resolution by the following commands:
```
python3 gravity_sim -r <width> <height>
```
### Image scale
All the images of the gravitational objects are scaled by some factors for better visibility. You can change the scale by the following commands:
```
python3 gravity_sim -i <sun_img_scale> <obj_img_scale>
```
### Distance scale
In some cases, celestial bodies may be too far out of the screen to be seen. You can address this by changing the distance scale. Simply run:
```
python3 gravity_sim -d <distance_scale>
```
You may change the parameters altogether. The following command can be used to set a resolution of 2560 x 1440, sun image scale of 10, object image scale of 30 and distance scale of 0.5:
```
python3 gravity_sim -r 2560 1440 -i 10 30 -d 0.5
```


## Functions in development:
1. Simulation function (Main focus)
2. `Left click`: Create new gravitational object
3. `Arrow up and down`: Change the rate of mass change when creating new object
4. `P`: Pause


## Bibliography
1. Roa, Javier, et al. *Moving Planets Around: An Introduction to N-Body Simulations Applied to Exoplanetary Systems*, MIT Press, 2020
2. Matthes, Eric. *Python Crash Course, 3rd Edition: A Hands-On, Project-Based Introduction to Programming*, No Starch Press, 2023