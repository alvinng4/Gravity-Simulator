# grav_sim

Gravity simulation library written in C with Python API.

Target Audience: Physics students, or anyone who is interested in N-body simulations.

Note: This project is currently in active development.
There could be many incoming changes.

Features:
* Ten integrators including WHFast and IAS15
* Barnes-Hut algorithm for large $N$
* Particle-Mesh method for cosmological structure formation
* Multiple sample projects

<img src="./examples/media/solar_plus_3d.png" alt="Image" width="300">
<img src="./examples/media/rel_energy.png" alt="Image" width="300">

## Documentation
The documentations are currently outdated. I will update them as soon as possible.
Check `tutorial.ipynb`, `kirkword_gaps/`, `galaxy_collision` and `cosmic_structure` in the `examples/` folder for sample usage.

## Sample projects

Some projects are done with the API. The scripts are stored at the `examples` folder.

#### Evolving the solar system for 1 million years

<img src="./examples/media/eccentricity.png" alt="Image" width="300">
<img src="./examples/media/inclination.png" alt="Image" width="300">

#### Asteroid belt simulation

<img src="./examples/media/asteroid_belt.png" alt="Image" width="300">

#### Formation of Kirkwood gaps

<img src="./examples/media/Kirkwood_gap_semi_major_axes.png" alt="Image" width="300">
<img src="./examples/media/Kirkwood_gap_visualization.png" alt="Image" width="225">

Videos:
* https://www.youtube.com/watch?v=AEyjIF-8zT0
* https://www.youtube.com/watch?v=jHLLr7ACvDQ

#### Galaxy collision
Note: Initial condition is taken from Gadget-2 by Volker Springel.
Visualization is done with gadgetviewer.

<img src="./examples/media/galaxy_collision.png" alt="Image" width="400">

## Feedback and Bugs
If you found any bugs or want to leave some feedback, please feel free to let me know by opening an issue or sending an email to alvinng324@gmail.com.

## Acknowledgement
The following book helped me a lot in learning the basics of N-body simulations and also
the implementation of the integrators.
* J. Roa, et al. *Moving Planets Around: An Introduction to N-Body Simulations Applied to Exoplanetary Systems*, MIT Press, 2020