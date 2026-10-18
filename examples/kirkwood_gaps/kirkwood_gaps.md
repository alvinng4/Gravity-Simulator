# Formation of Kirkwood gaps

??? Note "Source code (Python) (Click to expand)"
    ```python linenums="1" title="examples/kirkwood_gaps/python/kirkwood_gaps.py"
    --8<-- "examples/kirkwood_gaps/python/kirkwood_gaps.py"
    ```
??? Note "Source code (C) (Click to expand)"
    ```c linenums="1" title="examples/kirkwood_gaps/c/kirkwood_gaps.c"
    --8<-- "examples/kirkwood_gaps/c/kirkwood_gaps.c"
    ```

??? Note "Source code (Animation) (Click to expand)"
    ```python linenums="1" title="examples/kirkwood_gaps/animate_kirkwood_gaps.py"
    --8<-- "examples/kirkwood_gaps/animate_kirkwood_gaps.py"
    ```

Due to gravitational resonance of Jupiter and other planets, gaps are formed in the asteroid belt
at certain ratio of the semi-major axis $a$. These gaps are called Kirkwood gaps.
To simulate its formation, we use the same initial conditions as the
[Asteroid belt animation](../asteroid_belt_animation/asteroid_belt.md) example, where
the initial value of $a$ are sampled with $a \sim \text{Uniform}(2.0, 3.5)$.

For the simulation, we chose the WHFast integrator with $\Delta t = 180$ days and
$t_f = 5$ million years. The large time step is chosen after a convergence test done without
the asteroids, as the secular evolutions seems to be fairly accurate even with this time step.
Note that this time step is only possible for the WHFast integrator, but not other integrators like
RK4 or LeapFrog.

Although the time scales are quite long, we assume the asteroids to be massless for performance reasons.
So, just keep in mind that this is not a very accurate simulation. Furthermore, even if we tries to include
the gravity due to the asteroids, their contribution would likely be smaller than the round off error,
so we can't really include them even if we wanted to.

## Results
The results seems to be quite successful. The simulation is done on Macbook Air M1 and it only took 24 hours.
About 7500 asteroids are removed from the simulation due to failure to converge in the Kepler's solver.

<iframe width="560" height="315" src="https://www.youtube.com/embed/AEyjIF-8zT0?si=M5G8cS0i3D71PGLg" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>

<iframe width="560" height="315" src="https://www.youtube.com/embed/jHLLr7ACvDQ?si=O_YzbUGx_grNV2ct" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>


## Running the example (C)
To build the example, run
```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_HDF5=ON -DUSE_OPENMP=ON --fresh
cmake --build build
```
Then, you can run it by
```
cd build
./kirkwood_gaps
```

## Running the example (Python)
Simply run
```
python kirkwood_gaps.py
```
