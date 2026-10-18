# Step 4: Higher-order algorithms

In the last step, we implemented a simple Euler method to integrate
the solar system for 200 years. The poor accuracy is expected, because
we used a simple first order algorithm with global error $\mathcal{O}(\Delta t)$.
In this step, we will implement 3 new algorithms to improve our simulation.

## Relative energy error

Before implementing new algorithms, we need a way to measure the accuracy
of the simulation. By conservation of energy, we can assume that better
algorithms will have a better conservation on the total energy of the system.
Computing the total energy is trivial. All we need is the solution array.

$$
    \text{Total Energy} 
    = \overset{\text{KE}}{\overbrace{\sum_{i=1}^{N} \frac{1}{2} m_i v_i^2}} 
    - \overset{\text{PE}}{\overbrace{\sum_{i=1}^{N} \sum_{j = i + 1}^{N} \frac{G m_i m_j}{r_{ij}}}}.
$$

Then, we can compute the relative energy error as

$$
    \text{Relative Energy Error} = \frac{|\text{Energy} - \text{Initial Energy}|}{\text{Initial Energy}}.
$$

This gives the following code:

```python title="common.py"
def compute_rel_energy_error(
    sol_x: np.ndarray, sol_v: np.ndarray, system: System
) -> np.ndarray:
    """
    Compute the relative energy error of the simulation.

    Parameters
    ----------
    sol_x : np.ndarray
        Solution position array with shape (N_steps, num_particles, 3).
    sol_v : np.ndarray
        Solution velocity array with shape (N_steps, num_particles, 3).
    system : System
        System object.

    Returns
    -------
    energy_error : np.ndarray
        Relative energy error of the simulation, with shape (N_steps,).
    """
    # Allocate memory and initialize arrays
    n_steps = sol_x.shape[0]
    num_particles = system.num_particles
    m = system.m
    G = system.G
    rel_energy_error = np.zeros(n_steps)

    # Compute the total energy (KE + PE)
    for count in range(n_steps):
        x = sol_x[count]
        v = sol_v[count]
        for i in range(num_particles):
            # KE
            rel_energy_error[count] += 0.5 * m[i] * np.linalg.norm(v[i]) ** 2
            # PE
            for j in range(i + 1, num_particles):
                rel_energy_error[count] -= G * m[i] * m[j] / np.linalg.norm(x[i] - x[j])

    # Compute the relative energy error
    initial_energy = rel_energy_error[0]
    rel_energy_error = (rel_energy_error - initial_energy) / initial_energy
    rel_energy_error = np.abs(rel_energy_error)

    return rel_energy_error
```

!!! tip
    If you found it too slow, you are encouraged to vectorize it, just like
    what we did in step 2

Now, we can plot the relative energy error with the following code, where
y-axis is in log scale.

```python title="common.py"
def plot_rel_energy_error(rel_energy_error: np.ndarray, sol_t: np.ndarray) -> None:
    """
    Plot the relative energy error.

    Parameters
    ----------
    rel_energy_error : np.ndarray
        Relative energy error of the simulation, with shape (N_steps,).
    sol_t : np.ndarray
        Solution time array with shape (N_steps,).
    """
    plt.figure()
    plt.plot(sol_t, rel_energy_error)
    plt.yscale("log")
    plt.xlabel("Time step")
    plt.ylabel("Relative Energy Error")
    plt.title("Relative Energy Error vs Time Step")
    plt.show()
```

We copy `step3.py` to `step4.py` and add the following code to the
end of the `main` function:

```python title="step4.py"
def main() -> None:
    ...
    # Compute and plot relative energy error
    rel_energy_error = common.compute_rel_energy_error(sol_x, sol_v, system)
    print(f"Relative energy error: {rel_energy_error[-1]:.3g}")
    common.plot_rel_energy_error(rel_energy_error, sol_t / 365.24)
```

## Euler method

Let us run the simulation again. We have the following plots:

![Trajectory](../figures/step3/trajectory.png)
![Relative Energy Error](../figures/step4/euler_rel_energy_error.png)

The final relative energy error $\sim 10^{-1}$, which is not very good.

## Euler-Cromer method

Now, we begin implementing our first new algorithm.
The Euler-Cromer method, also known as semi-implicit Euler method, is a simple
modification of the Euler method,

$$
    \mathbf{v}_{n+1} = \mathbf{v}_n + \mathbf{a}(\mathbf{r}_n) \Delta t,
$$

$$
    \mathbf{r}_{n+1} = \mathbf{r}_n + \mathbf{v}_{n+1} \Delta t.
$$

Notice that the position update is done using the updated velocity $\mathbf{v}_{n+1}$ instead
of $\mathbf{v}_n$. Therefore, it is an first order semi-implicit method. Although the global
error is still $\mathcal{O}(\Delta t)$, it is a *symplectic* method, which implies that
the energy error over time is bounded. This is a very useful property for long time-scale
simulations. We have the following code:

```python title="common.py"
def euler_cromer(a: np.ndarray, system: System, dt: float) -> None:
    """
    Advance one step with the Euler-Cromer method.

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array with shape (N, 3).
    system : System
        System object.
    dt : float
        Time step.
    """
    acceleration(a, system)
    system.v += a * dt
    system.x += system.v * dt
```

Running the simulation again, we have the following plots:

![Trajectory](../figures/step4/euler_cromer_trajectory.png)
![Relative Energy Error](../figures/step4/euler_cromer_rel_energy_error.png)

The trajectory looks a lot better than the Euler method! Also, the relative energy error is
bounded at $\sim 10^{-4}$. This provides long term stability for the simulation.

!!! Note "Energy conservation $\iff$ Higher accuracy?"
    Although it seems like energy conservation implies higher accuracy,
    this is not neccessarily true. Recall that the global error for Euler-Cromer method
    is still $\mathcal{O}(\Delta t)$, even though its energy error is bounded over time. 
    Therefore, we can only state the opposite direction: Higher accuracy $\implies$ Energy conservation.

## Runge-Kutta method

The Runge-Kutta method is a family of ODE solvers. First, let us look at a 
second order variation called the midpoint method. 

* In the Euler method, we are using
the slope evaluated at the beginning of the interval to update the position and velocity.

* In the Euler-Cromer method, for the
position update, we are using the slope evaluated at the end instead.

What about using the slope evaluated at the center of the interval?
We can approximate the position and velocity at the midpoint using one Euler step,
then perform the updates using the midpoint values. This gives us the following updates:

$$
    x_{n+1} = x_n + f\left(t_n + \frac{1}{2} \Delta t, x + \frac{1}{2} f(t_n) \Delta t \right) \Delta t.
$$

Using a more general notation, we can write the midpoint method as

$$
\begin{aligned}
    k_1 &= f(t_n, x_n), \\
    k_2 &= f\left(t_n + \frac{1}{2} \Delta t, x_n + \frac{1}{2} k_1 \Delta t \right), \\
    x_{n+1} &= x_n + k_2 \Delta t + \mathcal{O}(\Delta t^3).
\end{aligned}
$$

A simple taylor expansion will shows that the local truncation error is
indeed $\mathcal{O}(\Delta t^3)$, which means that it is a second order method.
A more commonly used Runge-Kutta method is the 4th order Runge-Kutta method (RK4),
which provides a good balance between accuracy and computational cost. It is given by

$$
\begin{aligned}
    k_1 &= f(t_n, x_n), \\
    k_2 &= f\left(t_n + \frac{1}{2} \Delta t, x_n + \frac{1}{2} k_1 \Delta t \right), \\
    k_3 &= f\left(t_n + \frac{1}{2} \Delta t, x_n + \frac{1}{2} k_2 \Delta t \right), \\
    k_4 &= f(t_n + \Delta t, x_n + k_3 \Delta t), \\
    x_{n+1} &= x_n + \frac{1}{6} (k_1 + 2 k_2 + 2 k_3 + k_4) \Delta t + \mathcal{O}(\Delta t^5).
\end{aligned}
$$


In our code, we can write out the computation of each term explicitly. However, I prefer a more general
approach by defining the `coeff` and `weights` arrays instead. The `coeff` array is given as

$$
    \text{coeff} = \begin{bmatrix}
        1/2, 1/2, 1
    \end{bmatrix}
$$

for the computation of $k_2$, $k_3$ and $k_4$. The `weights` array is given as

$$
    \text{weights} = \begin{bmatrix}
        1/6, 1/3, 1/3, 1/6
    \end{bmatrix}
$$

for the final update. The code is given as follows:

```python title="common.py"
def rk4(a: np.ndarray, system: System, dt: float) -> None:
    """
    Advance one step with the RK4 method.

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array with shape (N, 3).
    system : System
        System object.
    dt : float
        Time step.
    """
    num_stages = 4
    coeff = np.array([0.5, 0.5, 1.0])
    weights = np.array([1.0, 2.0, 2.0, 1.0]) / 6.0

    # Allocate memory and initialize arrays
    x0 = system.x.copy()
    v0 = system.v.copy()
    xk = np.zeros((num_stages, system.num_particles, 3))
    vk = np.zeros((num_stages, system.num_particles, 3))

    # Initial stage
    acceleration(a, system)
    xk[0] = v0
    vk[0] = a

    # Compute the stages
    for stage in range(1, num_stages):
        # Compute acceleration
        system.x = x0 + dt * coeff[stage - 1] * xk[stage - 1]
        acceleration(a, system)

        # Compute xk and vk
        xk[stage] = v0 + dt * coeff[stage - 1] * vk[stage - 1]
        vk[stage] = a

    # Advance step
    dx = 0.0
    dv = 0.0
    for stage in range(num_stages):
        dx += weights[stage] * xk[stage]
        dv += weights[stage] * vk[stage]

    system.x = x0 + dt * dx
    system.v = v0 + dt * dv
```

!!! Tip
    The final loop can be vectorized to improve performance:
    ```python
    dx = np.einsum("i,ijk->jk", weights, xk)
    dv = np.einsum("i,ijk->jk", weights, vk)
    ```

Let's run the simulation again. We only show the relative energy error 
plot:

![Relative Energy Error](../figures/step4/rk4_rel_energy_error.png)

The final error is in the order of $10^{-6}$, which is quite nice. However,
the error is growing over time, so this may not be a good choice for long term simulations.

## Leapfrog method

Our final algorithm is the leapfrog method, which is a second-order symplectic method.
Similar to the Euler-Cromer method, it conserves energy over time. We will implement the 
Kick-Drift-Kick (KDK) variant, which is given by a velocity kick for half time step,

$$
    \mathbf{v}_{n+1/2} = \mathbf{v}_n + \frac{1}{2} \mathbf{a}(\mathbf{r}_n) \Delta t,
$$

a position drift for a full time step,

$$
    \mathbf{r}_{n+1} = \mathbf{r}_n + \mathbf{v}_{n+1/2} \Delta t,
$$

and a final velocity kick for half time step,

$$
    \mathbf{v}_{n+1} = \mathbf{v}_{n+1/2} + \frac{1}{2} \mathbf{a}(\mathbf{r}_{n+1}) \Delta t.
$$

!!! Tip
    For optimization, you can combine the final velocity kick from the last time step with the first velocity kick.
    However, you need to be careful because now the velocity and position is not synchronized.
    There is also a synchronized version of the leapfrog method called the velocity Verlet method.

The implementation is simple:

```python
def leapfrog(a: np.ndarray, system: System, dt: float) -> None:
    """
    Advance one step with the LeapFrog method.

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array with shape (N, 3).
    system : System
        System object.
    dt : float
        Time step.
    """
    # Velocity kick (v_1/2)
    acceleration(a, system)
    system.v += a * 0.5 * dt

    # Position drift (x_1)
    system.x += system.v * dt

    # Velocity kick (v_1)
    acceleration(a, system)
    system.v += a * 0.5 * dt
```

Running the simulation again, we have the relative energy error plot:

![Relative Energy Error](../figures/step4/leapfrog_rel_energy_error.png)

The relative energy error is bounded at $\sim 10^{-6}$, which is better than the Euler-Cromer method!

## Which algorithm should I choose?

To choose between RK4 and LeapFrog, below are some
factors that you may consider:

* **Time scale**: If you are simulating a long time scale, LeapFrog should be better as it conserves energy.
* **Accuracy**: If accuracy matters to you, then RK4 is a better choice since it is higher order.
* **Computational cost**:
    RK4 is very expensive as it requires 4 acceleration evaluations per step.
    In contrast, LeapFrog only requires 1 acceleration evaluation although it is second order.
* **Experiment**:
    If you are not sure which one to choose, you can always experiment and compare the results.
    Below is a relative energy error plot I made for the solar system simulation using both RK4 
    and LeapFrog. I chose two very small dt while ensuring that the runtime is the same for both methods.
    In terms of the results, seems like RK4 is better in this case, but note that I used some techniques
    to remove the rounding error. (See [Reducing round off error](../docs/documentations/reducing_round_off_error.md))

    <img src="../../examples/media/rel_energy.png" alt="Relative energy error comparison" width="400"/>

## Summary

In this step, we have implemented 3 new algorithms: Euler-Cromer, RK4 and Leapfrog.
RK4 and Leapfrog are both very popular algorithms for N-body simulations.
All algorithms we implemented so far are fixed time step methods, which may not be very flexible
for chaotic systems with close encounters. It is also a headache to tune the time step.
In the next step, we will implement an adaptive time-stepping method which uses a tolerance
parameter to control the time step instead.

## Full scripts
The full scripts are available at `5_steps_to_n_body_simulation/python/`,
or https://github.com/alvinng4/grav_sim/blob/main/5_steps_to_n_body_simulation/python/

??? note "step4.py (Click to expand)"
    ```python linenums="1" title="5_steps_to_n_body_simulation/python/step4.py"
    --8<-- "5_steps_to_n_body_simulation/python/step4.py"
    ```

??? note "common.py (Click to expand)"
    ```python linenums="1" title="5_steps_to_n_body_simulation/python/common.py"
    --8<-- "5_steps_to_n_body_simulation/python/common.py"
    ```
