# Step 2: Gravity

Welcome to step 2. This is the most important step:
computing the gravitational acceleration.
Turns out this is also the most expensive part in N-body simulation, 
so we will spend some time on optimization.

## Newton's law of gravitation
I believe most of you are familiar with Newton's law of gravitation

$$
    \mathbf{F}_{ij} = m_{i} \mathbf{a}_{ij} = \frac{G m_i m_j}{r_{ij}^2} \hat{r}_{ij},
$$

where $\mathbf{F}_{ij}$ is the force on particle $i$ due to particle $j$, and
$\hat{r}_{ij}$ is the unit vector pointing from particle $i$ to particle $j$. That is,

$$
    \mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i, \quad
    \hat{r}_{ij} = \frac{\mathbf{r}_{ij}}{r_{ij}}.
$$

In practice, we are only interested in the acceleration. To compute the acceleration of
particle $i \in \{1, \dots, N\}$, we have

$$
    \mathbf{a}_{i} = \sum_{j \neq i} \frac{G m_j}{r_{ij}^3} \mathbf{r}_{ij},
$$

which can be easily done as this only involves simple vector operations.

## Implementation 1
Below shows our first naive implementation of the acceleration function.
The outer loop $i$ iterates over all particles, and the inner loop $j$ iterates
over all particles again to compute the acceleration between all pairs of particles.
However, this implementation is very slow.

```python linenums="1"
def acceleration_1(
    a: np.ndarray,
    system: common.System,
) -> None:
    """
    Compute the gravitational acceleration

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array to be modified in-place,
        with shape (N, 3)
    system : System
        System object.
    """
    # Empty acceleration array
    a.fill(0.0)

    # Declare variables
    num_particles = system.num_particles
    x = system.x
    m = system.m
    G = system.G

    # Calculations
    for i in range(num_particles):
        for j in range(num_particles):
            if i == j:
                continue

            R = x[j] - x[i]
            a[i] += G * m[j] * R / (np.linalg.norm(R) ** 3)
```

!!! tip "Where is the return statement?"
    Actually, there is no need to return the acceleration array `a` because
    we are modifying the memory in-place.

## Implementation 2

To optimize the code, we utilize the fact that the distance between particles $i$ and $j$ is the same as that between $j$ and $i$:

$$
    \mathbf{r}_{ij} = - \mathbf{r}_{ji} \implies r_{ij} = r_{ji}.
$$

!!! Note
    In our notation, the lowercase, non-bold $r_{ij}$ is the vector norm, which is always non-negative.

This allows us to cut the number of distance calculations in half.
(Calculating the distance is quite expensive as it involves the computation of `sqrt`.)
The outer loop $i$ still iterates over all particles,
but the inner loop $j$ only iterates from $i + 1$ to $N$. (Why? because all combinations of $i$ 
and $j \leq i$ has already been computed in the previous iterations. Therefore, we only need to
care about $j > i$.)

!!! tip
    If you want to rewrite this in C, this is the implementation you should use.

```python linenums="1"
def acceleration_2(
    a: np.ndarray,
    system: common.System,
) -> None:
    """
    Compute the gravitational acceleration

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array to be modified in-place,
        with shape (N, 3)
    system : System
        System object.
    """
    # Empty acceleration array
    a.fill(0.0)

    # Declare variables
    num_particles = system.num_particles
    x = system.x
    m = system.m
    G = system.G

    # Calculations
    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            R = x[j] - x[i]
            temp_value = G * R / (np.linalg.norm(R) ** 3)
            a[i] += temp_value * m[j]
            a[j] -= temp_value * m[i]
```

## Implementation 3* (Advanced)

The above implementation is still quite slow because we are using Python loops
to iterate over the particles. NumPy is (partly) implemented in C, which
makes it much faster than Python operations. If we can avoid 
Python loops completely, we can achieve a significant speedup.
This can be done by vectorizing the code. Note that this may be difficult
for beginners.

1. We first compute a displacement matrix $\mathbf{R}$, where $\mathbf{R}_{ij} = \mathbf{r}_{ij}$.
    It is a 3D array of shape $(N, N, 3)$, and the diagonal elements are all
    zero. This is computed by broadcasting: 
    
    * `x[np.newaxis, :, :]` has shape $(1, N, 3)$ and
    * `x[:, np.newaxis, :]` has shape $(N, 1, 3)$,
    
    so `r_ij = x[np.newaxis, :, :] - x[:, np.newaxis, :]` has entry $(i, j)$ equal to
    $\mathbf{r}_{ij}$.

    $$
        \mathbf{R} =
            \begin{bmatrix}
                \mathbf{r}_1 - \mathbf{r}_1 & \mathbf{r}_2 - \mathbf{r}_1 & \cdots & \mathbf{r}_N - \mathbf{r}_1 \\
                \mathbf{r}_1 - \mathbf{r}_2 & \mathbf{r}_2 - \mathbf{r}_2 & \cdots & \mathbf{r}_N - \mathbf{r}_2 \\
                \vdots & \vdots & \ddots & \vdots \\
                \mathbf{r}_1 - \mathbf{r}_N & \mathbf{r}_2 - \mathbf{r}_N & \cdots & \mathbf{r}_N - \mathbf{r}_N
            \end{bmatrix}
    $$

2. We compute the distance matrix $\mathbf{R}_{\mathrm{norm}}$, which is a 2D array of 
    shape $(N, N)$. It is computed by taking the norm along the last axis with length 3
    (`r_norm = np.linalg.norm(r_ij, axis=2)`).

3. We compute $1 / \mathbf{R}_\text{norm}^3$ (element-wise).
    Because the diagonal elements are all zero, division by zero will produce
    undefined values along the diagonal. Therefore, we only divide for positions where `r_norm` $> 0$.

    ```python
    inv_r_cubed = np.divide(
        1.0, r_norm * r_norm * r_norm, out=np.zeros_like(r_norm), where=r_norm > 0.0
    )
    ```

4. We compute the acceleration by 

$$
    \mathbf{a}_{i} = \sum_{j \neq i} \frac{G m_j}{r_{ij}^3} \mathbf{r}_{ij},
$$

The last step can again be done by using broadcasting. The resulting
acceleration array will be of shape $(N, 3)$.
```python
a[:] = G * np.sum(
    r_ij * inv_r_cubed[:, :, np.newaxis] * m[np.newaxis, :, np.newaxis], axis=1
)
```

This line of code is a bit complicated. Let us break it down:
 
* G is a constant, so it can be factored out of the summation.
* Ignore the last dimension with length 3 for now. We have
```python
a[:, 0] = G * np.sum(
    r_ij[:, :, 0] * inv_r_cubed[:, :] * m[np.newaxis, :], axis=1
)
```
We are summing over the axis 1 (column), so
the mass vector $\mathbf{m}$ needs to be broadcasted along the axis 0 (rows)
to $\mathbf{M}$ with the shape of $(N, N)$. We have

$$
    \mathbf{M} = 
    \begin{bmatrix}
        m_1 & m_2 & \cdots & m_N \\
        m_1 & m_2 & \cdots & m_N \\
        \vdots & \vdots & \ddots & \vdots \\
        m_1 & m_2 & \cdots & m_N
    \end{bmatrix}.
$$

The element-wise multiplication of $\mathbf{M}$ with $\mathbf{R}$ divided by $\mathbf{R}_\text{norm}^3$ gives

$$
    \begin{bmatrix}
        0                          & m_2 x_{12} / r_{12}^3 & \cdots & m_N x_{1N} / r_{1N}^3 \\
        m_1 x_{21} / r_{21}^3      & 0                     & \cdots & m_N x_{2N} / r_{2N}^3 \\
        \vdots                     & \vdots                & \ddots & \vdots \\
        m_1 x_{N1} / r_{N1}^3      & m_2 x_{N2} / r_{N2}^3 & \cdots & 0
    \end{bmatrix}
$$

where $x_{ij}$ denotes the $x$-component of $\mathbf{r}_{ij}$.
We are summing along the axis 1 (column). For particle $i$, we have

$$
    a_{i, x} = G \left[m_1 \frac{x_{i1}}{r_{i1}^3} + m_2 \frac{x_{i2}}{r_{i2}^3}
    + \cdots + 0 + \cdots + m_N \frac{x_{iN}}{r_{iN}^3} \right]
    = \sum_{j \neq i} \frac{G m_j}{r_{ij}^3} x_{ij}
$$

This is exactly what we want! Now, to also include the last dimension with length 3,
we simply add `np.newaxis` to the last axis. This gives us the vectorized version of the
full acceleration function. Later, we will see in the benchmark that this is much faster
than the previous implementations.

```python linenums="1"
def acceleration_3(
    a: np.ndarray,
    system: common.System,
) -> None:
    """
    Compute the gravitational acceleration

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array to be modified in-place,
        with shape (N, 3)
    system : System
        System object.
    """
    # Declare variables
    x = system.x
    m = system.m
    G = system.G

    # Compute the displacement vector
    r_ij = x[np.newaxis, :, :] - x[:, np.newaxis, :]

    # Compute the distance
    r_norm = np.linalg.norm(r_ij, axis=2)

    # Compute 1 / r^3
    inv_r_cubed = np.divide(
        1.0, r_norm * r_norm * r_norm, out=np.zeros_like(r_norm), where=r_norm > 0.0
    )

    # Compute the acceleration
    a[:] = G * np.sum(
        r_ij * inv_r_cubed[:, :, np.newaxis] * m[np.newaxis, :, np.newaxis], axis=1
    )
```

## Implementation 4* (Advanced)

In the last implementation, we are using `np.sum` and broadcasting to compute the
acceleration. In NumPy, there is a faster method `np.einsum`. Therefore, in this
implementation, we will replace the `np.sum` with `np.einsum`.

In our original implementation, notice how the broadcasting is done:
```python
a[:] = G * np.sum(
    r_ij * inv_r_cubed[:, :, np.newaxis] * m[np.newaxis, :, np.newaxis], axis=1
)
```

Denote the axis 0, 1, and 2 as $i$, $j$, and $k$ respectively.

* `r_ij` is a 3D array multiplied without broadcasting $\implies ijk$.
* `inv_r_cubed` is a 2D array multiplied with broadcasting along axis 2 $\implies ij$.
* `m` is a 1D array multiplied with broadcasting along axis 0 and 2 $\implies j$.

The final sum is done along axis 1 $\implies ijk \to ik$.

Using `np.einsum`, we can specify the indices to be summed over.
The following line of code is equivalent to the original implementation:   
```python
a[:] = G * np.einsum("ijk,ij,j->ik", r_ij, inv_r_cubed, m)
```

Full implementation:
```python linenums="1" hl_lines="33"
def acceleration_4(
    a: np.ndarray,
    system: common.System,
) -> None:
    """
    Compute the gravitational acceleration

    Parameters
    ----------
    a : np.ndarray
        Gravitational accelerations array to be modified in-place,
        with shape (N, 3)
    system : System
        System object.
    """
    # Declare variables
    x = system.x
    m = system.m
    G = system.G

    # Compute the displacement vector
    r_ij = x[np.newaxis, :, :] - x[:, np.newaxis, :]

    # Compute the distance
    r_norm = np.linalg.norm(r_ij, axis=2)

    # Compute 1 / r^3
    inv_r_cubed = np.divide(
        1.0, r_norm * r_norm * r_norm, out=np.zeros_like(r_norm), where=r_norm > 0.0
    )

    # Compute the acceleration
    a[:] = G * np.einsum("ijk,ij,j->ik", r_ij, inv_r_cubed, m)
```

## Benchmark

To benchmark the performance, we will use the `timeit` module and 
repeat each function 10000 times. We will take the mean with standard
error = $\sigma / \sqrt{N_\text{repeats}}$.

```python title="step2.py"
import math
import timeit

import numpy as np

import common

INITIAL_CONDITION = "solar_system"
NUM_REPEATS = 10000


def main() -> None:
    # Get initial conditions
    system, _, _, _ = common.get_initial_conditions(INITIAL_CONDITION)

    ### Benchmark ###
    print("Benchmarking with 10000 repetitions")
    print()

    # Allocate memory
    a = np.zeros((system.num_particles, 3))

    # Acceleration 1
    run_time_1 = np.zeros(NUM_REPEATS)
    for i in range(NUM_REPEATS):
        start = timeit.default_timer()
        acceleration_1(a, system)
        end = timeit.default_timer()
        run_time_1[i] = end - start
    print(
        f"acceleration_1: {run_time_1.mean():.6f} +- {run_time_1.std(ddof=1) / math.sqrt(NUM_REPEATS):.3g} seconds"
    )

    ...  # (Repeat for acceleration_2, 3, and 4)
```

Finally, we do an error check by comparing the results with those from the first naive implementation.
```python
def main() -> None:
    ...
    # Check for relative errors
    ### Error check ###
    acceleration_1(a, system)
    a_1 = a.copy()
    acceleration_2(a, system)
    a_2 = a.copy()
    acceleration_3(a, system)
    a_3 = a.copy()
    acceleration_4(a, system)
    a_4 = a.copy()

    rel_error_2 = np.sum(np.abs(a_1 - a_2)) / np.sum(np.abs(a_1))
    rel_error_3 = np.sum(np.abs(a_1 - a_3)) / np.sum(np.abs(a_1))
    rel_error_4 = np.sum(np.abs(a_1 - a_4)) / np.sum(np.abs(a_1))

    print()
    print("Error check: (relative difference from acceleration_1)")
    print(f"acceleration_2: {rel_error_2:.3g}")
    print(f"acceleration_3: {rel_error_3:.3g}")
    print(f"acceleration_4: {rel_error_4:.3g}")
```

The results are as follows:
```
Benchmarking with 10000 repetitions

acceleration_1: 0.000203 +- 8.08e-08 seconds
acceleration_2: 0.000164 +- 1.25e-06 seconds
acceleration_3: 0.000013 +- 2.21e-08 seconds
acceleration_4: 0.000012 +- 1.38e-08 seconds

Error check: (relative difference from acceleration_1)
acceleration_2: 0
acceleration_3: 1.9e-16
acceleration_4: 1.9e-16
```
The vectorized implementation is about 10 - 20 times faster than the naive implementation!
As for the error check, the small relative difference is likely due to rounding errors,
which could be ignored. (For 64-bit floating point numbers, the machine epsilon is about $10^{-16}$.)
By the way, since `acceleration_4` is the fastest, we put it into `common.py`.

!!! Tip "Performance in C"
    If you are interested in the performance in C,
    below is a benchmark using our `grav_sim` package:
    ```
    Test 0:    Method: Pairwise
        Number of times: 10000000
        Avg time: 2.06e-07 (+- 1.30e-10) s
    ```
    This is about 58 times faster than the vectorized NumPy implementation. But beware that
    this may not be totally accurate as the run time for each run is too short.

## Full scripts
The full scripts are available at `5_steps_to_n_body_simulation/python/`,
or https://github.com/alvinng4/grav_sim/blob/main/5_steps_to_n_body_simulation/python/

??? note "step2.py (Click to expand)"
    ```python linenums="1" title="5_steps_to_n_body_simulation/python/step2.py"
    --8<-- "5_steps_to_n_body_simulation/python/step2.py"
    ```

??? note "common.py (Click to expand)"
    ```python linenums="1" title="5_steps_to_n_body_simulation/python/common.py"
    --8<-- "5_steps_to_n_body_simulation/python/common.py"
    ```
