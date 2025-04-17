# Reducing rounding error with compensated summation

A method known as compensated summation [@CompensatedSummation] is implemented for
all integrators in grav_sim except WHFast:

When we advance our system by $\Delta t$, we have 

$x_{n+1} = x_n + \delta x$

Since $\delta x$ is very small compared to $x_n$, many digits of precision will be lost.
By compensated summation, we keep track of the losing digits using another variable, which
allows us to effectively eliminates round off error with very little cost.

The algorithm is as follows:

1. Calculate $\delta x$
2. $x_0 \leftarrow x$
3. $e \leftarrow  e + \delta x$
4. $x \leftarrow  x_0 + e$
5. $e \leftarrow  e + (x_0 - x)$

!!! example

    Consider 

    $x_n$      = <span style="color:green">1.312412512512412</span><span style="color:red">2</span>

    $\delta x$ = <span style="color:green">0.000000012412512</span><span style="color:red">412512</span>

    $x_{n+1} = x_n + \delta x$ = <span style="color:green">1.312412524924924</span><span style="color:red">612512</span>

    Let us try to add them together in Python:
    ```python
    >>> x = 1.3124125125124122
    >>> dx = 0.000000012412512412512
    >>> x = x + dx
    >>> print(f"{x:.19f}")
    1.3124125249249245506
    ```

    The result is <span style="color:green">1.312412524924924</span><span style="color:red">5506</span>.

    Round off error $\approx 6.1912 \times 10^{-17}$

    Now we try compensated summation:
    ```python
    >>> x = 1.3124125125124122
    >>> dx = 0.000000012412512412512
    >>> e = 0.0
    >>> x_0 = x
    >>> e = e + dx
    >>> x = x_0 + e
    >>> e = e + (x_0 - x)
    >>> e
    6.12222940391771e-17
    ```
    As we can see, the value of $e$ is very close to the round off error!
    Therefore, the new variable $e$ effectively keeps track of the lost digits.
