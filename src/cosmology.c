#include <math.h>

#include "common.h"

WIN32DLL_API double compute_omega_k(
    const double omega_m,
    const double omega_lambda
)
{
    return 1.0 - omega_m - omega_lambda;
}

WIN32DLL_API double compute_G(
    const double omega_m,
    const double h0
)
{
    return 3.0 / 2.0 * omega_m * h0 * h0;
}

WIN32DLL_API double compute_da(
    const double a,
    const double h0,
    const double omega_m,
    const double omega_lambda
)
{
    const double omega_k = compute_omega_k(omega_m, omega_lambda);
    return h0 * a * sqrt(
        omega_lambda
        + omega_m / (a * a * a)
        + omega_k / (a * a)
    );
}
