#ifndef COSMOLOGY_H
#define COSMOLOGY_H


double compute_omega_k(
    const double omega_m,
    const double omega_lambda
);

double compute_da(
    const double a,
    const double h0,
    const double omega_m,
    const double omega_lambda
);

#endif
