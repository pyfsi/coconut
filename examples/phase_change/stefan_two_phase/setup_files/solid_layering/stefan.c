#include "udf.h"
#include <math.h>

real f(real x, real St_l, real St_s, real nu) {
    return St_l / (exp(x*x) * erf(x)) - St_s / (nu * exp(nu*nu*x*x) * erfc(nu*x)) - x * sqrt(M_PI);
}

real df(real x, real St_l, real St_s, real nu) {
    return -St_l * (2 * exp(-x*x) * x / erf(x) + 2 * exp(-2*x*x) / (sqrt(M_PI) * pow(erf(x), 2.))) - St_s * (2 * exp(-2*nu*nu*x*x) / (sqrt(M_PI) * pow(erfc(nu*x), 2.)) - 2 * nu * x * exp(-nu*nu*x*x) / erfc(nu*x)) - sqrt(M_PI);
}

real newton_raphson(real St_l, real St_s, real nu, int tol) {
    // first guess is an approximation
    real x_prev = 0.5 * (-St_s / (nu * sqrt(M_PI)) + sqrt(2 * St_l + pow(St_s / nu * sqrt(M_PI), 2.)));
    real x_current = x_prev - f(x_prev, St_l, St_s, nu) / df(x_prev, St_l, St_s, nu);
    real diff = fabs(x_prev - x_current);
    x_prev = x_current;

    while(diff > tol) {
        x_current = x_prev - f(x_prev, St_l, St_s, nu) / df(x_prev, St_l, St_s, nu);
        diff = fabs(x_prev - x_current);
        x_prev = x_current;
    }

    return(x_current);
}

DEFINE_INIT (temp_ini, d)
{
    real cp, k_s, k_l, rho, T_L, T_m, T_S, L, x_target, alpha_s, alpha_l, St_s, St_l, nu, temperature;

    // PCM properties
    cp = 2500.0;
    k_l = 1.5;
    k_s = 0.024;
    rho = 870.0;
    T_L = 309.15;
    T_m = 299.15;
    T_S = 289.15;
    L = 179000.0;
    x_target = 0.050;

    alpha_s = k_s / (rho * cp);
    alpha_l = k_l / (rho * cp);
    St_s = cp * (T_m - T_S) / L;
    St_l = cp * (T_L - T_m) / L;
    nu = sqrt(alpha_l / alpha_s);

    cell_t c;
    Thread *t;
    real lambda, t_sol;
    real xc[ND_ND];

    lambda = newton_raphson(St_l, St_s, nu, 0.000000001);
    t_sol = pow(x_target, 2.) / (4 * pow(lambda, 2.) * alpha_l);

    /* loop over all cell threads in the domain */
    thread_loop_c(t,d) {
        /* loop over all cells */
        begin_c_loop_all(c,t) {
            C_CENTROID(xc,c,t);
            if (xc[0] <= x_target) {
                temperature = T_L - (T_L - T_m) * erf(xc[0] / (2 * sqrt(alpha_l * t_sol))) / erf(lambda);
            } else {
                temperature = T_S + (T_m - T_S) * erfc(xc[0] / (2 * sqrt(alpha_s * t_sol))) / erfc(lambda * nu);
            }
            C_T(c,t) = temperature;
            C_H(c,t) = cp * (temperature - 298.15);
        } end_c_loop_all(c,t)
    }
}