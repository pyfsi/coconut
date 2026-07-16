#include "udf.h"
#include <math.h>

real f(real x, real St_l, real St_s, real munu) {
    return St_l / (exp(x*x) * erf(x)) - St_s / (munu * exp(munu*munu*x*x) * erfc(munu*x)) - x * sqrt(M_PI);
}

real df(real x, real St_l, real St_s, real munu) {
    return -St_l * (2 * exp(-x*x) * x / erf(x) + 2 * exp(-2*x*x) / (sqrt(M_PI) * pow(erf(x), 2.))) - St_s * (2 * exp(-2*munu*munu*x*x) / (sqrt(M_PI) * pow(erfc(munu*x), 2.)) - 2 * munu * x * exp(-munu*munu*x*x) / erfc(munu*x)) - sqrt(M_PI);
}

real newton_raphson(real St_l, real St_s, real munu, int tol) {
    // first guess is an approximation
    real x_prev = 0.5 * (-St_s / (munu * sqrt(M_PI)) + sqrt(2 * St_l + pow(St_s / munu * sqrt(M_PI), 2.)));
    real x_current = x_prev - f(x_prev, St_l, St_s, munu) / df(x_prev, St_l, St_s, munu);
    real diff = fabs(x_prev - x_current);
    x_prev = x_current;

    while(diff > tol) {
        x_current = x_prev - f(x_prev, St_l, St_s, munu) / df(x_prev, St_l, St_s, munu);
        diff = fabs(x_prev - x_current);
        x_prev = x_current;
    }

    return(x_current);
}

DEFINE_INIT (temp_ini, d)
{
    real cp_s, cp_l, cp, k_s, k_l, rho_s, rho_l, T_L, T_m, T_S, L, x_target, alpha_s, alpha_l, St_s, St_l, nu, mu, munu, temperature;

    // PCM properties
    cp_s = 1900.0;
    cp_l = 2240.0;
    k_s = 0.334;
    k_l = 0.15;
    rho_s = 867.914;
    rho_l = 775.0;
    T_L = 311.13;
    T_m = 301.13;
    T_S = 291.13;
    L = 236980.0;
    x_target = 0.0008;

    alpha_s = k_s / (rho_s * cp_s);
    alpha_l = k_l / (rho_l * cp_l);
    St_s = cp_s * (T_m - T_S) / L;
    St_l = cp_l * (T_L - T_m) / L;
    nu = sqrt(alpha_l / alpha_s);
    mu = rho_l / rho_s;
    munu = mu * nu;

    cell_t c;
    Thread *t;
    real lambda, t_sol;
    real xc[ND_ND];

    lambda = newton_raphson(St_l, St_s, munu, 0.000000001);
    t_sol = pow(x_target, 2.) / (4 * pow(lambda, 2.) * alpha_l);

    /* loop over all cell threads in the domain */
    thread_loop_c(t,d) {
        /* loop over all cells */
        begin_c_loop_all(c,t) {
            C_CENTROID(xc,c,t);
            if (xc[0] <= x_target) {
                temperature = T_L - (T_L - T_m) * erf(xc[0] / (2 * sqrt(alpha_l * t_sol))) / erf(lambda);
            } else {
                temperature = T_S + (T_m - T_S) * erfc((xc[0] / (2 * sqrt(alpha_s * t_sol))) - (1 - mu) * nu * lambda) / erfc(lambda * munu);
            }
            C_T(c,t) = temperature;
            cp = -1029.0 + 9.797 * temperature;
            C_H(c,t) = cp * (temperature - 298.15);
        } end_c_loop_all(c,t)
    }
}