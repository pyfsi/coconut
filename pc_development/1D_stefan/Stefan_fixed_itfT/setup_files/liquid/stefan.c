#include "udf.h"
#include <math.h>

real f(real x, real St) {
    return x * exp(x*x) * erf(x) - St / sqrt(M_PI);
}

real df(real x) {
    return (2*x*x + 1) * exp(x*x) * erf(x) + (2*x) / sqrt(M_PI);
}

real newton_raphson(real x0, real St, int tol) {
    real x_prev = x0;
    real x_current = x_prev - f(x_prev, St) / df(x_prev);
    real diff = fabs(x_prev - x_current);
    x_prev = x_current;

    while(diff > tol) {
        x_current = x_prev - f(x_prev, St) / df(x_prev);
        diff = fabs(x_prev - x_current);
        x_prev = x_current;
    }

    return(x_current);
}

DEFINE_INIT (temp_ini, d)
{
    real cp, k, rho, T_L, T_m, dT, L, x_target, alpha, St;
    cp = 2500.0;
    k = 1.5;
    rho = 870.0;
    T_L = 309.15;
    T_m = 299.15;
    dT = T_L - T_m;
    L = 179000.0;
    x_target = 0.05;
    alpha = k / (rho * cp);
    St = (cp * dT) / L;

    cell_t c;
    Thread *t;
    real lambda, t_sol;
    real xc[ND_ND];

    lambda = newton_raphson(sqrt(St / 2), St, 0.00000001);
    t_sol = pow(x_target, 2.) / (4 * pow(lambda, 2.) * alpha);

    /* loop over all cell threads in the domain */
    thread_loop_c(t,d) {
        /* loop over all cells */
        begin_c_loop_all(c,t) {
            C_CENTROID(xc,c,t);
            C_T(c,t) = T_L - dT * erf(xc[0] / (2 * sqrt(alpha * t_sol))) / erf(lambda);
        } end_c_loop_all(c,t)
    }
}