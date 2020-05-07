# Predictors

This documentation describes the available predictors.
A predictor is used to determine the initial guess in each time step by extrapolating the solution from previous time steps.
The different predictors differ in the number of previous iterations they take into account 
and the polynomial degree that is used.

Only the `type` has to be provided, no `settings` dictionary is required.

Specification of a predictor is mandatory, also for a steady simulation. In that case however, it does not matter which
predictor is chosen as only 1 "time step" is performed.

## Linear
This predictor uses the results from the last two time steps to determine the initial guess in the current time step 
using a linear extrapolation as follows
$$
x^{n+1}=2x^{n}-1x^{n-1}.
$$
If no two previous solutions are available, the solution from the last time step is used.

## Legacy
This predictor uses the results from the last three time steps to determine the initial guess in the current time step
using a linear extrapolation as follows
$$
x^{n+1}=\frac{5}{2}x^{n}-2x^{n-1}+\frac{1}{2}x^{n-2}.
$$
If no three previous solutions are available, the predictor `linear` is used.
This predictor is called `legacy` as it corresponds to the 2nd order extrapolator in the coupling code Tango.

## Quadratic
This predictor uses the results from the last three time steps to determine the initial guess in the current time step
using a quadratic extrapolation as follows
$$
x^{n+1}=3x^{n}-3x^{n-1}+1x^{n-2}.
$$
If no three previous solutions are available, the predictor `linear` is used.

## Cubic
This predictor uses the results from the last three time steps to determine the initial guess in the current time step
using a cubic extrapolation as follows
$$
x^{n+1}=4x^{n}-6x^{n-1}+4x^{n-2}-1x^{n-3}.
$$
If no four previous solutions are available, the predictor `quadratic` is used.
