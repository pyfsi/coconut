# Predictors

This documentation describes the available predictors.
A predictor is used to determine the initial guess in each time step by extrapolating the solution from previous time steps.
The predictors differ in the number of previous time steps they take into account 
and the polynomial degree that is used.
The formulas in this document, use the index $n$ to refer to the time step, where $n+1$ is the current time step, i.e. the time step for which the initial guess is made.
The vector $x$ is the input for the first solver, conform the [coupled solvers documentation](../coupled_solvers/coupled_solvers.md).

A predictor is intialized in the coupled solver using an initial solution as determined by the coupled solver.
As such, there is at least one previous solution available.

Only the `type` has to be provided (e.g. `predictors.linear`), no `settings` dictionary is required.

Specification of a predictor is mandatory, also for a steady simulation. In that case, however, it does not matter which
predictor is chosen as only one "time step" is performed.

## Constant
The `type` for this predictor is `predictors.constant`.
This predictor uses the result from the previous solution as the initial guess in the current time step:
$$
x^{n+1}=x^{n}.
$$

## Linear
The `type` for this predictor is `predictors.linear`.
This predictor uses the results from the last two time steps to determine the initial guess in the current time step 
using a linear extrapolation as follows
$$
x^{n+1}=2x^{n}-1x^{n-1}.
$$
If no two previous solutions are available, the solution from the last time step is used, i.e. the predictor `PredictorConstant`.

## Legacy
The `type` for this predictor is `predictors.legacy`.
This predictor uses the results from the last three time steps to determine the initial guess in the current time step
using a linear extrapolation as follows
$$
x^{n+1}=\frac{5}{2}x^{n}-2x^{n-1}+\frac{1}{2}x^{n-2}.
$$
If no three previous solutions are available, the predictor `PredictorLinear` is used.
This predictor is called `PredictorLegacy` as it corresponds to the 2^{nd} order extrapolator in the coupling code _Tango_.

## Quadratic
The `type` for this predictor is `predictors.quadratic`.
This predictor uses the results from the last three time steps to determine the initial guess in the current time step
using a quadratic extrapolation as follows
$$
x^{n+1}=3x^{n}-3x^{n-1}+1x^{n-2}.
$$
If no three previous solutions are available, the predictor `PredictorLinear` is used.

## Cubic
The `type` for this predictor is `predictors.cubic`.
This predictor uses the results from the last three time steps to determine the initial guess in the current time step
using a cubic extrapolation as follows
$$
x^{n+1}=4x^{n}-6x^{n-1}+4x^{n-2}-1x^{n-3}.
$$
If no four previous solutions are available, the predictor `PredictorQuadratic` is used.
