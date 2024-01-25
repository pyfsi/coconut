# Predictors

This documentation describes the available predictors.
A predictor is used to determine the initial guess in each time step by extrapolating the solution from previous time steps.
The predictors differ in the number of previous time steps they take into account and the polynomial degree that is used.
Additionally, there is one special type that uses a surrogate model.
The formulas in this document, use the index $n$ to refer to the time step, where $n+1$ is the current time step, i.e. the time step for which the initial guess is made.
The vector $x$ is the input for the first solver, conform the [coupled solvers' documentation](../coupled_solvers/coupled_solvers.md).

A predictor is intialized in the coupled solver using an initial solution as determined by the coupled solver.
As such, there is at least one previous solution available.

For the polynomial predictors, the `predictor` dictionary only requires a `type` (e.g. `predictors.linear`) in the JSON file, no `settings` dictionary has to be provided.

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
This predictor is called `PredictorLegacy` as it corresponds to the second order extrapolator in the coupling code _Tango_, the predecessor of _CoCoNuT_.

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

## Surrogate
The `type` for this predictor is `predictors.surrogate`.

The following parameters may be included in the `settings` dictionary.

This predictor requires that a `surrogate` model is used in the coupled solver [`CoupledSolverIQNISM`](../coupled_solvers/coupled_solvers.md#iqnism) (defined in the settings of the coupled solver).
Moreover, this `surrogate` model must provide a surrogate solution, which is used to update the values in this predictor.

There are two options available:

- The surrogate solution can be used directly: $$ x^{n+1}=x_s^{n+1}, $$ where $x_s^{n+1}$ is the surrogate solution of the current time step.
- Or the change in surrogate solution can be used: $$ x^{n+1}=x^{n} + (x_s^{n+1} - x_s^{n}), $$ where $x_s$ is the surrogate solution, and the superscript $n+1$ and $n$ indicate the current and previous time steps, respectively.

The following parameters need to be included in the `settings` dictionary.

|                     parameter | type | description                                                                                                                                                 |
|------------------------------:|:----:|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| <nobr>`predict_change`</nobr> | dict | (optional) Default: `true`. Indicates it the change in surrogate solution should be used. If `false`, the surrogate solution serves as prediction directly. |

## Restart

Upon restart, the predictor may be changed.
When changing an extrapolator (Constant, Linear, Quadratic, Legacy and Cubic) to a higher order, the first extrapolation will still be of the original order, but the order will increase with each time step until the new order is reached,
just like at the start of a simulation.
Changing to a lower order has direct effect.

No information is transferred when switching from or to a [surrogate predictor](#surrogate).

## Dummy predictor

This dummy predictor can be used in the [one-way](../coupled_solvers/coupled_solvers.md#one-way) coupled solvers, which doesn't require a predictor.

If the use of a dummy predictor is allowed, the `type` (`convergence_criteria.dummy_convergence_criterion`) can be written explicitly or omitted. No `settings` are required.
Note that the key `predictor` is still required.
