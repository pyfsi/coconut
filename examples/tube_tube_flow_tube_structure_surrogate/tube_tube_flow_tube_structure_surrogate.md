# Tube case with TubeFlow and TubeStructure and with Surrogate model

This example is identical to Python TubeFlow - Python TubeStructure, except for the use of a surrogate model.
This is reflected in a different coupled solver and corresponding models, and the type of predictor.

## Coupling algorithm

The coupled solver used is `coupled_solvers.iqnism`, which has a `model` and `surrogate` setting, both of which are a type of model component.
The secant model used is `coupled_solvers.models.mvmf` where 100 time steps are reused: the previous time steps are used as surrogate model.
In addition, a `surrogate` model is defined, which uses similar solvers, but with `50` degrees of freedom (`m`), instead of `100`.
Note, that the settings in the JSON file for these solvers overwrite the setting given in the `input_file`.
Because of this different discretization, mapping is required. This is achieved using the `mapped` model.
The `surrogate` model has its own coupled solver with its own predictor, convergence criterion, solver wrappers and models.
The model is once more `coupled_solvers.models.mvmf` that reuses 100 time steps as surrogate model.
Note that the surrogate solvers have their own working directory.

The `settings` of `coupled_solvers.iqnism` dictate that `surrogate_synchronize` is enabled.
This means that at the end of each time step, the surrogate solvers are solved once more with final solution of the actual solvers, such that they are synchronized.
This is also shown in the printing of the residuals.

Note that the relaxation factor $\omega$ is omitted in `coupled_solvers.iqnism`, such that its default value, equal to 1, is used.
This is possible because the surrogate model provides a Jacobian, when the secant model is not yet able to.
The coupled solver `coupled_solvers.iqni` inside of the `surrogate` does require a relaxation factor.

## Predictor

The predictor of the `coupled_solvers.iqnism` is `predictors.surrogate`, which means the surrogate solution is used for the prediction of the actual solver.
As the `predict_change` is not disabled, it is the change in surrogate solution with respect to the previous time step that will be used to determine the prediction for the current time step starting from the previous one.

## Mappers

Mappers are used in three locations:

- To map the input and output of the actual solvers
- To map the input and output of the surrogate solvers
- To map the discretization of the surrogate and actual solvers

Here, `mappers.linear` are used each time.