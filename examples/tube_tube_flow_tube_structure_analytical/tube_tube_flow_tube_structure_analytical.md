# Tube case with TubeFlow and TubeStructure and with Analytical1D model

This example is identical to Python TubeFlow - Python TubeStructure, with the only difference that a different model is used.

## Coupling algorithm

The coupled solver used is `coupled_solvers.iqni`, just as in the example Python TubeFlow - Python TubeStructure.
Only, now, the model is not `coupled_solvers.models.ls`, but `coupled_solvers.models.analytical_1d`, which uses two solver models (that should be the same as the actual solvers) to calculate an analytical Jacobian.

Note that the relaxation factor $omega$ may be equal to 1, as the analytical model is ready from the start, in contradiction to a secant model, which needs to collect input-output-pairs first.