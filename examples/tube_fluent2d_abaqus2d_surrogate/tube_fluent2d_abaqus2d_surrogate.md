# Tube case with TubeFlow and TubeStructure and Surrogate model

This example is identical to the example Fluent2D and Abaqus2D, except for the use of a surrogate model.
This is reflected in a different coupled solver and corresponding models, and the type of predictor.
The used surrogate model is the combination of the 1D Python solvers TubeFlow and TubeStructure.

## Coupling algorithm

The coupled solver used is `coupled_solvers.iqnism`, which has a `model` and `surrogate` setting, both of which are a type of model component.
The secant model used is `coupled_solvers.models.mvmf` where 100 time steps are reused: the previous time steps are used as surrogate model.
In addition, a `surrogate` model is defined, which uses the 1D Python solvers TubeFlow and TubeStructure.
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
 
## Solvers

The flow solver is Fluent, used to solve an axisymmetric representation of the tube,
with 100 cells on the fluid-structure interface. 
When setting up the case, the mesh is build based on the file *`mesh.jou`* using Gambit.
The displacements are applied in the nodes, of which there are 101. 
In contrast, the loads (pressure and traction) are calculated in the cell centers, of which there are 100.
The axial direction is along the x-axis, the radial direction along the y-axis. 
The setup script runs Fluent with the *`case.jou`* journal file to setup the case parameters, starting from the mesh file *`mesh_tube2d.msh`*.
This case is written to the *`case_tube2d.cas`* file, which serves as input for CoCoNuT. 
Additionally, a folder *`create_mesh`* is provided containing a script to create the mesh in Gambit using a journal file.
The mesh can be created by running the script *`create_mesh.sh`*, given that Gambit v2.4.6 is available.

The structure solver is Abaqus, used to solve an axisymmetric representation of the tube,
with 50 elements on the fluid-structure interface.
The Abaqus case is built when setting up the case starting from the file *`mesh_tube2d.inp`* containing nodes and elements. 
This is done by running Abaqus with the *`make_inp.py`* Python script to set all parameters, such as surface definitions, material parameters, boundary conditions and time step information.
The result of the setup is a completed input file *`case_tube2d.inp`*.
The Abaqus element type used is CAX8RH. These are continuum elements for axisymmetric calculations, for stress and displacement without twist. 
They are: 8-node biquadratic, reduced integration, hybrid with linear pressure.
See the [Abaqus documentation](http://130.149.89.49:2080/v6.14/books/usb/default.htm?startat=book01.html#usb) for more information. 
The loads are applied on the faces in three points per element, which means on 150 load points in total. 
The displacement is calculated in the nodes. There are 101 nodes on the fluid-structure interface.
The axial direction is along the y-axis, the radial direction along the x-axis.

The difference in reference frames and number of cells on the fluid-structure interface requires the use of mappers.
In the structure solver wrapper, a permutation mapper is introduced to match the coordinate frames, flipping the x- and y-axis of the input.
Thereafter, a linear interpolation mapper is used to interpolate in the x- and y-direction from the 100 cell centers of Fluent to the 150 load points in Abaqus.
For the output the same is done in the opposite order: first interpolating and then flipping the axes.