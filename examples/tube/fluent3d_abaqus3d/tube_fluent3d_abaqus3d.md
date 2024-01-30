# Tube case with Fluent3D and Abaqus3D

This example calculates the flow inside and the deformation and stresses of a straight flexible tube, where a pressure pulse is applied at the inlet.
This done by using Fluent and Abaqus, both with a fully 3D case.

## Coupling algorithm

The coupling technique used is the *interface quasi-Newton algorithm with an approximation for the inverse of the Jacobian from a least-squares model* (IQNI-LS).

## Predictor

The initial guess in every time step is done using the linear predictor.

## Convergence criterion

Two convergence criteria have been specified:

- The number of iterations in every time step is larger than 20.
- The residual norm on the displacement is a factor $10^{-3}$ lower than the initial value.
 
When either criterion is satisfied the simulation stops.
 
## Solvers

The flow solver is Fluent, used to solve a fully 3D tube,
with 48 cells on the fluid-structure interface in the length direction and 8 in the circumferential direction.
When setting up the case, the mesh is build based on the file *`mesh.jou`* using Gambit.
The displacements are applied in the nodes. In contrast, the loads (pressure and traction) are calculated in the cell centers.
The axial direction is along the x-axis.
The setup script runs Fluent with the *`case.jou`* journal file to set up the case parameters, starting from the mesh file *`mesh_tube3d.msh`*.
This case is written to the *`case_tube3d.cas`* file, which serves as input for CoCoNuT. 
Additionally, a folder *`create_mesh`* is provided containing a script to create the mesh in Gambit using a journal file.
The mesh can be created by running the script *`create_mesh.sh`*, given that Gambit v2.4.6 is available.

The structure solver is Abaqus, used to solve a fully 3D tube,
with 12 elements on the fluid-structure interface in the length direction and 8 in the circumferential direction.
The Abaqus case is built when setting up the case starting from the file *`mesh_tube3d.inp`* containing nodes and elements. 
This is done by running Abaqus with the *`make_inp.py`* Python script to set all parameters, such as surface definitions, material parameters, boundary conditions and time step information.
The result of the setup is a completed input file *`case_tube3d.inp`*.
The Abaqus element type used is S8R. These are 8-node shell elements. The R refers to reduced integration.
See the [Abaqus documentation](http://130.149.89.49:2080/v6.14/books/usb/default.htm?startat=book01.html#usb) for more information. 
The loads are applied on the faces in 9 points per element, which means on 864 load points in total. 
The displacement is calculated in the nodes. There are 304 nodes on the fluid-structure interface.
The axial direction is along the x-axis.

The coordinate frames of both solvers are the same so there is no need for a permutation mapper.
In contrast, the difference of the points where loads and displacements are applied or calculated, require the use of interpolation mappers.
Therefore, a radial basis mapper is introduced in the structure solver to interpolate in the x-, y- and z-direction.