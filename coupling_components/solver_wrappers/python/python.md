# Python

This is the documentation for all Python solver wrappers.
Currently, only solvers exist for the one-dimensional (1D) calculation of a straight flexible tube.

## Tube

There are three tube solvers for a straight tube with circular cross-section.
The axial direction is along the z-axis.
All of them are 1D solvers.
This means the tube is divided in `m` intervals in the axial z-direction
and only variations in that direction are taken into account.
In other words, all variables are uniform over a cross-section.
They are calculated in the center of each of the `m` cells.

There is one flow solver `SolverWrapperTubeFlow` and two structure solvers, one with inertia `SolverWrapperTubeStructure` and one without `SolverWrapperTubeRingmodel`.
These solvers are very simple and provide insight in the physics, especially in the stability of fluid-structure interaction simulation.
Nonetheless, they are not meant to provide an accurate representation of reality.
Especially the pressure stabilization term in the flow solver smooths the pressure distribution.

### Settings

The following parameters, listed in alphabetical order, need to be provided in the main JSON parameter file as `settings` of the solver wrapper.

|                        parameter |  type  | description                                                                                                                                                                                                                                                                                                                                                                                                      |
|---------------------------------:|:------:|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|                          `debug` |  bool  | (optional) Default: `false`. For every iteration, text files are saved with the input and output data of the solver.                                                                                                                                                                                                                                                                                             |
|                        `delta_t` | double | Fixed time step size. This parameter is usually specified in a higher component.                                                                                                                                                                                                                                                                                                                                 |
|                     `input_file` |  str   | (optional) Name of the input file, which may be present in the folder given in `working_directory`. The file contains parameters required for the solver, in JSON-format. The parameters specified in the main parameter JSON file have priority over the parameters defined in this file.                                                                                                                       |
|                `interface_input` |  list  | List of dictionaries; each dictionary requires two keys: `model_part` and `variables`. The former contains the name of the `ModelPart` as a string. The value of the latter is a list of variables. Even if there is only one variable, a list is required. For the Python solver wrappers these variables are fixed: `['displacement']` for a flow solver and `['pressure','traction']` for a structure solver. |
|               `interface_output` |  list  | Analogous to `interface_input`. However, the variables are different: `['pressure','traction']` for a flow solver and `['displacement']` for a structure solver.                                                                                                                                                                                                                                                 |
|                       `unsteady` |  bool  | (optional) Default: `true`. Indicates if case is steady or unsteady.                                                                                                                                                                                                                                                                                                                                             |
|                   `save_restart` |  int   | (optional) Default: 0. Determines the time step interval upon which a pickle file *`case_timestep<time step>.pickle`* is written, to be used for restart purposes. A minus sign indicates only the file from the last interval is retained.                                                                                                                                                                      |
|                `time_step_start` |  int   | (optional) Default: 0. Time step number to (re)start a transient FSI calculation. If `0` is given, the simulation starts from scratch. Otherwise, the code looks for the pickle file *`case_timestep<timestep_start>.pickle`* to start from the corresponding time step. For a steady simulation, the value should be `0`.                                                                                       |
| <nobr>`working_directory`</nobr> |  str   | Absolute path to the working directory or relative path with respect to the current directory.                                                                                                                                                                                                                                                                                                                   |

`delta_t` is a necessary parameter, while `timestep_start` and `save_restart` are optional, but all are usually defined in a higher component. However, they can also be given directly as parameter of the solver wrapper (e.g. for standalone testing). If they are defined both in higher component and in the solver wrapper, then the former value is used and a warning is printed.

Because these solvers are simple, it is possible to omit the use of interpolation between two solver wrappers.
In that case, the names of the `ModelParts` of both flow and structure solver need to be the same.

### Tube flow solver
This flow solver calculates the flow inside a 1D straight and flexible tube.
The `type` for this solver wrapper is `solver_wrappers.python.tube_flow_solver`.
The required input is the radial displacement of the tube wall.
The other components of displacement are ignored.
The resulting output is the pressure on the tube wall.
Although a required variable in `interface_output`, `traction` is always returned as being zero.

From the radial displacements the area of each cross-section is calculated.
The flow is governed by the continuity and momentum equation
$$
\frac{\partial a}{\partial t}+\frac{\partial av}{\partial z}=0
$$
$$
\frac{\partial av}{\partial t}+\frac{\partial av^2}{\partial z}+\frac{1}{\rho_f}\left(\frac{\partial ap}{\partial z}-p\frac{\partial a}{\partial z}\right)=0
$$
with $a=\pi r^2$ the cross-sectional area of the tube and $r$ the inner radius.
Furthermore, $t$ is the time, $v$ the velocity along the axis of the tube, $p$ the pressure and $\rho_f$ the density of the fluid.
At the start and at the end of the tube the boundary conditions have to be applied.
Discretizing these equations results in a system of linear algebraic equations which can be solved for the pressure and velocity
by performing Newton-Raphson iterations.
Inside the solver the kinematic pressure is used, which is the pressure divided by the density of the fluid.

The boundary conditions are implemented by four additional equations: at the in- and outlet for both pressure and velocity.
At the inlet, the pressure, velocity or total pressure can be specified.
At the outlet, the pressure can be set to a fixed value or a non-reflecting boundary conditions.
These settings are specified in the `input_file` in the `working directory`.

For more information about the implementation of this solver refer to [[1](#1)].

Finally, this solver also provides a Jacobian of the change in radius with pressure, which can be used as surrogate Jacobian [[4](#4)].

#### Solver parameters

The following parameters, listed in alphabetical order, need to be specified in the main JSON file or in a file with name *`input_file`*, located in the `working_directory`.
Care should be taken that the values of `d`, `e`, `h`, `l` and `rhof` match the corresponding values of the structural solver.

parameter|type|description
---:|:---:|---
`axial_offset`|double|(optional) Default: `0`. Distance over which tube is displaced axially in the coordinate system.
`d`|double|Nominal diameter of the tube.
`e`|double|Modulus of elasticity of the tube wall.
`h`|double|Thickness of the tube wall.
`inlet_boundary`|dict|Dictionary containing all information with respect to the inlet boundary condition.
`l`|double|Length of the tube.
`m`|int|Number of cells for discretization. The values are calculated in the cell centers.
`newtonmax`|double|Maximum number of Newton-Raphson iterations for the flow solver calculation.
`newtontol`|double|Relative tolerance for the Newton-Raphson iterations for the flow solver calculation.
`outlet_boundary`|dict|Dictionary containing all information with respect to the outlet boundary condition.
`preference`|double|(optional) Default: `0`. Reference pressure and initial pressure.
`u0`|double|(optional) Default: `ureference`. Initial velocity throughout the tube.
`ureference`|double|Reference velocity used for determination of pressure stabilization term.
`rhof`|double|Density of the fluid.

##### Inlet Boundary

This section describes all parameters that need to be specified in the dictionary `inlet_boundary`, listed in alphabetical order.

parameter|type|description
---:|:---:|---
`amplitude`|double|Amplitude of the inlet boundary condition.
`period`|double|Period of the inlet boundary condition. Period of oscillation for a periodic boundary condition, duration for a non-periodic boundary condition. Not used for a fixed value boundary condition (type `4`).
`reference`|double|(optional) Reference value of inlet boundary condition. If not provided, the value of this parameter is the corresponding reference value provided above, i.e. `ureference`, `preference` or `preference` + `rhof` * `ureference`^2 / 2.
`type`|int|Type of inlet boundary condition. <br>If `1`, a sine wave is used with amplitude, reference and period as specified. <br>If `2`, a pulse is used with amplitude as specified and a duration equal to the parameter period. After the pulse the variable is equal to the reference value. <br>If `3`, a quadratic sine wave is used with amplitude, reference and period as specified. <br>If `4`, a fixed value equal to the sum of the reference value and the amplitude. Used for steady cases. <br>If other, a steady increase of the value at the inlet with slope of amplitude divided by period is used.
`variable`|str|Variable upon which the inlet boundary condition is defined, either `'pressure'`, `'velocity'` or `'total pressure'`.

##### Outlet Boundary

This section describes all parameters that need to be specified in the dictionary `outlet_boundary`, listed in alphabetical order.

parameter|type|description
---:|:---:|---
`type`|int|Type of outlet boundary condition. <br>If `1`, a non-reflecting boundary condition is applied. This type cannot be used for a steady calculation. <br>If other, a fixed value equal to the reference pressure is applied.

### Tube Ringmodel solver

This structure solver calculates the deformation of the wall of a straight and flexible tube.
The `type` for this solver wrapper is `solver_wrappers.python.ring_model_solver`.
The tube is regarded as made up of independent rings and no inertia is taken into account.
Therefore, there is no dependence on previous time steps and the parameters `delta_t` and `timestep_start` are not used.
The required input is the pressure on the tube wall.
Traction is not taken into account, even though it is a required variable in `interface_input`. 
The resulting output is the radial displacement.
For the other components of displacement, zero is returned.
This solver is not suited to calculate the propagation of a pressure pulse.

The behaviour of the elastic tube wall is described by a Hookean relation,
which results in the following equation
$$
a=a_0\left(\frac{p_0-2c^2_{MK}}{p_0-2c^2_{MK}}\right)^2
$$
with $a=\pi r^2$ the cross-sectional area of the tube and $r$ the inner radius.
Furthermore, $p$ the pressure, $p_0$ the reference pressure, $\rho_f$ the density of the fluid
and $c^2_{MK}$ the Moens-Korteweg wave speed given by
$$
c^2_{MK}=\sqrt{\frac{Eh}{2\rho_f r_0}}
$$
Here, $E$ is the modulus of elasticity, $h$ the thickness of the tube wall and $r_0$ the reference radius.
No boundary conditions are required.
Inside the solver the kinematic pressure is used, which is the pressure divided by the density of the fluid.

More information about the implementation of this solver can be found in [[2](#2)].

#### Solver parameters

The following parameters, listed in alphabetical order, need to be specified in the main JSON file or in a file with name *`input_file`*, located in the `working_directory`.
Care should be taken that the values of `d`, `e`, `h`, `l` and `rhof` match the corresponding values of the flow solver.

parameter|type|description
---:|:---:|---
<nobr>`axial_offset`</nobr>|double|(optional) Default: `0`. Distance over which tube is displaced axially in the coordinate system.
`d`|double|Nominal diameter of the tube.
`e`|double|Modulus of elasticity of the tube wall.
`h`|double|Thickness of tube wall.
`l`|double|Length of the tube.
`m`|int|Number of cells for discretization. The values are calculated in the cell centers.
`preference`|double|(optional) Default: `0`. Reference pressure.
`rhof`|double|Density of the fluid.

### Tube structure solver

This structure solver calculates the deformation of the wall of a straight and flexible tube.
The `type` for this solver wrapper is `solver_wrappers.python.tube_structure_solver`.
In this model inertia is taken into account, but still only radial displacement is considered.
The required input is the pressure on the tube wall.
Traction is not taken into account, even though it is a required variable in `interface_input`. 
The resulting output is the radial displacement.
For the other components of displacement, zero is returned.

The deformation of the tube in the radial direction is determined by
$$
	\rho_s h\frac{\partial^2 r}{\partial t^2}+b_1\frac{\partial^4 r}{\partial z^4}-b_2\frac{\partial^2 r}{\partial z^2}+b_3(r-r_o)=p-p_o
$$
with $\rho_s$ the solid density and $h$ the thickness of the tube wall.
Further, $r$ is the inner radius, $p$ pressure and $t$ time. The variables $r_0$ and $p_0$ are a reference radius and pressure, respectively.
The parameters $b_1$ and $b_2$ ($b_1, b_2 \ge 0$) account for the inner action of the bending and the axial tension in the wall, while the parameter $b_3$ accounts for the circumferential stress in the tube wall.
For a thin-walled tube that is clamped in the axial direction, they are defined as
$$
	b_1=\frac{hE}{1-\nu^2}\frac{h^2}{12} 
	\textrm{, }
	b_2=\frac{hE}{1-\nu^2}\frac{h^2}{12}\frac{2\nu}{r_o^2}
	\textrm{ and }
	b_3=\frac{hE}{1-\nu^2}\frac{1}{r_o^2}+\frac{hE}{1-\nu^2}\frac{h^2}{12}\frac{1}{r_o^4}
$$
with $E$ the Young's modulus and $\nu$ Poisson's coefficient.
The second term of $b_3$ is considered small compared to the first one because $h\ll r_o$ and is thus neglected.
Discretizing these equations results in a system of linear algebraic equations with a Jacobian that does not depend on the pressure, nor the radius.
This system is solved for the radius without requiring iterations.

There are two solving options available, specified with the parameter `solver`. The fastest and most memory efficient is the "solve_banded" option, which writes the matrix in a diagonal ordered form.
In this way, memory use and the number of operations is reduced.
However, when the number of discretization intervals increases, the condition number of the coefficient matrix quickly increases.
For high condition number, the "direct" method is more accurate and allows deeper convergence.
This option calculates the inverse of the coefficient matrix directly, at the start of the calculation.
The inverse is only calculated once and thereafter stored. When the number of intervals is high, it is clear that this option will use a high amount of memory.

The tube is considered clamped at both ends. This boundary condition is imposed by adding four equations: two at the inlet and two at the outlet.

For more information about the implementation of this solver refer to [[1](#1)].
In this work the Newmark-beta time discretization is used. Here, however, backward Euler time discretization is used as well.

For the Newmark-beta time discretization, two Newmark parameters $\beta$ and $\gamma$ are required, which result in an unconditionally stable integration scheme if
$$
\gamma\geq\frac{1}{2} 
\textrm{ and }
\beta\geq\frac{1}{4}\left(\frac{1}{2}+\gamma\right)^2.
$$
Typical values are $\gamma$ equal to 1/2 and $\beta$ equal to 1/4.

A different time discretization for flow and structure can lead to difficulties for strongly coupled problems, especially looking at the resulting pressure distributions.
As most flow solvers are discretized using the backward Euler method, it is advised to chose the same method for the structure solver (the time discretization of `SolverWrapperTubeFlowSolver` is also backward Euler).
This avoids the occurrence of spurious oscillations of the pressure in time [[3](#3)].

Finally, this solver also provides a Jacobian of the change in radius with pressure, which can be used as surrogate Jacobian [[4](#4)].

#### Solver parameters

The following parameters, listed in alphabetical order, need to be specified in the main JSON file or in a file with name *`input_file`*, located in the `working_directory`.
Care should be taken that the values of `d`, `e`, `h`, `l` and `rhof` match the corresponding values of the flow solver.

parameter|type|description
---:|:---:|---
`axial_offset`|double|(optional) Default: `0`. Distance over which tube is displaced axially in the coordinate system.
`beta`|double|(optional) Newmark parameter $\beta$. Only required when the Newmark-beta time discretization is used.
`d`|double|Nominal diameter of the tube.
`e`|double|Modulus of elasticity of the tube wall.
`h`|double|Thickness of the tube wall.
`gamma`|double|(optional) Newmark parameter $\gamma$. Only required when the Newmark-beta time discretization is used.
`l`|double|Length of the tube.
`m`|int|Number of cells for discretization. The values are calculated in the cell centers.
`nu`|double|Poisson's ratio.
`preference`|double|(optional) Default: `0`. Reference pressure.
`rhof`|double|Density of the fluid.
`rhos`|double|Density of the tube wall.
`solver`|str|(optional) Default: `solve_banded`. Either `solve_banded` or `direct`. Specifies the solution method for the linear system of equations.
`time_disretization`|str|(optional) Default: `backward Euler`. Specifies the time discretization: either `Newmark` or `backward Euler`. Not case sensitive.

## References
<a id="1">[1]</a> 
[Degroote J., Annerel S. and Vierendeels J., "Stability analysis of Gauss-Seidel iterations in a partitioned simulation of fluid-structure interaction", Computers & Structures, vol. 88, no. 5-6, pp. 263, 2010.](http://hdl.handle.net/1854/LU-940283)

<a id="2">[2]</a> 
[Degroote J., Bruggeman P., Haelterman R. and Vierendeels J., "Stability of a coupling technique for partitioned solvers in FSI applications", Computers & Structures, vol. 86, no. 23–24, pp. 2224–2234, 2008.](http://hdl.handle.net/1854/LU-533350)

<a id="3">[3]</a> 
[Vierendeels J., Dumont K., Dick E. and Verdonck P., "Analysis and stabilization of fluid-structure interaction algorithm for rigid-body motion", American Institute of Aeronautics and Astronautics Journal", vol. 43, no. 12, pp. 2549–2557, 2005.](http://hdl.handle.net/1854/LU-325786)

<a id="4">[4]</a>
[Delaissé N., Demeester T., Fauconnier D. and Degroote J., "Surrogate-based acceleration of quasi-Newton techniques for fluid-structure interaction simulations", Computers & Structures, vol. 260, pp. 106720, 2022.](http://hdl.handle.net/1854/LU-8728347)