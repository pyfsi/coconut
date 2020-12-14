# SolverWrapperPython

This is the documentation for all Python SolverWrappers.
Currently only solvers exist for the 1D calculation of a straight flexible tube.

## Tube

There are three tube solvers for a straight tube with circular cross-section.
The axial direction is along the z-axis.
All of them are 1D solvers.
This means the tube is divided in `m` intervals in the axial z-direction
and only variations in that direction are taken into account.
In other words, all variables are uniform over a cross-section.
They are calculated in the center of each of the `m` cells.

There is one flow solver (TubeFlow) and two structure solvers, one with inertia (TubeStructure) and one without (TubeRingmodel).

### Settings

The following parameters, listed in alphabetical order, need to be provided in the JSON parameter file 
as `setttings` of the solver wrapper.

parameter|type|description
---:|:---:|---
`delta_t`|double|Fixed time step size in flow solver. This parameter is usually specified in a higher Component.
`input_file`|string|Name of the input file, which must be present in the folder `working_directory`. The file contains all parameters required for the solver, in JSON-format.
`interface_input`|dict|Keys are names of ModelParts. The values are (lists of) names of Variables. <br>For a flow solver:`"DISPLACEMENT"`. <br>For a structure solver:`["PRESSURE","TRACTION"]`.
`interface_output`|dict|Analogous to `interface_input`, but for the flow solver output. <br>For a flow solver:`["PRESSURE","TRACTION"]`. <br>For a structure solver:`"DISPLACEMENT"`.
`unsteady`|bool|(optional) Indicates if case is steady or unsteady. If omitted, `true` is assumed.
<nobr>`working_directory`</nobr>|string|Absolute path to the working directory or relative path w.r.t the current directory.

`delta_t` is a necessary parameter, but is usually defined in a higher Component. However, it can also be given directly as parameter of the solver wrapper (e.g. for standalone testing). If it is defined both in higher Component and in the solver wrapper, then the former value is used and a warning is printed.

If no intepolation is applied, the names of the ModelParts of both flow and structure solver need to be the same.

There is no parameter `timestep_start`, as currently, restart is not implemented in this solver wrapper.

### TubeFlow

This flow solver calculates the flow inside a 1D straight and flexible tube.
The required input is the radial displacement of the tube wall.
The resulting output is the pressure on the tube wall.

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
Add the start and end of the tube the boundary conditions have to be applied.
Discretizing these equations results in a system of linear algebraic equations which can be solved for the pressure and velocity
by performing Newton-Raphson iterations.
Inside the solver the kinematic pressure is used, which is the pressure divided by the density of the fluid.

For more information about the implementation this solver refer to .

#### Solver parameters

The following parameters, listed in alphabetical order, need to be specified in a file with name `input_file`.
This file should be located in the `working_directory`.
Care should be taken that the values of `d`, `e`, `h`, `l` and `rhof` match the corresponding values of the structural solver.

parameter|type|description
---:|:---:|---
`axial_offset`|double|(optional) Distance over which tube is displaced axially in the coordinate system. If not provided, the value of this parameter is 0.
`d`|double|Nominal diameter of the tube.
`e`|double|Modulus of elasticity of the tube wall.
`h`|double|Thickness of tube wall.
`inlet_boundary`|dict|Dictionary containing all information with respect to the inlet boundary condition.
`l`|double|Length of the tube.
`m`|int|Number of cells for discretization. The values are calculated in the cell centers.
`newtonmax`|double|Maximum number of Newton-Raphson iterations for the flow solver calculation.
`newtontol`|double|Relative tolerance for the Newton-Raphson iterations for the flow solver calculation.
<nobr>`outlet_boundary`</nobr>|dict|Dictionary containing all information with respect to the outlet boundary condition.
`preference`|double|(optional) Reference pressure and initial pressure. If not provided, the value of this parameter is 0.
`u0`|double|(optional) Initial velocity. If omitted, `ureference` is used.
`ureference`|double|Reference velocity used for determination of pressure stabilization term.
`rhof`|double|Density of the fluid.

##### Inlet Boundary

This section describes all parameters that need to be specified in the dictionary `inlet_boundary`, listed in alphabetical order.

parameter|type|description
---:|:---:|---
<nobr>`amplitude`</nobr>|double|Amplitude of the inlet boundary condition.
`period`|double|Period of the inlet boundary condition. Period of oscillation for a periodic boundary condition, duration for a non-periodic boundary condition. Not used for a fixed value boundary condition (type 4).
`reference`|double|(optional) Reference value of inlet boundary condition. If not provided, the value of this parameter is the corresponding reference value provided above, i.e. `ureference`, `preference` or `preference` + `ureference`^2.
`type`|int|Type of inlet boundary condition. <br>If 1, a sine wave is used with amplitude, reference and period as specified. <br>If 2, a pulse is used with amplitude as specified and a duration equal to the parameter period. After the pulse the variable is equal to the reference value. <br>If 3, a quadratic sine wave is used with amplitude, reference and period as specified. <br>If 4, a fixed value equal to the sum of the reference value and the amplitude. Used for steady cases. <br>If other, a steady increase of the value at the inlet with slope of amplitude divided by period is used.
`variable`|string|Variable upon which the inlet boundary condition is definded, either 'pressure', 'velocity' or 'total pressure'.

##### Outlet Boundary

This section describes all parameters that need to be specified in the dictionary `outlet_boundary`, listed in alphabetical order.

parameter|type|description
---:|:---:|---
`type`|int|Type of outlet boundary condition. <br>If 1, a non-reflecting boundary condition is applied. This type can not be used for a steady calculation. <br>If other, a fixed value equal to the reference pressure is applied.

### TubeRingModel

This structure solver calculates the deformation of the wall of a straight and flexible tube.
The tube is regarded as made up of independent rings.
The required input is the pressure on the tube wall.
The resulting output is the radial displacement.

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

For more information about the implementation this solver refer to .

#### Solver parameters

The following parameters, listed in alphabetical order, need to be specified in a file with name `input_file`.
This file should be located in the `working_directory`.
Care should be taken that the values of `d`, `e`, `h`, `l` and `rhof` match the corresponding values of the flow solver.

parameter|type|description
---:|:---:|---
<nobr>`axial_offset`</nobr>|double|(optional) Distance over which tube is displaced axially in the coordinate system. If not provided, the value of this parameter is 0.
`d`|double|Nominal diameter of the tube.
`e`|double|Modulus of elasticity of the tube wall.
`h`|double|Thickness of tube wall.
`l`|double|Length of the tube.
`m`|int|Number of cells for discretization. The values are calculated in the cell centers.
`preference`|double|(optional) Reference pressure. If not provided, the value of this parameter is 0.
`rhof`|double|Density of the fluid.

### TubeStructure

This structure solver calculates the deformation of the wall of a straight and flexible tube.
In this model inertia is taken into account, but still only radial displacements are considered.
The required input is the pressure on the tube wall.
The resulting output is the radial displacement.

The deformation of the tube in the radial direction is determined by
$$
	\rho_s h\frac{\partial^2 r}{\partial t^2}+b_1\frac{\partial^4 r}{\partial z^4}-b_2\frac{\partial^2 r}{\partial z^2}+b_3(r-r_o)=p-p_o
$$
with $\rho_s$ the solid density and $h$ the thickness of the tube's wall.
Further, $r$ is the inner radius, $p$ pressure and $t$ time. The variables $r_0$ and $p_0$ are a reference radius and pressure, respectively.
The tube is considered clamped at both ends.
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
Discretizing these equations results in a system of linear algebraic equations whit a Jacobian that does not depend on the pressure.
The system can be solved for the radius in one step.

For more information about the implementation this solver refer to .

#### Solver parameters

The following parameters, listed in alphabetical order, need to be specified in a file with name `input_file`.This file should be located in the `working_directory`.
Care should be taken that the values of `d`, `e`, `h`, `l` and `rhof` match the corresponding values of the flow solver.

parameter|type|description
---:|:---:|---
`axial_offset`|double|(optional) Distance over which tube is displaced axially in the coordinate system. If not provided, the value of this parameter is 0.
`beta`|double|(optional) Newmark parameter $\beta$.  Only required when the Newmark-beta time discretization is used.
`d`|double|Nominal diameter of the tube.
`e`|double|Modulus of elasticity of the tube wall.
`h`|double|Thickness of tube wall.
`gamma`|double|(optional) Newmark parameter $\gamma$. Only required when the Newmark-beta time discretization is used.
`l`|double|Length of the tube.
`m`|int|Number of cells for discretization. The values are calculated in the cell centers.
`nu`|double|Poisson's ratio.
`preference`|double|(optional) Reference pressure. If not provided, the value of this parameter is 0.
`rhof`|double|Density of the fluid.
`rhos`|double|Density of the tube wall.
<nobr>`time_disretization`</nobr>|string|(optional) Either 'Newmark' or 'backward Euler'. If not provided, backward Euler is used.

The equations are discretized in time with the Newmark-beta method or backward Euler method.
In the former case, two Newmark parameters $\beta$ and $\gamma$ are required, which result in an un-conditionally stable integration scheme if
$$
\gamma\geq\frac{1}{2} 
\textrm{ and }
\beta\geq\frac{1}{4}\left(\frac{1}{2}+\gamma\right)^2.
$$
Typical values are $\gamma$ equal to 1/2 and $\beta$ equal to 1/4.

A different time discretization for flow and structure can lead to difficulties for strongly coupled problems, especially looking at the resulting pressure distributions.
As most flow solvers are discretized using the backward Euler method, it is advised to chose the same method for the structure solver.