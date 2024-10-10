# Conjugate heat transfer in Fluent

This is the documentation for the Fluent solver wrapper adapted for conjugate heat transfer simulations within CoCoNuT.
The functioning of this solver wrapper is entirely based on the original [Fluent solver wrapper](../fluent/fluent.md).


## Coupling strategy for conjugate heat transfer

Conjugate heat transfer (CHT) considers the simultaneous solution of heat transfer in adjacent solid and fluid domains [[1](#1)].
This solver wrapper allows to solve conjugate heat transfer problems using a partitioned approach in CoCoNuT.
Both the solid and fluid domains are solved separately and are coupled during coupling iterations by exchanging interface data until convergence.
Clearly, the coupling strategy for a partitioned fluid-structure interaction (FSI) problem is very similar to the coupling strategy used here for CHT problems.
Instead of exchanging forces and interface displacements, temperature and heat flux are exchanged.
As a result, it makes sense to use the existing CoCoNuT framework for the coupling of the solvers, and adapt the solver wrapper to handle the thermal variables.


## Limitations and untested conditions

* Only cases in 2D have been tested so far. Even though the code is written with axisymmetric or 3D cases in mind and no issues are expected, trouble-free operation is not guaranteed.
* The restart functionality in CoCoNuT has not been tested yet and is not guaranteed to work.
* Only temperature - heat flux exchanging coupling schemes for conjugate heat transfer are possible. Coupling schemes with virtual heat transfer coefficients, as in [[2](#2)], are not provided.
* FSI functionality of the solver wrapper is still intact, but combined CHT - FSI problems will need additional programming to run.


## Parameters

All parameters used in the orginal [Fluent solver wrapper](../fluent/fluent.md) remain available.
A new subdictionary with keyword `CHT` should be provided, however, containing the following keywords:

|         parameter | type  | description                                                                                                                                                                                             |
|------------------:|:-----:|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `moving_boundary` | bool  | (optional) Default: `false`. `false` in case of a CHT problem with fixed boundary, `true` if motion of the coupling interface is allowed, e.g. in case of a combined CHT - FSI problem.                 |
|   `ini_condition` | float | Scalar value as initial condition for the output thermal boundary condition: temperature (in K) or heat flux (in W/m$\cdot$K). Not used when `bc_from_case` is `true`.                                  |
|    `bc_from_case` | bool  | (optional) Default: `false`. If `true`, `ini_condition` value is not used and the existing temperature or heat flux profile in the provided cases is used as initial condition for the output variable. |


## Overview of operation

The solver wrapper consists of 3 files (with X the Fluent version, e.g. "v2023R1"):

-   *`X.py`*: defines the `SolverWrapperCHTFluentX` class, 
-   *`X.jou`*: Fluent journal file to interactively run the CHT simulation, written in Scheme, 
-   *`udf_thermal.c`*: Fluent UDF file that implements additional functionality (I/O) used in Fluent, written in C.


### Added functionality

The most important feature of the CHT Fluent solver wrapper compared to the original is the ability to exchange temperature and heat flux as interface variables between the solvers.
To that end, the necessary UDFs have been added.
*`store_temperature`* and *`store_heat_flux`* write the temperature and heat flux at the cell faces on the coupling interface, very similar to *`store_pressure_traction`* for pressure and traction.
Furthermore, *`set_temperature`* and *`set_heat_flux`* are *DEFINE_PROFILE* UDFs used to read the input files and set either the temperature or heat flux profiles as boundary conditions at the coupling interface.
To accommodate these new variables, the solver wrapper journal file and python script are adapted as well.
The python script is now suited for input variables stored at the faces, such as temperature and heat flux, and not only at the face nodes, such as displacement.
Furthermore, the thermal interface variables can be initialized with a constant value or profile.


### Files created during simulation

Only the files created for pure CHT problems are listed below, for FSI related problems, the reader is referred to the [Fluent solver wrapper](../fluent/fluent.md).
In these file conventions, A is the time step number and B the Fluent thread ID.

-   Fluent case and data files are saved as files of the form *`case_timestepA.cas`* and *`case_timestepA.dat`*.
-   Current node and face coordinates are passed from Fluent to CoCoNuT with files of the form *`nodes_timestepA_threadB.dat`* and *`faces_timestepA_threadB.dat`*. 
-   Temperature is passed from CoCoNuT to Fluent with files of the form *`temperature_timestepA_threadB.dat`*.
-   Heat flux is passed from Fluent to CoCoNuT with files of the form *`heat flux_timestepA_threadB.dat`*.
-   Files with extension *`.coco`* are used to exchange messages between CoCoNuT and Fluent.


## Setting up a new case

Following items should be set up and saved in the Fluent case file (this list may be non-exhaustive):

-   additional UDFs must be configured, 
-   steady/unsteady (should match with the `unsteady` parameter),
-   2D, 3D or axisymmetric (should match with the `dimensions` parameter),
-   initial interface conditions if a profile is required,
-   boundary conditions, material properties, numerical models, discretization schemes, operating conditions, turbulence modeling, convergence criteria.

A data file should also be present with the fields either initialized or containing the results of a previous calculation.

Following items are taken care of by CoCoNuT, and must therefore not be included in the Fluent case file:

-   initial interface conditions if a constant value is sufficient,
-   the time step (`delta_t`).


## Version specific documentation

### v2023R1 (23.1.0)

Base version.

### v2024R1 (24.1.0)

No changes.


## References

<a id="1">[1]</a>
[Van Riet V., Beyne W., De Paepe M., Degroote J., "Convergence behaviour of partitioned methods for conjugate heat transfer problems", in ECCOMAS 2024, Proceedings, Lisbon, Portugal, 2024.]

<a id="2">[2]</a>
[Verstraete T. and Scholl S., "Stability analysis of partitioned methods for predicting conjugate heat transfer", International Journal of Heat and Mass Transfer, vol. 101, pp. 852–869, 2016.](https://doi.org/10.1016/J.IJHEATMASSTRANSFER.2016.05.041)
