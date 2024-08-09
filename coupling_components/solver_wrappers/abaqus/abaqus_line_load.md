# Abaqus

This is the preliminary documentation for the *`abaqus_line_load`* solver wrapper. It is now possible to apply fluid forces on a 1D beam structure.
Only the traction variable is used (pressure does not make sense).

## Changes with respect to the standard Abaqus solver wrapper

The solver wrapper itself is in principle the same. This can thus easily be merged with the standard Abaqus solver wrapper (in a later stage).

The DLOAD-subroutine is changed in the following way: add a new integer variable `AXIS`, that is related to `JLTYP`.
The latter variable is an input from Abaqus and denotes which load type calls the DLOAD subroutine. This ranges between 41 and 43 for the PXNU, PYNU and PZNU line loads.

Most importantly, since the line loads are given on element sets rather than surfaces (surface line loads are not possible), the variable `SNAME` is blank in the DLOAD subroutine.
This implies that there can be only ONE model part for the line loads. This is therefore hardcoded in the solver wrapper.

Future work, for e.g. three bladed wind turbines: write some sort of 'combine'-mapper, that has as function to merge three ALM model parts in the flow solver 
(you usually want to separate these per blade) to one beam elements model part in the structural solver (as it is required to have only one model part here).

The UTRACLOAD subroutine is removed as it is unnecessary.

## Notes to include *`abaqus_line_load`* in *`abaqus`*

In the USRInit.f, write the faces file only if `JLTYP` equals 0 (pressure via the `*dsload` key in the input file, as is the case now) or 41 (PXNU via `*dload` for line loads in the x-direction).
Reconsider pressure in all cases to keep the same formatting of the *`CSM_TimeBSurfaceACpu0Input.dat`*-file, but set it to zero. So set `AXIS = JLTYP - 39` instead.

Add a new keyword to the parameter file so that the user can indicate whether the model is a beam elements model (and takes line loads) or a general surface based model.
How to check the correct setting of this parameter?