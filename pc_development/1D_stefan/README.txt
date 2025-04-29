WARNING: the cases in the 1D_stefan directory will not work out of the box.
They need to be adapted to the new energy source term in the solid solver and the new face-to-node mapper.
At least the following changes will need to be made:
- Correct settings for the face-to-node mapper in parameters.json file
- Changes to the case.jou files in the steup files for the solid solver
    -> Disable Fluent solidification-melting model and flow equations
    -> Only the energy equation is solved with the new Darcy source term for temperature

EXCEPT: Stefan_subcooled(_expl) are two new cases compatible with the latest solver wrappers.