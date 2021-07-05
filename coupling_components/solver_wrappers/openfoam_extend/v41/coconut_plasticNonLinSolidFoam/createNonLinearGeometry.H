// Add nonlinear geometry entry to solution dict
// This can then be accessed by boundary conditions, mechanical model, etc to
// know what type of solver is being used.
{
    Info<< "Setting nonLinear to updatedLagrangianKirchhoff" << nl << endl;
    mesh.solutionDict().subDict("solidMechanics").set
    (
        "nonLinear",
        nonLinearGeometry::nonLinearNames_
        [
            nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF
        ]
    );
}
