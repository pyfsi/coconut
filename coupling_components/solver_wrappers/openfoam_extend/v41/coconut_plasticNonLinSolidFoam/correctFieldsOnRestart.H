// On restart, the values on symmetryPlanes of surface fields are all zero for
// some reason. So we will reset them here from the equivalent volume fields
// Also, at the start of a new pass, fields such as F are set to zero in the new
// roller; these cells must be reset to I.

// Check if case has been restart: if the F field exists
const bool restarted
(
    IOobject
    (
        "F",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    ).headerOk()
);

if (restarted)
{
    Info<< "Correcting fields on restart" << endl;

    volScalarField newRho("newRho", mechanical.rho());
    volScalarField& rhoO = rho.oldTime();
    volScalarField& rhoOO = rhoO.oldTime();
    volScalarField& rhoOOO = rhoOO.oldTime();

    int nResetCells = 0;
    forAll(rhoO, cellI)
    {
        if (rhoO[cellI] < SMALL)
        {
            rho[cellI] = newRho[cellI];
            rhoO[cellI] = newRho[cellI];
            rhoOO[cellI] = newRho[cellI];
            rhoOOO[cellI] = newRho[cellI];
            nResetCells++;
        }
    }

    int nResetFaces = 0;
    forAll(rhoO.boundaryField(), patchI)
    {
        forAll(rhoO.boundaryField()[patchI], faceI)
        {
            if (rhoO.boundaryField()[patchI][faceI] < SMALL)
            {
                scalar nRho = newRho.boundaryField()[patchI][faceI];
                rho.boundaryField()[patchI][faceI] = nRho;
                rhoO.boundaryField()[patchI][faceI] = nRho;
                rhoOO.boundaryField()[patchI][faceI] = nRho;
                rhoOOO.boundaryField()[patchI][faceI] = nRho;

                nResetFaces++;
            }
        }
    }

    rho.correctBoundaryConditions();
    rhoO.correctBoundaryConditions();
    rhoOO.correctBoundaryConditions();
    rhoOOO.correctBoundaryConditions();

    Info<< "Resetting rho" << nl
        << "    cells: " << nResetCells << nl
        << "    boundary faces: " << nResetFaces << endl;

    // Check if F or J are zero anywhere
    tensorField& FI = F.internalField();
    scalarField& JI = J.internalField();

    nResetCells = 0;

    forAll(FI, cellI)
    {
        if (JI[cellI] < SMALL)
        {
            FI[cellI] = I;
            JI[cellI] = 1.0;

            nResetCells++;
        }
    }

    // Check all boundary values

    nResetFaces = 0;

    forAll(F.boundaryField(), patchI)
    {
        if
        (
            !F.boundaryField()[patchI].coupled()
         && mesh.boundaryMesh()[patchI].type() != "empty"
        )
        {
            tensorField& pF = F.boundaryField()[patchI];
            scalarField& pJ = J.boundaryField()[patchI];

            forAll(pF, faceI)
            {
                if (pJ[faceI] < SMALL)
                {
                    pF[faceI] = I;
                    pJ[faceI] = 1.0;

                    nResetFaces++;
                }
            }
        }
    }

    F.correctBoundaryConditions();
    J.correctBoundaryConditions();

    Info<< "Resetting F and J" << nl
        << "    cells: " << returnReduce(nResetCells, sumOp<int>()) << nl
        << "    boundary faces: " << returnReduce(nResetFaces, sumOp<int>())
        << endl;
}
