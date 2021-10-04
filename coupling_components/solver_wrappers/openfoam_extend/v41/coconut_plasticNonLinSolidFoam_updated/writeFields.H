if (runTime.outputTime())
{
    // Equivalent strain
    // volScalarField epsilonEq
    // (
    //     IOobject
    //     (
    //         "epsilonEq",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     sqrt((2.0/3.0)*magSqr(dev(epsilon)))
    // );

    // Info<< "Max epsilonEq = " << max(epsilonEq).value()
    //     << endl;


    // Cauchy von Mises stress
    volScalarField sigmaCauchyEq
    (
        IOobject
        (
            "sigmaCauchyEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigmaCauchy)))
    );

    Info<< "Max sigmaCauchyEq = " << gMax(sigmaCauchyEq) << endl;


    // True/log/Hencky/natural strain tensor
    volSymmTensorField epsilonTrue
    (
        IOobject
        (
            "epsilonTrue",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5*log(symm(F.T() & F))
    );

    volScalarField epsilonTrueEq
    (
        IOobject
        (
            "epsilonTrueEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((2.0/3.0)*magSqr(dev(epsilonTrue)))
    );

    Info<< "Max epsilonTrueEq = " << gMax(epsilonTrueEq) << nl << endl;

    // Write out friction coefficient from the solidContact boundary conditions
    // This can be stored as a field within the frictionModel: to-do
    // volScalarField friCoeff
    // (
    //     IOobject
    //     (
    //         "frictionCoeff",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar("zero", dimless, 0)
    // );

    // forAll(DU.boundaryField(), patchI)
    // {
    //     if
    //     (
    //         DU.boundaryField()[patchI].type()
    //      == solidContactFvPatchVectorField::typeName
    //     )
    //     {
    //         const solidContactFvPatchVectorField& contactDU =
    //             refCast<const solidContactFvPatchVectorField>
    //             (
    //                 DU.boundaryField()[patchI]
    //             );

    //         const unallocLabelList& faceCells = contactDU.patch().faceCells();

    //         vectorField n = mesh.boundary()[patchI].nf();

    //         if (!contactDU.master())
    //         {
    //             forAll(contactDU.contact(),faceI)
    //             {
    //                 scalar contactP =
    //                     mag(contactDU.traction()[faceI] & n[faceI]);
    //                 scalar contactTau =
    //                     mag((I-sqr(n[faceI])) & contactDU.traction()[faceI]);
    //                 if (contactDU.contact()[faceI])
    //                 {
    //                     friCoeff[faceCells[faceI]] = contactTau/contactP;
    //                     friCoeff.boundaryField()[patchI][faceI] =
    //                         contactTau/contactP;
    //                 }
    //             }
    //         }
    //     }
    // }


    // Write out contact field from the thermalContact boundary conditions
    // if (false)
    // {
    //     volScalarField contact
    //         (
    //             IOobject
    //             (
    //                 "contact",
    //                 runTime.timeName(),
    //                 mesh,
    //                 IOobject::NO_READ,
    //                 IOobject::AUTO_WRITE
    //             ),
    //             mesh,
    //             dimensionedScalar("zero", dimless, -1)
    //         );

    //     volScalarField& T = *TPtr;

    //     forAll(T.boundaryField(), patchI)
    //     {
    //         if
    //             (
    //                 T.boundaryField()[patchI].type()
    //                 == thermalContactFvPatchScalarField::typeName
    //             )
    //         {
    //             const thermalContactFvPatchScalarField& contactT =
    //                 refCast<const thermalContactFvPatchScalarField>
    //                 (
    //                     T.boundaryField()[patchI]
    //                 );

    //             const unallocLabelList& faceCells =
    //                 contactT.patch().faceCells();

    //             forAll(contactT.contact(), faceI)
    //             {
    //                 contact.boundaryField()[patchI][faceI] =
    //                     contactT.contact()[faceI];

    //                 contact[faceCells[faceI]] = contactT.contact()[faceI];
    //             }
    //         }
    //         else if
    //             (
    //                 T.boundaryField()[patchI].type()
    //                 == thermalGeneralContactFvPatchScalarField::typeName
    //             )
    //         {
    //             const thermalGeneralContactFvPatchScalarField& contactT =
    //                 refCast<const thermalGeneralContactFvPatchScalarField>
    //                 (
    //                     T.boundaryField()[patchI]
    //                 );

    //             const unallocLabelList& faceCells =
    //                 contactT.patch().faceCells();

    //             forAll(contactT.contact(), faceI)
    //             {
    //                 const label faceCellI = faceCells[faceI];

    //                 if (contact[faceCellI] < 1.0)
    //                 {
    //                     contact[faceCellI] = contactT.contact()[faceI];

    //                     contact.boundaryField()[patchI][faceI] =
    //                         contactT.contact()[faceI];
    //                 }
    //             }
    //         }
    //     }
    // }

    // Write out Rhie-Chie fields
    // const volVectorField RhieChowDU
    // (
    //     "RhieChowDU",
    //     RhieChowScaleFactor
    //    *(
    //        fvc::laplacian(twoMuLambdaf, DU, "laplacian(DDU,DU)")
    //      - fvc::div(twoMuLambdaf*(mesh.Sf() & gradDUf))
    //     )
    // );
    // const surfaceTensorField gradUf = fvc::interpolate(gradU);
    // const volVectorField RhieChowU
    // (
    //     "RhieChowU",
    //     RhieChowScaleFactor
    //    *(
    //        fvc::laplacian(twoMuLambdaf, U, "laplacian(DDU,DU)")
    //      - fvc::div(twoMuLambdaf*(mesh.Sf() & gradUf))
    //     )
    // );
    // RhieChowDU.write();
    // RhieChowU.write();

    runTime.write();
}
