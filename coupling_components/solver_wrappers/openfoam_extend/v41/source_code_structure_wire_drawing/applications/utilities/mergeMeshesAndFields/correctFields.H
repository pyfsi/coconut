    // IOobject Jheader
    // (
    //     "J",
    //     masterMesh.time().timeName(),
    //     masterMesh,
    //     IOobject::MUST_READ,
    //     IOobject::AUTO_WRITE
    // );

    // IOobject Jfheader
    // (
    //     "Jf",
    //     masterMesh.time().timeName(),
    //     masterMesh,
    //     IOobject::MUST_READ,
    //     IOobject::AUTO_WRITE
    // );

    IOobject Fheader
    (
        "F",
        masterMesh.time().timeName(),
        masterMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    IOobject Ffheader
    (
        "Ff",
        masterMesh.time().timeName(),
        masterMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    Warning
        << "Fields for new materialModels may have to be correct here!!!"
        << nl << "Or corrected within the constructor or the model"
        << endl;

    // IOobject bEheader
    // (
    //     "bE",
    //     masterMesh.time().timeName(),
    //     masterMesh,
    //     IOobject::MUST_READ,
    //     IOobject::AUTO_WRITE
    // );

    // IOobject bEfheader
    // (
    //     "bEf",
    //     masterMesh.time().timeName(),
    //     masterMesh,
    //     IOobject::MUST_READ,
    //     IOobject::AUTO_WRITE
    // );

    // IOobject bEbarheader
    // (
    //     "bEbar",
    //     masterMesh.time().timeName(),
    //     masterMesh,
    //     IOobject::MUST_READ,
    //     IOobject::AUTO_WRITE
    // );

    // IOobject bEbarfheader
    // (
    //     "bEbarf",
    //     masterMesh.time().timeName(),
    //     masterMesh,
    //     IOobject::MUST_READ,
    //     IOobject::AUTO_WRITE
    // );

    // vol fields

    // No longer needed as solver checked this in createFields and also
    // mechanicalLaws will include checks

    // if
    // (
    //     Fheader.headerOk()
    //     //&& bEheader.headerOk()
    //     //&& bEbarheader.headerOk()
    // )
    // {
    //     // Info<< "Correcting F, bE, and bEbar for roller"
    //     Info<< "Correcting F for roller"
    //         << " cells/faces" << endl;

    //     volTensorField F(Fheader, mesh);
    //     //volSymmTensorField bE(bEheader, mesh);
    //     //volSymmTensorField bEbar(bEbarheader, mesh);

    //     volScalarField J = det(F);

    //     forAll(J.internalField(), cellI)
    //     {
    //         if (J.internalField()[cellI] < SMALL)
    //         {
    //             F.internalField()[cellI] = 1.0*I;
    //             // bE.internalField()[cellI] = 1.0*I;
    //             // bEbar.internalField()[cellI] = 1.0*I;
    //         }
    //     }

    //     forAll(J.boundaryField(), patchI)
    //     {
    //         forAll(J.boundaryField()[patchI], faceI)
    //         {
    //             if (J.boundaryField()[patchI][faceI] < SMALL)
    //             {
    //                 F.boundaryField()[patchI][faceI] = 1.0*I;
    //                 // bE.boundaryField()[patchI][faceI] = 1.0;
    //                 // bEbar.boundaryField()[patchI][faceI] = 1.0*I;
    //             }
    //         }
    //     }

    //     F.write();
    //     // bE.write();
    //     // bEbar.write();
    // }

    // // surface fields

    // if
    // (
    //     Ffheader.headerOk()
    //  // && bEfheader.headerOk()
    //  // && bEbarfheader.headerOk()
    // )
    // {
    //     // Info<< "Correcting Ff, bEf, and bEbarf for roller"
    //     Info<< "Correcting Ff for roller"
    //         << " cells/faces" << endl;

    //     surfaceTensorField Ff(Ffheader, mesh);
    //     // surfaceSymmTensorField bEf(bEfheader, mesh);
    //     // surfaceSymmTensorField bEbarf(bEbarfheader, mesh);

    //     surfaceScalarField Jf = det(Ff);

    //     forAll(Jf.internalField(), faceI)
    //     {
    //         if (Jf.internalField()[faceI] < SMALL)
    //         {
    //             Ff.internalField()[faceI] = 1.0*I;
    //             // bEf.internalField()[faceI] = 1.0;
    //             // bEbarf.internalField()[faceI] = 1.0*I;
    //         }
    //     }

    //     forAll(Jf.boundaryField(), patchI)
    //     {
    //         forAll(Jf.boundaryField()[patchI], faceI)
    //         {
    //             if (Jf.boundaryField()[patchI][faceI] < SMALL)
    //             {
    //                 Ff.boundaryField()[patchI][faceI] = 1.0*I;
    //                 // bEf.boundaryField()[patchI][faceI] = 1.0;
    //                 // bEbarf.boundaryField()[patchI][faceI] = 1.0*I;
    //             }
    //         }
    //     }

    //     Ff.write();
    //     // bEf.write();
    //     // bEbarf.write();
    // }
