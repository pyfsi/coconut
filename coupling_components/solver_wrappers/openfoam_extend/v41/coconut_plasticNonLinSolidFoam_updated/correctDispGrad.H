forAll(DU.boundaryField(), patchI)
{
    // if
    // (
    //     isA<solidContactImplicitFvPatchVectorField>
    //     (
    //         DU.boundaryField()[patchI]
    //     )
    // )
    // {
    //     const vectorField patchN = mesh.boundary()[patchI].nf();

    //     const solidContactImplicitFvPatchVectorField& miDU =
    //         refCast<const solidContactImplicitFvPatchVectorField>
    //         (
    //             DU.boundaryField()[patchI]
    //         );

    //     // gradDU.boundaryField()[patchI] -=
    //     //     (sqr(n) & gradDU.boundaryField()[patchI]);

    //     gradDU.boundaryField()[patchI] +=
    //          patchN*miDU.localSnGrad()
    //        - (sqr(patchN) & gradDU.boundaryField()[patchI]);

    //     //Info<< "sol: patch " << patchI
    //     //     << ", gradDU " << (n & gradDU.boundaryField()[patchI]) << endl;

    //     // gradDUf.boundaryField()[patchI] +=
    //     //     n*miDU.localSnGrad()
    //     //   - (sqr(n) & gradDUf.boundaryField()[patchI]);
    // }
    if
    (
        isA<materialGgiFvPatchVectorField>
        (
            DU.boundaryField()[patchI]
        )
    )
    {
        const vectorField patchN = mesh.boundary()[patchI].nf();

        const materialGgiFvPatchVectorField& miDU =
            refCast<const materialGgiFvPatchVectorField>
            (
                DU.boundaryField()[patchI]
            );

        // gradDU.boundaryField()[patchI] -=
        //     (sqr(n) & gradDU.boundaryField()[patchI]);

        // First extrapolate gradient from internal cells
        gradDU.boundaryField()[patchI] =
            gradDU.boundaryField()[patchI].patchInternalField();

        // Now correct normal gradient
        gradDU.boundaryField()[patchI] +=
             patchN*miDU.localSnGrad()
           - (sqr(patchN) & gradDU.boundaryField()[patchI]);

        //Info<< "sol: patch " << patchI
        //     << ", gradDU " << (n & gradDU.boundaryField()[patchI]) << endl;

        // gradDUf.boundaryField()[patchI] +=
        //     patchN*miDU.localSnGrad()
        //   - (sqr(patchN) & gradDUf.boundaryField()[patchI]);
    }
}
