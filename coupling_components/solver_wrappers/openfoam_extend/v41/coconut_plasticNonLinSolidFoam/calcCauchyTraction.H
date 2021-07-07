{
    // Update current area vectors surface field
    curS = (mesh.Sf() & fvc::interpolate(relFinv*relJ));

    // Update Cauchy traction surface field
    cauchyTraction = (curS & fvc::interpolate(tau/J))/mag(curS);

    // Update gradient of displacement increment surface field
    gradDUf = fvc::interpolate(gradDU);

    // Correct Cauchy traction at materialGgi patches
    forAll(DU.boundaryField(), patchI)
    {
        // If the patch is a materialGgi then update curS,
        // cauchyTraction and gradDUf
        // if (materialGgi || regionCoupled) Commented out by ZT (parallel run)
        if
        (
            isA<materialGgiFvPatchVectorField>(DU.boundaryField()[patchI])
         || isA<materialCouplingFvPatchVectorField>(DU.boundaryField()[patchI])
        )
        {
            curS.boundaryField()[patchI] =
            (
                mesh.Sf().boundaryField()[patchI]
              & (
                    relFinv.boundaryField()[patchI]
                   *relJ.boundaryField()[patchI]
                )
            );

            cauchyTraction.boundaryField()[patchI] =
            (
                curS.boundaryField()[patchI]
              & (
                    tau.boundaryField()[patchI]
                   /J.boundaryField()[patchI]
                )
            )/mag(curS.boundaryField()[patchI]);

            const vectorField patchN = mesh.boundary()[patchI].nf();

            gradDUf.boundaryField()[patchI] +=
                patchN*DU.boundaryField()[patchI].snGrad()
              - (sqr(patchN) & gradDUf.boundaryField()[patchI]);
        }
    }
}
