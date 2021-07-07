{
    // scalarField magDU = mag(DU.internalField());
    // scalarField magDU = mag(U.internalField() - U.oldTime().internalField());

    // forAll(magDU, cellI)
    // {
    //     if (magDU[cellI] < SMALL)
    //     {
    //         magDU[cellI] = SMALL;
    //     }
    // }

    relativeResidual =
        gMax
        (
            mag
            (
                DU.internalField()
              - DU.prevIter().internalField()
              //   U.internalField()
              // - U.prevIter().internalField()
            )/(gMax(mag(DU.internalField())) + SMALL)
           // /magDU
        );

    // Info<< "Writing relResField" << endl;
    // relResField = (DU - DU.prevIter())/gMax(magDU);
    // relResField.write();
}
