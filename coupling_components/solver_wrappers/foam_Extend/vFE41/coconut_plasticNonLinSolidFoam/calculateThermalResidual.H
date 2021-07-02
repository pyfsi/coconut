{
    // scalarField magT = mag(T.internalField());

    // forAll(magT, cellI)
    // {
    //     if (magT[cellI] < SMALL)
    //     {
    //         magT[cellI] = SMALL;
    //     }
    // }

    scalar maxMagT = gMax(mag(T.internalField()));

    const volScalarField& TprevIter = *TprevIterPtr;

    thermalResidual =
        // gMax(mag(T.internalField() - T.prevIter().internalField()))/magT;
        gMax(mag(T.internalField() - TprevIter.internalField()))/maxMagT;
}
