{
    scalarField magT = max(mag(T.internalField()), SMALL);
    //scalar maxMagT = gMax(mag(T.internalField()));

    relativeResidualT =
        gMax(mag(T.internalField() - T.prevIter().internalField())/magT);
    //gMax(mag(T.internalField() - T.prevIter().internalField()))/maxMagT;
}
