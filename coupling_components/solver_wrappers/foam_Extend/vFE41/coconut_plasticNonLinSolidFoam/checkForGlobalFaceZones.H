// GlobalFaceZones are no longer used for the contact, so we will throw an error
// if we find a globalFaceZone
if (!materialGgi && !regionCoupled)
{
    IOobject decomposeParDictHeader
    (
        "decomposeParDict",
        runTime.system(),
        runTime,
        IOobject::MUST_READ
    );

    if (decomposeParDictHeader.headerOk())
    {
        // Read the dict
        IOdictionary decomposeParDict
        (
            IOobject
            (
                "decomposeParDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        if (decomposeParDict.found("globalFaceZones"))
        {
            FatalError
                << "globalFaceZones should no be specified in the "
                << "decomposeParDict! solidContact not longer uses "
                << "globalFaceZones."
                << abort(FatalError);
        }
    }
}
