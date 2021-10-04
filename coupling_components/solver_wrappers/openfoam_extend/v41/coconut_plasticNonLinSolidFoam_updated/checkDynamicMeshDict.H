// Check if the dynamicMeshDict is found, if not, then we will add it setting
// the mesh to be a solidDynamicFvMesh
{
    IOobject dynamicMeshDictHeader
    (
        "dynamicMeshDict",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ
    );

    if (!dynamicMeshDictHeader.headerOk())
    {
        WarningIn("checkDynamicMeshDict.H")
            << "dynamicMeshDict not found!" << nl
            << "    Adding a dynamicMeshDict and setting dynamicFvMesh to "
            << "solidDynamicFvMesh" << endl;

        // Create the dict
        IOdictionary dynamicMeshDict
        (
            IOobject
            (
                "dynamicMeshDict",
                runTime.constant(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );

        dynamicMeshDict.add("dynamicFvMesh", "solidDynamicFvMesh");

        dynamicMeshDict.regIOobject::write();
    }
    else
    {
        // Read the dict
        IOdictionary dynamicMeshDict
        (
            IOobject
            (
                "dynamicMeshDict",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        if (word(dynamicMeshDict.lookup("dynamicFvMesh")) == "staticFvMesh")
        {
            FatalError
                << "dynamicFvMesh should never be staticFvMesh!"
                << abort(FatalError);
        }
    }
}
