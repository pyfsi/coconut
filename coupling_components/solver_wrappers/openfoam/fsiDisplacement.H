#include "fixedValuePointPatchField.H"

void ApplyFSIPointDisplacement(dynamicFvMesh& mesh, word boundary_name)
{

    label patchWallID = mesh.boundaryMesh().findPatchID(boundary_name);

    // *** Set patch displacement for motion solver *** //
    // Find the reference to the pointDisplacement field
    pointVectorField& pointDisp = const_cast<pointVectorField&>
    (
        mesh.objectRegistry::lookupObject<pointVectorField>
        (
        "pointDisplacement"
        )
    );

    Info<< "Reading pointDisplacement\n" << endl;

    pointVectorField pointDispTmp_
    (
        IOobject
        (
            "pointDisplacementTmp",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        pointMesh::New(mesh)
    );

    pointVectorField& pointDispTmp = const_cast<pointVectorField&>(pointDispTmp_);

    vectorField& patchpointDispTmp = refCast<vectorField>(pointDispTmp.boundaryFieldRef()[patchWallID]);

    // Assign the new boundary displacements
    pointDisp.boundaryFieldRef()[patchWallID] == patchpointDispTmp;
}
