#include "fixedValuePointPatchField.H"

void ApplyFSIPointDisplacement(dynamicFvMesh& mesh, word boundary_name)
{
    // Read input pointDisplacement
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

    // Get the displacement on the patch
    Foam::pointVectorField* pointDisplacement_;
    pointDisplacement_ = const_cast<pointVectorField*>(&mesh.lookupObject<pointVectorField>("pointDisplacement"));

    pointVectorField& pointDispTmp = const_cast<pointVectorField&>(pointDispTmp_);

    label patchWallID = mesh.boundaryMesh().findPatchID(boundary_name);
    vectorField& patchpointDispTmp = refCast<vectorField>(pointDispTmp.boundaryFieldRef()[patchWallID]);

    // Overwrite the nodes on the interface directly
    pointDisplacement_->boundaryFieldRef()[patchWallID] == patchpointDispTmp;
}
