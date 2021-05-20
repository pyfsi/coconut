#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "fixedValuePointPatchField.H"
#include "fvOptions.H"
#include "dimensionedScalar.H"

void ApplyFSIPointDisplacement(dynamicFvMesh& mesh, word boundary_name, word prev_runTime )
    {

        label patchWallID = mesh.boundaryMesh().findPatchID(boundary_name);

        // *** Set patch displacement for motion solver.*** //
        // Find the reference to the pointDisplacement field
        pointVectorField& pointDisp = const_cast<pointVectorField&>
        (
            mesh.objectRegistry::lookupObject<pointVectorField>
            (
            "pointDisplacement"
            )
        );

        Info<< "Reading pointDisplacement\n" << endl;

        pointVectorField pointDispNext_
        (
            IOobject
            (
                "pointDisplacement_Next",
                prev_runTime,
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            pointMesh::New(mesh)
        );

        pointVectorField& pointDispNext = const_cast<pointVectorField&>(pointDispNext_);

        vectorField& patchpointDispNext = refCast<vectorField>(pointDispNext.boundaryFieldRef()[patchWallID]);

        // Assign the new boundary displacements
        pointDisp.boundaryFieldRef()[patchWallID] ==  patchpointDispNext;

    }
