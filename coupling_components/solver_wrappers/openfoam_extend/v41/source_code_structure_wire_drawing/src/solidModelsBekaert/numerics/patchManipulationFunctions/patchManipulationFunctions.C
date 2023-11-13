/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "patchManipulationFunctions.H"
#include "directTopoChange.H"
#include "polyModifyFace.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// namespace patchManipulationFunctions
// {

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void patchManipulationFunctions::reorderPatchesConsistently
(
    const polyMesh& meshSource,
    polyMesh& meshTarget
)
{
    // Check both meshes have the same number of patches, excluding processor
    // patches

    // Count number of processor boundaries on each mesh
    int nMeshSourceProcPatches = 0;
    boolList meshSourceProcPatches(meshSource.boundaryMesh().size(), false);
    forAll(meshSource.boundaryMesh(), patchI)
    {
        if
        (
            meshSource.boundaryMesh()[patchI].type()
         == processorPolyPatch::typeName
        )
        {
            nMeshSourceProcPatches++;
            meshSourceProcPatches[patchI] = true;
        }
    }

    int nMeshTargetProcPatches = 0;
    boolList meshTargetProcPatches(meshTarget.boundaryMesh().size(), false);
    forAll(meshTarget.boundaryMesh(), patchI)
    {
        if
        (
            meshTarget.boundaryMesh()[patchI].type()
         == processorPolyPatch::typeName
        )
        {
            nMeshTargetProcPatches++;
            meshTargetProcPatches[patchI] = true;
        }
    }

    if
    (
        (meshSource.boundaryMesh().size() - nMeshSourceProcPatches)
     != (meshTarget.boundaryMesh().size() - nMeshTargetProcPatches)
    )
    {
        FatalErrorIn
        (
            "void reorderPatchesConsistently"
            "(polyMesh& meshSource, polyMesh& meshTarget)"
        )   << "The source mesh should have the same number of patches "
            << "(excluding processor patches) as the target mesh!"
            << abort(FatalError);
    }

    // Store the source mesh patch names and the map from the target patches
    // to the source patches
    wordList sourcePatchNames(meshSource.boundaryMesh().size());
    labelList targetToSourcePatchMap(meshSource.boundaryMesh().size(), -1);
    forAll(sourcePatchNames, patchI)
    {
        sourcePatchNames[patchI] = meshSource.boundaryMesh()[patchI].name();

        if (!meshSourceProcPatches[patchI])
        {
            // Find the patch in the target mesh with the same name

            const label targetPatchID =
                meshTarget.boundaryMesh().findPatchID(sourcePatchNames[patchI]);

            if (targetPatchID == -1)
            {
                FatalError
                    << "patch " << sourcePatchNames[patchI]
                    << " exists in the source "
                    << "mesh but does not exist in the target mesh!"
                    << abort(FatalError);
            }

            targetToSourcePatchMap[targetPatchID] = patchI;
        }
    }


    // Take a reference to the boundary patches
    const polyBoundaryMesh& targetPatches = meshTarget.boundaryMesh();
    const polyBoundaryMesh& sourcePatches = meshSource.boundaryMesh();


    // Repatch faces

    // Create mesh modifier
    directTopoChange meshMod(meshTarget);

    // Move the faces to the correct patches
    // Here we are making sure all the boundary faces are in the correct order,
    // and we are not worrying about the patch names
    forAll(targetPatches, patchI)
    {
        if (!meshTargetProcPatches[patchI])
        {
            const polyPatch& pp = targetPatches[patchI];
            const label targetPatchID = targetToSourcePatchMap[patchI];

            forAll(pp, faceI)
            {
                const label faceID = pp.start() + faceI;

                const label zoneID = meshTarget.faceZones().whichZone(faceID);

                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = meshTarget.faceZones()[zoneID];

                    zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
                }

                // Move face to the correct patch
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        meshTarget.faces()[faceID],        // face
                        faceID,                            // face ID,
                        meshTarget.faceOwner()[faceID],    // owner
                        -1,                                // neighbour,
                        false,                             // flip flux
                        targetPatchID,                     // patchID,
                        false,                             // remove from zone
                        zoneID,                            // zone ID
                        zoneFlip                           // zone flip
                    )
                );
            }
        }
    }

    // Perform the topo changes on the target mesh to reorder the faces
    meshMod.changeMesh
    (
        meshTarget,
        false    // inflate points
    );

    // Re-create patches to ensure they have the correct type
    DynamicList<polyPatch*> allPatches(meshTarget.boundaryMesh().size());
    forAll(meshTarget.boundaryMesh(), patchI)
    {
        const polyPatch& pp = meshTarget.boundaryMesh()[patchI];

        if (!meshTargetProcPatches[patchI])
        {
            // Create patch dictionary and use sourcePatch type and name
            dictionary patchDict;
            patchDict.set("nFaces", pp.size());
            patchDict.set("startFace", pp.start());

            patchDict.set("type", sourcePatches[patchI].type());

            // Add patch
            allPatches.append
            (
                polyPatch::New
                (
                    sourcePatches[patchI].name(),
                    patchDict,
                    patchI,
                    meshTarget.boundaryMesh()
                ).ptr()
            );
        }
        else
        {
            // clone processor patches
            allPatches.append(pp.clone(meshTarget.boundaryMesh()).ptr());
        }
    }

    allPatches.shrink();
    meshTarget.removeBoundary();
    meshTarget.addPatches(allPatches);
}


void patchManipulationFunctions::removePatchesWithNoFaces(polyMesh& mesh)
{
    // Find patches with no faces
    DynamicList<word> patchesWithNoFaces(mesh.boundaryMesh().size());
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (mesh.boundaryMesh()[patchI].size() == 0)
        {
            patchesWithNoFaces.append(mesh.boundaryMesh()[patchI].name());
        }
    }

    patchesWithNoFaces.shrink();

    // Remove all patches with no faces
    forAll(patchesWithNoFaces, pI)
    {
        removePatchWithNoFaces(mesh, patchesWithNoFaces[pI]);
    }
}


void patchManipulationFunctions::removePatchWithNoFaces
(
    polyMesh& mesh,
    const word& patchName
)
{
    // Find the patchID
    const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    if (patchID == -1)
    {
        FatalErrorIn
        (
            "void removePatchWithNoFaces(polyMesh& mesh, const word& patchName)"
        )   << "Cannot find patch " << patchName << " in the mesh!"
            << abort(FatalError);
    }

    if (mesh.boundaryMesh()[patchID].size() > 0)
    {
        WarningIn
        (
            "void removePatchWithNoFaces(polyMesh& mesh, const word& patchName)"
        )   << patchName << " is not removed as it contains faces!"
            << endl;

        return;
    }

    // Make a copy of the patches, without the patch to be removed
    DynamicList<polyPatch*> allPatches(mesh.boundaryMesh().size() - 1);
    label newPatchID = 0;
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (patchI != patchID)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchI];

            // Create patch dictionary and use sourcePatch type and name
            dictionary patchDict;
            patchDict.set("type", pp.type());
            patchDict.set("nFaces", pp.size());
            patchDict.set("startFace", pp.start());

            // Add patch
            allPatches.append
            (
                polyPatch::New
                (
                    pp.name(),
                    patchDict,
                    newPatchID,
                    mesh.boundaryMesh()
                ).ptr()
            );

            newPatchID++;
        }
    }

    allPatches.shrink();
    mesh.removeBoundary();
    mesh.addPatches(allPatches);
}


void patchManipulationFunctions::mergeTwoPatches
(
    polyMesh& mesh,
    const word& patchName1,
    const word& patchName2
)
{
    // Create mesh modifier
    directTopoChange meshMod(mesh);

    // Find the patch1
    const label patch1ID = mesh.boundaryMesh().findPatchID(patchName1);
    if (patch1ID == -1)
    {
        FatalErrorIn
        (
            "void mergeTwoPatches\n"
            "(\n"
            "    polyMesh& mesh,\n"
            "    const word& patchName1,\n"
            "    const word& patchName2\n"
            ")"
        )   << "Cannot find patch " << patchName1 << " in the mesh!"
            << abort(FatalError);
    }

    // Find the patch2
    const label patch2ID = mesh.boundaryMesh().findPatchID(patchName2);
    if (patch2ID == -1)
    {
        FatalErrorIn
        (
            "void mergeTwoPatches\n"
            "(\n"
            "    polyMesh& mesh,\n"
            "    const word& patchName1,\n"
            "    const word& patchName2\n"
            ")"
        )   << "Cannot find patch " << patchName2 << " in the mesh!"
            << abort(FatalError);
    }

    // Move all the faces in the patch2 to patch1
    const polyPatch& patch2 = mesh.boundaryMesh()[patch2ID];
    forAll(patch2, faceI)
    {
        const label faceID = patch2.start() + faceI;

        const label zoneID = mesh.faceZones().whichZone(faceID);
        bool zoneFlip = false;
        if (zoneID >= 0)
        {
            const faceZone& fZone = mesh.faceZones()[zoneID];
            zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
        }

        // Move face to the correct patch
        meshMod.setAction
        (
            polyModifyFace
            (
                mesh.faces()[faceID],              // face
                faceID,                            // face ID,
                mesh.faceOwner()[faceID],          // owner
                -1,                                // neighbour,
                false,                             // flip flux
                patch1ID,                          // patchID,
                false,                             // remove from zone
                zoneID,                            // zone ID
                zoneFlip                           // zone flip
            )
        );
    }

    // Perform the topo changes on the mesh
    meshMod.changeMesh(mesh, false);


    // Remove patch2, which now has zero faces
    removePatchWithNoFaces(mesh, patchName2);
}


void patchManipulationFunctions::addFaceZoneFromPatch
(
    polyMesh& mesh,
    const word& patchName,
    const word& faceZoneName,
    const bool flipMap
)
{
    // Find patchID
    const label patchID = mesh.boundaryMesh().findPatchID(patchName);

    if (patchID == -1)
    {
        FatalErrorIn("void patchManipulationFunctions::addFaceZoneFromPatch()")
            << "Cannot find " << patchID << " patch in the mesh"
            << abort(FatalError);
    }

    // Add a new faceZone to the end of the faceZone list

    mesh.faceZones().setSize(mesh.faceZones().size() + 1);

    //
    const polyPatch& pp = mesh.boundaryMesh()[patchID];

    labelList faceIDs(pp.size(), 0);
    forAll(pp, faceI)
    {
        faceIDs[faceI] = pp.start() + faceI;
    }

    // Add face zone
    mesh.faceZones().set
    (
        mesh.faceZones().size() - 1,
        new faceZone
        (
            faceZoneName,
            faceIDs,
            boolList(faceIDs.size(), flipMap),
            mesh.faceZones().size() - 1, // faceZoneID
            mesh.faceZones()
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// } // End namespace patchManipulationFunctions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
