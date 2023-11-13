/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "polyMeshGenModifier.H"
#include "VRWGraphList.H"
#include "demandDrivenData.H"
#include "polyMesh.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::addMesh(const polyMeshGen& mesh)
{
    Info << "Merging meshes" << endl;

    pointFieldPMG& points = mesh_.points_;
    faceListPMG& faces = mesh_.faces_;
    cellListPMG& cells = mesh_.cells_;

    const label nOrigPoints = points.size();
    const label nOrigFaces = faces.size();
    const label nOrigIntFaces = mesh_.nInternalFaces();
    const label nOrigCells = cells.size();
    const label nOrigPatches = mesh_.boundaries_.size();

    DynList<word> patchNames;
    DynList<word> patchTypes;
    std::map<word, label> patchNameToIndex;
    DynList<label> patchSize;
    forAll(mesh_.boundaries_, patchI)
    {
        const boundaryPatch& patch = mesh_.boundaries_[patchI];

        patchNames.append(patch.patchName());
        patchTypes.append(patch.patchType());
        patchNameToIndex[patch.patchName()] = patchI;
        patchSize.append(patch.patchSize());
    }

    //- append vertices to the mesh
    points.setSize(nOrigPoints + mesh.points().size());
    faces.setSize(nOrigFaces + mesh.faces().size());
    cells.setSize(nOrigCells + mesh.cells().size());

    //- update patches and their size after the meshes are merged
    PtrList<boundaryPatch>& boundaries = mesh_.boundaries_;
    forAll(mesh.boundaries(), patchI)
    {
        const boundaryPatch& newPatch = mesh.boundaries()[patchI];

        const word pName = newPatch.patchName();

        if( patchNameToIndex.find(pName) == patchNameToIndex.end() )
        {
            //- patch with the same name does not exist in the mesh
            patchNameToIndex[pName] = patchNames.size();
            patchSize.append(newPatch.patchSize());
            patchTypes.append(newPatch.patchType());
            patchNames.append(pName);
        }
        else
        {
            //- add faces for the other mesh into the patch
            patchSize[patchNameToIndex[pName]] += newPatch.patchSize();
        }
    }

    //- calculate the starting face of each patch
    labelList patchStart(patchSize.size());
    patchStart[0] = nOrigIntFaces + mesh.nInternalFaces();
    for(label patchI=1;patchI<patchSize.size();++patchI)
    {
        patchStart[patchI] =
            patchStart[patchI-1] + patchSize[patchI-1];
    }

    labelLongList newFaceLabel(nOrigFaces);
    labelLongList newAddedFaceLabel(mesh.faces().size());

    faceList newFaces(newFaceLabel.size()+newAddedFaceLabel.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        //- copy mesh points
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1) nowait
        # endif
        forAll(mesh.points(), pointI)
            points[nOrigPoints+pointI] = mesh.points()[pointI];

        //- copy internal faces of the orig mesh
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1) nowait
        # endif
        for(label faceI=0;faceI<nOrigIntFaces;++faceI)
        {
            newFaceLabel[faceI] = faceI;
            newFaces[faceI].transfer(faces[faceI]);
        }

        //- copy internal faces of the second mesh
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1) nowait
        # endif
        for(label faceI=0;faceI<mesh.nInternalFaces();++faceI)
        {
            //- set the face label
            newAddedFaceLabel[faceI] = nOrigIntFaces + faceI;

            const face& origF = mesh.faces()[faceI];

            //- update the point labels in the face
            face newF(origF.size());
            forAll(origF, pI)
                newF[pI] = nOrigPoints + origF[pI];

            //- copy the faces into the common list of faces
            newFaces[nOrigIntFaces+faceI].transfer(newF);
        }

        //- copy the patches of the orig mesh
        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            //- copy the faces from the orig mesh
            forAll(boundaries, patchI)
            {
                # ifdef USE_OMP
                # pragma omp task default(shared) firstprivate(patchI)
                # endif
                {
                    const boundaryPatch& patch = boundaries[patchI];

                    const label start = patch.patchStart();
                    const label size = patch.patchSize();

                    for(label fI=0;fI<size;++fI)
                    {
                        const label newLabel = patchStart[patchI] + fI;

                        newFaces[newLabel].transfer(faces[start+fI]);
                        newFaceLabel[start+fI] = newLabel;
                    }
                }
            }

            //- copy boundary faces from the other mesh
            forAll(mesh.boundaries(), patchI)
            {
                # ifdef USE_OMP
                # pragma omp task default(shared) firstprivate(patchI)
                # endif
                {
                    const boundaryPatch& patch = mesh.boundaries()[patchI];

                    const label start = patch.patchStart();
                    const label size = patch.patchSize();

                    //- find the number of faces in the patch with the same
                    //- name located in the orig mesh
                    const label patchIndex =
                        patchNameToIndex[patch.patchName()];

                    label startFace = patchStart[patchIndex];
                    if( patchIndex < nOrigPatches )
                    {
                        startFace += boundaries[patchIndex].patchSize();
                    }

                    for(label fI=0;fI<size;++fI)
                    {
                        const label newLabel = startFace + fI;

                        const label faceI = start + fI;

                        //- set the face label
                        newAddedFaceLabel[faceI] = newLabel;

                        const face& origF = mesh.faces()[faceI];

                        //- update the point labels in the face
                        face newF(origF.size());
                        forAll(origF, pI)
                            newF[pI] = nOrigPoints + origF[pI];

                        //- copy the faces into the common list of faces
                        newFaces[newLabel].transfer(newF);
                    }
                }
            }
        }

        //- renumber the cells in the orig mesh
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1) nowait
        # endif
        for(label cellI=0;cellI<nOrigCells;++cellI)
        {
            cell& c = cells[cellI];

            forAll(c, fI)
                c[fI] = newFaceLabel[c[fI]];
        }

        //- copy the cells into the mesh
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(mesh.cells(), cellI)
        {
            const cell& c = mesh.cells()[cellI];

            cell newC(c.size());
            forAll(c, fI)
                newC[fI] = newAddedFaceLabel[c[fI]];

            cells[nOrigCells+cellI].transfer(newC);
        }

        //- tranfer faces back to
        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            # if USE_OMP
            const label nTasks = 5 * omp_get_num_threads();
            # else
            const label nTasks = 1;
            # endif

            faces.setSize(newFaces.size());

            for(label taskI=0;taskI<nTasks;++taskI)
            {
                # ifdef USE_OMP
                # pragma omp task default(shared) firstprivate(taskI)
                # endif
                {
                    for(label faceI=taskI;faceI<newFaces.size();faceI+=nTasks)
                    {
                        faces[faceI].transfer(newFaces[faceI]);
                    }
                }
            }
        }
    }

    newFaces.clear();

    //- update patches
    boundaries.clear();
    boundaries.setSize(patchNames.size());
    forAll(patchNames, patchI)
    {
        boundaries.set
        (
            patchI,
            new boundaryPatch
            (
                patchNames[patchI],
                patchTypes[patchI],
                patchSize[patchI],
                patchStart[patchI]
            )
        );
    }

    //- update face subset and zones
    mesh_.updateFaceSubsets(newFaceLabel);
    mesh_.updateFaceZones(newFaceLabel);

    //- copy point subsets
    DynList<label> ids;
    mesh.pointSubsetIndices(ids);
    forAll(ids, i)
    {
        const label origId = ids[i];
        const word name = mesh.pointSubsetName(origId);

        label newId = mesh_.pointSubsetIndex(name);
        if( newId < 0 )
            newId = mesh_.addPointSubset(name);

        labelLongList elmts;
        mesh.pointsInSubset(origId, elmts);

        forAll(elmts, i)
            mesh_.addPointToSubset(newId, nOrigPoints+elmts[i]);
    }

    //- copy face subsets
    ids.clear();
    mesh.faceSubsetIndices(ids);
    forAll(ids, i)
    {
        const label origId = ids[i];
        const word name = mesh.faceSubsetName(origId);

        label newId = mesh_.faceSubsetIndex(name);
        if( newId < 0 )
        {
            newId = mesh_.addFaceSubset(name);
        }

        labelLongList elmts;
        mesh.facesInSubset(origId, elmts);

        forAll(elmts, i)
        {
            mesh_.addFaceToSubset(newId, newAddedFaceLabel[elmts[i]]);
        }
    }

    //- copy cell subsets
    ids.clear();
    mesh.cellSubsetIndices(ids);
    forAll(ids, i)
    {
        const label origId = ids[i];
        const word name = mesh.cellSubsetName(origId);

        label newId = mesh_.cellSubsetIndex(name);
        if( newId < 0 )
            newId = mesh_.addCellSubset(name);

        labelLongList elmts;
        mesh.cellsInSubset(origId, elmts);

        forAll(elmts, i)
            mesh_.addCellToSubset(newId, nOrigCells+elmts[i]);
    }

    //- copy point zones
    std::map<label, label> zoneToZone;
    ids.clear();
    mesh.pointZoneIndices(ids);

    forAll(ids, i)
    {
        const word name = mesh.pointZoneName(ids[i]);

        label newId = mesh_.pointZoneIndex(name);
        if( newId < 0 )
            newId = mesh_.addPointZone(name);

        zoneToZone[ids[i]] = newId;
    }

    forAll(mesh.points(), pointI)
    {
        const label zoneI = mesh.pointZone(pointI);

        if( zoneI >= 0 )
            mesh_.addPointToZone(zoneToZone[zoneI], nOrigPoints+pointI);
    }

    //- copy face zones
    zoneToZone.clear();
    ids.clear();
    mesh.faceZoneIndices(ids);

    forAll(ids, i)
    {
        const word name = mesh.faceZoneName(ids[i]);

        label newId = mesh_.faceZoneIndex(name);
        if( newId < 0 )
            newId = mesh_.addFaceZone(name);

        zoneToZone[ids[i]] = newId;
    }

    forAll(mesh.faces(), faceI)
    {
        const label zoneI = mesh.faceZone(faceI);

        if( zoneI >= 0 )
            mesh_.addFaceToZone(zoneToZone[zoneI], newAddedFaceLabel[faceI]);
    }

    //- copy cell zones
    zoneToZone.clear();
    ids.clear();
    mesh.cellZoneIndices(ids);

    forAll(ids, i)
    {
        const word name = mesh.cellZoneName(ids[i]);

        label newId = mesh_.cellZoneIndex(name);
        if( newId < 0 )
            newId = mesh_.addCellZone(name);

        zoneToZone[ids[i]] = newId;
    }

    forAll(mesh.cells(), cellI)
    {
        const label zoneI = mesh.cellZone(cellI);

        if( zoneI >= 0 )
            mesh_.addCellToZone(zoneToZone[zoneI], nOrigCells+cellI);
    }

    //- delete all demand driven data
    mesh_.clearAddressingData();
    mesh_.clearOut();

    Info << "Finished merging meshes" << endl;
}

void polyMeshGenModifier::addMesh(const polyMesh& mesh)
{
    Info << "Merging meshes" << endl;

    pointFieldPMG& points = mesh_.points_;
    faceListPMG& faces = mesh_.faces_;
    cellListPMG& cells = mesh_.cells_;

    const label nOrigPoints = points.size();
    const label nOrigFaces = faces.size();
    const label nOrigIntFaces = mesh_.nInternalFaces();
    const label nOrigCells = cells.size();
    const label nOrigPatches = mesh_.boundaries_.size();

    DynList<word> patchNames;
    DynList<word> patchTypes;
    std::map<word, label> patchNameToIndex;
    DynList<label> patchSize;
    forAll(mesh_.boundaries_, patchI)
    {
        const boundaryPatch& patch = mesh_.boundaries_[patchI];

        patchNames.append(patch.patchName());
        patchTypes.append(patch.patchType());
        patchNameToIndex[patch.patchName()] = patchI;
        patchSize.append(patch.patchSize());
    }

    //- append vertices to the mesh
    points.setSize(nOrigPoints + mesh.nPoints());
    faces.setSize(nOrigFaces + mesh.nFaces());
    cells.setSize(nOrigCells + mesh.nCells());

    //- update patches and their size after the meshes are merged
    PtrList<boundaryPatch>& boundaries = mesh_.boundaries_;
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& newPatch = mesh.boundaryMesh()[patchI];

        const word pName = newPatch.name();

        if( patchNameToIndex.find(pName) == patchNameToIndex.end() )
        {
            //- patch with the same name does not exist in the mesh
            patchNameToIndex[pName] = patchNames.size();
            patchSize.append(newPatch.size());
            patchTypes.append(newPatch.type());
            patchNames.append(pName);
        }
        else
        {
            //- add faces for the other mesh into the patch
            patchSize[patchNameToIndex[pName]] += newPatch.size();
        }
    }

    //- calculate the starting face of each patch
    labelList patchStart(patchSize.size());
    patchStart[0] = nOrigIntFaces + mesh.nInternalFaces();
    for(label patchI=1;patchI<patchSize.size();++patchI)
    {
        patchStart[patchI] =
            patchStart[patchI-1] + patchSize[patchI-1];
    }

    labelLongList newFaceLabel(nOrigFaces);
    labelLongList newAddedFaceLabel(mesh.faces().size());

    faceList newFaces(newFaceLabel.size()+newAddedFaceLabel.size());

    const polyBoundaryMesh& bndMesh = mesh.boundaryMesh();

    //- copy mesh points
    forAll(mesh.points(), pointI)
        points[nOrigPoints+pointI] = mesh.points()[pointI];

    //- copy internal faces of the orig mesh
    for(label faceI=0;faceI<nOrigIntFaces;++faceI)
    {
        newFaceLabel[faceI] = faceI;
        newFaces[faceI].transfer(faces[faceI]);
    }

    //- copy internal faces of the second mesh
    for(label faceI=0;faceI<mesh.nInternalFaces();++faceI)
    {
        //- set the face label
        newAddedFaceLabel[faceI] = nOrigIntFaces + faceI;

        const face& origF = mesh.faces()[faceI];

        //- update the point labels in the face
        face newF(origF.size());
        forAll(origF, pI)
            newF[pI] = nOrigPoints + origF[pI];

        //- copy the faces into the common list of faces
        newFaces[nOrigIntFaces+faceI].transfer(newF);
    }

    //- copy the patches of the orig mesh
    //- copy the faces from the orig mesh
    forAll(boundaries, patchI)
    {
        # ifdef USE_OMP
        # pragma omp task default(shared) firstprivate(patchI)
        # endif
        {
            const boundaryPatch& patch = boundaries[patchI];

            const label start = patch.patchStart();
            const label size = patch.patchSize();

            for(label fI=0;fI<size;++fI)
            {
                const label newLabel = patchStart[patchI] + fI;

                newFaces[newLabel].transfer(faces[start+fI]);
                newFaceLabel[start+fI] = newLabel;
            }
        }
    }

    //- copy boundary faces from the other mesh
    forAll(bndMesh, patchI)
    {
        const polyPatch& patch = bndMesh[patchI];

        const label start = patch.start();
        const label size = patch.size();

        //- find the number of faces in the patch with the same
        //- name located in the orig mesh
        const label patchIndex =
            patchNameToIndex[patch.name()];

        label startFace = patchStart[patchIndex];
        if( patchIndex < nOrigPatches )
        {
            startFace += boundaries[patchIndex].patchSize();
        }

        for(label fI=0;fI<size;++fI)
        {
            const label newLabel = startFace + fI;

            const label faceI = start + fI;

            //- set the face label
            newAddedFaceLabel[faceI] = newLabel;

            const face& origF = mesh.faces()[faceI];

            //- update the point labels in the face
            face newF(origF.size());
            forAll(origF, pI)
                newF[pI] = nOrigPoints + origF[pI];

            //- copy the faces into the common list of faces
            newFaces[newLabel].transfer(newF);
        }
    }

    //- renumber the cells in the orig mesh
    for(label cellI=0;cellI<nOrigCells;++cellI)
    {
        cell& c = cells[cellI];

        forAll(c, fI)
            c[fI] = newFaceLabel[c[fI]];
    }

    //- copy the cells into the mesh
    forAll(mesh.cells(), cellI)
    {
        const cell& c = mesh.cells()[cellI];

        cell newC(c.size());
        forAll(c, fI)
            newC[fI] = newAddedFaceLabel[c[fI]];

        cells[nOrigCells+cellI].transfer(newC);
    }

    //- tranfer faces back to
    faces.setSize(newFaces.size());

    forAll(faces, faceI)
    {
        faces[faceI].transfer(newFaces[faceI]);
    }

    newFaces.clear();

    //- update patches
    boundaries.clear();
    boundaries.setSize(patchNames.size());
    forAll(patchNames, patchI)
    {
        boundaries.set
        (
            patchI,
            new boundaryPatch
            (
                patchNames[patchI],
                patchTypes[patchI],
                patchSize[patchI],
                patchStart[patchI]
            )
        );
    }

    //- update face subset and zones
    mesh_.updateFaceSubsets(newFaceLabel);
    mesh_.updateFaceZones(newFaceLabel);

    //- copy point zones
    std::map<label, label> zoneToZone;
    wordList zoneNames = mesh.pointZones().names();

    forAll(zoneNames, i)
    {
        const word& name = zoneNames[i];

        const label id = mesh.pointZones().findZoneID(name);

        label newId = mesh_.pointZoneIndex(name);
        if( newId < 0 )
            newId = mesh_.addPointZone(name);

        zoneToZone[id] = newId;
    }

    forAll(mesh.points(), pointI)
    {
        const label zoneI = mesh.pointZones().whichZone(pointI);

        if( zoneI >= 0 )
            mesh_.addPointToZone(zoneToZone[zoneI], nOrigPoints+pointI);
    }

    //- copy face zones
    zoneToZone.clear();
    zoneNames = mesh.faceZones().names();

    forAll(zoneNames, i)
    {
        const word name = zoneNames[i];

        const label id = mesh.faceZones().findZoneID(name);

        label newId = mesh_.faceZoneIndex(name);
        if( newId < 0 )
            newId = mesh_.addFaceZone(name);

        zoneToZone[id] = newId;
    }

    forAll(mesh.faces(), faceI)
    {
        const label zoneI = mesh.faceZones().whichZone(faceI);

        if( zoneI >= 0 )
            mesh_.addFaceToZone(zoneToZone[zoneI], newAddedFaceLabel[faceI]);
    }

    //- copy cell zones
    zoneToZone.clear();
    zoneNames = mesh.cellZones().names();

    forAll(zoneNames, i)
    {
        const word name = zoneNames[i];

        const label id = mesh.cellZones().findZoneID(name);

        label newId = mesh_.cellZoneIndex(name);
        if( newId < 0 )
            newId = mesh_.addCellZone(name);

        zoneToZone[id] = newId;
    }

    forAll(mesh.cells(), cellI)
    {
        const label zoneI = mesh.cellZones().whichZone(cellI);

        if( zoneI >= 0 )
            mesh_.addCellToZone(zoneToZone[zoneI], nOrigCells+cellI);
    }

    //- delete all demand driven data
    mesh_.clearAddressingData();
    mesh_.clearOut();

    Info << "Finished merging meshes" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
