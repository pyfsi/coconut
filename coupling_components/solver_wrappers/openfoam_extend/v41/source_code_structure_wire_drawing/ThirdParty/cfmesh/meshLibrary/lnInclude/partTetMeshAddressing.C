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

#include "demandDrivenData.H"
#include "polyMeshGenModifier.H"
#include "partTetMesh.H"
#include "polyMeshGenAddressing.H"
#include "helperFunctions.H"
#include "VRWGraphSMPModifier.H"

#include "partTetMeshSimplex.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void partTetMesh::createPointsAndTets
(
    const List<direction>& useCell,
    const boolList& lockedPoints
)
{
    const pointFieldPMG& points = origMesh_.points();
    const faceListPMG& faces = origMesh_.faces();
    const PtrList<boundaryPatch>& boundaries = origMesh_.boundaries();
    const PtrList<processorBoundaryPatch>& procBoundaries =
        origMesh_.procBoundaries();
    const label nInternalFaces = origMesh_.nInternalFaces();

    //- check how many neighbours of a face are marked for smoothing
    labelList usedFace(faces.size());
    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        //- initialise he array to zero
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(usedFace, faceI)
            usedFace[faceI] = 0;

        //- mark faces
        # ifdef USE_OMP
        # pragma omp for schedule(guided, 100)
        # endif
        forAll(faces, faceI)
        {
            if( useCell[owner_[faceI]] )
                ++usedFace[faceI];

            if( neighbour_[faceI] < 0 )
                continue;

            if( useCell[neighbour_[faceI]] )
                ++usedFace[faceI];
        }
    }

    //- send data at processor boundaries
    forAll(procBoundaries, patchI)
    {
        const label start = procBoundaries[patchI].patchStart();
        const label size = procBoundaries[patchI].patchSize();

        labelLongList dataToSend;
        for(label faceI=0;faceI<size;++faceI)
        {
            if( usedFace[start+faceI] )
                dataToSend.append(faceI);
        }

        OPstream toOtherProc
        (
            Pstream::blocking,
            procBoundaries[patchI].neiProcNo(),
            dataToSend.byteSize()
        );

        toOtherProc << dataToSend;
    }

    //- receive data at proc boundaries
    forAll(procBoundaries, patchI)
    {
        labelList receivedData;

        IPstream fromOtherProc
        (
            Pstream::blocking,
            procBoundaries[patchI].neiProcNo()
        );

        fromOtherProc >> receivedData;

        const label start = procBoundaries[patchI].patchStart();
        forAll(receivedData, faceI)
            ++usedFace[start+receivedData[faceI]];
    }

    smoothVertex_.setSize(points.size());

    List<direction> faceType(faces.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        //- initialize smoothVertex_
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(smoothVertex_, pointI)
        {
            if( !lockedPoints[pointI] )
            {
                smoothVertex_[pointI] = NONE;
            }
            else
            {
                smoothVertex_[pointI] = LOCKED;
            }
        }

        //- mark face types
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(usedFace, faceI)
        {
            direction fType(NONE);
            if( usedFace[faceI] == 2 )
            {
                fType |= SMOOTH;
            }
            else if( usedFace[faceI] == 1 )
            {
                if( faceI < nInternalFaces )
                {
                    fType |= INTERNALBOUNDARY;
                }
                else if( !Pstream::parRun() )
                {
                    fType |= BOUNDARY;
                }
                else if( faceI < procBoundaries[0].patchStart() )
                {
                    fType |= BOUNDARY;
                }
                else
                {
                    fType |= PARALLELBOUNDARY;
                }
            }

            faceType[faceI] = fType;
        }

        //- mark point types
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(pointFaces_, pointI)
        {
            forAllRow(pointFaces_, pointI, pfI)
                smoothVertex_[pointI] |= faceType[pointFaces_(pointI, pfI)];
        }
    }

    //- create face-to-patch addressing
    label nBndFaces(0);
    forAll(boundaries, patchI)
        nBndFaces += boundaries[patchI].patchSize();

    bndFaceToPatch_.setSize(nBndFaces);

    nBndFaces = 0;

    forAll(boundaries, patchI)
    {
        const label size = boundaries[patchI].patchSize();
        for(label fI=0;fI<size;++fI)
            bndFaceToPatch_[nBndFaces++] = patchI;
    }

    //- create addressing for parallel runs
    if( Pstream::parRun() )
    {
        createParallelAddressing();

        createBufferLayers();
    }

    //- create lists of movable points
    movableInternal_.clear();
    movableBoundary_.clear();

    forAll(smoothVertex_, pI)
    {
        if( smoothVertex_[pI] & SMOOTH )
        {
            if( smoothVertex_[pI] == SMOOTH )
            {
                movableInternal_.append(pI);
            }
            else
            {
                smoothVertex_[pI] ^= SMOOTH;
            }
        }
        if( smoothVertex_[pI] & BOUNDARY )
        {
            if( smoothVertex_[pI] == BOUNDARY )
            {
                movableBoundary_.append(pI);
            }
            else
            {
                smoothVertex_[pI] ^= BOUNDARY;
            }
        }
    }

    # ifdef DEBUGSmooth
    Info << "Finished constructing tet mesh" << endl;
    returnReduce(1, sumOp<label>());

    forAll(smoothVertex_, pI)
    {
        const label origPointI = nodeLabelInOrigMesh_[pI];

        if( origPointI < 0 && (smoothVertex_[pI] & SMOOTH) )
            Pout << "Strange crap " << endl;
        if( origPointI < 0&& (smoothVertex_[pI] & BOUNDARY) )
            Pout << "1.Strange crap " << endl;
    }

    forAll(nodeLabelInOrigMesh_, pI)
        if(
            (nodeLabelInOrigMesh_[pI] != -1) &&
            (mag(points_[pI] - points[nodeLabelInOrigMesh_[pI]]) > SMALL)
        )
            FatalErrorIn
            (
                "void partTetMesh::createPointsAndTets"
                "(const boolList& useCell)"
            ) << "Node " << pI << " is dislocated" << abort(FatalError);
    # endif
}

void partTetMesh::createTetsAtFace
(
    const label faceI,
    const label pointI,
    std::map<label, label>& pToPts,
    std::map<label, label>& fcToPts,
    std::map<label, label>& ccToPts,
    DynList<point, 128>& pts,
    DynList<partTet, 256>& tets,
    DynList<labelledTri, 32>& bndTriangles
) const
{
    const face& f = faces_[faceI];

    if( pToPts.find(pointI) == pToPts.end() )
    {
        pToPts[pointI] = pts.size();
        pts.append(points_[pointI]);
    }

    //- insert face centre
    if( f.size() > 3 )
    {
        fcToPts[faceI] = pts.size();
        pts.append(faceCentres_[faceI]);
    }

    //- insert centre of owner cell
    const label cOwn = owner_[faceI];
    if( ccToPts.find(cOwn) == ccToPts.end() )
    {
        ccToPts[cOwn] = pts.size();
        pts.append(cellCentres_[cOwn]);
    }

    //- insert centre of neighbour cell
    const label cNei = neighbour_[faceI];
    if( cNei >= 0 && ccToPts.find(cNei) == ccToPts.end() )
    {
        ccToPts[cNei] = pts.size();
        pts.append(cellCentres_[cNei]);
    }

    //- insert face vertices
    const label pos = f.which(pointI);

    forAll(f, pI)
    {
        if( pToPts.find(f[pI]) == pToPts.end() )
        {
            pToPts[f[pI]] = pts.size();
            pts.append(points_[f[pI]]);
        }
    }

    //- check if the face is at the boundary
    const label bfI = faceI - origMesh_.nInternalFaces();
    if( bfI >= 0 && (bfI < bndFaceToPatch_.size()) )
    {
        const label facePatch = bndFaceToPatch_[bfI];

        if( f.size() > 3 )
        {
            bndTriangles.append
            (
                labelledTri
                (
                    0,
                    pToPts[f.nextLabel(pos)],
                    fcToPts[faceI],
                    facePatch
                )
            );

            bndTriangles.append
            (
                labelledTri
                (
                    0,
                    fcToPts[faceI],
                    pToPts[f.prevLabel(pos)],
                    facePatch
                )
            );
        }
        else
        {
            bndTriangles.append
            (
                labelledTri
                (
                    0,
                    pToPts[f.nextLabel(pos)],
                    pToPts[f.prevLabel(pos)],
                    facePatch
                )
            );
        }
    }

    //- create tets on the owner side
    if( f.size() > 3 )
    {
        //- first tet attached to the face centre
        tets.append
        (
            partTet
            (
                pToPts[f.nextLabel(pos)],
                fcToPts[faceI],
                ccToPts[cOwn],
                0
            )
        );

        //- second tet attached to the face centre
        tets.append
        (
            partTet
            (
                pToPts[f.prevLabel(pos)],
                ccToPts[cOwn],
                fcToPts[faceI],
                0
            )
        );
    }
    else
    {
        //- check if the cell is a tet
        const cell& c = cells_[cOwn];

        if( c.size() == 4 )
        {
            //- cell is a tet
            forAll(c, fI)
            {
                const face& tf = faces_[c[fI]];

                //- skip faces containing pointI
                if( tf.which(pointI) >= 0 )
                    continue;

                //- check if all points are available in the simplex
                forAll(tf, pI)
                {
                    if( pToPts.find(tf[pI]) == pToPts.end() )
                    {
                        pToPts[tf[pI]] = pts.size();
                        pts.append(points_[tf[pI]]);
                    }
                }

                if( owner_[c[fI]] == cOwn )
                {
                    //- tet owns the face
                    tets.append
                    (
                        partTet
                        (
                            pToPts[tf[0]],
                            pToPts[tf[2]],
                            pToPts[tf[1]],
                            0
                        )
                    );
                }
                else
                {
                    //- tet is a neighbour of a face
                    tets.append
                    (
                        partTet
                        (
                            pToPts[tf[0]],
                            pToPts[tf[1]],
                            pToPts[tf[2]],
                            0
                        )
                    );
                }

                break;
            }
        }
        else
        {
            //- cell is not a tet
            tets.append
            (
                partTet
                (
                    pToPts[f.nextLabel(pos)],
                    pToPts[f.prevLabel(pos)],
                    ccToPts[cOwn],
                    0
                )
            );
        }
    }

    //- check if the neighbour exists
    if( cNei < 0 )
        return;

    //- generate tets on the neighbour side
    if( f.size() > 3 )
    {
        //- first tet attached to the face centre
        tets.append
        (
            partTet
            (
                pToPts[f.nextLabel(pos)],
                ccToPts[cNei],
                fcToPts[faceI],
                0
            )
        );

        //- second tet attached to the face centre
        tets.append
        (
            partTet
            (
                pToPts[f.prevLabel(pos)],
                fcToPts[faceI],
                ccToPts[cNei],
                0
            )
        );
    }
    else
    {
        //- check if the cell is a tet
        const cell& c = cells_[cNei];

        if( c.size() == 4 )
        {
            //- cell is a tet
            forAll(c, fI)
            {
                const face& tf = faces_[c[fI]];

                //- skip faces containing pointI
                if( tf.which(pointI) >= 0 )
                    continue;

                //- check if all points are available in the simplex
                forAll(tf, pI)
                {
                    if( pToPts.find(tf[pI]) == pToPts.end() )
                    {
                        pToPts[tf[pI]] = pts.size();
                        pts.append(points_[tf[pI]]);
                    }
                }

                if( owner_[c[fI]] == cNei )
                {
                    //- tet owns the face
                    tets.append
                    (
                        partTet
                        (
                            pToPts[tf[0]],
                            pToPts[tf[2]],
                            pToPts[tf[1]],
                            0
                        )
                    );
                }
                else
                {
                    //- tet is a neighbour
                    tets.append
                    (
                        partTet
                        (
                            pToPts[tf[0]],
                            pToPts[tf[1]],
                            pToPts[tf[2]],
                            0
                        )
                    );
                }

                break;
            }
        }
        else
        {
            //- cell is not a tet
            tets.append
            (
                partTet
                (
                    pToPts[f.prevLabel(pos)],
                    pToPts[f.nextLabel(pos)],
                    ccToPts[cNei],
                    0
                )
            );
        }
    }
}

void partTetMesh::createSimplex
(
    const label pointI,
    partTetMeshSimplex& simplex
) const
{
    std::map<label, DynList<partTet, 256> >::const_iterator tetIter =
        tetsAtPoint_.find(pointI);

    if( tetIter != tetsAtPoint_.end() )
    {
        //- generate the simplex at an inter-processor boundary
        const DynList<partTet, 256>& tets = tetIter->second;

        if( tets.size() == 0 )
            FatalErrorIn
            (
                "void partTetMesh::createSimplex"
                "(const label, partTetMeshSimplex&) const"
            ) << "Tets are not available at point "
              << pointI << abort(FatalError);

        simplex.pts_.clear();
        simplex.tets_.setSize(tets.size());

        std::map<label, label> parPointToPts;

        simplex.pts_.append(points_[pointI]);
        parPointToPts[tets[0].d()] = 0;

        forAll(tets, tetI)
        {
            const partTet& origTet = tets[tetI];

            if( parPointToPts.find(origTet.a()) == parPointToPts.end() )
            {
                parPointToPts[origTet.a()] = simplex.pts_.size();
                simplex.pts_.append(parPoints_[origTet.a()]);
            }

            if( parPointToPts.find(origTet.b()) == parPointToPts.end() )
            {
                parPointToPts[origTet.b()] = simplex.pts_.size();
                simplex.pts_.append(parPoints_[origTet.b()]);
            }

            if( parPointToPts.find(origTet.c()) == parPointToPts.end() )
            {
                parPointToPts[origTet.c()] = simplex.pts_.size();
                simplex.pts_.append(parPoints_[origTet.c()]);
            }

            //- copy tet
            simplex.tets_[tetI] =
                partTet
                (
                    parPointToPts[origTet.a()],
                    parPointToPts[origTet.b()],
                    parPointToPts[origTet.c()],
                    0
                );
        }

        std::map<label, DynList<labelledTri, 32> >::const_iterator bIt =
            bndTrianglesAtPoint_.find(pointI);
        if( bIt != bndTrianglesAtPoint_.end() )
        {
            const DynList<labelledTri, 32>& bndTriangles = bIt->second;

            simplex.bndTriangles_.setSize(bndTriangles.size());

            forAll(bndTriangles, triI)
            {
                const labelledTri& origTri = bndTriangles[triI];

                if( parPointToPts.find(origTri[0]) == parPointToPts.end() )
                {
                    parPointToPts[origTri[0]] = simplex.pts_.size();
                    simplex.pts_.append(parPoints_[origTri[0]]);
                }

                if( parPointToPts.find(origTri[1]) == parPointToPts.end() )
                {
                    parPointToPts[origTri[1]] = simplex.pts_.size();
                    simplex.pts_.append(parPoints_[origTri[1]]);
                }

                if( parPointToPts.find(origTri[2]) == parPointToPts.end() )
                {
                    parPointToPts[origTri[2]] = simplex.pts_.size();
                    simplex.pts_.append(parPoints_[origTri[2]]);
                }

                //- copy tet
                simplex.bndTriangles_[triI] =
                    labelledTri
                    (
                        parPointToPts[origTri[0]],
                        parPointToPts[origTri[1]],
                        parPointToPts[origTri[2]],
                        origTri.region()
                    );
            }
        }
    }
    else
    {
        //- create vertices
        std::map<label, label> ccToPts, fcToPts, pToPts;

        //- create tetrahedra for all faces
        forAllRow(pointFaces_, pointI, pfI)
        {
            createTetsAtFace
            (
                pointFaces_(pointI, pfI),
                pointI,
                pToPts,
                fcToPts,
                ccToPts,
                simplex.pts_,
                simplex.tets_,
                simplex.bndTriangles_
            );
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
