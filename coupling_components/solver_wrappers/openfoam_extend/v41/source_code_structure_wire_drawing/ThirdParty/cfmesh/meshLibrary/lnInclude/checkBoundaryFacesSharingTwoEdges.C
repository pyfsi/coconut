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

#include "checkBoundaryFacesSharingTwoEdges.H"
#include "demandDrivenData.H"
#include "helperFunctions.H"
#include "meshSurfaceEngine.H"
#include "decomposeFaces.H"
#include "decomposeCells.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGCheck

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

void checkBoundaryFacesSharingTwoEdges::createMeshSurface() const
{
    meshSurfacePtr_ = new meshSurfaceEngine(mesh_);
}

void checkBoundaryFacesSharingTwoEdges::clearMeshSurface()
{
    deleteDemandDrivenData(meshSurfacePtr_);
}

void checkBoundaryFacesSharingTwoEdges::findFacesAtBndEdge()
{
    const meshSurfaceEngine& mse = meshSurface();

    const labelLongList& bp = mse.bp();
    const edgeLongList& edges = mse.edges();
    const VRWGraph& pointEdges = mse.boundaryPointEdges();

    const label nIntFaces = mesh_.nInternalFaces();
    const faceListPMG& faces = mesh_.faces();

    //- find the internal faces attached to the boundary points
    removeBndPoint_.setSize(pointEdges.size());
    removeBndPoint_ = true;

    # ifdef USE_OMP
    # pragma omp parallel for if( nIntFaces > 100 ) schedule(guided, 100)
    # endif
    for(label fI=0;fI<nIntFaces;++fI)
    {
        const face& f = faces[fI];

        forAll(f, pI)
        {
            const label bpI = bp[f[pI]];

            if( bpI < 0 )
                continue;

            if( nBndFacesAtBndPoint_[bpI] == 2 )
            {
                const edge ePrev = f.faceEdge(f.rcIndex(pI));
                const edge eNext = f.faceEdge(pI);

                //- check if both edges attached to this point
                //- are boundary edges
                bool foundNext(false), foundPrev(false);
                forAllRow(pointEdges, bpI, peI)
                {
                    const label beI = pointEdges(bpI, peI);

                    if( edges[beI] == ePrev )
                    {
                        foundPrev = true;
                    }
                    else if( edges[beI] == eNext )
                    {
                        foundNext = true;
                    }
                }

                //- point can be removed only if both edges are at the boundary
                if( !(foundPrev && foundNext) )
                {
                    removeBndPoint_[bpI] = false;
                }
            }
            else
            {
                removeBndPoint_[bpI] = false;
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- check processor faces
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();

        const label start = procBoundaries[0].patchStart();
        const label end = faces.size();

        # ifdef USE_OMP
        # pragma omp parallel for schedule(guided, 100)
        # endif
        for(label faceI=start;faceI<end;++faceI)
        {
            const face& f = faces[faceI];

            forAll(f, pI)
            {
                const label bpI = bp[f[pI]];

                if( bpI < 0 )
                    continue;

                if( nBndFacesAtBndPoint_[bpI] == 2 )
                {
                    const edge ePrev = f.faceEdge(f.rcIndex(pI));
                    const edge eNext = f.faceEdge(pI);

                    bool foundNext(false), foundPrev(false);
                    forAllRow(pointEdges, bpI, peI)
                    {
                        const label beI = pointEdges(bpI, peI);

                        if( edges[beI] == ePrev )
                        {
                            foundPrev = true;
                        }
                        else if( edges[beI] == eNext )
                        {
                            foundNext = true;
                        }
                    }

                    if( !(foundPrev && foundNext) )
                        removeBndPoint_[bpI] = false;
                }
                else
                {
                    removeBndPoint_[bpI] = false;
                }
            }
        }

        //- make sure that all processors have the same information
        const DynList<label>& bpNei = mse.bpNeiProcs();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();

        std::map<label, labelLongList> exchangeData;
        forAll(bpNei, i)
            exchangeData.insert(std::make_pair(bpNei[i], labelLongList()));

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( removeBndPoint_[bpI] )
                continue;

            //- the point shall not be removed
            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);
                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append(it.key());
            }
        }

        //- exchange data
        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        //- set remove flag to false
        forAll(receivedData, i)
            removeBndPoint_[globalToLocal[receivedData[i]]] = false;
    }
}

void checkBoundaryFacesSharingTwoEdges::findBndFacesAtBndVertex()
{
    const meshSurfaceEngine& mse = meshSurface();
    const VRWGraph& pointFaces = mse.pointFaces();

    nBndFacesAtBndPoint_.setSize(pointFaces.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(nBndFacesAtBndPoint_, bpI)
        nBndFacesAtBndPoint_[bpI] = pointFaces.sizeOfRow(bpI);

    if( Pstream::parRun() )
    {
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();
        const DynList<label>& neiProcs = mse.bpNeiProcs();

        //- create data that shall be exhcnaged
        std::map<label, labelLongList> exchangeData;
        forAll(neiProcs, i)
            exchangeData.insert(std::make_pair(neiProcs[i], labelLongList()));

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);
                if( neiProc == Pstream::myProcNo() )
                    continue;

                labelLongList& data = exchangeData[neiProc];
                data.append(it.key());
                data.append(nBndFacesAtBndPoint_[bpI]);
            }
        }

        //- exchange data with other processors
        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        label counter(0);
        while( counter < receivedData.size() )
        {
            const label bpI = globalToLocal[receivedData[counter++]];
            nBndFacesAtBndPoint_[bpI] += receivedData[counter++];
        }
    }
}

bool checkBoundaryFacesSharingTwoEdges::removeExcessiveVertices()
{
    const labelLongList& bp = meshSurface().bp();

    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();

    //- remove points causing faces to share two edges
    //- if all faces at those edges share the same two edges
    //- then these edges can be merged into a single edge

    //- start processing internal faces
    const label nIntFaces = mesh_.nInternalFaces();

    bool changed(false);

    # ifdef USE_OMP
    # pragma omp parallel for if( nIntFaces > 100 ) schedule(guided, 100)
    # endif
    for(label faceI=0;faceI<nIntFaces;++faceI)
    {
        const face& f = faces[faceI];

        DynList<label> newF;
        forAll(f, pI)
        {
            const label bpI = bp[f[pI]];

            if(
                (bpI >= 0) && removeBndPoint_[bpI] &&
                (nBndFacesAtBndPoint_[bpI] == 2)
            )
                continue;

            newF.append(f[pI]);
        }

        if( newF.size() < f.size() )
        {
            changed = true;
            face& mf = faces[faceI];
            mf.setSize(newF.size());
            forAll(mf, i)
                mf[i] = newF[i];
        }
    }

    //- boundary faces
    const faceList::subList& bFaces = meshSurface().boundaryFaces();
    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100)
    # endif
    forAll(bFaces, bfI)
    {
        const label faceI = nIntFaces + bfI;

        const face& f = faces[faceI];

        DynList<label> newF;
        forAll(f, pI)
        {
            const label bpI = bp[f[pI]];

            if( removeBndPoint_[bpI] && (nBndFacesAtBndPoint_[bpI] == 2) )
                continue;

            newF.append(f[pI]);
        }

        if( newF.size() < f.size() )
        {
            changed = true;
            face& mf = faces[faceI];
            mf.setSize(newF.size());
            forAll(mf, i)
                mf[i] = newF[i];
        }
    }

    //- processor boundaries
    forAll(mesh_.procBoundaries(), patchI)
    {
        const processorBoundaryPatch& patch = mesh_.procBoundaries()[patchI];
        const label start = patch.patchStart();
        const label end = start + patch.patchSize();

        # ifdef USE_OMP
        # pragma omp parallel for schedule(guided, 100)
        # endif
        for(label faceI=start;faceI<end;++faceI)
        {
            const face& f = faces[faceI];

            DynList<label> newF;
            forAll(f, pI)
            {
                const label bpI = bp[f[pI]];

                if(
                    (bpI >= 0) && removeBndPoint_[bpI] &&
                    (nBndFacesAtBndPoint_[bpI] == 2)
                )
                    continue;

                newF.append(f[pI]);
            }

            if( newF.size() < f.size() )
            {
                changed = true;
                face& mf = faces[faceI];
                mf.setSize(newF.size());

                if( !patch.owner() && (newF[0] != f[0]) )
                {
                    forAll(mf, i)
                        mf[i] = newF[mf.rcIndex(i)];
                }
                else
                {
                    forAll(mf, i)
                        mf[i] = newF[i];
                }
            }
        }
    }

    //- delete mesh surface
    clearMeshSurface();

    return returnReduce(changed, maxOp<bool>());
}

label checkBoundaryFacesSharingTwoEdges::findBndFacesForDecomposition
(
    boolList& decomposeFace
)
{
    const meshSurfaceEngine& mse = meshSurface();
    const labelLongList& bp = mse.bp();
    const faceList::subList& bFaces = mse.boundaryFaces();

    label nDecomposed(0);

    //- find boundary sharing two common edges
    //- these faces have to be decomposed to end up with a valid mesh
    # ifdef USE_OMP
    # pragma omp parallel for if( bFaces.size() > 100 ) \
    schedule(guided, 100) reduction(+ : nDecomposed)
    # endif
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        forAll(bf, pI)
        {
            const label bpI = bp[bf[pI]];

            if( nBndFacesAtBndPoint_[bpI] == 2 )
            {
                ++nDecomposed;
                decomposeFace[bfI] = true;
            }
        }
    }

    reduce(nDecomposed, sumOp<label>());

    return nDecomposed;
}

void checkBoundaryFacesSharingTwoEdges::decomposeBoundaryFaces
(
    const boolList& decomposeFace
)
{
    const pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = meshSurface();
    const labelLongList& bp = mse.bp();
    const VRWGraph& faceEdges = mse.faceEdges();
    const edgeLongList& edges = mse.edges();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelLongList& faceOwner = mse.faceOwners();
    const labelLongList& facePatch = mse.boundaryFacePatches();

    //- create points at split edges
    labelLongList newEdgeLabel(edges.size(), -1);

    VRWGraph newBoundaryFaces;
    labelLongList newFaceOwners;
    labelLongList newFacePatches;

    # ifdef DEBUGCheck
    const label cId = mesh_.addCellSubset("decBndFacesCells");
    # endif

    forAll(decomposeFace, bfI)
    {
        if( decomposeFace[bfI] )
        {
            const face& bf = bFaces[bfI];

            # ifdef DEBUGCheck
            mesh_.addCellToSubset(cId, faceOwner[bfI]);
            # endif

            //- find the vertices that were generated by splitting an edge
            DynList<bool> edgeCentreVertex(bf.size(), false);

            forAll(bf, pI)
            {
                const point& p = points[bf[pI]];

                vector prev = points[bf.prevLabel(pI)] - p;
                prev /= (mag(prev) + VSMALL);
                vector next = points[bf.nextLabel(pI)] - p;
                next /= (mag(next) + VSMALL);

                if( mag(prev & next) > 0.5 )
                    edgeCentreVertex[pI] = true;
            }

            const label fLabel = points.size();
            mesh_.appendVertex(help::faceCentre(points, bf));

            //- split remaining edges
            forAllRow(faceEdges, bfI, feI)
            {
                const label beI = faceEdges(bfI, feI);

                const edge& e = edges[beI];

                if( edgeCentreVertex[feI] || edgeCentreVertex[bf.fcIndex(feI)] )
                    continue;

                if( newEdgeLabel[beI] == -1 )
                {
                    newEdgeLabel[beI] = points.size();
                    const point newP = e.centre(points);
                    mesh_.appendVertex(newP);
                }
            }

            //- split the face into 4 new faces
            forAll(bf, pI)
            {
                if( edgeCentreVertex[pI] )
                    continue;

                FixedList<label, 4> newF;
                newF[0] = bf[pI];
                if( edgeCentreVertex[bf.fcIndex(pI)] )
                {
                    newF[1] = bf.nextLabel(pI);
                }
                else
                {
                    newF[1] = newEdgeLabel[faceEdges(bfI, pI)];
                }
                newF[2] = fLabel;
                if( edgeCentreVertex[bf.rcIndex(pI)] )
                {
                    newF[3] = bf.prevLabel(pI);
                }
                else
                {
                    newF[3] = newEdgeLabel[faceEdges(bfI, bf.rcIndex(pI))];
                }

                newBoundaryFaces.appendList(newF);
                newFaceOwners.append(faceOwner[bfI]);
                newFacePatches.append(facePatch[bfI]);
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- reduce the information over all processor sharing each split edge
        const Map<label>& globalToLocal = mse.globalToLocalBndEdgeAddressing();
        const VRWGraph& beAtProcs = mse.beAtProcs();

        std::map<label, LongList<label> > exchangeData;
        forAll(mse.beNeiProcs(), i)
            exchangeData[mse.beNeiProcs()[i]].clear();

        //- assemble the labels of split edges that
        //- shall be sent to other processors
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label beI = it();

            if( newEdgeLabel[beI] < 0 )
                continue;

            forAllRow(beAtProcs, beI, i)
            {
                const label neiProc = beAtProcs(beI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append(it.key());
            }
        }

        //- exchange data
        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        //- check if the edge is marked for decomposition
        forAll(receivedData, i)
        {
            const label beI = globalToLocal[receivedData[i]];

            if( newEdgeLabel[beI] == -1 )
            {
                newEdgeLabel[beI] = points.size();
                const point newP = edges[beI].centre(points);
                mesh_.appendVertex(newP);
            }
        }
    }

    //- update remaining boundary faces
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        if( !decomposeFace[bfI] )
        {
            DynList<label> newF;

            forAllRow(faceEdges, bfI, feI)
            {
                newF.append(bf[feI]);

                const label eLabel = newEdgeLabel[faceEdges(bfI, feI)];
                if( eLabel != -1 )
                    newF.append(eLabel);
            }

            newBoundaryFaces.appendList(newF);
            newFaceOwners.append(faceOwner[bfI]);
            newFacePatches.append(facePatch[bfI]);
        }
    }

    //- modify the mesh
    polyMeshGenModifier meshModifier(mesh_);

    //- replace the boundary with te new faces
    wordList patchNames(mesh_.boundaries().size());
    forAll(patchNames, patchI)
        patchNames[patchI] = mesh_.boundaries()[patchI].patchName();
    meshModifier.replaceBoundary
    (
        patchNames,
        newBoundaryFaces,
        newFaceOwners,
        newFacePatches
    );

    //- update internal faces with newly create vertices
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const label nIntFaces = mesh_.nInternalFaces();
    faceListPMG& faces = meshModifier.facesAccess();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    for(label faceI=0;faceI<nIntFaces;++faceI)
    {
        face& f = faces[faceI];

        DynList<label> newF;

        forAll(f, eI)
        {
            newF.append(f[eI]);

            const label bpI = bp[f[eI]];

            if( bpI < 0 )
                continue;

            const edge e = f.faceEdge(eI);

            forAllRow(bpEdges, bpI, peI)
            {
                const label beI = bpEdges(bpI, peI);

                if( e == edges[beI] )
                {
                    if( newEdgeLabel[beI] != -1 )
                        newF.append(newEdgeLabel[beI]);
                }
            }
        }

        if( newF.size() > f.size() )
        {
            //- face has been altered
            f.setSize(newF.size());
            forAll(f, pI)
                f[pI] = newF[pI];
        }
    }

    if( Pstream::parRun() )
    {
        //- update faces at processor boundaries
        const label start = mesh_.procBoundaries()[0].patchStart();
        const label end = faces.size();
        # ifdef USE_OMP
        # pragma omp parallel for schedule(static, 1)
        # endif
        for(label faceI=start;faceI<end;++faceI)
        {
            face& f = faces[faceI];

            DynList<label> newF;

            forAll(f, eI)
            {
                newF.append(f[eI]);

                const label bpI = bp[f[eI]];

                if( bpI < 0 )
                    continue;

                const edge e = f.faceEdge(eI);

                forAllRow(bpEdges, bpI, peI)
                {
                    const label beI = bpEdges(bpI, peI);

                    if( e == edges[beI] )
                    {
                        if( newEdgeLabel[beI] != -1 )
                            newF.append(newEdgeLabel[beI]);
                    }
                }
            }

            if( newF.size() > f.size() )
            {
                //- face has been altered
                f.setSize(newF.size());
                forAll(f, pI)
                    f[pI] = newF[pI];
            }
        }
    }

    meshModifier.clearAll();

    //- delete meshSurfaceEngine because it is not valid any more
    clearMeshSurface();
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //
// Constructors

checkBoundaryFacesSharingTwoEdges::checkBoundaryFacesSharingTwoEdges
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    meshSurfacePtr_(NULL),
    nBndFacesAtBndPoint_(),
    removeBndPoint_()
{}

// * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * * * * //

checkBoundaryFacesSharingTwoEdges::~checkBoundaryFacesSharingTwoEdges()
{
    clearMeshSurface();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkBoundaryFacesSharingTwoEdges::findPoints(labelHashSet& badPoints)
{
    badPoints.clear();

    findBndFacesAtBndVertex();

    const labelLongList& bPoints = meshSurface().boundaryPoints();
    forAll(nBndFacesAtBndPoint_, bpI)
    {
        if( nBndFacesAtBndPoint_[bpI] != 2 )
            continue;

        badPoints.insert(bPoints[bpI]);
    }
}

bool checkBoundaryFacesSharingTwoEdges::improveTopology()
{
    bool changed(false);

    findBndFacesAtBndVertex();

    findFacesAtBndEdge();

    if( removeExcessiveVertices() )
    {
        findBndFacesAtBndVertex();

        findFacesAtBndEdge();
    }

    boolList decomposeFace(meshSurface().boundaryFaces().size(), false);
    const label nDecomposed = findBndFacesForDecomposition(decomposeFace);

    Info << "Marked " << nDecomposed << " faces for decomposition" << endl;

    if( nDecomposed != 0 )
    {
        decomposeBoundaryFaces(decomposeFace);

        changed = true;
    }

    polyMeshGenModifier(mesh_).removeUnusedVertices();
    polyMeshGenModifier(mesh_).clearAll();

    return changed;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
