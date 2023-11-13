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
#include "boundaryLayerOptimisation.H"
#include "meshSurfacePartitioner.H"
#include "meshSurfaceEngine.H"
#include "detectBoundaryLayers.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngineModifier.H"
#include "OFstream.H"

//#define DEBUGLayer

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayerOptimisation::writeVTK
(
    const fileName& fName,
    const pointField& origin,
    const vectorField& vecs
)
{
    if( origin.size() != vecs.size() )
        FatalErrorIn
        (
            "void boundaryLayerOptimisation::writeVTK(const fileName&,"
            " const pointField&, const vectorField&)"
        ) << "Sizes do not match" << abort(FatalError);

    OFstream file(fName);

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "POINTS " << 2*origin.size() << " float\n";
    forAll(origin, pI)
    {
        const point& p = origin[pI];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;

        const point op = p + vecs[pI];

        file << op.x() << ' ' << op.y() << ' ' << op.z() << nl;
    }

    //- write lines
    file << "\nLINES " << vecs.size()
         << " " << 3*vecs.size() << nl;
    forAll(vecs, eI)
    {
        file << 2 << " " << 2*eI << " " << (2*eI+1) << nl;
    }

    file << "\n";
}

void boundaryLayerOptimisation::writeHairEdges
(
    const fileName& fName,
    const direction eType,
    const vectorField& vecs
) const
{
    if( vecs.size() != hairEdges_.size() )
        FatalErrorIn
        (
            "void boundaryLayerOptimisation::writeHairEdges"
            "(const fileName&, const direction, const vectorField&) const"
        ) << "Sizes do not match" << abort(FatalError);

    //- count the number of hair edges matching this criteria
    label counter(0);

    forAll(hairEdgeType_, heI)
        if( hairEdgeType_[heI] & eType )
            ++counter;

    //- copy edge vector
    vectorField copyVecs(counter);
    pointField pts(counter);

    counter = 0;

    const pointFieldPMG& points = mesh_.points();

    forAll(hairEdgeType_, heI)
    {
        if( hairEdgeType_[heI] & eType )
        {
            const edge& he = hairEdges_[heI];

            pts[counter] = points[he.start()];
            copyVecs[counter] = vecs[heI] * he.mag(points);

            ++counter;
        }
    }

    //- write data to file
    writeVTK(fName, pts, copyVecs);
}

void boundaryLayerOptimisation::writeHairEdges
(
    const fileName& fName,
    const direction eType
) const
{
    //- count the number of hair edges matching this criteria
    label counter(0);

    forAll(hairEdgeType_, heI)
        if( hairEdgeType_[heI] & eType )
            ++counter;

    //- copy edge vector
    vectorField vecs(counter);
    pointField pts(counter);

    counter = 0;

    const pointFieldPMG& points = mesh_.points();

    forAll(hairEdgeType_, heI)
    {
        if( hairEdgeType_[heI] & eType )
        {
            const edge& he = hairEdges_[heI];

            pts[counter] = points[he.start()];
            vecs[counter] = he.vec(points);

            ++counter;
        }
    }

    //- write data to file
    writeVTK(fName,pts, vecs);
}

const meshSurfaceEngine& boundaryLayerOptimisation::meshSurface() const
{
    if( !meshSurfacePtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const meshSurfaceEngine&"
                " boundaryLayerOptimisation::meshSurface()"
            ) << "Cannot generate meshSurfaceEngine" << abort(FatalError);
        # endif

        meshSurfacePtr_ = new meshSurfaceEngine(mesh_);
    }

    return *meshSurfacePtr_;
}

const meshSurfacePartitioner&
boundaryLayerOptimisation::surfacePartitioner() const
{
    if( !partitionerPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const meshSurfacePartitioner& "
                "boundaryLayerOptimisation::surfacePartitioner()"
            ) << "Cannot generate meshSurfacePartitioner" << abort(FatalError);
        # endif

        partitionerPtr_ = new meshSurfacePartitioner(meshSurface());
    }

    return *partitionerPtr_;
}

void boundaryLayerOptimisation::calculateHairEdges()
{
    const meshSurfaceEngine& mse = meshSurface();
    const edgeLongList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const labelLongList& faceOwner = mse.faceOwners();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelLongList& bp = mse.bp();

    const meshSurfacePartitioner& mPart = surfacePartitioner();

    //- detect layers in the mesh
    const detectBoundaryLayers detectLayers(mPart);

    hairEdges_ = detectLayers.hairEdges();
    hairEdgesAtBndPoint_ = detectLayers.hairEdgesAtBndPoint();

    //- mark boundary faces which are base face for the boundary layer
    const labelLongList& layerAtBndFace = detectLayers.faceInLayer();
    isBndLayerBase_.setSize(bFaces.size());
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(layerAtBndFace, bfI)
    {
        if( layerAtBndFace[bfI] < 0 )
        {
            isBndLayerBase_[bfI] = false;
        }
        else
        {
            isBndLayerBase_[bfI] = true;
        }
    }

    # ifdef DEBUGLayer
    const label bndLayerFaceId = mesh_.addFaceSubset("bndLayerFaces");
    const label startBndFaceI = mesh_.boundaries()[0].patchStart();
    forAll(isBndLayerBase_, bfI)
        if( isBndLayerBase_[bfI] )
            mesh_.addFaceToSubset(bndLayerFaceId, startBndFaceI+bfI);
    # endif

    //- check if a face is an exiting face for a bnd layer
    isExitFace_.setSize(isBndLayerBase_.size());
    isExitFace_ = false;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(edgeFaces, edgeI)
    {
        //- avoid edges at inter-processor boundaries
        if( edgeFaces.sizeOfRow(edgeI) != 2 )
            continue;

        const label f0 = edgeFaces(edgeI, 0);
        const label f1 = edgeFaces(edgeI, 1);

        //- both faces have to be part of the same cell
        if( faceOwner[f0] != faceOwner[f1] )
            continue;

        //- check if the feature edge is attachd to exittin faces
        if
        (
            (isBndLayerBase_[f0] && (bFaces[f1].size() == 4)) &&
            (isBndLayerBase_[f1] && (bFaces[f0].size() == 4))
        )
        {
            isExitFace_[f0] = true;
            isExitFace_[f1] = true;
        }
    }

    # ifdef DEBUGLayer
    const label exittingFaceId = mesh_.addFaceSubset("exittingFaces");
    forAll(isExitFace_, bfI)
        if( isExitFace_[bfI] )
            mesh_.addFaceToSubset(exittingFaceId, startBndFaceI+bfI);
    # endif

    //- classify hair edges as they require different tretment
    //- in the smoothing procedure
    hairEdgeType_.setSize(hairEdges_.size());

    const labelHashSet& corners = mPart.corners();
    const labelHashSet& edgePoints = mPart.edgePoints();
    const labelHashSet& featureEdges = mPart.featureEdges();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        hairEdgeType_[hairEdgeI] = NONE;

        const edge& e = hairEdges_[hairEdgeI];
        const label bpI = bp[e.start()];

        //- check if it is a boundary edge
        forAllRow(bpEdges, bpI, peI)
        {
            const label beI = bpEdges(bpI, peI);

            const edge& be = edges[bpEdges(bpI, peI)];

            if( be == e )
            {
                hairEdgeType_[hairEdgeI] |= BOUNDARY;

                if( featureEdges.found(beI) )
                    hairEdgeType_[hairEdgeI] |= FEATUREEDGE;
            }

            if( corners.found(bpI) )
            {
                hairEdgeType_[hairEdgeI] |= ATCORNER;
            }
            else if( edgePoints.found(bpI) )
            {
                hairEdgeType_[hairEdgeI] |= ATEDGE;
            }
        }

        if( !(hairEdgeType_[hairEdgeI] & BOUNDARY) )
            hairEdgeType_[hairEdgeI] |= INSIDE;
    }

    thinnedHairEdge_.setSize(hairEdges_.size());

    //- calculate which other hair edges influence a hair edges
    //- and store it in a graph
    hairEdgesNearHairEdge_.setSize(hairEdges_.size());
    hairEndPointsAtPoint_.clear();
    nEndPointsAtPoint_.clear();

    const cellListPMG& cells = mesh_.cells();
    const faceList& faces = mesh_.faces();

    VRWGraph bpFacesHelper(bpEdges.size());
    forAll(faceOwner, bfI)
    {
        const label cellI = faceOwner[bfI];

        const cell& c = cells[cellI];

        forAll(c, fI)
        {
            const face& f = faces[c[fI]];

            forAll(f, pI)
            {
                const label bpI = bp[f[pI]];
                if( bpI < 0 )
                    continue;

                bpFacesHelper.appendIfNotIn(bpI, c[fI]);
            }
        }
    }

    forAll(hairEdges_, hairEdgeI)
    {
        const edge& e = hairEdges_[hairEdgeI];
        const label bpI = bp[e.start()];

        hairEndPointsAtPoint_[e.end()].append(hairEdgeI);

        DynList<label> neiHairEdges;

        //- find mesh faces comprising of the current hair edge
        forAllRow(bpFacesHelper, bpI, pfI)
        {
            const face& f = faces[bpFacesHelper(bpI, pfI)];

            //- face must be a quad
            if( f.size() != 4 )
                continue;

            //- check if the current face comprises of the hair edge
            label faceEdge(-1);
            forAll(f, eI)
                if( f.faceEdge(eI) == e )
                {
                    faceEdge = eI;
                    break;
                }

            if( faceEdge != -1 )
            {
                //- check if the opposite edge is also a hair edge
                const label eJ = (faceEdge+2) % 4;

                const edge fe = f.faceEdge(eJ);

                for(label i=0;i<2;++i)
                {
                    const label bpJ = bp[fe[i]];

                    if( bpJ >= 0 )
                    {
                        forAllRow(hairEdgesAtBndPoint_, bpJ, pI)
                        {
                            const label heJ = hairEdgesAtBndPoint_(bpJ, pI);
                            if( hairEdges_[heJ] == fe )
                                neiHairEdges.append(heJ);
                        }
                    }
                }
            }
        }

        hairEdgesNearHairEdge_.setRow(hairEdgeI, neiHairEdges);
    }

    //- calculate the number of end points attached to a point
    for
    (
        std::map<label, DynList<label, 3> >::const_iterator it =
            hairEndPointsAtPoint_.begin();
        it!=hairEndPointsAtPoint_.end();
        ++it
    )
        nEndPointsAtPoint_[it->first] = it->second.size();

    //- calculate hair edge to procs addressing
    hairEdgeAtProcs_.clear();
    hairEdgeNeiProcs_.clear();
    hairEndPointOwnedByProc_.clear();

    if( Pstream::parRun() )
    {
        const polyMeshGenAddressing& addr = mesh_.addressingData();
        const Map<label>& globalToLocal =
            addr.globalToLocalPointAddressing();
        const VRWGraph& pAtProcs = addr.pointAtProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(addr.pointNeiProcs(), i)
            exchangeData[addr.pointNeiProcs()[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label pointI = it();

            std::map<label, label>::const_iterator endIt =
                nEndPointsAtPoint_.find(pointI);
            if( endIt != nEndPointsAtPoint_.end() )
            {
                //- create hair edge to processors addressing
                const DynList<label, 3>& pHairs = hairEndPointsAtPoint_[pointI];
                forAll(pHairs, i)
                {
                    const edge& he = hairEdges_[pHairs[i]];
                    DynList<label, 3>& neiProcs = hairEdgeAtProcs_[pHairs[i]];
                    neiProcs.clear();

                    forAllRow(pAtProcs, he.end(), j)
                    {
                        const label procI = pAtProcs(he.end(), j);

                        if( pAtProcs.contains(he.start(), procI) )
                            neiProcs.append(procI);
                    }

                    //- this hair edge is located at the following processors
                    forAll(neiProcs, j)
                        hairEdgeNeiProcs_.appendIfNotIn(neiProcs[j]);

                    //- the end point of this hair edge is moved by the proc
                    if( (endIt->second == 1) && (neiProcs.size() > 1) )
                        hairEndPointOwnedByProc_.insert(pointI);
                }

                //- prepare data that will be exchanged among processors
                forAllRow(pAtProcs, pointI, i)
                {
                    const label neiProc = pAtProcs(pointI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    labelLongList& dts = exchangeData[neiProc];
                    dts.append(it.key());
                    dts.append(nEndPointsAtPoint_[pointI]);
                }
            }
        }

        labelLongList receiveData;
        help::exchangeMap(exchangeData, receiveData);

        for(label i=0;i<receiveData.size();)
        {
            const label pointI = globalToLocal[receiveData[i++]];
            const label nEndPoints = receiveData[i++];

            if( nEndPointsAtPoint_.find(pointI) == nEndPointsAtPoint_.end() )
                nEndPointsAtPoint_[pointI] = 0;

            label& nEnds = nEndPointsAtPoint_[pointI];
            nEnds = max(nEnds, nEndPoints);
        }
    }

    # ifdef DEBUGLayer
    const label hairEdgesId = mesh_.addPointSubset("hairEdgePoints");
    const label bndHairEdgeId = mesh_.addPointSubset("bndHairEdgePoints");
    const label featureHairEdgeId = mesh_.addPointSubset("featureEdgePoints");
    const label cornerHairEdgeId = mesh_.addPointSubset("cornerHairEdgePoints");
    const label hairEdgeAtEdgeId = mesh_.addPointSubset("hairEdgeAtEdgePoints");

    forAll(hairEdgeType_, heI)
    {
        const edge& e = hairEdges_[heI];

        mesh_.addPointToSubset(hairEdgesId, e.start());
        mesh_.addPointToSubset(hairEdgesId, e.end());
        if( hairEdgeType_[heI] & FEATUREEDGE)
        {
            mesh_.addPointToSubset(featureHairEdgeId, e.start());
            mesh_.addPointToSubset(featureHairEdgeId, e.end());
        }

        if( hairEdgeType_[heI] & BOUNDARY)
        {
            mesh_.addPointToSubset(bndHairEdgeId, e.start());
            mesh_.addPointToSubset(bndHairEdgeId, e.end());
        }

        if( hairEdgeType_[heI] & ATCORNER)
        {
            mesh_.addPointToSubset(cornerHairEdgeId, e.start());
            mesh_.addPointToSubset(cornerHairEdgeId, e.end());
        }

        if( hairEdgeType_[heI] & ATEDGE)
        {
            mesh_.addPointToSubset(hairEdgeAtEdgeId, e.start());
            mesh_.addPointToSubset(hairEdgeAtEdgeId, e.end());
        }
    }

    mesh_.write();
    # endif
}

void boundaryLayerOptimisation::classifyBoundaryPoints()
{
    const meshSurfaceEngine& mse = meshSurface();
    const labelLongList& bPoints = mse.boundaryPoints();
    const labelLongList& bp = mse.bp();

    const labelHashSet& edgePoints = surfacePartitioner().edgePoints();
    const labelHashSet& corners = surfacePartitioner().corners();

    boundaryPointType_.setSize(bPoints.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        //- classify corners and edges
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(bPoints, bpI)
        {
            direction bType(NONE);

            if( edgePoints.found(bpI) )
                bType |= ATEDGE;
            if( corners.found(bpI) )
                bType |= ATCORNER;

            boundaryPointType_[bpI] = bType;
        }

        //- mark points and ends of hair edges
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(hairEdges_, heI)
        {
            if( hairEdgeType_[heI] & BOUNDARY )
            {
                boundaryPointType_[bp[hairEdges_[heI].end()]] |= ENDOFHAIR;
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- reduce classification across all processors
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(mse.bpNeiProcs(), i)
            exchangeData[mse.bpNeiProcs()[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                labelLongList& dts = exchangeData[neiProc];

                //- send the global point index and its type
                dts.append(it.key());
                dts.append(boundaryPointType_[bpI]);
            }
        }

        labelLongList receiveData;
        help::exchangeMap(exchangeData, receiveData);

        for(label i=0;i<receiveData.size();i+=2)
        {
            const label bpI = globalToLocal[receiveData[i]];
            boundaryPointType_[bpI] |= receiveData[i+1];
        }
    }

    //- find and mark floating points at feature edges
    const edgeLongList& edges = mse.edges();
    const labelHashSet& featureEdges = surfacePartitioner().featureEdges();

    forAllConstIter(labelHashSet, featureEdges, it)
    {
        const edge& e = edges[it.key()];

        bool hasCorner(false), hasEndOfHair(false);
        forAll(e, pI)
        {
            if( boundaryPointType_[bp[e[pI]]] & ATCORNER )
                hasCorner = true;
            if( boundaryPointType_[bp[e[pI]]] & ENDOFHAIR )
                hasEndOfHair = true;
        }

        if( hasEndOfHair && !hasCorner )
        {
            forAll(e, pI)
            {
                if( boundaryPointType_[bp[e[pI]]] & ENDOFHAIR )
                    continue;

                boundaryPointType_[bp[e[pI]]] |= FLOATINGATEDGE;
            }
        }
    }

    # ifdef DEBUGLayer
    const label cornerId = mesh_.addPointSubset("cornerPoints");
    const label edgeId = mesh_.addPointSubset("pointsAtEdge");
    const label endOfHairId = mesh_.addPointSubset("endOfHair");
    const label floatingId = mesh_.addPointSubset("floatingPtsAtFeatureEdges");

    forAll(boundaryPointType_, bpI)
    {
        if( boundaryPointType_[bpI] & ATEDGE )
            mesh_.addPointToSubset(edgeId, bPoints[bpI]);
        if( boundaryPointType_[bpI] & ATCORNER )
            mesh_.addPointToSubset(cornerId, bPoints[bpI]);
        if( boundaryPointType_[bpI] & ENDOFHAIR )
            mesh_.addPointToSubset(endOfHairId, bPoints[bpI]);
        if( boundaryPointType_[bpI] & FLOATINGATEDGE )
            mesh_.addPointToSubset(floatingId, bPoints[bpI]);
    }

    mesh_.write();
    throw std::string("finished classification");
    # endif
}

bool boundaryLayerOptimisation::optimiseLayersAtExittingFaces()
{
    bool modified(false);

    //- find the points with more than one hair edge which was modified
    //- in the previous procedure
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(thinnedHairEdge_, heI)
    {
        std::map<label, label>::const_iterator it =
            nEndPointsAtPoint_.find(hairEdges_[heI].end());
        if
        (
            thinnedHairEdge_[heI] &&
            (it != nEndPointsAtPoint_.end()) &&
            (it->second > 2)
        )
        {
            modified = true;
        }
    }

    reduce(modified, maxOp<bool>());

    if( !modified )
        return false;

    Info << "Hair edges at exitting faces shall "
         << "be modified due to inner constraints" << endl;

    return true;
}

void boundaryLayerOptimisation::unifyCoordinatesParallel
(
    const boolList& modifiedEdge
)
{
    if( !Pstream::parRun() )
        return;

    //- this function updates coordinates of end points on processes that do
    //- not modify it
    const pointFieldPMG& points = mesh_.points();
    polyMeshGenModifier meshModifier(mesh_);

    const polyMeshGenAddressing& addr = mesh_.addressingData();
    const Map<label>& globalToLocal = addr.globalToLocalPointAddressing();
    const VRWGraph& pAtProcs = addr.pointAtProcs();
    const DynList<label>& pNeiProcs = addr.pointNeiProcs();

    //- create a map for exchanging of data
    std::map<label, LongList<labelledPoint> > exchangeData;
    forAll(pNeiProcs, i)
        exchangeData[pNeiProcs[i]].clear();

    //- fill the map with data
    forAllConstIter(Map<label>, globalToLocal, it)
    {
        const label pointI = it();

        //- do not continue if this processor does not contain that edge
        std::map<label, DynList<label, 3> >::const_iterator mIt =
            hairEndPointsAtPoint_.find(pointI);

        if( mIt == hairEndPointsAtPoint_.end() )
            continue;

        if
        (
            hairEndPointOwnedByProc_.find(pointI) ==
            hairEndPointOwnedByProc_.end()
        )
            continue;

        //- check if the point has been moved
        const DynList<label, 3>& nHairs = mIt->second;

        bool isMoved(false);
        forAll(nHairs, i)
            isMoved |= modifiedEdge[nHairs[i]];

        if( isMoved )
        {
            //- prepare data for sending
            forAllRow(pAtProcs, pointI, i)
            {
                const label neiProc = pAtProcs(pointI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append
                (
                    labelledPoint(it.key(), points[pointI])
                );
            }
        }
    }

    //- exchange data between processors
    LongList<labelledPoint> receiveData;
    help::exchangeMap(exchangeData, receiveData);

    //- update point coordinates at processors that do not control
    //- the position of the points
    forAll(receiveData, i)
    {
        const labelledPoint& lp = receiveData[i];
        const label pointI = globalToLocal[lp.pointLabel()];

        std::map<label, DynList<label, 3> >::const_iterator it =
            hairEndPointsAtPoint_.find(pointI);

        if( it == hairEndPointsAtPoint_.end() )
            continue;

        //- check if any of the hairs attached to this point
        //- are controlled by this processor
        if
        (
            hairEndPointOwnedByProc_.find(pointI) ==
            hairEndPointOwnedByProc_.end()
        )
        {
            //- update point coordinates
            meshModifier.movePoint(pointI, lp.coordinates());
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayerOptimisation::optimiseLayer()
{
    //- optimise hairs near convex corners
    optimiseHairsNearCorners();

    //- create surface smoother
    meshSurfaceOptimizer surfOpt(meshSurface());

    //- lock exitting faces and feature edges
    {
        labelLongList lockedFaces;
        forAll(isExitFace_, bfI)
            if( isExitFace_[bfI] )
                lockedFaces.append(bfI);
        surfOpt.lockBoundaryFaces(lockedFaces);
        surfOpt.lockFeatureEdges();
    }

    # ifdef DEBUGLayer
    const label lockedFacesId = mesh_.addFaceSubset("lockedBlFaces");
    const label lockedEdgePointsId =
        mesh_.addPointSubset("lockedBoundaryPoints");

    forAll(lockedFaces, i)
        mesh_.addFaceToSubset(lockedFacesId, mesh_.nInternalFaces()+lockedFaces[i]);
    const labelHashSet& edgePts = surfacePartitioner().edgePoints();
    const labelLongList& bPoints = meshSurface().boundaryPoints();
    forAllConstIter(labelHashSet, edgePts, it)
        mesh_.addPointToSubset(lockedEdgePointsId, bPoints[it.key()]);
    # endif

    label nIter(0);
    do
    {
        thinnedHairEdge_ = false;

        //- optimise positions of constrained points
        optimiseConstrainedEdges();

        //- calculate normals at the boundary
        optimiseHairNormalsAtTheBoundary();

        //- optimise positions of constrained points
        optimiseConstrainedEdges();

        //- smoothing thickness variation of boundary hairs
        optimiseThicknessVariation(BOUNDARY);

        //- optimise positions of constrained points
        optimiseConstrainedEdges();

        if( true )
        {
            meshSurfaceEngineModifier bMod(meshSurface());
            bMod.updateGeometry();

            surfOpt.optimizeSurface(2);
            bMod.updateGeometry();
        }

        # ifdef DEBUGLayer
        label counter(0);
        forAll(thinnedHairEdge_, heI)
            if( thinnedHairEdge_[heI] )
                ++counter;
        reduce(counter, sumOp<label>());
        Info << "Thinned " << counter << " bnd hair edges" << endl;
        # endif

        //- optimise normals inside the mesh
        optimiseHairNormalsInside();

        //- optimise thickness variation inside the mesh
        optimiseThicknessVariation(INSIDE);

        # ifdef DEBUGLayer
        label intCounter = 0;
        forAll(thinnedHairEdge_, heI)
            if( thinnedHairEdge_[heI] )
                ++intCounter;
        Info << "Thinned " << (intCounter - counter)
             << " inner hair edges" << endl;
        # endif
    } while( optimiseLayersAtExittingFaces() && (++nIter < maxNumIterations_) );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
