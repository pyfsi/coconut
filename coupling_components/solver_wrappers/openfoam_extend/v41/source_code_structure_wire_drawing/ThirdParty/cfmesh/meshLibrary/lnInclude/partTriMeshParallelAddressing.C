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
#include "meshSurfacePartitioner.H"
#include "partTriMesh.H"
#include "triSurfModifier.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "parTriFace.H"

#include <map>

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void partTriMesh::createParallelAddressing()
{
    const meshSurfaceEngine& mse = mPart_.surfaceEngine();

    const pointField& pts = surf_.points();

    //- allocate global point labels
    if( !globalPointLabelPtr_ )
        globalPointLabelPtr_ = new labelLongList();
    labelLongList& globalPointLabel = *globalPointLabelPtr_;
    globalPointLabel.setSize(pts.size());
    globalPointLabel = -1;

    //- allocated point-processors addressing
    if( !pAtProcsPtr_ )
        pAtProcsPtr_ = new VRWGraph();
    VRWGraph& pProcs = *pAtProcsPtr_;
    pProcs.setSize(0);
    pProcs.setSize(pts.size());

    //- allocate global-to-local point addressing
    if( !globalToLocalPointAddressingPtr_ )
        globalToLocalPointAddressingPtr_ = new Map<label>();
    Map<label>& globalToLocal = *globalToLocalPointAddressingPtr_;
    globalToLocal.clear();

    //- allocate storage for points at parallel boundaries
    if( !pAtParallelBoundariesPtr_ )
        pAtParallelBoundariesPtr_ = new labelLongList();
    labelLongList& pAtParallelBoundaries = *pAtParallelBoundariesPtr_;
    pAtParallelBoundaries.clear();

    //- create point-processors addressing
    std::map<label, labelLongList> exchangeData;
    std::map<label, labelLongList>::iterator iter;

    const Map<label>& globalToLocalPointAddressing =
        mse.globalToLocalBndPointAddressing();
    const VRWGraph& pAtProcs = mse.bpAtProcs();
    const DynList<label>& pNeiProcs = mse.bpNeiProcs();

    forAll(pNeiProcs, procI)
        exchangeData.insert(std::make_pair(pNeiProcs[procI], labelLongList()));

    //- exchange data such that all points at inter-processor boundaries
    //- have the same classification
    forAllConstIter(Map<label>, globalToLocalPointAddressing, it)
    {
        const label bpI = it();

        const label tpI = meshSurfacePointLabelInTriMesh_[bpI];

        if( tpI < 0 )
            continue;

        forAllRow(pAtProcs, bpI, procI)
        {
            const label neiProc = pAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            //- send the global label in the mesh surface
            //- and the classification in the partTriMesh
            labelLongList& dts = exchangeData[neiProc];
            dts.append(it.key());
            dts.append(pointType_[tpI]);
        }
    }

    //- exchange data with other processors
    labelLongList receivedData;
    help::exchangeMap(exchangeData, receivedData);

    //- set the values according to other processors
    std::map<label, label> addPoints;
    const label nPoints = pts.size();
    for(label i=0;i<receivedData.size();)
    {
        const label bpI = globalToLocalPointAddressing[receivedData[i++]];
        const direction pType = receivedData[i++];

        if( meshSurfacePointLabelInTriMesh_[bpI] == -1 )
        {
            //- points does not yet exist in the partTriMesh
            addPoints[pointType_.size()] = bpI;
            meshSurfacePointLabelInTriMesh_[bpI] = pointType_.size();
            pointLabelInMeshSurface_.append(bpI);
            pProcs.appendList(DynList<label>());
            globalPointLabel.append(-1);

            if( pType & BOUNDARY )
            {
                pointType_.append(BOUNDARY);
            }
            else
            {
                pointType_.append(pType);
            }
        }
        else
        {
            if( pType & BOUNDARY )
            {
                pointType_[meshSurfacePointLabelInTriMesh_[bpI]] |= BOUNDARY;
            }
            else
            {
                pointType_[meshSurfacePointLabelInTriMesh_[bpI]] |= pType;
            }
        }
    }

    if( addPoints.size() )
    {
        //- add new points into the partTriMesh
        const pointFieldPMG& points = mse.points();
        const labelLongList& bPoints = mse.boundaryPoints();

        pointField& triPoints = triSurfModifier(surf_).pointsAccess();
        triPoints.setSize(nPoints+addPoints.size());

        for
        (
            std::map<label, label>::const_iterator it=addPoints.begin();
            it!=addPoints.end();
            ++it
        )
            triPoints[it->first] = points[bPoints[it->second]];
    }

    for(iter=exchangeData.begin();iter!=exchangeData.end();++iter)
        iter->second.clear();

    //- start creating global-to-local addressing
    //- find the starting point labels
    label startPoint(0), nLocalPoints(0), nSharedPoints(0);

    //- count the number of points at processor boundaries
    forAllConstIter(Map<label>, globalToLocalPointAddressing, it)
    {
        const label bpI = it();

        const label tpI = meshSurfacePointLabelInTriMesh_[bpI];

        if( tpI == -1 )
            continue;

        ++nSharedPoints;

        label pMin(Pstream::myProcNo());
        forAllRow(pAtProcs, bpI, procI)
            pMin = Foam::min(pMin, pAtProcs(bpI, procI));

        if( pMin == Pstream::myProcNo() )
            ++nLocalPoints;
    }

    labelList nPointsAtProc(Pstream::nProcs());
    nSharedPoints -= nLocalPoints;
    nPointsAtProc[Pstream::myProcNo()] = pts.size() - nSharedPoints;
    Pstream::gatherList(nPointsAtProc);
    Pstream::scatterList(nPointsAtProc);

    for(label i=0;i<Pstream::myProcNo();++i)
        startPoint += nPointsAtProc[i];

    //- create global labels for points at processor boundaries
    forAllConstIter(Map<label>, globalToLocalPointAddressing, it)
    {
        const label bpI = it();

        const label tpI = meshSurfacePointLabelInTriMesh_[bpI];

        if( tpI == -1 )
            continue;

        label pMin(Pstream::myProcNo());
        forAllRow(pAtProcs, bpI, procI)
        {
            const label neiProc = pAtProcs(bpI, procI);
            pProcs.append(tpI, neiProc);
            pMin = Foam::min(pMin, neiProc);
        }

        if( pMin != Pstream::myProcNo() )
            continue;

        globalPointLabel[tpI] = startPoint++;

        forAllRow(pAtProcs, bpI, procI)
        {
            const label neiProc = pAtProcs(bpI, procI);

            if( neiProc == Pstream::myProcNo() )
                continue;

            //- the following information is sent to other processor
            //- 1. global point label in the original mesh
            //- 2. global point label in the tet mesh
            labelLongList& dts = exchangeData[neiProc];
            dts.append(it.key());
            dts.append(globalPointLabel[tpI]);
        }
    }

    //- exchange data with other processors
    receivedData.clear();
    help::exchangeMap(exchangeData, receivedData);

    label counter(0);
    while( counter < receivedData.size() )
    {
        const label gpI = receivedData[counter++];
        const label tgI = receivedData[counter++];
        const label pLabel =
            meshSurfacePointLabelInTriMesh_[globalToLocalPointAddressing[gpI]];

        globalPointLabel[pLabel] = tgI;
    }

    //- set global labels for remaining points
    forAll(globalPointLabel, pI)
    {
        if( globalPointLabel[pI] == -1 )
            globalPointLabel[pI] = startPoint++;
    }

    //- create global to local mapping
    forAll(globalPointLabel, pI)
    {
        if( pProcs.sizeOfRow(pI) != 0 )
        {
            pAtParallelBoundaries.append(pI);
            globalToLocal.insert(globalPointLabel[pI], pI);
        }
    }

    //- create neighbour processors addressing
    if( !neiProcsPtr_ )
        neiProcsPtr_ = new DynList<label>();
    DynList<label>& neiProcs = *neiProcsPtr_;

    for(iter=exchangeData.begin();iter!=exchangeData.end();++iter)
        neiProcs.append(iter->first);
}

void partTriMesh::createBufferLayers()
{
    const VRWGraph& pTriangles = surf_.pointFacets();
    pointField& pts = triSurfModifier(surf_).pointsAccess();

    //- point at processors addressing
    const VRWGraph& pProcs = *pAtProcsPtr_;
    labelLongList& globalPointLabel = *globalPointLabelPtr_;
    Map<label>& globalToLocal = *globalToLocalPointAddressingPtr_;
    const DynList<label>& neiProcs = *this->neiProcsPtr_;

    //- mark points in buffer layers
    if( !pAtBufferLayersPtr_ )
        pAtBufferLayersPtr_ = new labelLongList();
    labelLongList& pAtBufferLayers = *pAtBufferLayersPtr_;
    pAtBufferLayers.clear();

    //- buffer points at processors addressing
    if( !updateBufferLayerAtProcsPtr_ )
        updateBufferLayerAtProcsPtr_ =
            new Map<DynList<label, 4> >(globalToLocal.size());
    Map<DynList<label, 4> >& updateBufferPoints = *updateBufferLayerAtProcsPtr_;
    updateBufferPoints.clear();

    //- create the maps used for exchanging data between processors
    std::map<label, std::set<label> > exchangePointsHelper;
    std::map<label, LongList<labelledPoint> > exchangePoints;
    std::map<label, std::set<label> > exchangeTriasHelper;
    std::map<label, labelLongList> exchangeTrias;
    forAll(neiProcs, procI)
    {
        exchangePointsHelper[neiProcs[procI]].clear();
        exchangeTriasHelper[neiProcs[procI]].clear();
        exchangeTrias[neiProcs[procI]].clear();
        exchangePoints[neiProcs[procI]].clear();
    }

    //- loop over triangles and add the ones having vertices at parallel
    //- boundaries for sending
    forAllConstIter(Map<label>, globalToLocal, it)
    {
        const label pI = it();

        forAllRow(pProcs, pI, i)
        {
            const label neiProc = pProcs(pI, i);

            if( neiProc == Pstream::myProcNo() )
                continue;

            //- check if the index is within the scope
            if( pI >= pTriangles.size() )
                continue;

            forAllRow(pTriangles, pI, ptI)
            {
                const label triI = pTriangles(pI, ptI);

                exchangeTriasHelper[neiProc].insert(triI);

                const labelledTri& tri = surf_[triI];

                forAll(tri, j)
                {
                    if( !updateBufferPoints.found(tri[j]) )
                    {
                        //- add point to the list of buffer points
                        pAtBufferLayers.append(tri[j]);
                        updateBufferPoints.insert(tri[j], DynList<label, 4>());
                    }

                    if( !pProcs.contains(tri[j], neiProc) )
                    {
                        exchangePointsHelper[neiProc].insert(tri[j]);
                        updateBufferPoints[tri[j]].appendIfNotIn(neiProc);
                    }
                }
            }
        }
    }

    //- prepare points for sending
    for
    (
        std::map<label, std::set<label> >::iterator it =
            exchangePointsHelper.begin();
        it!=exchangePointsHelper.end();
        ++it
    )
    {
        LongList<labelledPoint>& ptsToSend = exchangePoints[it->first];

        ptsToSend.setSize(it->second.size());

        label counter(0);
        forAllConstIter(std::set<label>, it->second, pIt)
            ptsToSend[counter++] =
                labelledPoint(globalPointLabel[*pIt], pts[*pIt]);

        it->second.clear();
    }

    LongList<labelledPoint> receivedPoints;
    help::exchangeMap(exchangePoints, receivedPoints);

    Map<label> newGlobalToLocal(globalToLocal.size());
    std::map<label, point> addCoordinates;
    label nPoints = pts.size();
    forAll(receivedPoints, i)
    {
        const labelledPoint& lp = receivedPoints[i];

        if( !newGlobalToLocal.found(lp.pointLabel()) )
        {
            pointLabelInMeshSurface_.append(-1);
            pointType_.append(BUFFER);

            globalPointLabel.append(lp.pointLabel());
            newGlobalToLocal.insert(lp.pointLabel(), nPoints);
            addCoordinates[nPoints++] = lp.coordinates();
        }
    }
    exchangePoints.clear();

    //- store newly added points
    pts.setSize(nPoints);
    for
    (
        std::map<label, point>::const_iterator it=addCoordinates.begin();
        it!=addCoordinates.end();
        ++it
    )
        pts[it->first] = it->second;

    addCoordinates.clear();

    //- insert the global labels of the buffer points
    //- into the globalToLocal map
    forAllConstIter(Map<label>, newGlobalToLocal, it)
        globalToLocal.insert(it.key(), it());

    //- Finally, exchange triangles such that all processors
    //- have the same simplices at inter-processor boundaries
    //- prepare triangles for sending
    for
    (
        std::map<label, std::set<label> >::iterator it =
            exchangeTriasHelper.begin();
        it!=exchangeTriasHelper.end();
        ++it
    )
    {
        labelLongList& triasToSend = exchangeTrias[it->first];

        triasToSend.setSize(3*it->second.size());

        label counter(0);
        forAllConstIter(std::set<label>, it->second, pIt)
        {
            const labelledTri& tri = surf_[*pIt];
            triasToSend[counter++] = globalPointLabel[tri[0]];
            triasToSend[counter++] = globalPointLabel[tri[1]];
            triasToSend[counter++] = globalPointLabel[tri[2]];
        }

        it->second.clear();
    }

    //- receive triangles sent to this processor
    labelLongList receivedTriangles;
    help::exchangeMap(exchangeTrias, receivedTriangles);
    exchangeTrias.clear();

    //- add triangles into the mesh and update the addressing
    label i(0);
    while( i < receivedTriangles.size() )
    {
        const label p0 = globalToLocal[receivedTriangles[i++]];
        const label p1 = globalToLocal[receivedTriangles[i++]];
        const label p2 = globalToLocal[receivedTriangles[i++]];

        //- append tet
        surf_.appendTriangle(labelledTri(p0, p1, p2, -1));
    }

    //- update addressing of the surface mesh
    surf_.clearAddressing();
}

void partTriMesh::updateBufferLayers()
{
    const labelLongList& bufferLayerPoints = this->bufferLayerPoints();

    updateBufferLayers(bufferLayerPoints);
}

void partTriMesh::updateBufferLayers(const labelLongList& movedPoints)
{
    if( returnReduce(movedPoints.size(), sumOp<label>()) == 0 )
        return;

    const pointField& points = surf_.points();
    const Map<DynList<label, 4> >& bufferLayerPointsAtProcs =
        this->updateBufferLayerPointsAtProcs();
    const labelLongList& globalPointLabel = this->globalPointLabel();
    const Map<label>& globalToLocal = this->globalToLocalPointAddressing();
    const DynList<label>& neiProcs = this->neiProcs();

    //- create the map
    std::map<label, LongList<labelledPoint> > exchangeData;
    forAll(neiProcs, i)
        exchangeData.insert
        (
            std::make_pair(neiProcs[i], LongList<labelledPoint>())
        );

    //- add points into the map
    forAll(movedPoints, pI)
    {
        const label pointI = movedPoints[pI];

        const DynList<label, 4>& updateAtProcs =
            bufferLayerPointsAtProcs[pointI];

        forAll(updateAtProcs, i)
        {
            const label neiProc = updateAtProcs[i];

            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeData[neiProc].append
            (
                labelledPoint(globalPointLabel[pointI], points[pointI])
            );
        }
    }

    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
        const labelledPoint& lp = receivedData[i];

        this->updateVertex
        (
            globalToLocal[lp.pointLabel()],
            lp.coordinates()
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
