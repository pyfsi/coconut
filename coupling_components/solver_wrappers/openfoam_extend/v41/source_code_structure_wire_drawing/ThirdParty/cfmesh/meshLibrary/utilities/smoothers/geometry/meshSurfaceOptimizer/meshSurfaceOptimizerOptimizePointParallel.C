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
#include "meshSurfaceOptimizer.H"
#include "plane.H"
#include "meshSurfaceMapper.H"
#include "surfaceOptimizer.H"
#include "refLabelledPoint.H"
#include "labelledPointScalar.H"
#include "labelledPoint.H"
#include "helperFunctions.H"

//#define DEBUGSmoothing

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceOptimizer::nodeDisplacementFaceFlatnessParallel
(
    const labelLongList& nodesToSmooth,
    pointField& newPoints,
    const bool weighted
)
{
    if( returnReduce(nodesToSmooth.size(), sumOp<label>()) == 0 )
        return;

    //- create storage for data
    typedef std::map<label, DynList<std::pair<point, scalar> > > mapType;
    mapType localData;

    //- exchange data with other processors
    std::map<label, LongList<labelledPointScalar> > exchangeData;

    forAll(surfaceEngine_.bpNeiProcs(), i)
        exchangeData[surfaceEngine_.bpNeiProcs()[i]].clear();

    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const vectorLongList& fNormals = surfaceEngine_.faceNormals();
    const vectorLongList& fCentres = surfaceEngine_.faceCentres();

    const labelLongList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();

    //- perform smoothing
    labelLongList nodesAtParBnd;
    forAll(nodesToSmooth, pI)
    {
        const label bpI = nodesToSmooth[pI];

        if( !(vertexType_[bpI] & PROCBND) )
            continue;

        nodesAtParBnd.append(pI);

        const point& p = points[bPoints[bpI]];

        //- collect local points
        DynList<std::pair<point, scalar> >& data = localData[bpI];
        data.clear();

        scalar w = 1.0;
        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces(bpI, pfI);

            const point& c = fCentres[bfI];
            const vector& fn = fNormals[bfI];
            const scalar fnSq = magSqr(fn) + VSMALL;

            const vector d = ((c - p) & fn) * fn / fnSq;

            if( weighted )
                w = 1.0 / fnSq;

            const vector disp = w * d;

            data.append(std::make_pair(disp, w));
        }

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            //- add data to the list which will be sent to other processor
            LongList<labelledPointScalar>& dts = exchangeData[neiProc];
            forAll(data, i)
            {
                dts.append
                (
                    labelledPointScalar
                    (
                        globalPointLabel[bpI],
                        data[i].first,
                        data[i].second
                    )
                );
            }
        }
    }

    //- exchange data with other processors
    LongList<labelledPointScalar> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
        const labelledPointScalar& lp = receivedData[i];
        const label bpI = globalToLocal[lp.pointLabel()];

        DynList<std::pair<point, scalar> >& data = localData[bpI];

        data.append(std::make_pair(lp.coordinates(), lp.scalarValue()));
    }

    //- calculate new positions and move the vertices
    forAll(nodesAtParBnd, i)
    {
        const label bpI = nodesToSmooth[nodesAtParBnd[i]];

        if( localData.find(bpI) == localData.end() )
            continue;

        //- create new point position
        const point& p = points[bPoints[bpI]];
        const DynList<std::pair<point, scalar> >& data = localData[bpI];

        point disp(vector::zero);
        scalar sumW(0.0);

        forAll(data, j)
        {
            disp += data[j].first;
            sumW += data[j].second;
        }

        disp /= sumW;

        //- calculate new position of the vertex
        newPoints[nodesAtParBnd[i]] = p + disp;
    }
}

void meshSurfaceOptimizer::nodeDisplacementLaplacianParallel
(
    const labelLongList& nodesToSmooth,
    pointField& newPoints,
    const bool transformIntoPlane,
    const bool weighted
)
{
    if( returnReduce(nodesToSmooth.size(), sumOp<label>()) == 0 )
        return;

    //- create storage for data
    std::map<label, std::map<label, point> > localData;

    //- exchange data with other processors
    std::map<label, LongList<refLabelledPoint> > exchangeData;

    forAll(surfaceEngine_.bpNeiProcs(), i)
        exchangeData[surfaceEngine_.bpNeiProcs()[i]].clear();

    const pointField& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const VRWGraph& pPoints = surfaceEngine_.pointPoints();
    const vectorLongList& pNormals = surfaceEngine_.pointNormals();
    const labelLongList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();

    //- perform smoothing
    labelLongList nodesAtParBnd;

    forAll(nodesToSmooth, pI)
    {
        const label bpI = nodesToSmooth[pI];

        if( !(vertexType_[bpI] & PROCBND) )
            continue;

        nodesAtParBnd.append(pI);

        if( magSqr(pNormals[bpI]) < VSMALL )
            continue;

        //- collect local points
        std::map<label, point>& data = localData[bpI];
        data.clear();
        forAllRow(pPoints, bpI, ppI)
        {
            const label bpJ = pPoints(bpI, ppI);

            data[globalPointLabel[bpJ]] = points[bPoints[bpJ]];
        }

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            //- add data to the list which will be sent to other processor
            LongList<refLabelledPoint>& dts = exchangeData[neiProc];
            for
            (
                std::map<label, point>::const_iterator it=data.begin();
                it!=data.end();
                ++it
            )
            {
                dts.append
                (
                    refLabelledPoint
                    (
                        globalPointLabel[bpI],
                        labelledPoint(globalPointLabel[it->first], it->second)
                    )
                );
            }
        }
    }

    //- exchange data with other processors
    LongList<refLabelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
        const refLabelledPoint& lp = receivedData[i];
        const label bpI = globalToLocal[lp.objectLabel()];

        std::map<label, point>& data = localData[bpI];

        data[lp.lPoint().pointLabel()] = lp.lPoint().coordinates();
    }

    //- calculate new positions and move the vertices
    forAll(nodesAtParBnd, pI)
    {
        const label bpI = nodesToSmooth[nodesAtParBnd[pI]];

        if( localData.find(bpI) == localData.end() )
            continue;

        //- create new point position
        const std::map<label, point>& data = localData[bpI];

        DynList<point> neiPoints;
        for
        (
            std::map<label, point>::const_iterator it=data.begin();
            it!=data.end();
            ++it
        )
            neiPoints.append(it->second);

        const point newP =
            newPositionLaplacian(bpI, neiPoints, transformIntoPlane, weighted);

        newPoints[nodesAtParBnd[pI]] = newP;
    }
}

void meshSurfaceOptimizer::nodeDisplacementLaplacianFCParallel
(
    const labelLongList& nodesToSmooth,
    pointField& newPoints,
    const bool transformIntoPlane,
    const bool weighted
)
{
    if( returnReduce(nodesToSmooth.size(), sumOp<label>()) == 0 )
        return;

    //- create storage for data
    std::map<label, DynList<point> > localData;

    //- exchange data with other processors
    std::map<label, LongList<labelledPoint> > exchangeData;
    forAll(surfaceEngine_.bpNeiProcs(), i)
        exchangeData[surfaceEngine_.bpNeiProcs()[i]].clear();

    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const vectorLongList& faceCentres = surfaceEngine_.faceCentres();
    const vectorLongList& pNormals = surfaceEngine_.pointNormals();

    const labelLongList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();

    //- perform smoothing
    labelLongList nodesAtParBnd;

    forAll(nodesToSmooth, pI)
    {
        const label bpI = nodesToSmooth[pI];

        if( !(vertexType_[bpI] & PROCBND) )
            continue;

        if( magSqr(pNormals[bpI]) < VSMALL )
            continue;

        nodesAtParBnd.append(pI);

        //- collect centres of neighbouring faces
        DynList<point>& data = localData[bpI];
        data.clear();

        forAllRow(pFaces, bpI, pfI)
            data.append(faceCentres[pFaces(bpI, pfI)]);

        //- append the data to the message
        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);
            if( neiProc == Pstream::myProcNo() )
                continue;

            //- add data to the list which will be sent to other processor
            LongList<labelledPoint>& dts = exchangeData[neiProc];
            forAll(data, i)
                dts.append(labelledPoint(globalPointLabel[bpI], data[i]));
        }
    }

    //- exchange data with other processors
    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    forAll(receivedData, i)
    {
        const labelledPoint& lp = receivedData[i];
        const label bpI = globalToLocal[lp.pointLabel()];

        DynList<point>& lpd = localData[bpI];

        lpd.append(lp.coordinates());
    }

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(nodesAtParBnd, pI)
    {
        const label bpI = nodesToSmooth[nodesAtParBnd[pI]];

        if( localData.find(bpI) == localData.end() )
        {
            continue;
        }

        //- create new point position
        const DynList<point>& neiPoints = localData[bpI];
        const point newP =
            newPositionLaplacian(bpI, neiPoints, transformIntoPlane, weighted);

        newPoints[nodesAtParBnd[pI]] = newP;
    }
}

void meshSurfaceOptimizer::edgeNodeDisplacementParallel
(
    const labelLongList& nodesToSmooth,
    pointField& newPoints
)
{
    if( returnReduce(nodesToSmooth.size(), sumOp<label>()) == 0 )
        return;

    # ifdef DEBUGSmoothing
    Info << "Num edge nodes at parallel boundaries "
         << returnReduce(nodesToSmooth.size(), sumOp<label>()) << endl;
    const scalar startEdge = omp_get_wtime();
    # endif

    std::map<label, DynList<labelledPoint, 2> > mPts;

    //- insert labels into the map
    const labelLongList& globalPointLabel =
        surfaceEngine_.globalBoundaryPointLabel();

    labelLongList nodesAtParBnd;

    forAll(nodesToSmooth, nI)
    {
        if( vertexType_[nodesToSmooth[nI]] & PROCBND )
        {
            nodesAtParBnd.append(nI);

            mPts.insert
            (
                std::make_pair
                (
                    globalPointLabel[nodesToSmooth[nI]],
                    DynList<labelledPoint, 2>()
                )
            );
        }
    }

    //- add local edge points
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const edgeLongList& edges = surfaceEngine_.edges();
    const VRWGraph& bpEdges = surfaceEngine_.boundaryPointEdges();
    const  labelLongList& bp = surfaceEngine_.bp();

    const labelHashSet& featureEdges = partitionerPtr_->featureEdges();

    forAll(nodesAtParBnd, nI)
    {
        const label bpI = nodesToSmooth[nodesAtParBnd[nI]];

        DynList<labelledPoint, 2>& neiPoints = mPts[globalPointLabel[bpI]];

        forAllRow(bpEdges, bpI, epI)
        {
            const label beI = bpEdges(bpI, epI);

            if( !featureEdges.found(beI) )
                continue;

            const edge& e = edges[beI];
            const label pI = bp[e.otherVertex(bPoints[bpI])];
            if( vertexType_[pI] & (EDGE+CORNER) )
            {
                neiPoints.append
                (
                    labelledPoint(globalPointLabel[pI], points[bPoints[pI]])
                );
            }
        }
    }

    # ifdef DEBUGSmoothing
    const scalar endPointGathering = omp_get_wtime();
    Info << "MPI nei points " << (endPointGathering - startEdge) << endl;
    # endif

    //- start preparing data which can be exchanged with other processors
    const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();

    std::map<label, LongList<refLabelledPoint> > mProcs;
    forAll(neiProcs, procI)
    {
        const label neiProc = neiProcs[procI];
        mProcs.insert(std::make_pair(neiProc, LongList<refLabelledPoint>()));
    }

    forAll(nodesAtParBnd, i)
    {
        const label bpI = nodesToSmooth[nodesAtParBnd[i]];

        const DynList<labelledPoint, 2>& neiPoints =
            mPts[globalPointLabel[bpI]];

        forAllRow(bpAtProcs, bpI, procI)
        {
            const label neiProc = bpAtProcs(bpI, procI);

            if( neiProc == Pstream::myProcNo() )
                continue;

            forAll(neiPoints, npI)
            {
                LongList<refLabelledPoint>& neiProcPts = mProcs[neiProc];
                neiProcPts.append
                (
                    refLabelledPoint(globalPointLabel[bpI], neiPoints[npI])
                );
            }
        }
    }

    //- exchange data with other processors
    LongList<refLabelledPoint> receivedData;
    help::exchangeMap(mProcs, receivedData);

    forAll(receivedData, prI)
    {
        const refLabelledPoint& lp = receivedData[prI];
        DynList<labelledPoint, 2>& lPts = mPts[lp.objectLabel()];
        lPts.appendIfNotIn(receivedData[prI].lPoint());
    }

    # ifdef DEBUGSmoothing
    const scalar exchangeTime = omp_get_wtime();
    Info << "MPI exchange time " << (exchangeTime - endPointGathering) << endl;
    # endif

    //- Finally, the data is ready to start smoothing
    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(nodesAtParBnd, i)
    {
        const label bpI = nodesToSmooth[nodesAtParBnd[i]];

        const DynList<labelledPoint, 2>& nPts = mPts[globalPointLabel[bpI]];

        const point& p = points[bPoints[bpI]];

        if( nPts.size() != 2 )
        {
            newPoints[nodesAtParBnd[i]] = p;
            continue;
        }

        //- calculate the position of the new point
        const point newP =
            newEdgePositionLaplacian
            (
                nPts[0].coordinates(),
                p,
                nPts[1].coordinates()
            );

        newPoints[nodesAtParBnd[i]] = newP;
    }

    # ifdef DEBUGSmoothing
    const scalar smoothTime = omp_get_wtime();
    Info << "MPI smoothing time " << (smoothTime - exchangeTime) << endl;
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
