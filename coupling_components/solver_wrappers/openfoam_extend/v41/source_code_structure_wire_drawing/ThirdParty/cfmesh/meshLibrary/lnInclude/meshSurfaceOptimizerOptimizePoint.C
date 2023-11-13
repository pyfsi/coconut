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
#include "refLabelledPoint.H"
#include "labelledPointScalar.H"

//#define DEBUGSmoothing

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceOptimizer::smoothEdgePoints
(
    const labelLongList& edgePoints,
    const labelLongList& procEdgePoints,
    const scalar relaxationFactor
)
{
    # ifdef DEBUGSmoothing
    const scalar edgeStart = omp_get_wtime();
    # endif

    pointField newPositions(edgePoints.size());

    const pointFieldPMG& pts = surfaceEngine_.points();
    const labelLongList& bp = surfaceEngine_.bp();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    const edgeLongList& edges = surfaceEngine_.edges();
    const VRWGraph& bpEdges = surfaceEngine_.boundaryPointEdges();
    const pointFieldPMG& points = surfaceEngine_.points();

    const labelHashSet& featureEdges = partitionerPtr_->featureEdges();

    # ifdef DEBUGSmoothing
    const scalar endAddressing = omp_get_wtime();
    Info << "Edge smoothing addressing " << (endAddressing - edgeStart) << endl;
    # endif

    //- smooth edge vertices
    labelLongList procPositions;
    std::map<label, label> bpToPos;

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 40)
        # endif
        forAll(edgePoints, i)
        {
            const label bpI = edgePoints[i];

            if( vertexType_[bpI] & PROCBND )
            {
                # ifdef USE_OMP
                # pragma omp critical(procPositions)
                # endif
                {
                    bpToPos[bpI] = i;
                    procPositions.append(i);
                }

                continue;
            }

            newPositions[i] = newEdgePositionLaplacian(bpI);
        }

        if( Pstream::parRun() )
        {
            # ifdef USE_OMP
            # pragma omp single
            # endif
            {
                const labelLongList& globalPointLabel =
                    surfaceEngine_.globalBoundaryPointLabel();
                const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
                const Map<label>& globalToLocal =
                    surfaceEngine_.globalToLocalBndPointAddressing();

                std::map<label, LongList<refLabelledPoint> > exchangeData;
                forAll(surfaceEngine_.bpNeiProcs(), i)
                    exchangeData[surfaceEngine_.bpNeiProcs()[i]].clear();

                typedef std::map<label, DynList<labelledPoint, 2> > localMap;
                localMap ePts;

                //- prepare a message
                forAll(procPositions, i)
                {
                    const label bpI = edgePoints[procPositions[i]];

                    ePts[bpI].clear();

                    forAllRow(bpEdges, bpI, j)
                    {
                        const label beI = bpEdges(bpI, j);

                        if( !featureEdges.found(beI) )
                            continue;

                        const edge& e = edges[bpEdges(bpI, j)];

                        const label otherPointI = e.otherVertex(bPoints[bpI]);
                        const label otherBpI = bp[otherPointI];

                        ePts[bpI].appendIfNotIn
                        (
                            labelledPoint
                            (
                                globalPointLabel[otherBpI],
                                points[otherPointI]
                            )
                        );

                        forAllRow(bpAtProcs, bpI, nProcI)
                        {
                            const label neiProc = bpAtProcs(bpI, nProcI);

                            if( neiProc == Pstream::myProcNo() )
                                continue;

                            exchangeData[neiProc].append
                            (
                                refLabelledPoint
                                (
                                    globalPointLabel[bpI],
                                    labelledPoint
                                    (
                                        globalPointLabel[otherBpI],
                                        points[otherPointI]
                                    )
                                )
                            );
                        }
                    }
                }

                //- exchange data among processors
                LongList<refLabelledPoint> receiveData;
                help::exchangeMap(exchangeData, receiveData);

                //- combine the data into a local map
                forAll(receiveData, i)
                {
                    const refLabelledPoint& rlp = receiveData[i];

                    const label bpI = globalToLocal[rlp.objectLabel()];
                    const labelledPoint& lp = rlp.lPoint();

                    ePts[bpI].appendIfNotIn(lp);
                }

                //- move the vertices
                forAllConstIter(localMap, ePts, it)
                {
                    const label bpI = it->first;
                    const DynList<labelledPoint, 2>& nPts = it->second;

                    if( bpToPos.find(bpI) == bpToPos.end() )
                        FatalErrorIn
                        (
                            "void meshSurfaceOptimizer::smoothEdgePoints("
                            "const labelLongList&, const labelLongList&,"
                            " const scalar)"
                        ) << "Map entry " << bpI << " does not exist!"
                          << abort(FatalError);

                    if( nPts.size() == 2 )
                    {
                        //- move the point to a new location
                        newPositions[bpToPos[bpI]] =
                            newEdgePositionLaplacian
                            (
                                nPts[0].coordinates(),
                                points[bPoints[bpI]],
                                nPts[1].coordinates()
                            );
                    }
                    else
                    {
                        //- do not move the point
                        newPositions[bpToPos[bpI]] = points[bPoints[bpI]];
                    }
                }

                # ifdef DEBUGSmoothing
                for(label procI=0;procI<Pstream::nProcs();++procI)
                {
                    if( Pstream::myProcNo() == procI )
                    {
                        forAll(newPositions, i)
                        {
                            const label bpI = edgePoints[i];
                            const bool isProc = vertexType_[bpI] & PROCBND;
                            Pout << "Edge point " << bpI << " coordinates "
                                 << pts[bPoints[bpI]]
                                 << " new coordinates " << newPositions[i]
                                 << " proc bnd " << isProc << endl;
                        }
                    }

                    returnReduce(1, sumOp<label>());
                }
                # endif
            }
        }

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(newPositions, i)
        {
            const label bpI = edgePoints[i];

            if( vertexType_[bpI] & LOCKED )
                continue;

            const point& np = newPositions[i];

            if( help::isnan(np) )
                continue;
            if( help::isinf(np) )
                continue;

            const point& orig = pts[bPoints[bpI]];
            const point newP = orig + relaxationFactor * (np - orig);

            surfaceModifier_.moveBoundaryVertexNoUpdate(bpI, newP);
        }
    }

    # ifdef DEBUGSmoothing
    const scalar endSmoothing = omp_get_wtime();
    Info << "Edge smoothing: End smoothing "
         << (endSmoothing - endAddressing) << endl;
    # endif

    //- sync vertices at inter-processor boundaries
    surfaceModifier_.syncVerticesAtParallelBoundaries(procEdgePoints);

    //- update geometric data
    surfaceModifier_.updateGeometry(edgePoints);

    # ifdef DEBUGSmoothing
    Info << "Edge smoothing : Finished updating geometry "
         << (omp_get_wtime() - endSmoothing) << endl;
    # endif
}

void meshSurfaceOptimizer::smoothFaceFlatness
(
    const labelLongList& selectedPoints,
    const labelLongList& selectedProcPoints,
    const scalar relaxationFactor,
    const bool weighted
)
{
    # ifdef DEBUGSmoothing
    const scalar startFlatness = omp_get_wtime();
    Info << "Number of selected points for the flatness smoother "
         << returnReduce(selectedPoints.size(), sumOp<label>()) << endl;
    # endif

    pointField newPositions(selectedPoints.size());

    const pointFieldPMG& pts = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    labelLongList procPositions;
    std::map<label, label> bpToPos;

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 40)
        # endif
        forAll(selectedPoints, i)
        {
            const label bpI = selectedPoints[i];

            if( vertexType_[bpI] & PROCBND )
            {
                # ifdef USE_OMP
                # pragma omp critical(procPositions)
                # endif
                {
                    procPositions.append(i);
                    bpToPos[bpI] = i;
                }

                continue;
            }

            newPositions[i] = newPositionFaceFlatness(bpI, weighted);
        }

        if( Pstream::parRun() )
        {
            # ifdef USE_OMP
            # pragma omp single
            # endif
            {
                const labelLongList& globalPointLabel =
                    surfaceEngine_.globalBoundaryPointLabel();
                const Map<label>& globalToLocal =
                    surfaceEngine_.globalToLocalBndPointAddressing();
                const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();

                std::map<label, LongList<labelledPointScalar> > exchangeData;
                forAll(surfaceEngine_.bpNeiProcs(), i)
                    exchangeData[surfaceEngine_.bpNeiProcs()[i]].clear();

                typedef std::map<label, std::pair<point, scalar> > localMap;
                localMap pFcs;

                //- prepare the message
                forAll(procPositions, i)
                {
                    const label bpI = selectedPoints[procPositions[i]];

                    //- calculate locally available information
                    vector disp(vector::zero);
                    scalar sumW(0.0);

                    localDataFaceFlatness(bpI, disp, sumW, weighted);

                    pFcs[bpI].first = disp;
                    pFcs[bpI].second = sumW;

                    forAllRow(bpAtProcs, bpI, j)
                    {
                        const label neiProc = bpAtProcs(bpI, j);

                        if( neiProc == Pstream::myProcNo() )
                            continue;

                        exchangeData[neiProc].append
                        (
                            labelledPointScalar
                            (
                                globalPointLabel[bpI],
                                disp,
                                sumW
                            )
                        );
                    }
                }

                //- exchange data among processors
                LongList<labelledPointScalar> receiveData;
                help::exchangeMap(exchangeData, receiveData);

                //- take the received data into account
                forAll(receiveData, i)
                {
                    const label bpI =
                        globalToLocal[receiveData[i].pointLabel()];

                    pFcs[bpI].first += receiveData[i].coordinates();
                    pFcs[bpI].second += receiveData[i].scalarValue();
                }

                //- calculate new positions of points
                forAllConstIter(localMap, pFcs, it)
                {
                    const label bpI = it->first;
                    const std::pair<vector, scalar>& vs = it->second;

                    newPositions[bpToPos[bpI]] =
                        newPositionFaceFlatness(bpI, vs.first, vs.second);
                }

                # ifdef DEBUGSmoothing
                for(label procI=0;procI<Pstream::nProcs();++procI)
                {
                    if( Pstream::myProcNo() == procI )
                    {
                        forAll(newPositions, i)
                        {
                            const label bpI = selectedPoints[i];
                            const bool isProc = vertexType_[bpI] & PROCBND;
                            Pout << "Flatness: point " << bpI << " coordinates "
                                 << pts[bPoints[bpI]]
                                 << " new coordinates " << newPositions[i]
                                 << " proc bnd " << isProc << endl;
                        }
                    }

                    returnReduce(1, sumOp<label>());
                }
                # endif
            }
        }

        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            //- project points back at the original surface
            projectPointsBackOnTheSurface(selectedPoints, newPositions);
        }

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(newPositions, i)
        {
            const label bpI = selectedPoints[i];

            if( vertexType_[bpI] & LOCKED )
                continue;

            const point& np = newPositions[i];

            if( help::isnan(np) )
                continue;
            if( help::isinf(np) )
                continue;

            const point& orig = pts[bPoints[bpI]];
            const point newP = orig + relaxationFactor * (np - orig);

            surfaceModifier_.moveBoundaryVertexNoUpdate(bpI, newP);
        }
    }

    # ifdef DEBUGSmoothing
    const scalar endFlatness = omp_get_wtime();
    Info << "Flatness: smoothing time "
         << (endFlatness - startFlatness) << endl;
    # endif

    //- ensure that vertices at inter-processor boundaries are at the same
    //- location at all processors
    surfaceModifier_.syncVerticesAtParallelBoundaries(selectedProcPoints);

    # ifdef DEBUGSmoothing
    const scalar endSync = omp_get_wtime();
    Info << "Flatness: Syncing time " << (endSync - endFlatness) << endl;
    # endif

    //- update geometric data
    surfaceModifier_.updateGeometry(selectedPoints);

    # ifdef DEBUGSmoothing
    Info << "Flatness: finished updating geometry "
         << (omp_get_wtime() - endSync) << endl;
    # endif
}

void meshSurfaceOptimizer::smoothLaplacianFC
(
    const labelLongList& selectedPoints,
    const labelLongList& selectedProcPoints,
    const scalar relaxationFactor,
    const bool transform,
    const bool weighted
)
{
    # ifdef DEBUGSmoothing
    Pout << "Starting laplacian smoothing" << endl;
    const scalar startLaplacian = omp_get_wtime();
    # endif

    pointField newPositions(selectedPoints.size());

    const pointFieldPMG& pts = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();
    const vectorLongList& faceCentres = surfaceEngine_.faceCentres();

    labelLongList procPositions;
    std::map<label, label> bpToPos;

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 40)
        # endif
        forAll(selectedPoints, i)
        {
            const label bpI = selectedPoints[i];

            if( vertexType_[bpI] & PROCBND )
            {
                # ifdef USE_OMP
                # pragma omp critical(procPositions)
                # endif
                {
                    procPositions.append(i);
                    bpToPos[bpI] = i;
                }

                continue;
            }

            DynList<point> neiPoints;
            forAllRow(pointFaces, bpI, pfI)
                neiPoints.append(faceCentres[pointFaces(bpI, pfI)]);

            newPositions[i] =
                newPositionLaplacian
                (
                    bpI,
                    neiPoints,
                    transform,
                    weighted
                );
        }

        if( Pstream::parRun() )
        {
            # ifdef USE_OMP
            # pragma omp single
            # endif
            {
                const labelLongList& globalPointLabel =
                    surfaceEngine_.globalBoundaryPointLabel();
                const Map<label>& globalToLocal =
                    surfaceEngine_.globalToLocalBndPointAddressing();
                const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();

                std::map<label, LongList<labelledPoint> > exchangeData;
                forAll(surfaceEngine_.bpNeiProcs(), i)
                    exchangeData[surfaceEngine_.bpNeiProcs()[i]].clear();

                typedef std::map<label, DynList<point> > localMap;
                localMap pFcs;

                //- create messages
                forAll(procPositions, i)
                {
                    const label bpI = selectedPoints[procPositions[i]];

                    DynList<point>& neiPts = pFcs[bpI];
                    neiPts.clear();

                    forAllRow(pointFaces, bpI, pfI)
                        neiPts.append(faceCentres[pointFaces(bpI, pfI)]);

                    forAllRow(bpAtProcs, bpI, npI)
                    {
                        const label neiProc = bpAtProcs(bpI, npI);

                        if( neiProc == Pstream::myProcNo() )
                            continue;

                        forAll(neiPts, npI)
                            exchangeData[neiProc].append
                            (
                                labelledPoint
                                (
                                    globalPointLabel[bpI],
                                    neiPts[npI]
                                )
                            );
                    }
                }

                //- exchange data
                LongList<labelledPoint> receiveData;
                help::exchangeMap(exchangeData, receiveData);

                //- merge data
                forAll(receiveData, i)
                    pFcs[globalToLocal[receiveData[i].pointLabel()]].append
                    (
                        receiveData[i].coordinates()
                    );

                //- calculate new coordinates
                forAllConstIter(localMap, pFcs, it)
                {
                    const label bpI = it->first;

                    newPositions[bpToPos[bpI]] =
                        newPositionLaplacian
                        (
                            bpI,
                            it->second,
                            transform,
                            weighted
                        );
                }

                # ifdef DEBUGSmoothing
                for(label procI=0;procI<Pstream::nProcs();++procI)
                {
                    if( Pstream::myProcNo() == procI )
                    {
                        forAll(newPositions, i)
                        {
                            const label bpI = selectedPoints[i];
                            const bool isProc = vertexType_[bpI] & PROCBND;
                            Pout << "LaplacianFC point " << bpI << " coordinates "
                                 << pts[bPoints[bpI]]
                                 << " new coordinates " << newPositions[i]
                                 << " proc bnd " << isProc << endl;
                        }
                    }

                    returnReduce(1, sumOp<label>());
                }
                # endif
            }
        }

        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            //- project points back at the original surface
            projectPointsBackOnTheSurface(selectedPoints, newPositions);
        }

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(newPositions, i)
        {
            const label bpI = selectedPoints[i];

            if( vertexType_[bpI] & LOCKED )
                continue;

            const point& np = newPositions[i];

            if( help::isnan(np) )
                continue;
            if( help::isinf(np) )
                continue;

            const point& orig = pts[bPoints[bpI]];
            const point newP = orig + relaxationFactor * (np - orig);

            surfaceModifier_.moveBoundaryVertexNoUpdate(bpI, newP);
        }
    }

    # ifdef DEBUGSmoothing
    const scalar endLaplacian = omp_get_wtime();
    Info << "Laplacian: smoothing time "
         << (endLaplacian - startLaplacian) << endl;
    # endif

    //- ensure that vertices at inter-processor boundaries are at the same
    //- location at all processors
    surfaceModifier_.syncVerticesAtParallelBoundaries(selectedProcPoints);

    # ifdef DEBUGSmoothing
    const scalar endSyncing = omp_get_wtime();
    Info << "Laplacian: Syncing time " << (endSyncing - endLaplacian) << endl;
    # endif

    //- update geometric data
    surfaceModifier_.updateGeometry(selectedPoints);

    # ifdef DEBUGSmoothing
    Info << "Laplacian: finished updating geometry "
         << (omp_get_wtime() - endSyncing) << endl;
    # endif
}

void meshSurfaceOptimizer::smoothSurfaceOptimizer
(
    const labelLongList& selectedPoints,
    const labelLongList& selectedProcPoints,
    const scalar relaxationFactor
)
{
    # ifdef DEBUGSmoothing
    Pout << "Starting surface optimizer" << endl;
    const scalar startOpt = omp_get_wtime();
    # endif

    //- create partTriMesh is it is not yet present
    this->triMesh();

    //- update coordinates of the triangulation
    updateTriMesh(selectedPoints);

    # ifdef DEBUGSmoothing
    const scalar startSmoothing = omp_get_wtime();
    Info << "Surf opt : tri mesh preparation "
         << (startSmoothing - startOpt) << endl;
    # endif

    pointField newPositions(selectedPoints.size());

    const pointFieldPMG& pts = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 20)
        # endif
        forAll(selectedPoints, i)
        {
            const label bpI = selectedPoints[i];

            newPositions[i] = newPositionSurfaceOptimizer(bpI);
        }

        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            //- project points back at the original surface
            projectPointsBackOnTheSurface(selectedPoints, newPositions);
        }

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(newPositions, i)
        {
            const label bpI = selectedPoints[i];

            const point& np = newPositions[i];

            if( help::isnan(np) || help::isinf(np) )
                continue;

            const point& orig = pts[bPoints[selectedPoints[i]]];
            const point newP = orig + relaxationFactor * (np - orig);

            surfaceModifier_.moveBoundaryVertexNoUpdate(bpI, newP);
        }
    }

    # ifdef DEBUGSmoothing
    const scalar endSurfOpt = omp_get_wtime();
    Info << "Surf opt: Smoothing time "
         << (endSurfOpt - startSmoothing) << endl;
    # endif

    //- ensure that vertices at inter-processor boundaries are at the same
    //- location at all processors
    surfaceModifier_.syncVerticesAtParallelBoundaries(selectedProcPoints);

    # ifdef DEBUGSmoothing
    const scalar endSync = omp_get_wtime();
    Info << "Surf optSyncing time " << (endSync - endSurfOpt) << endl;
    # endif

    //- update geometry addressing for moved points
    surfaceModifier_.updateGeometry(selectedPoints);

    # ifdef DEBUGSmoothing
    Info << "Surf opt: finished updating geometry "
         << (omp_get_wtime() - endSync) << endl;
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
