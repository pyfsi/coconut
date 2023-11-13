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
#include "meshSurfaceCheckInvertedVertices.H"
#include "meshOctree.H"
#include "helperFunctions.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceMapper2D.H"
#include "polyMeshGen2DEngine.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGenChecks.H"
#include "labelledPoint.H"
#include "FIFOStack.H"

#include <map>
#include <stdexcept>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSmooth

//#define smoothingTiming

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label meshSurfaceOptimizer::findInvertedVertices
(
    boolList& smoothVertex,
    const label nAdditionalLayers
) const
{
    # ifdef smoothingTiming
    const scalar startFindingInverted = omp_get_wtime();
    # endif

    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
    const labelLongList& bp = surfaceEngine_.bp();

    const VRWGraph* bpAtProcsPtr = NULL;
    if( Pstream::parRun() )
        bpAtProcsPtr = &surfaceEngine_.bpAtProcs();

    if( smoothVertex.size() != bPoints.size() )
    {
        smoothVertex.setSize(bPoints.size());
        # ifdef USE_OMP
        # pragma omp parallel for schedule(static, 1)
        # endif
        forAll(smoothVertex, i)
            smoothVertex[i] = true;
    }

    label nInvertedTria(0);

    # ifdef smoothingTiming
    const scalar startChecking = omp_get_wtime();
    Info << "Finding inverted: preparation time "
         << (startChecking - startFindingInverted) << endl;
    # endif

    //- check the vertices at the surface
    //- mark the ones where the mesh is tangled
    meshSurfaceCheckInvertedVertices vrtCheck(*partitionerPtr_, smoothVertex);
    const labelHashSet& inverted = vrtCheck.invertedVertices();

    # ifdef smoothingTiming
    const scalar startMarking = omp_get_wtime();
    Info << "Finding inverted: checking time "
         << (startMarking - startChecking) << endl;
    # endif

    boolList originallySelected;
    if( nAdditionalLayers )
        originallySelected.setSize(smoothVertex.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100) reduction(+:nInvertedTria)
        # endif
        forAll(smoothVertex, bpI)
        {
            smoothVertex[bpI] = false;

            if( vertexType_[bpI] & LOCKED )
                continue;

            const label pointI = bPoints[bpI];

            //- check if the point is inverted
            if( !inverted.found(pointI) )
                continue;

            smoothVertex[bpI] = true;

            if( bpAtProcsPtr )
            {
                label minProc = Pstream::myProcNo();

                forAllRow(*bpAtProcsPtr, bpI, i)
                    minProc = min(minProc, bpAtProcsPtr->operator()(bpI, i));

                if( minProc == Pstream::myProcNo() )
                    ++nInvertedTria;
            }
            else
            {
                ++nInvertedTria;
            }
        }

        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            reduce(nInvertedTria, sumOp<label>());
            Info << "Number of inverted boundary faces is "
                 << nInvertedTria << endl;

            if( nInvertedTria != 0 )
            {
                //- add additional layers around inverted points
                for(label i=0;i<nAdditionalLayers;++i)
                {
                    //- copy originally selected points
                    # ifdef USE_OMP
                    label nTasks = 2 * omp_get_num_threads();

                    for(label taskI=0;taskI<nTasks;++taskI)
                    {
                        # pragma omp task default(shared) firstprivate(taskI)
                        {
                            for
                            (
                                label bpI=taskI;
                                bpI<smoothVertex.size();
                                bpI+=nTasks
                            )
                                originallySelected[bpI] = smoothVertex[bpI];
                        }
                    }

                    # pragma omp taskwait
                    # pragma omp flush(originallySelected)

                    nTasks = 5 * omp_get_num_threads();
                    # else
                    boolList originallySelected = smoothVertex;
                    const label nTasks = 1;
                    # endif

                    for(label taskI=0;taskI<nTasks;++taskI)
                    {
                        # ifdef USE_OMP
                        # pragma omp task default(shared) firstprivate(taskI)
                        # endif
                        {
                            for
                            (
                                label bpI=taskI;
                                bpI<smoothVertex.size();
                                bpI+=nTasks
                            )
                            {
                                if( originallySelected[bpI] )
                                {
                                    forAllRow(pFaces, bpI, pfI)
                                    {
                                        const face& bf =
                                            bFaces[pFaces(bpI, pfI)];

                                        forAll(bf, pI)
                                        {
                                            const label bpJ = bp[bf[pI]];

                                            if( vertexType_[bpJ] & LOCKED )
                                                continue;

                                            smoothVertex[bpJ] = true;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    # ifdef USE_OMP
                    # pragma omp taskwait
                    # pragma omp flush(smoothVertex)
                    # endif

                    if( Pstream::parRun() )
                    {
                        //- exchange global labels of inverted points
                        const labelLongList& globalPointLabel =
                            surfaceEngine_.globalBoundaryPointLabel();
                        const Map<label>& globalToLocal =
                            surfaceEngine_.globalToLocalBndPointAddressing();
                        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
                        const DynList<label>& neiProcs =
                            surfaceEngine_.bpNeiProcs();

                        std::map<label, labelLongList> shareData;
                        forAll(neiProcs, procI)
                            shareData.insert
                            (
                                std::make_pair(neiProcs[procI], labelLongList())
                            );

                        forAllConstIter(Map<label>, globalToLocal, iter)
                        {
                            const label bpI = iter();

                            if( !smoothVertex[bpI] )
                                continue;

                            forAllRow(bpAtProcs, bpI, procI)
                            {
                                const label neiProc = bpAtProcs(bpI, procI);

                                if( neiProc == Pstream::myProcNo() )
                                    continue;

                                shareData[neiProc].append
                                (
                                    globalPointLabel[bpI]
                                );
                            }
                        }

                        //- exchange data with other processors
                        labelLongList receivedData;
                        help::exchangeMap(shareData, receivedData);

                        forAll(receivedData, j)
                        {
                            const label bpI = globalToLocal[receivedData[j]];

                            smoothVertex[bpI] = true;
                        }
                    }
                }
            }
        }
    }

    # ifdef smoothingTiming
    const scalar endMarking = omp_get_wtime();
    Info << "Finding inverted: marking time "
         << (endMarking - startMarking) << endl;
    Info << "Finding inverted: Total time "
         << (endMarking - startChecking) << endl;
    # endif

    return nInvertedTria;
}


bool meshSurfaceOptimizer::untangleSurfaceNonConstrainedLaplace
(
    const labelLongList& selectedBoundaryPoints,
    const label nAdditionalLayers,
    const label nGlobalIterations
)
{
    Info << "Starting non-constrained untangling of the mesh surface" << endl;

    bool changed(false);

    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    surfaceEngine_.boundaryPointEdges();

    if( Pstream::parRun() )
    {
        surfaceEngine_.bpAtProcs();
        surfaceEngine_.globalToLocalBndPointAddressing();
        surfaceEngine_.globalBoundaryPointLabel();
        surfaceEngine_.bpNeiProcs();
    }

    //- select vertices that may be moved
    boolList smoothVertex(bPoints.size(), false);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(selectedBoundaryPoints, i)
    {
        const label bpI = selectedBoundaryPoints[i];

        if( vertexType_[bpI] & LOCKED )
            continue;

        smoothVertex[bpI] = true;
    }

    label nInvertedTria;

    labelLongList procBndPoints, movedPoints;
    labelLongList procEdgePoints, movedEdgePoints;

    label minNumInverted(bPoints.size());
    pointField minInvertedPoints(bPoints.size());

    label nIter(0);

    minNumInverted = bPoints.size();

    do
    {
        nInvertedTria = findInvertedVertices(smoothVertex, nAdditionalLayers);

        if( nInvertedTria == 0 )
        {
            break;
        }
        else if( enforceConstraints_ )
        {
            polyMeshGen& mesh =
                const_cast<polyMeshGen&>(surfaceEngine_.mesh());

            const label subsetId =
                mesh.addPointSubset(badPointsSubsetName_);

            forAll(smoothVertex, bpI)
                if( smoothVertex[bpI] )
                    mesh.addPointToSubset(subsetId, bPoints[bpI]);

            WarningIn
            (
                "bool meshSurfaceOptimizer::"
                "untangleSurfaceNonConstrainedLaplace(const labelLongList&,"
                " const label, const label)"
            ) << "Writing mesh with " << badPointsSubsetName_
              << " subset. These points cannot be untangled"
              << " without sacrificing geometry constraints. Exitting.."
              << endl;

            returnReduce(1, sumOp<label>());

            throw std::logic_error
            (
                "bool meshSurfaceOptimizer::untangleSurface"
                "(const labelLongList&, const label)"
                "Cannot untangle mesh!!"
            );
        }

        //- find the min number of inverted points and
        //- add the last number to the stack
        if( nInvertedTria < minNumInverted )
        {
            minNumInverted = nInvertedTria;

            # ifdef USE_OMP
            # pragma omp parallel for schedule(dynamic, 100)
            # endif
            forAll(bPoints, bpI)
                minInvertedPoints[bpI] = points[bPoints[bpI]];
        }

        //- find points which will be handled by the smoothers
        changed = true;

        procBndPoints.clear();
        movedPoints.clear();
        procEdgePoints.clear();
        movedEdgePoints.clear();

        forAll(bPoints, bpI)
        {
            if( !smoothVertex[bpI] )
                continue;

            if( vertexType_[bpI] & PARTITION )
            {
                movedPoints.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procBndPoints.append(bpI);
            }
            else if( vertexType_[bpI] & EDGE )
            {
                movedEdgePoints.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procEdgePoints.append(bpI);
            }
        }

        //- smooth edge vertices
        smoothEdgePoints(movedEdgePoints, procEdgePoints);

        //- use laplacian smoothing
        smoothLaplacianFC(movedPoints, procBndPoints, false, false);

    } while( nInvertedTria && (++nIter < nGlobalIterations) );

    if( nInvertedTria > 0 )
    {
        //- use the combination with the minimum number of inverted points
        forAll(minInvertedPoints, bpI)
            surfaceModifier_.moveBoundaryVertexNoUpdate
            (
                bpI,
                minInvertedPoints[bpI]
            );

        surfaceModifier_.updateGeometry();
    }

    //- check the final number of inverted points
    nInvertedTria = findInvertedVertices(smoothVertex, nAdditionalLayers);

    if( returnReduce(nInvertedTria, sumOp<label>()) != 0 )
    {
        //- the procedure has given up without success
        //- there exist some remaining inverted faces in the mesh
        polyMeshGen& mesh =
            const_cast<polyMeshGen&>(surfaceEngine_.mesh());

        label subsetId = mesh.pointSubsetIndex(badPointsSubsetName_);
        if( subsetId >= 0 )
            mesh.removePointSubset(subsetId);
        subsetId = mesh.addPointSubset(badPointsSubsetName_);

        forAll(smoothVertex, bpI)
            if( smoothVertex[bpI] )
                mesh.addPointToSubset(subsetId, bPoints[bpI]);
    }

    Info << "Finished non-constrained untangling of the mesh surface" << endl;

    return changed;
}

bool meshSurfaceOptimizer::untangleSurfaceNonConstrainedLaplace
(
    const label nAdditionalLayers,
    const label nGlobalIterations
)
{
    labelLongList activeBoundaryPoints(surfaceEngine_.boundaryPoints().size());
    forAll(activeBoundaryPoints, i)
        activeBoundaryPoints[i] = i;

    return untangleSurfaceNonConstrainedLaplace
    (
        activeBoundaryPoints,
        nAdditionalLayers,
        nGlobalIterations
    );
}

bool meshSurfaceOptimizer::untangleSurface
(
    const labelLongList& selectedBoundaryPoints,
    const label nAdditionalLayers,
    const label nGlobalIterations,
    const label nLocalIterations,
    const label nLaplaceIterations
)
{
    Info << "Starting untangling the surface of the volume mesh" << endl;

    # ifdef smoothingTiming
    const scalar startUntangling = omp_get_wtime();
    # endif

    bool changed(false);

    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    surfaceEngine_.boundaryPointEdges();

    if( Pstream::parRun() )
    {
        surfaceEngine_.bpAtProcs();
        surfaceEngine_.globalToLocalBndPointAddressing();
        surfaceEngine_.globalBoundaryPointLabel();
        surfaceEngine_.bpNeiProcs();
    }

    boolList smoothVertex(bPoints.size(), false);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(selectedBoundaryPoints, i)
    {
        const label bpI = selectedBoundaryPoints[i];

        if( vertexType_[bpI] & LOCKED )
            continue;

        smoothVertex[bpI] = true;
    }

    meshSurfaceMapper* mapperPtr = NULL;
    if( octreePtr_ )
        mapperPtr =
            new meshSurfaceMapper
            (
                *partitionerPtr_,
                *octreePtr_,
                &surfaceModifier_
            );

    bool remapVertex(true);
    label nInvertedTria;
    label nGlobalIter(0);

    labelLongList procBndPoints, movedPoints;
    labelLongList procEdgePoints, movedEdgePoints;

    label minNumInverted(bPoints.size());
    FIFOStack<label> nInvertedHistory;
    pointField minInvertedPoints(bPoints.size());

    # ifdef smoothingTiming
    const scalar startOptimization = omp_get_wtime();
    Info << "untangleSurface: preparation time "
         << (startOptimization - startUntangling) << endl;
    # endif

    do
    {
        label nIter(0), nAfterRefresh(0);

        minNumInverted = bPoints.size();

        do
        {
            nInvertedTria =
                findInvertedVertices(smoothVertex, nAdditionalLayers);

            if( nInvertedTria == 0 )
            {
                break;
            }
            else if( enforceConstraints_ && !remapVertex && mapperPtr )
            {
                polyMeshGen& mesh =
                    const_cast<polyMeshGen&>(surfaceEngine_.mesh());

                const label subsetId =
                    mesh.addPointSubset(badPointsSubsetName_);

                forAll(smoothVertex, bpI)
                    if( smoothVertex[bpI] )
                        mesh.addPointToSubset(subsetId, bPoints[bpI]);

                WarningIn
                (
                    "bool meshSurfaceOptimizer::untangleSurface"
                    "(const labelLongList&, const label)"
                ) << "Writing mesh with " << badPointsSubsetName_
                  << " subset. These points cannot be untangled"
                  << " without sacrificing geometry constraints. Exitting.."
                  << endl;

                returnReduce(1, sumOp<label>());

                throw std::logic_error
                (
                    "bool meshSurfaceOptimizer::untangleSurface"
                    "(const labelLongList&, const label)"
                    "Cannot untangle mesh!!"
                );
            }

            # ifdef smoothingTiming
            const scalar startSmoothing = omp_get_wtime();
            # endif

            //- find the min number of inverted points and
            //- add the last number to the stack
            if( nInvertedTria < minNumInverted )
            {
                minNumInverted = nInvertedTria;
                nAfterRefresh = 0;

                # ifdef USE_OMP
                # pragma omp parallel for schedule(dynamic, 100)
                # endif
                forAll(bPoints, bpI)
                    minInvertedPoints[bpI] = points[bPoints[bpI]];
            }

            //- count the number of iteration after the last minimum occurence
            ++nAfterRefresh;

            //- update the stack
            nInvertedHistory.push(nInvertedTria);
            if( nInvertedHistory.size() > 2 )
                nInvertedHistory.pop();

            //- check if the number of inverted points reduces
            bool minimumInStack(false);
            forAllConstIter(FIFOStack<label>, nInvertedHistory, it)
                if( it() == minNumInverted )
                    minimumInStack = true;

            //- stop if the procedure does not minimise
            //- the number of inverted points
            if( !minimumInStack || (nAfterRefresh > 2) )
                break;

            //- find points which will be handled by the smoothers
            changed = true;

            procBndPoints.clear();
            movedPoints.clear();
            procEdgePoints.clear();
            movedEdgePoints.clear();

            # ifdef USE_OMP
            # pragma omp parallel for schedule(dynamic, 100)
            # endif
            forAll(selectedBoundaryPoints, i)
            {
                const label bpI = selectedBoundaryPoints[i];

                if( !smoothVertex[bpI] )
                    continue;

                if( vertexType_[bpI] & LOCKED )
                    continue;

                if( vertexType_[bpI] & PARTITION )
                {
                    # ifdef USE_OMP
                    # pragma omp critical(movedPoints)
                    # endif
                    {
                        movedPoints.append(bpI);
                    }

                    if( vertexType_[bpI] & PROCBND )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(procBndPoints)
                        # endif
                        {
                            procBndPoints.append(bpI);
                        }
                    }
                }
                else if( vertexType_[bpI] & EDGE )
                {
                    # ifdef USE_OMP
                    # pragma omp critical(movedEdgePoints)
                    # endif
                    {
                        movedEdgePoints.append(bpI);
                    }

                    if( vertexType_[bpI] & PROCBND )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(procEdgePoints)
                        # endif
                        {
                            procEdgePoints.append(bpI);
                        }
                    }
                }
            }

            # ifdef smoothingTiming
            const scalar afterSelectingPoints = omp_get_wtime();
            Info << "Untangling: time for selecting points "
                 << (afterSelectingPoints-startSmoothing) << endl;
            # endif

            for(label iterI=0;iterI<nLocalIterations;++iterI)
            {
                //- smooth edge vertices
                smoothEdgePoints(movedEdgePoints, procEdgePoints);

                # ifdef DEBUGSmooth
                surfaceEngine_.writeSurfaceVTK
                (
                    "untangleEdges_"+help::labelToText(nIter)+"_"+
                    help::labelToText(nGlobalIter)+".vtk"
                );
                # endif

                //- use laplacian smoothing
                smoothLaplacianFC(movedPoints, procBndPoints);

                # ifdef DEBUGSmooth
                surfaceEngine_.writeSurfaceVTK
                (
                    "untangleLaplacianFC_"+help::labelToText(nIter)+"_"+
                    help::labelToText(nGlobalIter)+".vtk"
                );
                # endif

                //- use surface optimizer
                smoothSurfaceOptimizer(movedPoints, procBndPoints);

                # ifdef DEBUGSmooth
                surfaceEngine_.writeSurfaceVTK
                (
                    "untangleSurfOpt_"+help::labelToText(nIter)+"_"+
                    help::labelToText(nGlobalIter)+".vtk"
                );
                # endif

                //- improve flatness of boundary faces
                smoothFaceFlatness(movedPoints, procBndPoints);

                # ifdef DEBUGSmooth
                surfaceEngine_.writeSurfaceVTK
                (
                    "untangleFlatness_"+help::labelToText(nIter)+"_"+
                    help::labelToText(nGlobalIter)+".vtk"
                );
                # endif
            }

            # ifdef smoothingTiming
            const scalar afterMovingPoints = omp_get_wtime();
            Info << "Untangling: time for moving "
                 << (afterMovingPoints-afterSelectingPoints) << endl;
            # endif

            if( remapVertex && mapperPtr )
            {
                mapperPtr->mapEdgeNodes(movedEdgePoints);

                mapperPtr->mapVerticesOntoSurface(movedPoints);

                # ifdef DEBUGSmooth
                surfaceEngine_.writeSurfaceVTK
                (
                    "untangleMapper_"+help::labelToText(nIter)+"_"+
                    help::labelToText(nGlobalIter)+".vtk"
                );
                # endif
            }

            # ifdef smoothingTiming
            Info << "Untangling: time for mapping "
                 << (omp_get_wtime()-afterMovingPoints) << endl;
            # endif

        } while( nInvertedTria && (++nIter < nLocalIterations) );

        # ifdef smoothingTiming
        const scalar startReverting = omp_get_wtime();
        # endif

        if( nInvertedTria > 0 )
        {
            //- use the combination with the minimum number of inverted points
            forAll(minInvertedPoints, bpI)
            {
                if( vertexType_[bpI] & LOCKED )
                    continue;

                surfaceModifier_.moveBoundaryVertexNoUpdate
                (
                    bpI,
                    minInvertedPoints[bpI]
                );
            }

            surfaceModifier_.updateGeometry();

            # ifdef DEBUGSmooth
            surfaceEngine_.writeSurfaceVTK
            (
                "untangleMinInverted_"+
                help::labelToText(nGlobalIter)+".vtk"
            );
            # endif
        }

        # ifdef smoothingTiming
        Info << "Untangling: Time for reverting "
             << (omp_get_wtime()-startReverting) << endl;
        # endif

        if( nInvertedTria )
        {
            Info << "Smoothing remaining inverted vertices " << endl;

            # ifdef smoothingTiming
            const scalar startSmoothingRemaining = omp_get_wtime();
            # endif

            movedPoints.clear();
            procBndPoints.clear();
            forAll(selectedBoundaryPoints, i)
            {
                const label bpI = selectedBoundaryPoints[i];

                if( smoothVertex[bpI] )
                {
                    movedPoints.append(bpI);

                    if( vertexType_[bpI] & PROCBND )
                        procBndPoints.append(bpI);
                }
            }

            # ifdef smoothingTiming
            const scalar startMovingRemaining = omp_get_wtime();
            Info << "Untangling: selecting remaining "
                 << (startMovingRemaining-startSmoothingRemaining) << endl;
            # endif

            for(label iter=0;iter<nLaplaceIterations;++iter)
            {
                smoothLaplacianFC
                (
                    movedPoints,
                    procBndPoints,
                    0.8,
                    false,
                    false
                );

                # ifdef DEBUGSmooth
                surfaceEngine_.writeSurfaceVTK
                (
                    "untangleRemainLaplacianFC_"+help::labelToText(iter)+"_"+
                    help::labelToText(nGlobalIter)+".vtk"
                );
                # endif

                if( remapVertex && mapperPtr )
                {
                    //- project the vertices back onto the surface mesh
                    mapperPtr->mapVerticesOntoSurface(movedPoints);

                    # ifdef DEBUGSmooth
                    surfaceEngine_.writeSurfaceVTK
                    (
                        "untangleRemainMapper_"+
                        help::labelToText(iter)+"_"+
                        help::labelToText(nGlobalIter)+".vtk"
                    );
                    # endif
                }

                if( nGlobalIter > 5 )
                    remapVertex = false;
            }

            # ifdef smoothingTiming
            Info << "Untangling: Time for moving and mapping "
                << (omp_get_wtime()-startMovingRemaining) << endl;
            Info << "Untangling: Total time for remaining "
                 << (omp_get_wtime()-startSmoothingRemaining) << endl;
            # endif
        }

    } while( nInvertedTria && (++nGlobalIter < nGlobalIterations) );

    deleteDemandDrivenData(mapperPtr);

    //- check the final number of inverted points
    nInvertedTria = findInvertedVertices(smoothVertex, nAdditionalLayers);

    if( returnReduce(nInvertedTria, sumOp<label>()) != 0 )
    {
        //- the procedure has given up without success
        //- there exist some remaining inverted faces in the mesh
        polyMeshGen& mesh =
            const_cast<polyMeshGen&>(surfaceEngine_.mesh());

        label subsetId = mesh.pointSubsetIndex(badPointsSubsetName_);
        if( subsetId >= 0 )
            mesh.removePointSubset(subsetId);
        subsetId = mesh.addPointSubset(badPointsSubsetName_);

        forAll(smoothVertex, bpI)
            if( smoothVertex[bpI] )
                mesh.addPointToSubset(subsetId, bPoints[bpI]);
    }

    # ifdef smoothingTiming
    Info << "untangleSurface: Total time "
         << (omp_get_wtime() - startUntangling) << endl;
    # endif

    Info << "Finished untangling the surface of the volume mesh" << endl;

    return changed;
}

bool meshSurfaceOptimizer::untangleSurface
(
    const label nAdditionalLayers,
    const label nGlobalIterations,
    const label nLocalIterations,
    const label nLaplaceIterations
)
{
    labelLongList selectedPts(surfaceEngine_.boundaryPoints().size());
    forAll(selectedPts, i)
        selectedPts[i] = i;

    return untangleSurface
    (
        selectedPts,
        nAdditionalLayers,
        nGlobalIterations,
        nLocalIterations,
        nLaplaceIterations
    );
}

void meshSurfaceOptimizer::untangleAndRevertSurface
(
    const labelLongList& activeBoundaryPoints,
    const label nAdditionalLayers,
    const label nGlobalIterations,
    const label nLocalIterations,
    const label nLaplaceIterations
)
{
    const polyMeshGen& mesh = surfaceEngine_.mesh();
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    //- check if there exist any inverted points at the surface
    //- and move them back to the original coordinates
    boolList smoothVertex;
    point pOrig(vector::zero);

    scalar alpha = 0.5;
    do
    {
        untangleSurface
        (
            activeBoundaryPoints,
            nAdditionalLayers,
            nGlobalIterations,
            nLocalIterations,
            nLaplaceIterations
        );

        # ifdef smoothingTiming
        const scalar startReverting = omp_get_wtime();
        # endif

        if( findInvertedVertices(smoothVertex, 0) )
        {
            if( !mesh.hasPointsBackup() )
            {
                Info << "Points have not been backed up" << endl;

                return;
            }

            Info << "Reverting points to their original positions " << endl;

            //- unlock any locked vertices
            removeUserConstraints();

            //- starting moving the inverted points towards their original
            //- positions
            labelLongList lockPoints;
            forAll(smoothVertex, bpI)
            {
                if( !smoothVertex[bpI] )
                    continue;

                const label pointI = bPoints[bpI];

                mesh.getOrigPoint(pointI, pOrig);

                surfaceModifier_.moveBoundaryVertexNoUpdate
                (
                    bpI,
                    pOrig + alpha * (points[pointI] - pOrig)
                );

                lockPoints.append(bpI);
            }

            surfaceModifier_.updateGeometry(lockPoints);

            lockBoundaryPoints(lockPoints);

            alpha *= 0.5;
        }
        else
        {
            //- exit the loop because there exist no inverted vertices
            Info << "Mesh surface Ok." << endl;
            break;
        }

        # ifdef smoothingTiming
        Info << "Time for reverting of vertices "
             << (omp_get_wtime()-startReverting) << endl;
        # endif

    } while( true );
}

void meshSurfaceOptimizer::untangleAndRevertSurface
(
    const label nAdditionalLayers,
    const label nGlobalIterations,
    const label nLocalIterations,
    const label nLaplaceIterations
)
{
    # ifdef smoothingTiming
    const scalar startSelectingPoints = omp_get_wtime();
    # endif

    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    labelLongList activePoints(bPoints.size());
    forAll(bPoints, bpI)
        activePoints[bpI] = bpI;

    # ifdef smoothingTiming
    Info << "Time for selecting points "
         << (omp_get_wtime()-startSelectingPoints) << endl;
    # endif

    untangleAndRevertSurface
    (
        activePoints,
        nAdditionalLayers,
        nGlobalIterations,
        nLocalIterations,
        nLaplaceIterations
    );
}

void meshSurfaceOptimizer::optimizeSurfaceLaplace
(
    const label nIterations,
    const bool allowShrinking
)
{
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    surfaceEngine_.boundaryPointEdges();

    meshSurfaceMapper* mapperPtr = NULL;
    if( octreePtr_ )
        mapperPtr =
            new meshSurfaceMapper
            (
                *partitionerPtr_,
                *octreePtr_,
                &surfaceModifier_
            );

    labelLongList procBndPoints, edgePoints, partitionPoints, procPoints;
    forAll(bPoints, bpI)
    {
        if( vertexType_[bpI] & LOCKED )
            continue;

        if( vertexType_[bpI] & EDGE )
        {
            edgePoints.append(bpI);

            if( vertexType_[bpI] & PROCBND )
                procBndPoints.append(bpI);
        }
        else if( vertexType_[bpI] & PARTITION )
        {
            partitionPoints.append(bpI);

            if( vertexType_[bpI] & PROCBND )
                procPoints.append(bpI);
        }
    }

    //- optimize edge vertices
    Info << "Optimizing edges. Iteration:" << flush;
    for(label i=0;i<nIterations;++i)
    {
        Info << "." << flush;

        smoothEdgePoints(edgePoints, procBndPoints);

        //- project vertices back onto the boundary
        if( mapperPtr )
            mapperPtr->mapEdgeNodes(edgePoints);
    }
    Info << endl;

    //- delete the mapper
    deleteDemandDrivenData(mapperPtr);

    //- optimize positions of surface vertices which are not on surface edges
    Info << "Optimizing surface vertices. Iteration:";
    for(label i=0;i<nIterations;++i)
    {
        smoothLaplacianFC(partitionPoints, procPoints, !allowShrinking);

        Info << "." << flush;
    }

    Info << endl;

    untangleSurface(0);
}

void meshSurfaceOptimizer::optimizeSurface(const label nIterations)
{
    # ifdef smoothingTiming
    const scalar startOptimizingSurface = omp_get_wtime();
    # endif

    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();
    surfaceEngine_.boundaryPointEdges();

    //- construct the mapper
    meshSurfaceMapper* mapperPtr = NULL;
    if( octreePtr_ )
        mapperPtr =
            new meshSurfaceMapper
            (
                *partitionerPtr_,
                *octreePtr_,
                &surfaceModifier_
            );

    //- classify vertices
    # ifdef smoothingTiming
    const scalar startSelecting = omp_get_wtime();
    Info << "optimizeSurface: Preparation time "
         << (startSelecting - startOptimizingSurface) << endl;
    # endif

    labelLongList procBndPoints, edgePoints, partitionPoints, procPoints;
    forAll(bPoints, bpI)
    {
        if( vertexType_[bpI] & LOCKED )
            continue;

        if( vertexType_[bpI] & EDGE )
        {
            edgePoints.append(bpI);

            if( vertexType_[bpI] & PROCBND )
                procBndPoints.append(bpI);
        }
        else if( vertexType_[bpI] & PARTITION )
        {
            partitionPoints.append(bpI);

            if( vertexType_[bpI] & PROCBND )
                procPoints.append(bpI);
        }
    }

    # ifdef smoothingTiming
    const scalar startSmoothingEdgePoints = omp_get_wtime();
    Info << "optimizeSurface: point selection time "
         << (startSmoothingEdgePoints - startSelecting) << endl;
    # endif

    //- optimize edge vertices
    Info << "Optimizing edges. Iteration:" << flush;
    for(label i=0;i<nIterations;++i)
    {
        Info << "." << flush;

        smoothEdgePoints(edgePoints, procBndPoints);

        # ifdef DEBUGSmooth
        surfaceEngine_.writeSurfaceVTK
        (
            "optimizeEdges_"+help::labelToText(i)+".vtk"
        );
        # endif

        //- project vertices back onto the boundary
        if( mapperPtr )
        {
            mapperPtr->mapEdgeNodes(edgePoints);

            # ifdef DEBUGSmooth
            surfaceEngine_.writeSurfaceVTK
            (
                "optimizeEdgesMapper_"+help::labelToText(i)+".vtk"
            );
            # endif
        }
    }
    Info << endl;

    # ifdef smoothingTiming
    const scalar startSmoothingPoints = omp_get_wtime();
    Info << "optimizeSurface: Edge smoothing "
         << (startSmoothingPoints - startSmoothingEdgePoints) << endl;
    # endif

    //- optimize positions of surface vertices which are not on surface edges
    Info << "Optimizing surface vertices. Iteration:";
    for(label i=0;i<nIterations;++i)
    {
        smoothFaceFlatness(partitionPoints, procPoints, 0.8, true);

        # ifdef DEBUGSmooth
        surfaceEngine_.writeSurfaceVTK
        (
            "optimizeFlatness_"+help::labelToText(i)+".vtk"
        );
        # endif

        smoothLaplacianFC(partitionPoints, procPoints, 0.8, true);

        # ifdef DEBUGSmooth
        surfaceEngine_.writeSurfaceVTK
        (
            "optimizeLaplace_"+help::labelToText(i)+".vtk"
        );
        # endif

        smoothSurfaceOptimizer(partitionPoints, procPoints);

        # ifdef DEBUGSmooth
        surfaceEngine_.writeSurfaceVTK
        (
            "optimizeSurfOpt_"+help::labelToText(i)+".vtk"
        );
        # endif

        Info << "." << flush;
    }

    Info << endl;

    //- delete the mapper
    deleteDemandDrivenData(mapperPtr);

    # ifdef smoothingTiming
    const scalar endSmoothing = omp_get_wtime();
    Info << "optimizeSurface: Smoothing time "
         << (endSmoothing - startOptimizingSurface) << endl;
    # endif

    //untangleAndRevertSurface();
}

void meshSurfaceOptimizer::optimizeSurface2D(const label nIterations)
{
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const edgeLongList& edges = surfaceEngine_.edges();
    const labelLongList& bp = surfaceEngine_.bp();

    polyMeshGen2DEngine mesh2DEngine
    (
        const_cast<polyMeshGen&>(surfaceEngine_.mesh())
    );
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();

    labelLongList procBndPoints, movedPoints, activeEdges, updatePoints;
    forAll(edges, beI)
    {
        const edge& e = edges[beI];

        if( zMinPoint[e.start()] ^ zMinPoint[e.end()] )
        {
            label bpI = bp[e.start()];
            if( !zMinPoint[e.start()] )
                bpI = bp[e.end()];

            if( vertexType_[bpI] & EDGE )
            {
                activeEdges.append(beI);

                updatePoints.append(bp[e.start()]);
                updatePoints.append(bp[e.end()]);

                movedPoints.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procBndPoints.append(bpI);
            }
        }
    }

    meshSurfaceMapper2D* mapperPtr = NULL;
    if( octreePtr_ )
        mapperPtr = new meshSurfaceMapper2D(surfaceEngine_, *octreePtr_);

    //- optimize edge vertices

    Info << "Optimizing edges. Iteration:" << flush;
    for(label i=0;i<nIterations;++i)
    {
        Info << "." << flush;

        smoothEdgePoints(movedPoints, procBndPoints);

        //- move points with maximum z coordinate
        mesh2DEngine.correctPoints();

        //- map boundary edges to the surface
        mapperPtr->mapVerticesOntoSurfacePatches(activeEdges);

        //- update normal, centres, etc, after the surface has been modified
        surfaceModifier_.updateGeometry(updatePoints);
    }
    Info << endl;

    //- optimize Pts of surface vertices which are not on surface edges
    procBndPoints.clear();
    movedPoints.clear();
    forAll(bPoints, bpI)
        if( zMinPoint[bPoints[bpI]] && (vertexType_[bpI] & PARTITION) )
        {
            movedPoints.append(bpI);

            if( vertexType_[bpI] & PROCBND )
                procBndPoints.append(bpI);
        }
    Info << "Optimizing surface vertices. Iteration:";
    for(label i=0;i<nIterations;++i)
    {
        Info << "." << flush;

        smoothSurfaceOptimizer(movedPoints, procBndPoints);

        //- move the points which are not at minimum z coordinate
        mesh2DEngine.correctPoints();

        //- update geometrical data due to movement of vertices
        surfaceModifier_.updateGeometry();
    }

    Info << endl;

    deleteDemandDrivenData(mapperPtr);
}

void meshSurfaceOptimizer::optimizeLowQualitySurface2D(const label nIterations)
{
    const polyMeshGen& mesh = surfaceEngine_.mesh();
    const faceListPMG& faces = mesh.faces();
    const VRWGraph& pointFaces = mesh.addressingData().pointFaces();

    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const labelLongList& bp = surfaceEngine_.bp();

    polyMeshGen2DEngine mesh2DEngine(const_cast<polyMeshGen&>(mesh));
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();
    const boolList& activeFace = mesh2DEngine.activeFace();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();

    boolList activeBoundaryPoint(bPoints.size());
    boolList changedFace(activeFace.size(), true);

    label iterationI(0);
    do
    {
        labelHashSet badFaces;
        const label nBadFaces =
            polyMeshGenChecks::findLowQualityFaces
            (
                mesh,
                badFaces,
                false,
                &changedFace
            );

        Info << "Iteration " << iterationI
             << ". Number of bad faces " << nBadFaces << endl;

        if( nBadFaces == 0 )
            break;

        //- update active points and faces affected by the movement
        //- of active points
        activeBoundaryPoint = false;
        changedFace = false;
        forAllConstIter(labelHashSet, badFaces, it)
        {
            const face& f = faces[it.key()];

            forAll(f, pI)
            {
                if( zMinPoint[f[pI]] )
                {
                    activeBoundaryPoint[bp[f[pI]]] = true;

                    forAllRow(pointFaces, f[pI], pfI)
                        changedFace[pointFaces(f[pI], pfI)] = true;
                }
            }
        }

        if( Pstream::parRun() )
        {
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();
            const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
            const VRWGraph& bpNeiProcs = surfaceEngine_.bpAtProcs();

            std::map<label, labelLongList> exchangeData;
            forAll(neiProcs, i)
                exchangeData[neiProcs[i]].clear();

            //- collect active points at inter-processor boundaries
            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                if( activeBoundaryPoint[bpI] )
                {
                    forAllRow(bpNeiProcs, bpI, i)
                    {
                        const label neiProc = bpNeiProcs(bpI, i);

                        if( neiProc == Pstream::myProcNo() )
                            continue;

                        exchangeData[neiProc].append(it.key());
                    }
                }
            }

            //- exchange active points among the processors
            labelLongList receivedData;
            help::exchangeMap(exchangeData, receivedData);

            //- ensure that all processors have the same Pts active
            forAll(receivedData, i)
            {
                const label bpI = globalToLocal[receivedData[i]];

                //- activate this boundary point
                activeBoundaryPoint[bpI] = true;

                //- set the changeFaces for the faces attached to this point
                forAllRow(pointFaces, bPoints[bpI], pfI)
                    changedFace[pointFaces(bPoints[bpI], pfI)] = true;
            }
        }

        //- apply smoothing to the activated points
        labelLongList movedPts, procBndPts, edgePts, procEdgePts;
        forAll(bPoints, bpI)
        {
            if( !activeBoundaryPoint[bpI] )
                continue;

            if( vertexType_[bpI] & EDGE )
            {
                edgePts.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procEdgePts.append(bpI);
            }
            else if( vertexType_[bpI] & PARTITION )
            {
                movedPts.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procBndPts.append(bpI);
            }
        }

        for(label i=0;i<5;++i)
        {
            smoothEdgePoints(edgePts, procEdgePts);

            smoothSurfaceOptimizer(movedPts, procBndPts);
        }

        //- move the points which are not at minimum z coordinate
        mesh2DEngine.correctPoints();

        //- update geometrical data due to movement of vertices
        surfaceModifier_.updateGeometry();

        //- update cell centres and face centres
        const_cast<polyMeshGenAddressing&>
        (
            mesh.addressingData()
        ).updateGeometry(changedFace);

    } while( ++iterationI < nIterations );

    //- delete invalid data
    mesh.clearAddressingData();
}

void meshSurfaceOptimizer::untangleSurface2D(const label maxNumIterations)
{
    const polyMeshGen& mesh = surfaceEngine_.mesh();
    const faceListPMG& faces = mesh.faces();
    const VRWGraph& pointFaces = mesh.addressingData().pointFaces();

    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const labelLongList& bp = surfaceEngine_.bp();

    polyMeshGen2DEngine mesh2DEngine(const_cast<polyMeshGen&>(mesh));
    const boolList& zMinPoint = mesh2DEngine.zMinPoints();
    const boolList& activeFace = mesh2DEngine.activeFace();

    //- needed for parallel execution
    surfaceEngine_.pointFaces();
    surfaceEngine_.faceCentres();
    surfaceEngine_.pointPoints();
    surfaceEngine_.boundaryPointEdges();
    surfaceEngine_.boundaryFacePatches();
    surfaceEngine_.pointNormals();

    boolList activeBoundaryPoint(bPoints.size());
    boolList changedFace(activeFace.size(), true);

    label iterationI(0);
    do
    {
        labelHashSet badFaces;
        const label nBadFaces =
            polyMeshGenChecks::findBadFaces
            (
                mesh,
                badFaces,
                false,
                &changedFace
            );

        Info << "Iteration " << iterationI
             << ". Number of bad faces " << nBadFaces << endl;

        if( nBadFaces == 0 )
            break;

        //- update active points and faces affected by the movement
        //- of active points
        activeBoundaryPoint = false;
        changedFace = false;
        forAllConstIter(labelHashSet, badFaces, it)
        {
            const face& f = faces[it.key()];

            forAll(f, pI)
            {
                if( zMinPoint[f[pI]] )
                {
                    activeBoundaryPoint[bp[f[pI]]] = true;

                    forAllRow(pointFaces, f[pI], pfI)
                        changedFace[pointFaces(f[pI], pfI)] = true;
                }
            }
        }

        if( Pstream::parRun() )
        {
            const Map<label>& globalToLocal =
                surfaceEngine_.globalToLocalBndPointAddressing();
            const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
            const VRWGraph& bpNeiProcs = surfaceEngine_.bpAtProcs();

            std::map<label, labelLongList> exchangeData;
            forAll(neiProcs, i)
                exchangeData[neiProcs[i]].clear();

            //- collect active points at inter-processor boundaries
            forAllConstIter(Map<label>, globalToLocal, it)
            {
                const label bpI = it();

                if( activeBoundaryPoint[bpI] )
                {
                    forAllRow(bpNeiProcs, bpI, i)
                    {
                        const label neiProc = bpNeiProcs(bpI, i);

                        if( neiProc == Pstream::myProcNo() )
                            continue;

                        exchangeData[neiProc].append(it.key());
                    }
                }
            }

            //- exchange active points among the processors
            labelLongList receivedData;
            help::exchangeMap(exchangeData, receivedData);

            //- ensure that all processors have the same Pts active
            forAll(receivedData, i)
            {
                const label bpI = globalToLocal[receivedData[i]];

                //- activate this boundary point
                activeBoundaryPoint[bpI] = true;

                //- set the changeFaces for the faces attached to this point
                forAllRow(pointFaces, bPoints[bpI], pfI)
                    changedFace[pointFaces(bPoints[bpI], pfI)] = true;
            }
        }

        //- apply smoothing to the activated points
        labelLongList movedPts, procBndPts, edgePts, procEdgePts;
        forAll(bPoints, bpI)
        {
            if( !activeBoundaryPoint[bpI] )
                continue;

            if( vertexType_[bpI] & EDGE )
            {
                edgePts.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procEdgePts.append(bpI);
            }
            else if( vertexType_[bpI] & PARTITION )
            {
                movedPts.append(bpI);

                if( vertexType_[bpI] & PROCBND )
                    procBndPts.append(bpI);
            }
        }

        for(label i=0;i<5;++i)
        {
            smoothEdgePoints(edgePts, procEdgePts);

            smoothSurfaceOptimizer(movedPts, procBndPts);
        }

        //- move the points which are not at minimum z coordinate
        mesh2DEngine.correctPoints();

        //- update geometrical data due to movement of vertices
        surfaceModifier_.updateGeometry();

        //- update cell centres and face centres
        const_cast<polyMeshGenAddressing&>
        (
            mesh.addressingData()
        ).updateGeometry(changedFace);

    } while( ++iterationI < maxNumIterations );

    //- delete invalid data
    mesh.clearAddressingData();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
