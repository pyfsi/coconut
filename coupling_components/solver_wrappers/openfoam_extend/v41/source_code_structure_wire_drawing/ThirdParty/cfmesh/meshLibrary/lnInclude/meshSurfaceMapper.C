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

#include "meshSurfaceMapper.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfacePartitioner.H"
#include "triSurf.H"
#include "triSurfacePartitioner.H"
#include "demandDrivenData.H"
#include "meshOctree.H"
#include "helperFunctions.H"
#include "findCellsIntersectingSurface.H"
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceCheckInvertedVertices.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper::createMeshSurfacePartitioner() const
{
    surfaceEnginePartitionerPtr_ = new meshSurfacePartitioner(surfaceEngine_);
}

void meshSurfaceMapper::createTriSurfacePartitioner() const
{
    surfPartitionerPtr_ = new triSurfacePartitioner(meshOctree_.surface());
}

meshSurfaceEngineModifier& meshSurfaceMapper::meshSurfaceModifier()
{
    if( !surfaceModifierPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
        {
            FatalErrorIn
            (
                "meshSurfaceEngineModifier&"
                " meshSurfaceMapper::meshSurfaceModifier()"
            ) << "Calculating addressing inside a parallel region"
              << exit(FatalError);
        }
        # endif

        surfaceModifierPtr_ = new meshSurfaceEngineModifier(surfaceEngine_);
    }

    return *surfaceModifierPtr_;
}

void meshSurfaceMapper::clearOut()
{
    if( deletePartitioner_ )
        deleteDemandDrivenData(surfaceEnginePartitionerPtr_);
    deleteDemandDrivenData(surfPartitionerPtr_);
    if( deleteModifier_ )
        deleteDemandDrivenData(surfaceModifierPtr_);
}

void meshSurfaceMapper::calculateMovementShrinkingSurfaceLaplaceFC
(
    const List<DynList<scalar> >& faceCentreWeights,
    vectorField& displacements
) const
{
    const labelLongList& boundaryPoints = surfaceEngine_.boundaryPoints();
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();
    const vectorLongList& faceCentres = surfaceEngine_.faceCentres();

    const pointFieldPMG& points = surfaceEngine_.points();

    List<labelledPointScalar> preMapPositions(boundaryPoints.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(boundaryPoints, bpI)
    {
        labelledPointScalar lp(bpI, vector::zero, 0.0);

        if( !surfaceEngine_.mesh().isLockedPoint(boundaryPoints[bpI]) )
        {
            forAllRow(pointFaces, bpI, pfI)
            {
                const point& fc = faceCentres[pointFaces(bpI, pfI)];

                const scalar w = faceCentreWeights[bpI][pfI];
                lp.coordinates() += w * fc;
                lp.scalarValue() += w;
            }
        }

        preMapPositions[bpI] = lp;
    }

    if( Pstream::parRun() )
    {
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
        const labelLongList& globalPointLabel =
            surfaceEngine_.globalBoundaryPointLabel();
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();

        //- collect data to be sent to other processors
        std::map<label, LongList<labelledPointScalar> > exchangeData;
        forAll(surfaceEngine_.bpNeiProcs(), i)
            exchangeData.insert
            (
                std::make_pair
                (
                    surfaceEngine_.bpNeiProcs()[i],
                    LongList<labelledPointScalar>()
                )
            );

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            forAllRow(bpAtProcs, bpI, procI)
            {
                const label neiProc = bpAtProcs(bpI, procI);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append
                (
                    labelledPointScalar
                    (
                        globalPointLabel[bpI],
                        preMapPositions[bpI].coordinates(),
                        preMapPositions[bpI].scalarValue()
                    )
                );
            }
        }

        //- exchange data with other processors
        LongList<labelledPointScalar> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        //- combine collected data with the available data
        forAll(receivedData, i)
        {
            const labelledPointScalar& lps = receivedData[i];

            const label bpI = globalToLocal[lps.pointLabel()];

            labelledPointScalar& lp = preMapPositions[bpI];
            lp.coordinates() += lps.coordinates();
            lp.scalarValue() += lps.scalarValue();
        }
    }

    //- calculate displacements of surface points
    displacements.setSize(boundaryPoints.size());
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(boundaryPoints, bpI)
    {
        if( !surfaceEngine_.mesh().isLockedPoint(boundaryPoints[bpI]) )
        {
            labelledPointScalar& lps = preMapPositions[bpI];

            lps.coordinates() /= lps.scalarValue();

            displacements[bpI] =
                lps.coordinates() - points[boundaryPoints[bpI]];
        }
        else
        {
            displacements[bpI] = vector::zero;
        }
    }
}

void meshSurfaceMapper::findNearestPointsToFaceCentres
(
    const labelLongList& activeFaces,
    List<labelledPointScalar>& nearestPoints
) const
{
    const vectorLongList& faceCentres = surfaceEngine_.faceCentres();

    nearestPoints.setSize(faceCentres.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100)
    # endif
    forAll(activeFaces, i)
    {
        const label bfI = activeFaces[i];

        labelledPointScalar& lps = nearestPoints[bfI];

        label patch;
        meshOctree_.findNearestSurfacePoint
        (
            lps.coordinates(),
            lps.scalarValue(),
            lps.pointLabel(),
            patch,
            faceCentres[bfI]
        );
    }
}

void meshSurfaceMapper::findNearestPointsToFaceCentres
(
    List<labelledPointScalar>& nearestPoints
)
{
    labelLongList activeFaces(surfaceEngine_.boundaryFaces().size());
    forAll(activeFaces, i)
        activeFaces[i] = i;

    findNearestPointsToFaceCentres(activeFaces, nearestPoints);
}

void meshSurfaceMapper::findNearestPoints
(
    const labelLongList& activePoints,
    List<labelledPointScalar>& nearestPoints
) const
{
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    nearestPoints.setSize(bPoints.size());

    const pointFieldPMG& points = surfaceEngine_.points();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100)
    # endif
    forAll(activePoints, i)
    {
        const label bpI = activePoints[i];

        const point& p = points[bPoints[bpI]];

        labelledPointScalar& lps = nearestPoints[bpI];

        label patch;
        meshOctree_.findNearestSurfacePoint
        (
            lps.coordinates(),
            lps.scalarValue(),
            lps.pointLabel(),
            patch,
            p
        );
    }
}

void meshSurfaceMapper::findNearestPoints
(
    List<labelledPointScalar>& nearestPoints
) const
{
    labelLongList activePoints(surfaceEngine_.boundaryPoints().size());
    forAll(activePoints, bpI)
        activePoints[bpI] = bpI;

    findNearestPoints(activePoints, nearestPoints);
}

void meshSurfaceMapper::findFacePatches
(
    const labelLongList& activeFaces,
    labelLongList& facePatch
) const
{
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

    const vectorLongList& faceCentres = surfaceEngine_.faceCentres();
    const pointFieldPMG& points = surfaceEngine_.points();

    const triSurf& surf = meshOctree_.surface();

    facePatch.setSize(bFaces.size());

    if( surf.patches().size() != surfaceEngine_.mesh().boundaries().size() )
    {
        facePatch = -1;

        //- reconstruct the patches in the volume mesh
        DynList<label> containedLeaves, ct;
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50) \
        private(containedLeaves, ct)
        # endif
        forAll(activeFaces, fI)
        {
            const label bfI = activeFaces[fI];

            const face& bf = bFaces[bfI];
            scalar range(VGREAT);
            forAll(bf, pI)
            {
                range =
                    Foam::min
                    (
                        range,
                        mag(faceCentres[bfI] - points[bf[pI]])
                    );
            }

            const boundBox bb
            (
                faceCentres[bfI] - vector(range, range, range),
                faceCentres[bfI] + vector(range, range, range)
            );

            containedLeaves.clear();
            meshOctree_.findLeavesContainedInBox(bb, containedLeaves);

            DynList<label> patches;
            forAll(containedLeaves, clI)
            {
                ct.clear();
                meshOctree_.containedTriangles(containedLeaves[clI], ct);

                forAll(ct, i)
                    patches.appendIfNotIn(surf[ct[i]].region());
            }

            scalar metric(VGREAT);
            label bestPatch(-1);
            forAll(patches, ptchI)
            {
                const scalar m = faceMetricInPatch(bfI, patches[ptchI]);

                if( m < metric )
                {
                    metric = m;
                    bestPatch = patches[ptchI];
                }
            }

            if( bestPatch < 0 )
            {
                facePatch[bfI] = -1;
                continue;
            }

            facePatch[bfI] = bestPatch;
        }
    }
    else
    {
        //- patches are already available
        facePatch = surfaceEngine_.boundaryFacePatches();
    }
}

void meshSurfaceMapper::findFacePatches(labelLongList& facePatch) const
{
    labelLongList activeFaces(surfaceEngine_.boundaryFaces().size());
    forAll(activeFaces, i)
        activeFaces[i] = i;

    findFacePatches(activeFaces, facePatch);
}

void meshSurfaceMapper::findPointPatches
(
    const labelLongList& activePoints,
    List<DynList<label> >& pointPatches
) const
{
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();

    //- find faces attached to active points
    boolList activeFace(surfaceEngine_.boundaryFaces().size(), false);
    forAll(activePoints, i)
    {
        const label bpI = activePoints[i];

        forAllRow(pointFaces, bpI, pfI)
            activeFace[pointFaces(bpI, pfI)] = true;
    }

    labelLongList activeFaces;
    forAll(activeFace, bfI)
        if( activeFace[bfI] )
            activeFaces.append(bfI);

    labelLongList facePatch;
    findFacePatches(activeFaces, facePatch);

    pointPatches.setSize(pointFaces.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(pointFaces, bpI)
    {
        DynList<label>& patches = pointPatches[bpI];
        patches.clear();

        forAllRow(pointFaces, bpI, pfI)
        {
            const label patchI = facePatch[pointFaces(bpI, pfI)];

            if( patchI >= 0 )
                patches.appendIfNotIn(patchI);
        }
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(surfaceEngine_.bpNeiProcs(), i)
            exchangeData[surfaceEngine_.bpNeiProcs()[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( pointPatches[bpI].size() == 0 )
                continue;

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                labelLongList& dts = exchangeData[neiProc];

                dts.append(it.key());
                dts.append(pointPatches[bpI].size());

                forAll(pointPatches[bpI], i)
                    dts.append(pointPatches[bpI][i]);
            }
        }

        labelLongList receiveData;
        help::exchangeMap(exchangeData, receiveData);

        label counter(0);
        while( counter < receiveData.size() )
        {
            const label bpI = globalToLocal[receiveData[counter++]];
            const label size = receiveData[counter++];

            for(label i=0;i<size;++i)
                pointPatches[bpI].appendIfNotIn(receiveData[counter++]);
        }
    }
}

void meshSurfaceMapper::findPointPatches
(
    List<DynList<label> >& pointPatches
) const
{
    labelLongList activePoints(surfaceEngine_.boundaryPoints().size());
    forAll(activePoints, bpI)
        activePoints[bpI] = bpI;

    findPointPatches(activePoints, pointPatches);
}

void meshSurfaceMapper::findMappingVertices
(
    const labelLongList& activePoints,
    pointField& newPoints,
    boolList& foundNearest
) const
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& boundaryPoints = surfaceEngine_.boundaryPoints();

    newPoints.setSize(boundaryPoints.size());
    foundNearest.setSize(boundaryPoints.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(activePoints, i)
    {
        const label bpI = activePoints[i];

        const point& p = points[boundaryPoints[bpI]];

        if( surfaceEngine_.mesh().isLockedPoint(boundaryPoints[bpI]) )
        {
            newPoints[bpI] = p;
            foundNearest[bpI] = true;
            continue;
        }

        //- find the nearest point at the surface contained inside the boxes
        //- located within the given range
        point nearest;
        scalar dSq;
        label nt, patch;
        meshOctree_.findNearestSurfacePoint(nearest, dSq, nt, patch, p);

        if( nt >= 0 )
        {
            foundNearest[bpI] = true;
            newPoints[bpI] = nearest;
        }
        else
        {
            foundNearest[bpI] = false;
            newPoints[bpI] = p;
        }
    }
}

void meshSurfaceMapper::findMappingVertices
(
    pointField& newPoints,
    boolList& foundNearest
) const
{
    labelLongList activePoints(surfaceEngine_.boundaryPoints().size());
    forAll(activePoints, bpI)
        activePoints[bpI] = bpI;

    findMappingVertices(activePoints, newPoints, foundNearest);
}

void meshSurfaceMapper::findMappingVertices
(
    const labelLongList& activePoints,
    const List<DynList<label> >& pointPatches,
    pointField& newPoints,
    boolList& foundNearest
) const
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& boundaryPoints = surfaceEngine_.boundaryPoints();

    newPoints.setSize(boundaryPoints.size());
    foundNearest.setSize(boundaryPoints.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(activePoints, i)
    {
        const label bpI = activePoints[i];

        const point& p = points[boundaryPoints[bpI]];

        //- do not move locked points
        if( surfaceEngine_.mesh().isLockedPoint(boundaryPoints[bpI]) )
        {
            newPoints[bpI] = p;
            foundNearest[bpI] = true;
            continue;
        }

        //- find the nearest point at the surface contained inside the boxes
        //- located within the given range
        if( pointPatches[bpI].size() > 1 )
        {
            point nearest;
            scalar dSq;
            const bool found =
                meshOctree_.findNearestPointToPatches
                (
                    nearest,
                    dSq,
                    p,
                    pointPatches[bpI],
                    1e-4
                );

            if( !found )
                Info << "Failed vertex " << bpI << " at patches "
                     << pointPatches[bpI] << endl;

            foundNearest[bpI] = found;
            newPoints[bpI] = nearest;
        }
        else if( pointPatches[bpI].size() )
        {
            point nearest;
            scalar dSq;
            label nt;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    nearest,
                    dSq,
                    nt,
                    pointPatches[bpI][0],
                    p
                );

            if( nt < 0 )
                Info << "Failed vertex " << bpI << " at patch "
                     << pointPatches[bpI][0] << endl;

            foundNearest[bpI] = (nt >= 0);
            newPoints[bpI] = nearest;
        }
        else
        {
            point nearest;
            scalar dSq;
            label nt, region;
                meshOctree_.findNearestSurfacePoint
                (
                    nearest,
                    dSq,
                    nt,
                    region,
                    p
                );

            if( nt < 0 )
                Info << "Failed vertex " << bpI << endl;

            foundNearest[bpI] = (nt >= 0);
            newPoints[bpI] = nearest;
        }
    }
}

void meshSurfaceMapper::findMappingVertices
(
    const List<DynList<label> >& pointPatches,
    pointField& newPoints,
    boolList& foundNearest
) const
{
    labelLongList activePoints(surfaceEngine_.boundaryPoints().size());
    forAll(activePoints, bpI)
        activePoints[bpI] = bpI;

    findMappingVertices(activePoints, pointPatches, newPoints, foundNearest);
}

void meshSurfaceMapper::analysePointsAndFaces()
{
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const labelLongList& faceOwner = surfaceEngine_.faceOwners();
    const edgeLongList& edges = surfaceEngine_.edges();
    const VRWGraph& bpEdges = surfaceEngine_.boundaryPointEdges();
    const labelLongList& bp = surfaceEngine_.bp();

    const label nIntFaces = surfaceEngine_.mesh().nInternalFaces();
    const cellListPMG& cells = surfaceEngine_.mesh().cells();
    const faceListPMG& faces = surfaceEngine_.mesh().faces();

    findCellsIntersectingSurface cis(surfaceEngine_.mesh(), meshOctree_);

    pointType_.setSize(bPoints.size());
    pointType_ = LOCKED;

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        LongList<edge> edgeCandidates;

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100) nowait
        # endif
        forAll(faceOwner, bfI)
        {
            const cell& c = cells[faceOwner[bfI]];

            forAll(c, fI)
            {
                const label faceI = c[fI];

                if
                (
                    (faceI < nIntFaces) ||
                    (faceI >= (nIntFaces+faceOwner.size()))
                )
                {
                    const face& f = faces[faceI];

                    //- face in an internal face of the mesh
                    forAll(f, eI)
                    {
                        const edge e = f.faceEdge(eI);

                        if( (bp[e.start()] != -1) || (bp[e.end()] != -1) )
                            edgeCandidates.append(e);
                    }
                }
            }
        }

        //- remove edges at the surface of the mesh from the list of candidates
        forAllReverse(edgeCandidates, eI)
        {
            const edge& e = edgeCandidates[eI];

            const label bps = bp[e.start()];

            if( bps < 0 )
                continue;

            forAllRow(bpEdges, bps, peI)
                if( edges[bpEdges(bps, peI)] == e )
                {
                    edgeCandidates.remove(eI);
                    continue;
                }
        }

        //- the remaining edges are internal edges
        forAll(edgeCandidates, eI)
        {
            const edge& e = edgeCandidates[eI];

            const label bps = bp[e.start()];
            const label bpe = bp[e.end()];

            if( (bps != -1) && (bpe != -1) )
            {
                //- this is an internal edge with both ends at the surface
                pointType_[bps] |= SPECIAL;
                pointType_[bpe] |= SPECIAL;
            }

            if( bps != -1 )
            {
                forAllRow(pFaces, bps, pfI)
                {
                    const face& bf = bFaces[pFaces(bps, pfI)];

                    forAll(bf, pI)
                        if( pointType_[bp[bf[pI]]] & LOCKED )
                            pointType_[bp[bf[pI]]] ^= LOCKED;
                }
            }

            if( bpe != -1 )
            {
                forAllRow(pFaces, bpe, pfI)
                {
                    const face& bf = bFaces[pFaces(bpe, pfI)];

                    forAll(bf, pI)
                        if( pointType_[bp[bf[pI]]] & LOCKED )
                            pointType_[bp[bf[pI]]] ^= LOCKED;
                }
            }
        }
    }

    //- classify faces
    faceType_.setSize(bFaces.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100)
    # endif
    forAll(faceType_, bfI)
    {
        const face& bf = bFaces[bfI];
        faceType_[bfI] = NONE;

        forAll(bf, pI)
        {
            if( pointType_[bp[bf[pI]]] & SPECIAL )
            {
                //- check if the owner cell is intersected by the surface
                if( cis.cellIntersectsSurface(faceOwner[bfI]) )
                {
                    faceType_[bfI] |= INNERATTRACTION;
                }
                else
                {
                    faceType_[bfI] |= OUTERATTRACTION;
                }

                break;
            }
        }
    }

    # ifdef DEBUGSearch
    polyMeshGen& mesh = const_cast<polyMeshGen&>(surfaceEngine_.mesh());
    const label iId = mesh.addFaceSubset("innerAttraction");
    const label oId = mesh.addFaceSubset("outerAttraction");
    forAll(faceType_, bfI)
    {
        if( faceType_[bfI] & INNERATTRACTION )
            mesh.addFaceToSubset(iId, nIntFaces+bfI);
        if( faceType_[bfI] & OUTERATTRACTION )
            mesh.addFaceToSubset(oId, nIntFaces+bfI);
    }

    const label spId = mesh.addPointSubset("specialPoints");
    const label lpId = mesh.addPointSubset("lockedPoints");
    forAll(pointType_, bpI)
    {
        if( pointType_[bpI] & LOCKED )
            mesh.addPointToSubset(lpId, bPoints[bpI]);
        if( pointType_[bpI] & SPECIAL )
            mesh.addPointToSubset(spId, bPoints[bpI]);
    }

    mesh.write();
    Info << "Exitting after surface analysis" << endl;
    ::exit(0);
    # endif
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceMapper::meshSurfaceMapper
(
    const meshSurfaceEngine& mse,
    const meshOctree& octree,
    meshSurfaceEngineModifier* sModPtr
)
:
    surfaceEngine_(mse),
    meshOctree_(octree),
    surfaceEnginePartitionerPtr_(NULL),
    deletePartitioner_(true),
    surfPartitionerPtr_(NULL),
    surfaceModifierPtr_(sModPtr),
    deleteModifier_(!sModPtr),
    pointType_(),
    faceType_()
{
    if( Pstream::parRun() )
    {
        //- allocate bpAtProcs and other addressing
        //- this is done here to prevent possible deadlocks
        surfaceEngine_.bpAtProcs();
    }
}

meshSurfaceMapper::meshSurfaceMapper
(
    const meshSurfacePartitioner& mPart,
    const meshOctree& octree,
    meshSurfaceEngineModifier* sModPtr
)
:
    surfaceEngine_(mPart.surfaceEngine()),
    meshOctree_(octree),
    surfaceEnginePartitionerPtr_(&mPart),
    deletePartitioner_(false),
    surfPartitionerPtr_(NULL),
    surfaceModifierPtr_(sModPtr),
    deleteModifier_(!sModPtr),
    pointType_(),
    faceType_()
{
    if( Pstream::parRun() )
    {
        //- allocate bpAtProcs and other addressing
        //- this is done here to prevent possible deadlocks
        surfaceEngine_.bpAtProcs();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceMapper::~meshSurfaceMapper()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
