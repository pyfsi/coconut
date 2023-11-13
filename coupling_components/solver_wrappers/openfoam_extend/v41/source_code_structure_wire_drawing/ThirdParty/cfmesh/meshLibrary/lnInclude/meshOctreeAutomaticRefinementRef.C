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

#include "meshOctreeAutomaticRefinement.H"
#include "triSurface.H"
#include "demandDrivenData.H"
#include "triSurfacePartitioner.H"
#include "triSurfaceCurvatureEstimator.H"
#include "meshOctreeAddressing.H"
#include "triSurf.H"
#include "helperFunctions.H"
#include "meshOctreeModifier.H"
#include "meshOctreeCreator.H"

#include "Map.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGAutoRef

# ifdef DEBUGAutoRef
#include "pointSet.H"
#include "IOdictionary.H"
#include "triSurfModifier.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeAutomaticRefinement::activateHexRefinement()
{
    hexRefinement_ = true;
}

void meshOctreeAutomaticRefinement::automaticRefinement()
{
    if( !maxRefLevel_ )
        return;

    Info << "Performing automatic refinement" << endl;

    curvatureRefinement();

    proximityRefinement();

    Info << "Finished with automatic refinement" << endl;
}

bool meshOctreeAutomaticRefinement::curvatureRefinement()
{
    Info << "Starting curvature refinement" << endl;
    labelList refineBox(octree_.numberOfLeaves(), 0);
    labelLongList refCandidates;
    List<direction> maxRefLevel;

    bool refine(false);

    //- refined based on the estimated surfac curvature
    findActiveBoxes(CURVATURE, refCandidates, maxRefLevel);

    if( returnReduce(refCandidates.size(), sumOp<label>()) )
    {
        while( refineBasedOnCurvature(refineBox, refCandidates, maxRefLevel) )
        {
            refineSelectedBoxes
            (
                refineBox,
                refCandidates
            );

            findMaxRefLevels(CURVATURE, refCandidates, maxRefLevel);

            refine = true;
        }
    }

    //- refine based on th shortest edge length in the surface mesh
    findActiveBoxes(SHORTESTEDGELENGTH, refCandidates, maxRefLevel);

    if( returnReduce(refCandidates.size(), sumOp<label>()) )
    {
        while( refineBasedOnEdgeLength(refineBox, refCandidates, maxRefLevel) )
        {
            refineSelectedBoxes
            (
                refineBox,
                refCandidates
            );

            findMaxRefLevels(SHORTESTEDGELENGTH, refCandidates, maxRefLevel);

            refine = true;
        }
    }

    Info << "Finished curvature refinement" << endl;

    return refine;
}

bool meshOctreeAutomaticRefinement::proximityRefinement()
{
    bool refine(false);

    Info << "Starting proximity refinement" << endl;

    labelList refineBox(octree_.numberOfLeaves(), 0);
    List<direction> maxRefLevel;
    labelLongList refCandidates;

    //- corner proximity
    findActiveBoxes(CORNERPROXIMITY, refCandidates, maxRefLevel);

    if( returnReduce(refCandidates.size(), sumOp<label>()) )
    {
        while
        (
            refineBasedOnContainedCorners(refineBox, refCandidates, maxRefLevel)
        )
        {
            refineSelectedBoxes
            (
                refineBox,
                refCandidates
            );

            findMaxRefLevels(CORNERPROXIMITY, refCandidates, maxRefLevel);

            refine = true;
        }
    }

    //- refinement based on distinct regions
    findActiveBoxes(DISTINCTREGIONS, refCandidates, maxRefLevel);

    if( returnReduce(refCandidates.size(), sumOp<label>()) )
    {
        while
        (
            refineBasedOnContainedPartitions
            (
                refineBox,
                refCandidates,
                maxRefLevel
            )
        )
        {
            refineSelectedBoxes
            (
                refineBox,
                refCandidates
            );

            findMaxRefLevels(CORNERPROXIMITY, refCandidates, maxRefLevel);

            refine = true;
        }
    }

    //- refinement based on edge proximity criteria
    findActiveBoxes(EDGEPROXIMITY, refCandidates, maxRefLevel);
    if( returnReduce(refCandidates.size(), sumOp<label>()) )
    {
        while
        (
            refineBasedOnEdgeProximityTests
            (
                refineBox,
                refCandidates,
                maxRefLevel
            )
        )
        {
            refineSelectedBoxes
            (
                refineBox,
                refCandidates
            );

            findMaxRefLevels(EDGEPROXIMITY, refCandidates, maxRefLevel);

            refine = true;
        }
    }

    //- refinement based on surface proximity criteria
    findActiveBoxes(SURFACEPROXIMITY, refCandidates, maxRefLevel);
    if( returnReduce(refCandidates.size(), sumOp<label>()) )
    {
        while
        (
            refineBasedOnSurfaceProximityTests
            (
                refineBox,
                refCandidates,
                maxRefLevel
            )
        )
        {
            refineSelectedBoxes
            (
                refineBox,
                refCandidates
            );

            findMaxRefLevels(SURFACEPROXIMITY, refCandidates, maxRefLevel);

            refine = true;
        }
    }

    //- refinement based on ray-casting criteria
    findActiveBoxes(RAYCASTING, refCandidates, maxRefLevel);
    if( returnReduce(refCandidates.size(), sumOp<label>()) )
    {
        while
        (
            refineBasedOnRayCasting
            (
                refineBox,
                refCandidates,
                maxRefLevel
            )
        )
        {
            refineSelectedBoxes
            (
                refineBox,
                refCandidates
            );

            findMaxRefLevels(RAYCASTING, refCandidates, maxRefLevel);

            refine = true;
        }
    }

    Info << "Finished proximity refinement" << endl;

    reduce(refine, maxOp<bool>());

    return refine;
}

bool meshOctreeAutomaticRefinement::refineBasedOnContainedCorners
(
    labelList& refineBox,
    const labelLongList& refCandidates,
    const List<direction>& maxRefLevel
)
{
    meshOctreeModifier octreeModifier(octree_);
    const LongList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();
    const boundBox& rootBox = octree_.rootBox();
    const triSurf& surface = octree_.surface();
    const pointField& points = surface.points();
    const triSurfacePartitioner& sPart = this->partitioner();

    //- find leaves which contains corner nodes
    std::map<label, DynList<label> > cornersInLeaf;
    const labelList& corners = sPart.corners();

    label nMarked(0);

    forAll(corners, cornerI)
    {
        const label cLabel =
        octree_.findLeafContainingVertex(points[corners[cornerI]]);

        if( cLabel < 0 )
            continue;

        cornersInLeaf[cLabel].append(corners[cornerI]);
    }

    DynList<label> leavesInBox;
    # ifdef USE_OMP
    # pragma omp parallel for if( refCandidates.size() > 1000 ) \
    private(leavesInBox) shared(cornersInLeaf) \
    reduction(+ : nMarked) schedule(dynamic, 20)
    # endif
    forAll(refCandidates, refI)
    {
        const label leafI = refCandidates[refI];
        if( leaves[leafI]->level() >= maxRefLevel[leafI] )
            continue;

        std::map<label, DynList<label> >::const_iterator currIt =
            cornersInLeaf.find(leafI);
        if( currIt == cornersInLeaf.end() )
            continue;

        if( currIt->second.size() > 1 )
        {
            refineBox[leafI] = 1;
            ++nMarked;
            continue;
        }

        // check if there exist some corners in the neighbour boxes
        // refine the box if some corners are found in the neighbouring boxes
        const point c = points[currIt->second[0]];
        const scalar r =
            nBoxesOverFeature_ * 1.732 * leaves[leafI]->size(rootBox);
        const scalar rSq = r * r;

        const boundBox bb(c - point(r, r, r), c + point(r, r, r));

        leavesInBox.clear();
        octree_.findLeavesContainedInBox(bb, leavesInBox);

        forAll(leavesInBox, i)
        {
            const label nei = leavesInBox[i];

            if( nei < 0 )
                continue;

            if( nei == leafI )
                continue;

            std::map<label, DynList<label> >::const_iterator it =
                cornersInLeaf.find(nei);
            if( it == cornersInLeaf.end() )
                continue;

            const DynList<label>& cornersInNei = it->second;

            forAll(cornersInNei, j)
            {
                if( magSqr(points[cornersInNei[j]] - c) < rSq )
                {
                    ++nMarked;
                    refineBox[nei] = 1;
                    refineBox[leafI] = 1;
                    break;
                }
            }
        }
    }

    nMarked += ensureCornerAndEdgeConsistency(refineBox, maxRefLevel);

    reduce(nMarked, sumOp<label>());
    Info << nMarked << " boxes marked by the corner criteria" << endl;

    if( nMarked )
        return true;

    return false;
}

bool meshOctreeAutomaticRefinement::refineBasedOnContainedPartitions
(
    labelList& refineBox,
    const labelLongList& refCandidates,
    const List<direction>& maxRefLevel
)
{
    const boundBox& rootBox = octree_.rootBox();
    const triSurfacePartitioner& sPart = this->partitioner();

    //- find leaves which contains corner nodes
    const List<labelHashSet>& pPatches = sPart.patchPatches();
    const labelList& edgeGroups = sPart.edgeGroups();
    const List<labelHashSet>& eNeiGroups = sPart.edgeGroupEdgeGroups();

    label nMarked(0);

    meshOctreeModifier octreeModifier(octree_);
    const triSurf& surf = octree_.surface();
    const LongList<meshOctreeCube*>& leaves = octreeModifier.leavesAccess();

    DynList<label> patches, eGroups, helper;
    # ifdef USE_OMP
    # pragma omp parallel for if( refCandidates.size() > 1000 ) \
    private(patches, eGroups, helper) \
    reduction(+ : nMarked) schedule(dynamic, 20)
    # endif
    forAll(refCandidates, refI)
    {
        const label leafI = refCandidates[refI];
        if( !leaves[leafI]->hasContainedElements() )
            continue;
        if( leaves[leafI]->level() >= maxRefLevel[leafI] )
            continue;

        const meshOctreeCubeBasic& oc = *leaves[leafI];
        const point c = oc.centre(rootBox);
        const scalar s = nBoxesOverFeature_ * 1.733 * oc.size(rootBox);
        const boundBox bb(c - point(s, s, s), c + point(s, s, s));

        //- find triangle patches contained in this box
        octree_.findTrianglesInBox(bb, helper);
        patches.clear();
        forAll(helper, i)
            patches.appendIfNotIn(surf[helper[i]].region());

        //- find edge partitions contained in this box
        helper.clear();
        octree_.findEdgesInBox(bb, helper);
        eGroups.clear();
        forAll(helper, i)
        {
            if( helper[i] < 0 || edgeGroups[helper[i]] < 0 )
                continue;

            eGroups.appendIfNotIn(edgeGroups[helper[i]]);
        }

        # ifdef DEBUGAutoRef
        Info << "patches for leaf " << leafI << " are " << patches << endl;
        # endif

        bool refine(false);
        forAll(patches, patchI)
        {
            for(label patchJ=(patchI+1);patchJ<patches.size();++patchJ)
                if( !pPatches[patches[patchI]].found(patches[patchJ]) )
                {
                    # ifdef DEBUGAutoRef
                    Info << "2.Here" << endl;
                    # endif

                    refine = true;
                    break;
                }
        }

        forAll(eGroups, egI)
        {
            for(label egJ=egI+1;egJ<eGroups.size();++egJ)
                if( !eNeiGroups[eGroups[egI]].found(eGroups[egJ]) )
                {
                    refine = true;
                    break;
                }
        }

        if( refine )
        {
            # ifdef DEBUGAutoRef
            Info << "Selecting leaf " << leafI
                << " for auto-refinement" << endl;
            # endif

            ++nMarked;
            refineBox[leafI] = 1;
        }
    }

    nMarked += ensureCornerAndEdgeConsistency(refineBox, maxRefLevel);

    reduce(nMarked, sumOp<label>());
    Info << nMarked << " boxed marked by partitioning criteria" << endl;

    if( nMarked )
        return true;

    return false;
}

bool meshOctreeAutomaticRefinement::refineBasedOnCurvature
(
    labelList& refineBox,
    const labelLongList& refCandidates,
    const List<direction>& maxRefLevel
)
{
    const triSurfaceCurvatureEstimator& curv = curvature();

    const boundBox& rootBox = octree_.rootBox();

    label nMarked(0);
    DynList<label> containedTrias;
    # ifdef USE_OMP
    # pragma omp parallel for private(containedTrias) \
    reduction(+ : nMarked) schedule(dynamic, 100)
    # endif
    forAll(refCandidates, refI)
    {
        const label leafI = refCandidates[refI];

        if( !octree_.hasContainedTriangles(leafI) )
            continue;

        const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);
        if( oc.level() >= maxRefLevel[leafI] )
            continue;

        //- search for the minimum curvature radius at surface triangles
        octree_.containedTriangles(leafI, containedTrias);

        scalar maxCurv(0.0);
        forAll(containedTrias, i)
        {
            maxCurv =
                Foam::max
                (
                    maxCurv,
                    max
                    (
                        mag(curv.maxCurvatureAtTriangle(containedTrias[i])),
                        mag(curv.minCurvatureAtTriangle(containedTrias[i]))
                    )
                );
        }

        //- check the edge curvature
        if( octree_.hasContainedEdges(leafI) )
        {
            octree_.containedEdges(leafI, containedTrias);

            forAll(containedTrias, i)
            {
                maxCurv =
                    Foam::max
                    (
                        maxCurv,
                        mag(curv.curvatureAtEdge(containedTrias[i]))
                    );
            }
        }

        if( oc.size(rootBox) > 0.2835 / (nBoxesOverFeature_*(maxCurv + SMALL)) )
        {
            refineBox[leafI] = 1;
            ++nMarked;
        }
    }

    nMarked += ensureCornerAndEdgeConsistency(refineBox, maxRefLevel);

    reduce(nMarked, sumOp<label>());
    Info << nMarked << " boxes marked by curvature criteria!" << endl;

    if( nMarked )
        return true;

    return false;
}

bool meshOctreeAutomaticRefinement::refineBasedOnEdgeLength
(
    labelList& refineBox,
    const labelLongList& refCandidates,
    const List<direction>& maxRefLevel
)
{
    const triSurf& surf = octree_.surface();
    const pointField& points = surf.points();
    const edgeLongList& edges = surf.edges();
    const VRWGraph& faceEdges  = surf.facetEdges();

    scalarList edgeLength(edges.size());
    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(edges, eI)
        edgeLength[eI] = edges[eI].mag(points);

    const boundBox& rootBox = octree_.rootBox();

    label nMarked(0);
    DynList<label> containedTrias;
    # ifdef USE_OMP
    # pragma omp parallel for private(containedTrias) \
    reduction(+ : nMarked) schedule(dynamic, 100)
    # endif
    forAll(refCandidates, refI)
    {
        const label leafI = refCandidates[refI];

        if( !octree_.hasContainedTriangles(leafI) )
            continue;

        const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);

        if( oc.level() >= maxRefLevel[leafI] )
            continue;

        const scalar cubeSize = oc.size(rootBox);

        //- search for the minimum edge length within the cube
        octree_.containedTriangles(leafI, containedTrias);

        scalar minEdgeLength(VGREAT);
        forAll(containedTrias, i)
        {
            const label triI = containedTrias[i];

            forAllRow(faceEdges, triI, teI)
                minEdgeLength =
                    Foam::min(minEdgeLength, edgeLength[faceEdges(triI, teI)]);
        }

        if( cubeSize > minEdgeLength )
        {
            refineBox[leafI] = 1;
            ++nMarked;
        }
    }

    nMarked += ensureCornerAndEdgeConsistency(refineBox, maxRefLevel);

    reduce(nMarked, sumOp<label>());
    Info << nMarked << " boxes marked by the edge length criteria!" << endl;

    if( nMarked )
        return true;

    return false;
}

bool meshOctreeAutomaticRefinement::refineBasedOnSurfaceProximityTests
(
    labelList& refineBox,
    const labelLongList& refCandidates,
    const List<direction>& maxRefLevel
)
{
    const boundBox& rootBox = octree_.rootBox();
    const triSurf& surf = octree_.surface();
    const pointField& pts = surf.points();
    const edgeLongList& edges = surf.edges();
    const VRWGraph& faceEdges = surf.facetEdges();
    const VRWGraph& edgeFaces = surf.edgeFacets();

    const List<DynList<label> >& triangleTriangles = trianglesInStickySpace();

    label nMarked(0);
    DynList<label> helper, ct;
    # ifdef USE_OMP
    # pragma omp parallel for private(helper, ct) reduction(+ : nMarked) \
    schedule(dynamic, 20)
    # endif
    forAll(refCandidates, refI)
    {
        const label leafI = refCandidates[refI];

        if( !octree_.hasContainedTriangles(leafI) )
            continue;

        const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);
        if( oc.level() >= maxRefLevel[leafI] )
            continue;

        const point c = oc.centre(rootBox);
        const scalar rCube = 0.866 * oc.size(rootBox);

        //- check if all leaves at vertex are intersected by the surface
        FixedList<point, 8> vertices;
        oc.vertices(rootBox, vertices);

        for(label vI=0;vI<8;++vI)
        {
            point np;
            label nt, patch;
            scalar dSq;
            octree_.findNearestSurfacePoint(np, dSq, nt, patch, vertices[vI]);

            if( dSq < SMALL * rCube )
                continue;

            FixedList<label, 8> leavesAtVertex;
            octree_.findLeavesForCubeVertex(leafI, vI, leavesAtVertex);

            bool allTriIntersected(true), allEdgeIntersected(true);
            ct.clear();

            forAll(leavesAtVertex, i)
            {
                if( leavesAtVertex[i] < 0 )
                    continue;

                if( !octree_.hasContainedTriangles(leavesAtVertex[i]) )
                {
                    allTriIntersected = false;
                }
                else
                {
                    helper.clear();
                    octree_.containedTriangles(leavesAtVertex[i], helper);

                    forAll(helper, i)
                        ct.appendIfNotIn(helper[i]);
                }

                if( !octree_.hasContainedEdges(leavesAtVertex[i]) )
                    allEdgeIntersected = false;
            }

            forAll(ct, i)
            {
                if( !allTriIntersected )
                    break;

                vector n = surf[ct[i]].normal(pts);
                const scalar magn = mag(n);
                if( magn < VSMALL )
                    continue;
                n /= (mag(n) + VSMALL);

                for(label j=i+1;j<ct.size();++j)
                {
                    vector nn = surf[ct[j]].normal(pts);
                    const scalar magnn = mag(nn);

                    if( magnn < VSMALL )
                        continue;

                    nn /= (mag(nn) + VSMALL);

                    if( mag(nn & n) < proximityAngleTol_ )
                    {
                        allTriIntersected = false;
                        break;
                    }
                }
            }

            if( allTriIntersected || allEdgeIntersected )
            {
                # ifdef DEBUGAutoRef
                triSurf boxSurf;
                triSurfModifier tMod(boxSurf);
                tMod.patchesAccess().setSize(2);
                tMod.patchesAccess()[0].name() = "box";
                tMod.patchesAccess()[1].name() = "surfTriangles";

                std::map<label, label> pointLabel;

                forAll(helper, ncI)
                {
                    const label leafI = helper[ncI];
                    const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);
                    oc.vertices(rootBox, vertices);
                    const label n = boxSurf.nPoints();
                    forAll(vertices, i)
                        boxSurf.appendVertex(vertices[i]);

                    boxSurf.appendTriangle(labelledTri(n, n+1, n+2, 0));
                    boxSurf.appendTriangle(labelledTri(n+1, n+3, n+2, 0));

                    boxSurf.appendTriangle(labelledTri(n+4, n+5, n+6, 0));
                    boxSurf.appendTriangle(labelledTri(n+5, n+7, n+6, 0));

                    boxSurf.appendTriangle(labelledTri(n, n+6, n+2, 0));
                    boxSurf.appendTriangle(labelledTri(n, n+4, n+6, 0));

                    boxSurf.appendTriangle(labelledTri(n+1, n+3, n+7, 0));
                    boxSurf.appendTriangle(labelledTri(n+1, n+7, n+5, 0));

                    boxSurf.appendTriangle(labelledTri(n, n+1, n+5, 0));
                    boxSurf.appendTriangle(labelledTri(n, n+5, n+4, 0));

                    boxSurf.appendTriangle(labelledTri(n+2, n+6, n+7, 0));
                    boxSurf.appendTriangle(labelledTri(n+2, n+7, n+3, 0));

                    DynList<label> containedTriangles;
                    octree_.containedTriangles(leafI, containedTriangles);
                    forAll(containedTriangles, ctI)
                    {
                        const labelledTri orig = surf[containedTriangles[ctI]];

                        labelledTri newTri;
                        newTri.region() = 1;
                        forAll(orig, i)
                        {
                            if( pointLabel.find(orig[i]) == pointLabel.end() )
                            {
                                pointLabel[orig[i]] = boxSurf.nPoints();
                                boxSurf.appendVertex(pts[orig[i]]);
                            }

                            newTri[i] = pointLabel[orig[i]];
                        }

                        boxSurf.appendTriangle(newTri);
                    }
                }

                const fileName fName
                (
                    "box_" +
                    help::labelToText(leafI) +
                    "_vertex_" +
                    help::labelToText(vI) +
                    ".fms"
                );

                boxSurf.writeSurface(fName);
                # endif

                ++nMarked;
                refineBox[leafI] = 1;
                break;
            }
        }

        if( refineBox[leafI] )
            continue;


        const scalar s = 2.0 * nBoxesOverFeature_ * rCube;
        const boundBox bb(c - point(s, s, s), c + point(s, s, s));

        //- find triangles in range
        helper.clear();
        octree_.findTrianglesInBox(bb, helper);

        labelHashSet triaInRange(helper.size());
        const scalar rangeSq = s * s;
        forAll(helper, i)
            triaInRange.insert(helper[i]);

        Map<label> elGroup(triaInRange.size());
        Map<label> testEdge(triaInRange.size());

        label nFaceGroups(0);

        DynList<label> front;

        forAllConstIter(labelHashSet, triaInRange, it)
        {
            if( !elGroup.found(it.key()) )
            {
                front.clear();
                front.append(it.key());
                elGroup.insert(it.key(), nFaceGroups);

                while( front.size() )
                {
                    const label fLabel = front.removeLastElement();

                    forAllRow(faceEdges, fLabel, feI)
                    {
                        const label edgeI = faceEdges(fLabel, feI);

                        //- check if the edge intersects the bounding box
                        if( testEdge.found(edgeI) )
                        {
                            if( !testEdge[edgeI] )
                                continue;
                        }
                        else
                        {
                            const point& s = pts[edges[edgeI][0]];
                            const point& e = pts[edges[edgeI][1]];
                            const point np =
                                help::nearestPointOnTheEdge(s, e, c);

                            if( magSqr(np - c) < rangeSq )
                            {
                                testEdge.insert(edgeI, 1);
                            }
                            else
                            {
                                testEdge.insert(edgeI, 0);
                                continue;
                            }
                        }

                        forAllRow(edgeFaces, edgeI, efI)
                        {
                            const label nei = edgeFaces(edgeI, efI);
                            if
                            (
                                triaInRange.found(nei) &&
                                !elGroup.found(nei)
                            )
                            {
                                elGroup.insert(nei, nFaceGroups);
                                front.append(nei);
                            }
                        }
                    }

                    const DynList<label>& neiTriangles =
                        triangleTriangles[fLabel];

                    forAll(neiTriangles, i)
                    {
                        const label nei = neiTriangles[i];

                        if
                        (
                            triaInRange.found(nei) &&
                            !elGroup.found(nei)
                        )
                        {
                            elGroup.insert(nei, nFaceGroups);
                            front.append(nei);
                        }
                    }
                }

                ++nFaceGroups;
            }
        }

        DynList<vector> avgNormal(nFaceGroups, vector::zero);
        DynList<label> nAppearances(nFaceGroups, 0);
        forAllConstIter(Map<label>, elGroup, it)
        {
            const label groupI = it();

            vector n = surf[it.key()].normal(pts);
            n /= (mag(n) + VSMALL);

            if( (avgNormal[groupI] & n) > -SMALL )
            {
                avgNormal[groupI] += n;
            }
            else
            {
                avgNormal[groupI] -= n;
            }

            ++nAppearances[groupI];
        }

        forAll(avgNormal, i)
            for(label j=0;j<avgNormal.size();++j)
                if( mag(avgNormal[i] & avgNormal[j]) < proximityAngleTol_ )
                {
                    --nFaceGroups;
                }

        //- refine boxes with more than two face groups
        if( nFaceGroups > 1 )
        {
            # ifdef DEBUGAutoRef
            triSurf boxSurf;
            triSurfModifier tMod(boxSurf);
            tMod.patchesAccess().setSize(nFaceGroups+1);
            tMod.patchesAccess()[0].name() = "box";
            for(label i=0;i<nFaceGroups;++i)
                tMod.patchesAccess()[i+1].name() =
                    "faceGroup_"+help::labelToText(i);
            tMod.pointsAccess().setSize(8);
            FixedList<point, 8> vertices;
            oc.vertices(rootBox, vertices);
            forAll(vertices, i)
                tMod.pointsAccess()[i] = vertices[i];
            boxSurf.appendTriangle(labelledTri(0, 1, 2, 0));
            boxSurf.appendTriangle(labelledTri(1, 3, 2, 0));

            boxSurf.appendTriangle(labelledTri(4, 5, 6, 0));
            boxSurf.appendTriangle(labelledTri(5, 7, 6, 0));

            boxSurf.appendTriangle(labelledTri(0, 6, 2, 0));
            boxSurf.appendTriangle(labelledTri(0, 4, 6, 0));

            boxSurf.appendTriangle(labelledTri(1, 3, 7, 0));
            boxSurf.appendTriangle(labelledTri(1, 7, 5, 0));

            boxSurf.appendTriangle(labelledTri(0, 1, 5, 0));
            boxSurf.appendTriangle(labelledTri(0, 5, 4, 0));

            boxSurf.appendTriangle(labelledTri(2, 6, 7, 0));
            boxSurf.appendTriangle(labelledTri(2, 7, 3, 0));

            std::map<label, label> pointLabel;
            forAllConstIter(Map<label>, elGroup, it)
            {
                const labelledTri orig = surf[it.key()];

                labelledTri newTri;
                newTri.region() = it() + 1;
                forAll(orig, i)
                {
                    if( pointLabel.find(orig[i]) == pointLabel.end() )
                    {
                        pointLabel[orig[i]] = boxSurf.nPoints();
                        boxSurf.appendVertex(pts[orig[i]]);
                    }

                    newTri[i] = pointLabel[orig[i]];
                }

                boxSurf.appendTriangle(newTri);
            }

            boxSurf.writeSurface("box_"+help::labelToText(leafI)+".fms");
            # endif

            ++nMarked;
            refineBox[leafI] = 1;
        }
    }

    nMarked += ensureCornerAndEdgeConsistency(refineBox, maxRefLevel);

    reduce(nMarked, sumOp<label>());
    Info << nMarked << " boxed marked by surface proximity criteria" << endl;

    if( nMarked != 0 )
        return true;

    return false;
}

bool meshOctreeAutomaticRefinement::refineBasedOnEdgeProximityTests
(
    labelList& refineBox,
    const labelLongList& refCandidates,
    const List<direction>& maxRefLevel
)
{
    const boundBox& rootBox = octree_.rootBox();
    const triSurf& surf = octree_.surface();

    const pointField& pts = surf.points();
    const edgeLongList& edges = surf.edges();
    const VRWGraph& pointEdges = surf.pointEdges();

    const List<DynList<label> >& edgeEdges = edgesInStickySpace();

    label nMarked(0);
    DynList<label> helper;
    # ifdef USE_OMP
    # pragma omp parallel for private(helper) reduction(+ : nMarked) \
    schedule(dynamic, 20)
    # endif
    forAll(refCandidates, refI)
    {
        const label leafI = refCandidates[refI];

        if( !octree_.hasContainedEdges(leafI) )
            continue;

        const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);
        if( oc.level() >= maxRefLevel[leafI] )
            continue;

        const point c = oc.centre(rootBox);
        const scalar rCube = 0.866 * oc.size(rootBox);

        const scalar s = 2.0 * nBoxesOverFeature_ * rCube;
        const boundBox bb(c - point(s, s, s), c + point(s, s, s));

        labelHashSet edgesInRange(100);

        //- find edges contained in the neighbourhood
        helper.clear();
        octree_.findEdgesInBox(bb, helper);
        forAll(helper, i)
            edgesInRange.insert(helper[i]);

        const scalar rangeSq = sqr(s);

        Map<label> elGroup(edgesInRange.size());
        Map<label> pointTest(edgesInRange.size());

        label nEdgeGroups(0);

        DynList<label> front;

        forAllConstIter(labelHashSet, edgesInRange, it)
        {
            if( !elGroup.found(it.key()) )
            {
                front.clear();
                front.append(it.key());
                elGroup.insert(it.key(), nEdgeGroups);

                while( front.size() )
                {
                    const label eLabel = front.removeLastElement();
                    const edge& e = edges[eLabel];

                    for(label i=0;i<2;++i)
                    {
                        if( pointTest.found(e[i]) )
                        {
                            if( !pointTest[e[i]] )
                                continue;
                        }
                        else
                        {
                            if( magSqr(pts[e[i]] - c) < rangeSq )
                            {
                                pointTest.insert(e[i], 1);
                            }
                            else
                            {
                                pointTest.insert(e[i], 0);
                                continue;
                            }
                        }

                        forAllRow(pointEdges, e[i], peI)
                        {
                            const label nei = pointEdges(e[i], peI);

                            if( edgesInRange.found(nei) && !elGroup.found(nei) )
                            {
                                elGroup.insert(nei, nEdgeGroups);
                                front.append(nei);
                            }
                        }
                    }

                    //- check connection over objects within the sticky space
                    const DynList<label>& neiEdges = edgeEdges[eLabel];

                    forAll(neiEdges, i)
                    {
                        const label nei = neiEdges[i];

                        if( edgesInRange.found(nei) && !elGroup.found(nei) )
                        {
                            elGroup.insert(nei, nEdgeGroups);
                            front.append(nei);
                        }
                    }
                }

                ++nEdgeGroups;
            }
        }

        if( nEdgeGroups > 1 )
        {
            ++nMarked;
            refineBox[leafI] = 1;
        }
    }

    nMarked += ensureCornerAndEdgeConsistency(refineBox, maxRefLevel);

    reduce(nMarked, sumOp<label>());
    Info << nMarked << " boxed marked by edge proximity criteria" << endl;

    if( nMarked != 0 )
        return true;

    return false;
}

bool meshOctreeAutomaticRefinement::refineBasedOnRayCasting
(
    labelList& refineBox,
    const labelLongList& refCandidates,
    const List<direction>& maxRefLevel
)
{
    const boundBox& rootBox = octree_.rootBox();
    const triSurf& surf = octree_.surface();
    const pointField& pts = surf.points();

    label nMarked(0);
    DynList<label> helper, ct;
    # ifdef USE_OMP
    # pragma omp parallel for private(helper,ct) reduction(+ : nMarked) \
    schedule(dynamic, 20)
    # endif
    forAll(refCandidates, refI)
    {
        const label leafI = refCandidates[refI];

        if( !octree_.hasContainedTriangles(leafI) )
            continue;

        const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);
        if( oc.level() >= maxRefLevel[leafI] )
            continue;

        const point c = oc.centre(rootBox);
        const scalar rCube = 0.866 * oc.size(rootBox);
        const scalar s = 2.0 * nBoxesOverFeature_ * rCube;

        const vector dvec = 0.55 * vector(s, s, s);

        //- find triangles in range
        helper.clear();
        octree_.containedTriangles(leafI, helper);

        point pMin, pMax;
        oc.cubeBox(rootBox, pMin, pMax);
        const boundBox bb(pMin-SMALL*(pMax-pMin), pMax+SMALL*(pMax-pMin));

        const scalar rangeSq = s * s;

        forAll(helper, i)
        {
            const label triI = helper[i];
            const labelledTri& tri = surf[triI];

            point sp = tri.centre(pts);
            vector n = tri.normal(pts);
            const scalar magn = mag(n);

            if( magn < VSMALL )
                continue;

            n /= (magn + VSMALL);

            if( !bb.contains(sp) )
            {
                sp = help::nearestPointOnTheTriangle(triI, surf, c);
            }

            const point opos = sp + s * n;
            const point oneg = sp - s * n;

            point cp = 0.5 * (opos + sp);

            const boundBox bbPos(cp - dvec, cp + dvec);
            ct.clear();
            octree_.findTrianglesInBox(bbPos, ct);

            forAll(ct, i)
            {
                if( ct[i] == triI )
                    continue;

                vector nn = surf[ct[i]].normal(pts);
                nn /= (mag(nn) + VSMALL);

                if( mag(nn & n) < proximityAngleTol_ )
                    continue;

                point intersection;

                if
                (
                    help::triLineIntersection
                    (
                        surf,
                        ct[i],
                        sp,
                        opos,
                        intersection
                    )
                )
                {
                    const scalar dSq = magSqr(intersection - sp);

                    if( dSq < sqr(stickyDistance_[triI]) )
                        continue;
                    if( dSq < sqr(SMALL * rCube) )
                        continue;

                    if( dSq < rangeSq )
                    {
                        ++nMarked;
                        refineBox[leafI] = 1;
                        break;
                    }
                }
            }

            if( refineBox[leafI] )
                continue;

            cp = 0.5 * (sp + oneg);
            const boundBox bbNeg(cp - dvec, cp + dvec);
            ct.clear();
            octree_.findTrianglesInBox(bbNeg, ct);

            forAll(ct, i)
            {
                if( ct[i] == triI )
                    continue;

                vector nn = surf[ct[i]].normal(pts);
                nn /= (mag(nn) + VSMALL);

                if( mag(nn & n) < proximityAngleTol_ )
                    continue;

                point intersection;

                if
                (
                    help::triLineIntersection
                    (
                        surf,
                        ct[i],
                        sp,
                        oneg,
                        intersection
                    )
                )
                {
                    const scalar dSq = magSqr(intersection - sp);

                    if( dSq < sqr(stickyDistance_[triI]) )
                        continue;
                    if( dSq < sqr(SMALL * rCube) )
                        continue;

                    if( dSq < rangeSq )
                    {
                        ++nMarked;
                        refineBox[leafI] = 1;
                        break;
                    }
                }
            }

            if( refineBox[leafI] )
                break;
        }
    }

    nMarked += ensureCornerAndEdgeConsistency(refineBox, maxRefLevel);

    reduce(nMarked, sumOp<label>());
    Info << nMarked << " boxed marked by ray casting criteria" << endl;

    if( nMarked != 0 )
        return true;

    return false;
}

label meshOctreeAutomaticRefinement::ensureCornerAndEdgeConsistency
(
    labelList& refineBox,
    const List<direction>& maxRefLevel
)
{
    const triSurf& surf = octree_.surface();
    DynList<label> ct;

    boolList containsFeature(refineBox.size()), refineFeature(refineBox.size());

    //- mark the boxes containing features
    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100) private(ct)
    # endif
    forAll(refineBox, leafI)
    {
        containsFeature[leafI] = false;
        refineFeature[leafI] = false;

        if( !octree_.hasContainedTriangles(leafI) )
            continue;

        if( refineBox[leafI] )
            continue;

        //- check the number of patches within the box
        ct.clear();
        octree_.containedTriangles(leafI, ct);

        //- check the number of patches within the box
        DynList<label> patches;
        forAll(ct, i)
            patches.appendIfNotIn(surf[ct[i]].region());

        if( patches.size() > 1 )
        {
            containsFeature[leafI] = true;
        }
    }

    //- check which boxes with features have a refine neighbour
    //- or a neighbour at the higher refinement level
    LongList<meshOctreeCubeCoordinates> transferLeaves;
    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100) private(ct)
    # endif
    forAll(refineBox, leafI)
    {
        if( !containsFeature[leafI] )
            continue;

        const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);
        if( oc.level() >= maxRefLevel[leafI] )
            continue;

        //- do not allow the box to be at the lower refinement level than
        //- its neighbour that contains less patches
        octree_.findAllLeafNeighbours(leafI, ct);

        bool transferLeaf(false);

        forAll(ct, i)
        {
            const label nei = ct[i];

            if( nei == -1 )
                continue;

            if( nei == meshOctreeCube::OTHERPROC )
            {
                transferLeaf = true;
                continue;
            }

            if
            (
                refineBox[nei] ||
                (oc.level() < octree_.returnLeaf(nei).level())
            )
            {
                refineFeature[leafI] |= 1;
            }
        }

        if( transferLeaf )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            transferLeaves.append(oc.coordinates());
        }
    }

    if( returnReduce(transferLeaves.size(), sumOp<label>()) )
    {
        LongList<meshOctreeCubeCoordinates> receivedLeaves;
        octree_.exchangeRequestsWithNeighbourProcessors
        (
            transferLeaves,
            receivedLeaves
        );

        //- find neighbours at other processors
        transferLeaves.clear();

        forAll(receivedLeaves, i)
        {
            const meshOctreeCubeCoordinates& cc = receivedLeaves[i];

            DynList<label> neighbours;
            octree_.findAllLeafNeighbours(cc, neighbours);

            forAll(neighbours, i)
            {
                const label nei = neighbours[i];

                if( nei < 0 )
                    continue;

                if
                (
                    refineBox[nei] ||
                    (cc.level() < octree_.returnLeaf(nei).level())
                )
                {
                    transferLeaves.append(cc);
                    break;
                }
            }
        }

        //- exchange refined leaves;
        receivedLeaves.clear();
        octree_.exchangeRequestsWithNeighbourProcessors
        (
            transferLeaves,
            receivedLeaves
        );

        //- mark received leaves for refinement
        forAll(receivedLeaves, i)
        {
            const label leafI =
                octree_.findLeafLabelForPosition(receivedLeaves[i]);

            if( leafI >= 0 )
                refineFeature[leafI] = true;
        }
    }

    //- set the refinement flag
    label nMarked(0);
    forAll(refineFeature, i)
    {
        if( refineFeature[i] )
            ++nMarked;

        refineBox[i] |= refineFeature[i];
    }

    return nMarked;
}

void meshOctreeAutomaticRefinement::refineSelectedBoxes
(
    labelList& refineBox,
    labelLongList& refCandidates
)
{
    deleteDemandDrivenData(octreeAddressingPtr_);

    meshOctreeModifier octreeModifier(octree_);
    LongList<meshOctreeCube*> leaves = octreeModifier.leavesAccess();

    octreeModifier.markAdditionalLayers(refineBox, nBoxesOverFeature_);
    octreeModifier.refineSelectedBoxes(refineBox, hexRefinement_);

    //- find the cubes which have been marked for refinement
    LongList<meshOctreeCubeCoordinates> refinedCubes;
    forAll(refineBox, i)
    {
        if( refineBox[i] )
            refinedCubes.append(leaves[i]->coordinates());
    }
    leaves.setSize(0);

    //- perform load distribution in case of parallel runs
    octreeModifier.loadDistribution();

    //- communicate the cubes selected for refinement with other processors
    LongList<meshOctreeCubeCoordinates> receivedCoordinates;
    if( Pstream::parRun() )
    {
        octree_.exchangeRequestsWithNeighbourProcessors
        (
            refinedCubes,
            receivedCoordinates
        );
    }

    forAll(refinedCubes, i)
        receivedCoordinates.append(refinedCubes[i]);
    refinedCubes.setSize(0);

    //- find the cubes which shall checked in the next iteration
    refCandidates.clear();
    forAll(receivedCoordinates, i)
    {
        const meshOctreeCubeCoordinates& cc = receivedCoordinates[i];

        for(label scI=0;scI<8;++scI)
        {
            const meshOctreeCubeCoordinates child = cc.refineForPosition(scI);

            meshOctreeCube* oc = octreeModifier.findCubeForPosition(child);

            if( !oc || !oc->isLeaf() )
                continue;

            refCandidates.append(oc->cubeLabel());
        }
    }

    refineBox.setSize(octree_.numberOfLeaves());
    refineBox = direction(0);
}

void meshOctreeAutomaticRefinement::findActiveBoxes
(
    const direction refType,
    labelLongList& activeLeaves,
    List<direction>& maxRefLevel
) const
{
    boolList activeLeaf(octree_.numberOfLeaves(), false);

    //- activate all leaves
    if( globalRefStrategies_ & refType )
    {
        activeLeaves.setSize(octree_.numberOfLeaves());

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(activeLeaves, leafI)
        {
            activeLeaf[leafI] = true;
        }
    }
    else
    {
        const triSurf& surf = octree_.surface();
        const geometricSurfacePatchList& patches = surf.patches();

        std::map<word, label> patchToIndex;
        forAll(patches, patchI)
            patchToIndex.insert(std::make_pair(patches[patchI].name(), patchI));

        boolList activeFacets(surf.size(), false);
        boolList activePatch(patches.size(), false);

        forAll(patches, patchI)
        {
            const word& pName = patches[patchI].name();

            std::map<word, direction>::const_iterator it =
                objectRefinementStrategies_.find(pName);
            if
            (
                (it != objectRefinementStrategies_.end()) &&
                (it->second & refType)
            )
            {
                activePatch[patchToIndex[pName]] = true;
            }
        }

        //- set active flag to facets belonging to active patches
        forAll(surf, triI)
        {
            if( activePatch[surf[triI].region()] )
                activeFacets[triI] = true;
        }

        DynList<label> subsetIndices;
        surf.facetSubsetIndices(subsetIndices);
        forAll(subsetIndices, i)
        {
            const label subsetI = subsetIndices[i];
            const word sName = surf.facetSubsetName(subsetI);

            std::map<word, direction>::const_iterator it =
                objectRefinementStrategies_.find(sName);
            if
            (
                (it != objectRefinementStrategies_.end()) &&
                (it->second & refType)
            )
            {
                labelLongList trianglesInSubset;
                surf.facetsInSubset(subsetI, trianglesInSubset);

                forAll(trianglesInSubset, i)
                    activeFacets[trianglesInSubset[i]] = true;
            }
        }

        //- select leaves containing active patches
        DynList<label> containedTriangles;
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100) \
        private(containedTriangles)
        # endif
        forAll(activeLeaf, leafI)
        {
            octree_.containedTriangles(leafI, containedTriangles);

            forAll(containedTriangles, i)
            {
                if( activeFacets[containedTriangles[i]] )
                {
                    activeLeaf[leafI] = true;
                    break;
                }
            }
        }
    }

    //- collect active leaves
    activeLeaves.clear();

    forAll(activeLeaf, leafI)
    {
        if( activeLeaf[leafI] )
            activeLeaves.append(leafI);
    }

    const label nCandidates = returnReduce(activeLeaves.size(), maxOp<label>());
    Info << "Found " << nCandidates << " candidates for refinement" << endl;

    activeLeaf.clear();

    //- find max refinement levels
    findMaxRefLevels(refType, activeLeaves, maxRefLevel);
}

void meshOctreeAutomaticRefinement::findMaxRefLevels
(
    const direction refType,
    const labelLongList& intersectedLeaves,
    List<direction>& maxRefLevel
) const
{
    //- set the size
    maxRefLevel.setSize(octree_.numberOfLeaves());
    maxRefLevel = NONE;

    //- find ref level for local objects
    const triSurf& surf = octree_.surface();
    const geometricSurfacePatchList& patches = surf.patches();

    std::map<word, label> patchToIndex;
    forAll(patches, patchI)
        patchToIndex.insert(std::make_pair(patches[patchI].name(), patchI));

    List<direction> activeFacets(surf.size(), NONE);
    List<direction> activePatch(patches.size(), NONE);

    forAll(patches, patchI)
    {
        const word& pName = patches[patchI].name();

        std::map<word, direction>::const_iterator it =
            objectMaxRefLevel_.find(pName);

        std::map<word, direction>::const_iterator sIt =
            objectRefinementStrategies_.find(pName);

        if
        (
            (it != objectMaxRefLevel_.end()) &&
            (sIt->second & refType) )
        {
            activePatch[patchToIndex[pName]] = it->second;
        }
    }

    //- set active flag to facets belonging to active patches
    forAll(surf, triI)
    {
        const direction refLevel = activePatch[surf[triI].region()];
        activeFacets[triI] = Foam::max(activeFacets[triI], refLevel);
    }

    DynList<label> subsetIndices;
    surf.facetSubsetIndices(subsetIndices);
    forAll(subsetIndices, i)
    {
        const label subsetI = subsetIndices[i];
        const word sName = surf.facetSubsetName(subsetI);

        std::map<word, direction>::const_iterator it =
            objectMaxRefLevel_.find(sName);

        std::map<word, direction>::const_iterator sIt =
            objectRefinementStrategies_.find(sName);

        if
        (
            (it != objectMaxRefLevel_.end()) &&
            (sIt->second & refType)
        )
        {
            labelLongList trianglesInSubset;
            surf.facetsInSubset(subsetI, trianglesInSubset);

            forAll(trianglesInSubset, i)
                activeFacets[trianglesInSubset[i]] =
                    Foam::max(activeFacets[trianglesInSubset[i]], it->second);
        }
    }

    //- select leaves containing active patches
    DynList<label> containedTriangles;
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100) \
    private(containedTriangles)
    # endif
    forAll(intersectedLeaves, rcI)
    {
        const label leafI = intersectedLeaves[rcI];

        octree_.containedTriangles(leafI, containedTriangles);

        forAll(containedTriangles, i)
        {
            const direction refLevel = activeFacets[containedTriangles[i]];
            maxRefLevel[leafI] = Foam::max(maxRefLevel[leafI], refLevel);
        }

        if( !maxRefLevel[leafI] )
            maxRefLevel[leafI] = maxRefLevel_;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
