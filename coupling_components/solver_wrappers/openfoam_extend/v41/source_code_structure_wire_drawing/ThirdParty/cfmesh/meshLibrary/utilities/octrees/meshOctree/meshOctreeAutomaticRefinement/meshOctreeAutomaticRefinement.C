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
#include "demandDrivenData.H"
#include "triSurfacePartitioner.H"
#include "triSurfaceCurvatureEstimator.H"
#include "meshOctreeAddressing.H"
#include "IOdictionary.H"
#include "triSurf.H"
#include "helperFunctions.H"

#include <set>

// #define DEBUGAutoRef

# ifdef DEBUGAutoRef
#include "OFstream.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void meshOctreeAutomaticRefinement::findTrianglesInStickySpace() const
{
    const triSurf& surf = octree_.surface();
    const pointField& points = surf.points();

    trianglesInStickySpacePtr_ = new List<DynList<label> >(surf.size());
    List<DynList<label> >& trianglesInStickySpace = *trianglesInStickySpacePtr_;

    if( !useStickySpace_ )
        return;

    DynList<label> containedTriangles;

    //- find connections
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 10) private(containedTriangles)
    # endif
    forAll(trianglesInStickySpace, triI)
    {
        scalar d = 2.0 * stickyDistance_[triI];

        boundBox bb;
        bb.min() = points[surf[triI][0]];
        bb.max() = bb.min();

        bb.min() = Foam::min(bb.min(), points[surf[triI][1]]);
        bb.max() = Foam::max(bb.max(), points[surf[triI][1]]);

        bb.min() = Foam::min(bb.min(), points[surf[triI][2]]);
        bb.max() = Foam::max(bb.max(), points[surf[triI][2]]);

        bb.min() -= vector(d, d, d);
        bb.max() += vector(d, d, d);

        octree_.findTrianglesInBox(bb, containedTriangles);

        const triangle<point, point> tri = help::surfaceTriangle(surf, triI);

        point ntr(vector::zero), notr(vector::zero);

        forAll(containedTriangles, i)
        {
            const label ntI = containedTriangles[i];

            if( ntI == triI )
                continue;

            const triangle<point, point> otri =
                help::surfaceTriangle(surf, ntI);

            help::nearestPointsOnTriangles(tri, otri, ntr, notr);

            if( magSqr(ntr - notr) < d * d )
            {
                trianglesInStickySpace[triI].append(ntI);
            }
        }
    }
}

const List<DynList<label> >&
meshOctreeAutomaticRefinement::trianglesInStickySpace() const
{
    if( !trianglesInStickySpacePtr_ )
    {
        findTrianglesInStickySpace();
    }

    return *trianglesInStickySpacePtr_;
}

void meshOctreeAutomaticRefinement::findEdgesInStickySpace() const
{
    const triSurf& surf = octree_.surface();
    const edgeLongList& edges = surf.edges();
    const VRWGraph& edgeFaces = surf.edgeFacets();

    const pointField& pts = surf.points();

    edgesInStickySpacePtr_ = new List<DynList<label> >(edges.size());
    List<DynList<label> >& edgesInStickySpace = *edgesInStickySpacePtr_;

    if( !useStickySpace_ )
        return;

    DynList<label> containedEdges;

    //- find connections
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 10) private(containedEdges)
    # endif
    forAll(edgesInStickySpace, edgeI)
    {
        scalar d(0.0);
        forAllRow(edgeFaces, edgeI, i)
            d = Foam::max(d, 2.0 * stickyDistance_[edgeFaces(edgeI, i)]);

        const edge& e = edges[edgeI];
        boundBox bb;
        bb.min() = pts[e.start()];
        bb.max() = bb.min();

        bb.min() = Foam::min(bb.min(), pts[e.end()]);
        bb.max() = Foam::max(bb.max(), pts[e.end()]);

        bb.min() -= vector(d, d, d);
        bb.max() += vector(d, d, d);

        octree_.findEdgesInBox(bb, containedEdges);

        point ne(vector::zero), noe(vector::zero);

        forAll(containedEdges, i)
        {
            const label neI = containedEdges[i];

            if( neI == edgeI )
                continue;

            const edge& oe = edges[neI];

            help::nearestPointsOnEdges
            (
                pts[e.start()], pts[e.end()],
                pts[oe.start()], pts[oe.end()],
                ne, noe
            );


            if( magSqr(ne - noe) < d * d )
            {
                edgesInStickySpace[edgeI].append(neI);
            }
        }
    }
}

const List<DynList<label> >&
meshOctreeAutomaticRefinement::edgesInStickySpace() const
{
    if( !edgesInStickySpacePtr_ )
    {
        findEdgesInStickySpace();
    }

    return *edgesInStickySpacePtr_;
}

void meshOctreeAutomaticRefinement::createOctreeAddressing() const
{
    octreeAddressingPtr_ =
        new meshOctreeAddressing(octree_, meshDict_, useDATABoxes_);
}

const meshOctreeAddressing& meshOctreeAutomaticRefinement::octreeAddressing()
const
{
    if( !octreeAddressingPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const meshOctreeAddressing& meshOctreeAutomaticRefinement"
                "::octreeAddressing() const"
            ) << "Cannot calculate addressing!" << abort(FatalError);
        # endif

        createOctreeAddressing();
    }

    return *octreeAddressingPtr_;
}

void meshOctreeAutomaticRefinement::createSurfacePartitioner() const
{
    partitionerPtr_ = new triSurfacePartitioner(octree_.surface());
}

const triSurfacePartitioner& meshOctreeAutomaticRefinement::partitioner() const
{
    if( !partitionerPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const triSurfacePartitioner& meshOctreeAutomaticRefinement"
                "::partitioner() const"
            ) << "Cannot calculate addressing!" << abort(FatalError);
        # endif

        createSurfacePartitioner();
    }

    return *partitionerPtr_;
}

void meshOctreeAutomaticRefinement::createCurvatureEstimator() const
{
    curvaturePtr_ = new triSurfaceCurvatureEstimator(octree_.surface());
}

const triSurfaceCurvatureEstimator&
meshOctreeAutomaticRefinement::curvature() const
{
    if( !curvaturePtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const triSurfaceCurvatureEstimator& "
                "meshOctreeAutomaticRefinement::curvature() const"
            ) << "Cannot calculate addressing!" << abort(FatalError);
        # endif

        createCurvatureEstimator();
    }

    return *curvaturePtr_;
}

void meshOctreeAutomaticRefinement::readParams()
{
    const triSurf& surf = octree_.surface();
    const geometricSurfacePatchList& patches = surf.patches();

    const boundBox& rootBox = octree_.rootBox();
    const scalar size = rootBox.max().x() - rootBox.min().x();
    const scalar maxSize(readScalar(meshDict_.lookup("maxCellSize")));

    maxRefLevel_ = globalRefLevel_;
    globalRefStrategies_ = NONE;
    nBoxesOverFeature_ = 1;
    stickyDistance_.setSize(surf.size());
    stickyDistance_ = -1.0;
    proximityAngleTol_ = 0.15;

    if( meshDict_.found("minCellSize") )
    {
        scalar cs(readScalar(meshDict_.lookup("minCellSize")));
        cs *= (1.0 + SMALL);

        if( cs > maxSize )
            return;

        bool finished;
        do
        {
            finished = false;

            const scalar lSize = size / pow(2, label(maxRefLevel_));

            if( lSize < cs )
            {
                finished = true;
            }
            else
            {
                ++maxRefLevel_;
            }
        } while( !finished );

        //useDATABoxes_ = true;

        Info << "Requested min cell size corresponds to octree level "
            << label(maxRefLevel_) << endl;

        globalRefStrategies_ |= CURVATURE;
        globalRefStrategies_ |= EDGEPROXIMITY;
        globalRefStrategies_ |= DISTINCTEDGEREGIONS;
        globalRefStrategies_ |= DISTINCTREGIONS;
        globalRefStrategies_ |= SURFACEPROXIMITY;
        globalRefStrategies_ |= CORNERPROXIMITY;
        globalRefStrategies_ |= RAYCASTING;
    }
    else if( meshDict_.found("maxAdditionalRefinementLevels") )
    {
        maxRefLevel_ +=
            readLabel(meshDict_.lookup("maxAdditionalRefinementLevels"));
        Info << "Requested min cell size corresponds to octree level "
            << label(maxRefLevel_) << endl;

        globalRefStrategies_ |= CURVATURE;
        globalRefStrategies_ |= EDGEPROXIMITY;
        globalRefStrategies_ |= DISTINCTEDGEREGIONS;
        globalRefStrategies_ |= DISTINCTREGIONS;
        globalRefStrategies_ |= SURFACEPROXIMITY;
        globalRefStrategies_ |= CORNERPROXIMITY;
        globalRefStrategies_ |= RAYCASTING;
    }

    if( meshDict_.found("automaticRefinement") )
    {
        const dictionary& autoRefDict =
            meshDict_.subDict("automaticRefinement");

        //- read global refinement level
        if( autoRefDict.found("minCellSize") )
        {
            scalar cs = readScalar(autoRefDict.lookup("minCellSize"));
            cs *= (1.0 + SMALL);

            if( cs > maxSize )
                return;

            bool finished;
            do
            {
                finished = false;

                const scalar lSize = size / pow(2, label(maxRefLevel_));

                if( lSize < cs )
                {
                    finished = true;
                }
                else
                {
                    ++maxRefLevel_;
                }
            } while( !finished );

            Info << "Requested min cell size corresponds to octree level "
                << label(maxRefLevel_) << endl;

            stickyDistance_ = cs / 4.0;
        }
        else if( autoRefDict.found("maxAdditionalRefinementLevels") )
        {
            maxRefLevel_ +=
                readLabel(autoRefDict.lookup("maxAdditionalRefinementLevels"));

            Info << "Requested min cell size corresponds to octree level "
                << label(maxRefLevel_) << endl;

            const scalar cs = maxSize / pow(2, label(maxRefLevel_));
            stickyDistance_ = cs / 4.0;
        }

        if( autoRefDict.found("stickyDistance") )
        {
            stickyDistance_ =
                readScalar(autoRefDict.lookup("stickyDistance"));
            useStickySpace_ = true;
        }

        if( autoRefDict.found("proximityAngle") )
        {
            const scalar angle =
                readScalar(autoRefDict.lookup("proximityAngle"));

            proximityAngleTol_ = Foam::cos(M_PI * angle / 180.0);
        }

        //- check the requested number of cells over feature
        if( autoRefDict.found("nCellsOverFeatureSize") )
        {
            nBoxesOverFeature_ =
                readLabel(autoRefDict.lookup("nCellsOverFeatureSize"));
        }

        //- check if global curvature shall be used
        if( autoRefDict.found("curvatureRefinement") )
        {
            const bool useCurvature
            (
                readBool(autoRefDict.lookup("curvatureRefinement"))
            );

            if( useCurvature )
                globalRefStrategies_ |= CURVATURE;
        }

        //- check if surface edge length shall be used globally
        if( autoRefDict.found("edgeLengthRefinement") )
        {
            const bool useLength
            (
                readBool(autoRefDict.lookup("edgeLengthRefinement"))
            );

            if( useLength )
                globalRefStrategies_ |= SHORTESTEDGELENGTH;
        }

        //- refine region with distinct parts
        if( autoRefDict.found("distinctPartsRefinement") )
        {
            const bool distinctParts
            (
                readBool(autoRefDict.lookup("distinctPartsRefinement"))
            );

            if( distinctParts )
            {
                globalRefStrategies_ |= DISTINCTEDGEREGIONS;
                globalRefStrategies_ |= DISTINCTREGIONS;
            }
        }

        //- refine based on proximity tests
        if( autoRefDict.found("proximityRefinement") )
        {
            const bool proximityRef
            (
                readBool(autoRefDict.lookup("proximityRefinement"))
            );

            if( proximityRef )
            {
                globalRefStrategies_ |= EDGEPROXIMITY;
                globalRefStrategies_ |= SURFACEPROXIMITY;
                globalRefStrategies_ |= CORNERPROXIMITY;
            }
        }

        //- refine based on ray-casting tests
        if( autoRefDict.found("rayCastingRefinement") )
        {
            const bool proximityRef
            (
                readBool(autoRefDict.lookup("rayCastingRefinement"))
            );

            if( proximityRef )
            {
                globalRefStrategies_ |= RAYCASTING;
            }
        }

        //- local refinement strategies
        if( autoRefDict.found("localRefinement") )
        {
            const dictionary& localDict =
                autoRefDict.subDict("localRefinement");

            const wordList objectNames = localDict.toc();

            forAll(objectNames, objI)
            {
                direction& localRefStrategies =
                    objectRefinementStrategies_[objectNames[objI]];
                localRefStrategies = NONE;

                const dictionary& dict = localDict.subDict(objectNames[objI]);

                //- read global refinement level
                if( dict.found("minCellSize") )
                {
                    scalar cs = readScalar(dict.lookup("minCellSize"));
                    cs *= (1.0 + SMALL);

                    if( cs > maxSize )
                        return;

                    direction level(0);
                    bool finished;
                    do
                    {
                        finished = false;

                        const scalar lSize = size / pow(2, label(level));

                        if( lSize < cs )
                        {
                            finished = true;
                        }
                        else
                        {
                            ++level;
                        }
                    } while( !finished );

                    objectMaxRefLevel_[objectNames[objI]] = level;
                }
                else if( dict.found("maxAdditionalRefinementLevels") )
                {
                    objectMaxRefLevel_[objectNames[objI]] =
                        globalRefLevel_ +
                        readLabel(dict.lookup("maxAdditionalRefinementLevels"));
                }
                else
                {
                    objectMaxRefLevel_[objectNames[objI]] = maxRefLevel_;
                }

                //- check if global curvature shall be used
                if( dict.found("curvatureRefinement") )
                {
                    const bool useCurvature
                    (
                        readBool(dict.lookup("curvatureRefinement"))
                    );

                    if( useCurvature )
                        localRefStrategies |= CURVATURE;
                }
                else
                {
                    localRefStrategies |= (globalRefStrategies_ & CURVATURE);
                }

                //- check if surface edge length shall be used
                if( dict.found("edgeLengthRefinement") )
                {
                    const bool useLength
                    (
                        readBool(dict.lookup("edgeLengthRefinement"))
                    );

                    if( useLength )
                        localRefStrategies |= SHORTESTEDGELENGTH;
                }
                else
                {
                    localRefStrategies |=
                        (globalRefStrategies_ & SHORTESTEDGELENGTH);
                }

                //- refine regions with distinct parts
                if( dict.found("distinctPartsRefinement") )
                {
                    const bool distinctParts
                    (
                        readBool(dict.lookup("distinctPartsRefinement"))
                    );

                    if( distinctParts )
                    {
                        localRefStrategies |= DISTINCTEDGEREGIONS;
                        localRefStrategies |= DISTINCTREGIONS;
                    }
                }
                else
                {
                    localRefStrategies |=
                        (globalRefStrategies_ & DISTINCTEDGEREGIONS);
                    localRefStrategies |=
                        (globalRefStrategies_ & DISTINCTREGIONS);
                }

                //- refine based on proximity tests
                if( dict.found("proximityRefinement") )
                {
                    const bool proximityRef
                    (
                        readBool(dict.lookup("proximityRefinement"))
                    );

                    if( proximityRef )
                    {
                        localRefStrategies |= EDGEPROXIMITY;
                        localRefStrategies |= SURFACEPROXIMITY;
                        localRefStrategies |= CORNERPROXIMITY;
                    }
                }
                else
                {
                    localRefStrategies |=
                        (globalRefStrategies_ & EDGEPROXIMITY);
                    localRefStrategies |=
                        (globalRefStrategies_ & SURFACEPROXIMITY);
                    localRefStrategies |=
                        (globalRefStrategies_ & CORNERPROXIMITY);
                }

                //- refine based on ray-casting tests
                if( dict.found("rayCastingRefinement") )
                {
                    const bool proximityRef
                    (
                        readBool(dict.lookup("rayCastingRefinement"))
                    );

                    if( proximityRef )
                    {
                        localRefStrategies |= RAYCASTING;
                    }
                }
                else
                {
                    localRefStrategies |=
                        (globalRefStrategies_ & RAYCASTING);
                }

                if( dict.found("stickyDistance") )
                {
                    const scalar dist =
                        readScalar(dict.lookup("stickyDistance"));
                    useStickySpace_ = true;

                    forAll(patches, patchI)
                    {
                        if( patches[patchI].name() == objectNames[objI] )
                        {
                            forAll(surf, triI)
                            {
                                if( surf[triI].region() == patchI )
                                {
                                    stickyDistance_[triI] = dist;
                                }
                            }
                        }
                    }

                    const label sId = surf.facetSubsetIndex(objectNames[objI]);
                    if( sId >= 0 )
                    {
                        labelLongList containedTriangles;
                        surf.facetsInSubset(sId, containedTriangles);

                        forAll(containedTriangles, i)
                            stickyDistance_[containedTriangles[i]] = dist;
                    }
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from surface and IOdictionary
meshOctreeAutomaticRefinement::meshOctreeAutomaticRefinement
(
    meshOctree& mo,
    const IOdictionary& dict,
    const direction globalRefLevel,
    bool useDATABoxes
)
:
    octree_(mo),
    meshDict_(dict),
    globalRefLevel_(globalRefLevel),
    useDATABoxes_(useDATABoxes),
    hexRefinement_(false),
    trianglesInStickySpacePtr_(NULL),
    edgesInStickySpacePtr_(NULL),
    octreeAddressingPtr_(NULL),
    partitionerPtr_(NULL),
    curvaturePtr_(NULL),
    maxRefLevel_(0),
    globalRefStrategies_(NONE),
    nBoxesOverFeature_(0),
    objectMaxRefLevel_(),
    objectRefinementStrategies_(),
    stickyDistance_(),
    proximityAngleTol_(0.15),
    useStickySpace_(false)
{
    if( !useDATABoxes_ && dict.found("keepCellsIntersectingBoundary") )
    {
        useDATABoxes_ = readBool(dict.lookup("keepCellsIntersectingBoundary"));
    }

    //- read parameters from meshDict and calculate ref levels
    readParams();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshOctreeAutomaticRefinement::~meshOctreeAutomaticRefinement()
{
    deleteDemandDrivenData(trianglesInStickySpacePtr_);
    deleteDemandDrivenData(edgesInStickySpacePtr_);
    deleteDemandDrivenData(octreeAddressingPtr_);
    deleteDemandDrivenData(curvaturePtr_);
    deleteDemandDrivenData(partitionerPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
