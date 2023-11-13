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

Class
    meshOctreeAutomaticRefinement

Description
    Automatic octree for mesh generation. This class checks if the current
    octree resolution is enough to resolve all surface features

SourceFiles
    meshOctreeAutomaticRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeAutomaticRefinement_H
#define meshOctreeAutomaticRefinement_H

#include "meshOctreeModifier.H"
#include "IOdictionary.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class triSurfaceCurvatureEstimator;
class triSurfacePartitioner;
class meshOctreeAddressing;

/*---------------------------------------------------------------------------*\
                Class meshOctreeAutomaticRefinement Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeAutomaticRefinement
{
    // Private data
        //- reference to meshOctree
        meshOctree& octree_;

        //- const reference to mesh dictionary
        const IOdictionary& meshDict_;

        //- refinement level correspondin to max cell size
        const direction globalRefLevel_;

        //- use DATA boxes
        bool useDATABoxes_;

        //- a flag for activating paired refinement
        bool hexRefinement_;

        //- a graph of triangles which are in the sticky space of a triangle
        mutable List<DynList<label> >* trianglesInStickySpacePtr_;

        //- a graph of edges within a sticky space of a edge
        mutable List<DynList<label> >* edgesInStickySpacePtr_;

        //- meshOctreeAddressing pointer
        mutable meshOctreeAddressing* octreeAddressingPtr_;

        //- surface partitioner pointer
        mutable triSurfacePartitioner* partitionerPtr_;

        //- curvature estimator pointer
        mutable triSurfaceCurvatureEstimator* curvaturePtr_;

    // Private member functions
        //- calculate triangles in sticky space
        void findTrianglesInStickySpace() const;
        const List<DynList<label> >& trianglesInStickySpace() const;

        //- calculate edges in sticky space
        void findEdgesInStickySpace() const;
        const List<DynList<label> >& edgesInStickySpace() const;

        //- calculate octree addresing
        void createOctreeAddressing() const;
        const meshOctreeAddressing& octreeAddressing() const;

        //- calculate partitioner
        void createSurfacePartitioner() const;
        const triSurfacePartitioner& partitioner() const;

        //- calculate curvature
        void createCurvatureEstimator() const;
        const triSurfaceCurvatureEstimator& curvature() const;

    // Functions and data performing automatic refinement
        //- maximum allowed refinement level
        direction maxRefLevel_;

        //- active refinement strategies
        direction globalRefStrategies_;

        //- global number of cells over gaps
        label nBoxesOverFeature_;

        //- active refinement strategies for each patch
        std::map<word, direction> objectMaxRefLevel_;
        std::map<word, direction> objectRefinementStrategies_;

        //- sticky distance. All disconnected triangles closer than
        //- the prescribed distances are considered to be connected
        scalarField stickyDistance_;

        //- proximity angle
        scalar proximityAngleTol_;

        //- shall the sticky space be used or not
        bool useStickySpace_;

        //- read parameters from meshDict
        void readParams();

        //- refine selected boxes
        void refineSelectedBoxes(labelList&, labelLongList&);

        //- find leaves intersecting local object
        void findActiveBoxes
        (
            const direction refinementType,
            labelLongList& intersectedLeaves,
            List<direction>& maxRefLevel
        ) const;

        //- find maximum refinement levels
        void findMaxRefLevels
        (
            const direction refinementType,
            const labelLongList& intersectedLeaves,
            List<direction>& maxRefLevel
        ) const;

        //- refine DATA boxes based on curvature
        bool curvatureRefinement();
        bool refineBasedOnCurvature
        (
            labelList&,
            const labelLongList&,
            const List<direction>&
        );

        //- refine octree based on local shortest edge length
        bool refineBasedOnEdgeLength
        (
            labelList&,
            const labelLongList&,
            const List<direction>&
        );

        //- refine DATA boxes which contain two or more independent parts
        //- of the surface in their vicinity
        bool proximityRefinement();
        bool refineBasedOnSurfaceProximityTests
        (
            labelList&,
            const labelLongList&,
            const List<direction>&
        );

        bool refineBasedOnEdgeProximityTests
        (
            labelList&,
            const labelLongList&,
            const List<direction>&
        );

        //- refine DATA boxes containing 2 or more corners, 2 or more distinct
        //- edges, etc.
        bool refineBasedOnContainedPartitions
        (
            labelList&,
            const labelLongList&,
            const List<direction>&
        );

        //- refine boxes based on the number of contained surface corners
        bool refineBasedOnContainedCorners
        (
            labelList&,
            const labelLongList&,
            const List<direction>&
        );

        //- refine boxes based on ray casting
        bool refineBasedOnRayCasting
        (
            labelList&,
            const labelLongList&,
            const List<direction>&
        );

        //- ensure that the boxes containing corner are never at lower
        //- refinement levels than any of its neighbouring leaves
        label ensureCornerAndEdgeConsistency
        (
            labelList&,
            const List<direction>&
        );

    // Private copy constructor
        //- Disallow default bitwise copy construct
        meshOctreeAutomaticRefinement(const meshOctreeAutomaticRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const meshOctreeAutomaticRefinement&);

    // Enumerators
        enum refinementStrategies_
        {
            NONE = 0,
            CURVATURE = 1,
            DISTINCTEDGEREGIONS = 2,
            EDGEPROXIMITY = 4,
            DISTINCTREGIONS = 8,
            SURFACEPROXIMITY = 16,
            CORNERPROXIMITY = 32,
            SHORTESTEDGELENGTH = 64,
            RAYCASTING = 128
        };

public:

    // Constructors

        //- Construct from meshOctree, dictionary and boolean
        meshOctreeAutomaticRefinement
        (
            meshOctree& mo,
            const IOdictionary& dict,
            const direction globalRefLevel,
            const bool useDATABoxes = false
        );

    // Destructor

        ~meshOctreeAutomaticRefinement();

    // Member Functions
        //- activate hex refinement
        void activateHexRefinement();

        //- refine octree boxes which are too coarse and cannot capture
        //- surface features
        void automaticRefinement();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
