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
    meshSurfaceMapper

Description
    Maps vertices to the nearest point on the geometry surface

SourceFiles
    meshSurfaceMapper.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceMapper_H
#define meshSurfaceMapper_H

#include "labelList.H"
#include "pointField.H"
#include "labelLongList.H"
#include "DynList.H"
#include "parMapperHelper.H"
#include "labelledPointScalar.H"
#include "boolList.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;
class meshSurfacePartitioner;
class triSurfacePartitioner;
class meshSurfaceEngineModifier;

/*---------------------------------------------------------------------------*\
                    Class meshSurfaceMapper Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceMapper
{
    // Private data
        //- mesh surface
        const meshSurfaceEngine& surfaceEngine_;

        //- reference to the octree
        const meshOctree& meshOctree_;

        //- mesh surface partitioner
        mutable const meshSurfacePartitioner* surfaceEnginePartitionerPtr_;
        const bool deletePartitioner_;

        //- triSurface partitioner
        mutable triSurfacePartitioner* surfPartitionerPtr_;

        //- mesh surface modifier
        mutable meshSurfaceEngineModifier* surfaceModifierPtr_;
        const bool deleteModifier_;

        //- classify points by types
        List<direction> pointType_;

        //- classify faces by type
        List<direction> faceType_;

    // Private member functions
        //- create and return mesh surface partitioner
        void createMeshSurfacePartitioner() const;
        inline const meshSurfacePartitioner& meshPartitioner() const
        {
            if( !surfaceEnginePartitionerPtr_ )
                createMeshSurfacePartitioner();

            return *surfaceEnginePartitionerPtr_;
        }

        //- create and return triSurfacePartitioner
        void createTriSurfacePartitioner() const;
        inline const triSurfacePartitioner& surfacePartitioner() const
        {
            if( !surfPartitionerPtr_ )
                createTriSurfacePartitioner();

            return *surfPartitionerPtr_;
        }

        //- create and return meshSurfaceEngineModifier
        meshSurfaceEngineModifier& meshSurfaceModifier();

        //- delete surfaceEnginePartitionerPtr_ and surfPartitionerPtr_
        void clearOut();

    // Private member functions
        //- surface laplace
        void calculateMovementShrinkingSurfaceLaplaceFC
        (
            const List<DynList<scalar> >& faceCentreWeights,
            vectorField& displacements
        ) const;

        //- find the nearest point at the surface mesh for every face centre
        void findNearestPointsToFaceCentres
        (
            const labelLongList& activeFaces,
            List<labelledPointScalar>& nearestPoints
        ) const;

        void findNearestPointsToFaceCentres
        (
            List<labelledPointScalar>& nearestPoints
        );

        //- find nearest points at the surface mesh for every
        void findNearestPoints
        (
            const labelLongList& activePoints,
            List<labelledPointScalar>& nearestPoints
        ) const;

        void findNearestPoints(List<labelledPointScalar>& nearestPoints) const;

        //- find the best patch for face
        void findFacePatches
        (
            const labelLongList& activePoints,
            labelLongList& facePatch
        ) const;
        void findFacePatches(labelLongList& facePatch) const;

        //- find patches at point
        void findPointPatches
        (
            const labelLongList& activePoints,
            List<DynList<label> >& pointPatches
        ) const;
        void findPointPatches(List<DynList<label> >& pointPatches) const;

        //- find vertices at the surface mesh
        void findMappingVertices
        (
            const labelLongList& activePoints,
            pointField& newPoints,
            boolList& foundNearest
        ) const;

        void findMappingVertices
        (
            pointField& newPoints,
            boolList& foundNearest
        ) const;

        //- find vertices at the surface mesh with respect to desired patches
        void findMappingVertices
        (
            const labelLongList& activePoints,
            const List<DynList<label> >& pointPatches,
            pointField& newPoints,
            boolList& findNearest
        ) const;

        void findMappingVertices
        (
            const List<DynList<label> >& pointPatches,
            pointField& newPoints,
            boolList& findNearest
        ) const;

        //- classify point and faces at the surface of the mesh
        void analysePointsAndFaces();

        //- map corner nodes to the boundary
        void mapCorners(const labelLongList& nodesToMap);

        //- find mapping distance for selected points
        void findMappingDistanceSquared
        (
            const labelLongList& nodesToMap,
            scalarList& mappingDistanceSq
        ) const;

        //- limit displacement
        void limitDisplacement
        (
            const label pointI,
            const scalar maxDistSq,
            point& newP,
            scalar& dSq
        ) const;

        //- calculate face metric
        scalar faceMetricInPatch(const label bfI, const label patchI) const;

    // Private member functions needed for parallel execution

        //- check if nodes at parallel boundaries are selected at all processors
        void selectNodesAtParallelBnd(const labelLongList&);

        //- map to the smallest distance
        //- makes sense for parallel calculations, only
        void mapToSmallestDistance(LongList<parMapperHelper>&);

        //- Disallow default bitwise copy construct
        meshSurfaceMapper(const meshSurfaceMapper&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceMapper&);

    // Enumerators

        //- point types
        enum pointTypes
        {
            NONE = 0,
            EDGE = 1,
            CORNER = 2,
            SPECIAL = 4,
            LOCKED = 8
        };

        //- face types
        enum faceTypes
        {
            OUTERATTRACTION = 1,
            INNERATTRACTION = 2
        };

public:

    // Constructors

        //- Construct from meshSurfaceEngine, octree
        //- and meshSurfaceEngineModifier
        meshSurfaceMapper
        (
            const meshSurfaceEngine&,
            const meshOctree&,
            meshSurfaceEngineModifier* surfaceModifierPtr = NULL
        );

        //- Construct from meshSurfacePartitioner, octree and surface modifier
        meshSurfaceMapper
        (
            const meshSurfacePartitioner&,
            const meshOctree&,
            meshSurfaceEngineModifier* surfaceModifierPtr = NULL
        );

    // Destructor

        ~meshSurfaceMapper();

    // Member Functions
        //- map the given surface node on the selected patch
        //- this does not make sense for vertices at parallel boundaries
        void mapNodeToPatch(const label bpI, const label patchI = -1);

        //- projects surface vertices onto their nearest location
        //- on the surface mesh
        void mapVerticesOntoSurface();

        //- projects selected surface vertices to their nearest location
        //- on the surface mesh
        void mapVerticesOntoSurface(const labelLongList& nodesToMap);

        //- projects corner and edge vertices onto their nearest location
        //- on the surface mesh
        void mapCornersAndEdges();

        //- projects selected edge vertices onto their nearest
        //- locations on the surface mesh
        void mapEdgeNodes(const labelLongList& nodesToMap);

        //- projects surface vertices onto the surface with respect
        //- to the surface patch they belong to. Edges and corner are respected
        void mapVerticesOntoSurfacePatches();

        //- projects selected surface vertices onto the surface with respect
        //- to the surface patch they belong to. Edges and corner are respected
        void mapVerticesOntoSurfacePatches(const labelLongList& nodesToMap);

        //- a combination of mapping and smoothing intended for better
        //- feature capturing
        void preMapVertices(const label nIterations = 3);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
