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
    meshSurfaceOptimizer

Description
    Surface smoothing without any topological changes

SourceFiles
    meshSurfaceOptimizer.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceOptimizer_H
#define meshSurfaceOptimizer_H

#include "polyMeshGen.H"
#include "boolList.H"
#include "labelLongList.H"
#include "DynList.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfacePartitioner.H"
#include "partTriMesh.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceMapper;
class plane;

/*---------------------------------------------------------------------------*\
                    Class meshSurfaceOptimizer Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceOptimizer
{
    // Private data
        //- const reference to the mesh surface
        const meshSurfaceEngine& surfaceEngine_;

        //- surface modifier
        meshSurfaceEngineModifier surfaceModifier_;

        //- type of surface vertex
        List<direction> vertexType_;

        //- locked faces which shall not be changed
        labelLongList lockedSurfaceFaces_;

        //- surface partitioner
        const meshSurfacePartitioner* partitionerPtr_;
        const bool deletePartitioner_;

        //- pointer to mesh octree
        const meshOctree* octreePtr_;

        //- mesh of surface triangles needed for some smoothers
        mutable partTriMesh* triMeshPtr_;

        //- enforce constraints
        bool enforceConstraints_;

        //- name of the subset contaning tangled points
        word badPointsSubsetName_;

    // Private member functions
        //- classify surface vertices as PARTITION, EDGE, CORNER
        void classifySurfaceVertices();

        //- calculate surface triangulation
        void calculateTrianglesAndAddressing() const;
        inline const partTriMesh& triMesh() const;
        inline void updateTriMesh(const labelLongList&);
        inline void updateTriMesh();

        //- find vertices which are inverted
        //- mark additional layers of vertices near the inverted ones
        label findInvertedVertices
        (
            boolList& smoothVertex,
            const label nAdditionalLayers = 2
        ) const;

        //- transform into a 2D space in plane
        inline bool transformIntoPlane
        (
            const label bpI,
            const plane& pl,
            vector& vecX,
            vector& vecY,
            scalar& minCosSq,
            DynList<point>& pts,
            DynList<triFace>& trias
        ) const;

        //- new position of a node after trying to make each face
        //- as flat as possible
        inline void localDataFaceFlatness
        (
            const label bpI,
            vector& disp,
            scalar& sumW,
            const bool weighted = true
        ) const;

        inline point newPositionFaceFlatness
        (
            const label bpI,
            const vector& disp,
            const scalar& sumW
        ) const;

        inline point newPositionFaceFlatness
        (
            const label bpI,
            const bool weighted = true
        ) const;

        //- new position of a point after laplacian smoothing
        //- vertices are provided
        inline point newPositionLaplacian
        (
            const label bpI,
            const DynList<point>& neiPoints,
            const bool transformIntoPlane = true,
            const bool weighted = true
        ) const;

        //- new position of a node after laplacian smoothing
        //- the position is the average of neighbouring vertex positions
        inline point newPositionLaplacian
        (
            const label bpI,
            const bool transformIntoPlane = true,
            const bool weighted = true
        ) const;

        //- new position of a node after laplacian smoothing
        //- the position is the average of neighbouring face centres
        inline point newPositionLaplacianFC
        (
            const label bpI,
            const bool transformIntoPlane = true,
            const bool weighted = true
        ) const;

        //- new position of a node after using surfaceOptimizer
        inline point newPositionSurfaceOptimizer
        (
            const label bpI,
            const scalar tol = 0.001
        ) const;

        //- new position of edge vertex
        //- the position is the average of neighbouring edge vertices
        static inline point newEdgePositionLaplacian
        (
            const point& p0,
            const point& p,
            const point& p1
        );

        inline point newEdgePositionLaplacian(const label bpI) const;

        //- smooth selected edge points
        void smoothEdgePoints
        (
            const labelLongList& edgePoints,
            const labelLongList& procEdgePoints,
            const scalar relaxationFactor = 0.8
        );

        //- smooth selected points to improve face flatness
        void smoothFaceFlatness
        (
            const labelLongList& selectedPoints,
            const labelLongList& selectedProcPoints,
            const scalar relaxationFactor = 0.8,
            const bool weighted = true
        );

        //- smooth selected points using laplace smoother and face centres
        void smoothLaplacianFC
        (
            const labelLongList& selectedPoints,
            const labelLongList& selectedProcPoints,
            const scalar relaxationFactor = 0.8,
            const bool transform = true,
            const bool weighted = true
        );

        //- smooth selected points using surface optimizer
        void smoothSurfaceOptimizer
        (
            const labelLongList& selectedPoints,
            const labelLongList& selectedProcPoints,
            const scalar relaxationFactor = 0.8
        );

    // Functions needed for parallel runs
        //- laplacian smoothing of points at processor boundaries
        //- new position is the average of the neighbouring vertices
        void nodeDisplacementFaceFlatnessParallel
        (
            const labelLongList& nodesToSmooth,
            pointField& newPoints,
            const bool weighted = true
        );

        //- laplacian smoothing of points at processor boundaries
        //- new position is the average of the neighbouring vertices
        void nodeDisplacementLaplacianParallel
        (
            const labelLongList& nodesToSmooth,
            pointField& newPoints,
            const bool transformIntoPlane = true,
            const bool weighted = true
        );

        //- laplacian smoothing of points at processor boundaries
        void nodeDisplacementLaplacianFCParallel
        (
            const labelLongList& nodesToSmooth,
            pointField& newPoints,
            const bool transformIntoPlane = true,
            const bool weighted = true
        );

        //- smooth edge nodes at processor boundaries
        void edgeNodeDisplacementParallel
        (
            const labelLongList& nodesToSmooth,
            pointField& newPoints
        );

        //- re-project vertices back onto the surface of the mesh
        void projectPointsBackOnTheSurface
        (
            const labelLongList& nodesToSmooth,
            pointField& newCoordinates
        ) const;

        //- Disallow default bitwise copy construct
        meshSurfaceOptimizer(const meshSurfaceOptimizer&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceOptimizer&);

    // Enumerators
        enum vertexTypes
        {
            PARTITION = 1,
            EDGE = 2,
            CORNER = 4,
            PROCBND = 8,
            LOCKED = 16
        };

public:

    // Constructors

        //- Construct from mesh surface
        meshSurfaceOptimizer(const meshSurfaceEngine&);

        //- Construct from partitioner
        meshSurfaceOptimizer(const meshSurfacePartitioner&);

        //- Construct from mesh surface and octree
        meshSurfaceOptimizer(const meshSurfaceEngine&, const meshOctree&);

        //- Construct from partitioner and octree
        meshSurfaceOptimizer
        (
            const meshSurfacePartitioner&,
            const meshOctree&
        );


    // Destructor

        ~meshSurfaceOptimizer();

    // Member Functions
        //- lock the boundary faces which shall not be modified
        template<class labelListType>
        inline void lockBoundaryFaces(const labelListType&);

        //- lock boundary points. They are not be moved.
        template<class labelListType>
        inline void lockBoundaryPoints(const labelListType&);

        //- lock edge points
        inline void lockFeatureEdges();

        //- reset to default constraints
        void removeUserConstraints();

        //- Flag stopping the meshing process if it is not possible
        //- to untangle the surface without sacrificing geometry constraints
        //- Points which must be moved away from the required position are
        //- stored into a point subset
        void enforceConstraints(const word subsetName="badPoints");

        //- agressive smoothing at the surface aimed at smoothing severly
        //- tangled parts
        bool untangleSurfaceNonConstrainedLaplace
        (
            const labelLongList& selectedBoundaryPoints,
            const label nAdditionalLayers = 1,
            const label nGlobalIterations = 20
        );

        bool untangleSurfaceNonConstrainedLaplace
        (
            const label nAdditionalLayers = 1,
            const label nGlobalIterations = 20
        );

        //- runs a surface smoother on the selected boundary points
        bool untangleSurface
        (
            const labelLongList& activeBoundaryPoints,
            const label nAdditionalLayers = 0,
            const label nGlobalIterations = 10,
            const label nLocalIterations = 10,
            const label nLaplaceIterations = 1
        );

        //- checks for inverted surface elements and tries to untangle them
        //- it tries to keep the points on the surface for a couple
        //- of iteration and gives up the final iterations
        //- by default, it smooths two additional layers of elements
        //- around the inverted ones
        bool untangleSurface
        (
            const label nAdditionalLayers = 0,
            const label nGlobalIterations = 10,
            const label nLocalIterations = 10,
            const label nLaplaceIterations = 1
        );

        //- checks for remaining inverted surface elements and tries to untangle
        //- them by reverting the surface vertices to their initial position
        //- when the state of the mesh was valid
        void untangleAndRevertSurface
        (
            const labelLongList& activeBoundaryPoints,
            const label nAdditionalLayers = 0,
            const label nGlobalIterations = 10,
            const label nLocalIterations = 10,
            const label nLaplaceIterations = 1
        );

        void untangleAndRevertSurface
        (
            const label nAdditionalLayers = 0,
            const label nGlobalIterations = 10,
            const label nLocalIterations = 10,
            const label nLaplaceIterations = 1
        );

        //- optimize boundary nodes using the surface laplace
        void optimizeSurfaceLaplace
        (
            const label nIterations = 5,
            const bool allowShrinking = false
        );

        //- optimize boundary nodes after boundary regions are created
        void optimizeSurface(const label nIterations = 5);

        //- optimize the surface of a 2D mesh
        void optimizeSurface2D(const label nIterations = 5);

        //- optimize low quality parts of a 2D mesh
        void optimizeLowQualitySurface2D(const label nIterations = 5);

        //- untangle the surface of a 2D mesh
        void untangleSurface2D(const label maxNumIterations = 10);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshSurfaceOptimizerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
