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
    partTetMesh

Description
    Mesh smoothing without any topological changes

SourceFiles
    partTetMesh.C

\*---------------------------------------------------------------------------*/

#ifndef partTetMesh_H
#define partTetMesh_H

#include "boolList.H"
#include "polyMeshGen.H"
#include "labelLongList.H"
#include "vectorLongList.H"
#include "VRWGraph.H"
#include "DynList.H"
#include "partTet.H"
#include "HashSet.H"
#include "labelledPoint.H"
#include "labelledTri.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMeshGen;
class VRWGraph;
class partTetMeshSimplex;

/*---------------------------------------------------------------------------*\
                        Class partTetMesh Declaration
\*---------------------------------------------------------------------------*/

class partTetMesh
{
    // Private data
        //- reference to the original mesh
        polyMeshGen& origMesh_;

        //- const reference to mesh points
        const pointFieldPMG& points_;

        //- const reference to mesh faces
        const faceListPMG& faces_;

        //- const reference to mesh cells
        const cellListPMG& cells_;

        //- const reference to owner cells
        const labelLongList& owner_;

        //- const reference to neighbour cells
        const labelLongList& neighbour_;

        //- const reference to face centres
        const vectorLongList& faceCentres_;

        //- const reference to cell centres
        const vectorLongList& cellCentres_;

        //- const reference to point-faces addressing
        const VRWGraph& pointFaces_;

        //- shall a node be used for smoothing or not
        LongList<direction> smoothVertex_;

        //- internal vertices that may be moved
        labelLongList movableInternal_;

        //- boundary vertices that may be moved
        labelLongList movableBoundary_;

        //- bnd face to patch
        labelLongList bndFaceToPatch_;

    // Private data for parallel runs
        //- vertices supporting tetrahedra attached to points
        //- at inter-processor boundaries
        LongList<point> parPoints_;

        //- tets at inter-processor vertices
        std::map<label, DynList<partTet, 256> > tetsAtPoint_;

        //- boundary triangles at inter-processor vertices
        std::map<label, DynList<labelledTri, 32> > bndTrianglesAtPoint_;

        //- orig points to parPoints mapping for points at this processor
        std::map<label, label> pointToParPoint_;

        //- face centre to parPoint mapping for face centres at this processor
        std::map<label, label> faceCentreToParPoint_;

        //- cell centre to parPoint mapping for cell centres at this processor
        std::map<label, label> cellCentreToParPoint_;

        //- global label of a parPoint
        labelLongList globalParPointLabel_;
        std::map<label, label> globalToLocalParPointLabel_;
        std::map<label, std::set<label> > parPointAtOtherProcs_;

    // Private member functions

        //- create points and tets
        void createPointsAndTets
        (
            const List<direction>& usedCells,
            const boolList& lockedPoints
        );

        //- create tetrahedral cells at a given face
        void createTetsAtFace
        (
            const label faceI,
            const label pointI,
            std::map<label, label>& ppToPts,
            std::map<label, label>& fcToPts,
            std::map<label, label>& ccToPts,
            DynList<point, 128>& pts,
            DynList<partTet, 256>& tets,
            DynList<labelledTri, 32>& bndTriangles
        ) const;

        //- create simplex for a given point
        void createSimplex(const label, partTetMeshSimplex&) const;

        //- create parallel addressing
        void createParallelAddressing();

        //- create buffer layers
        void createBufferLayers();

        //- ensure that vertices at inter-processor boundary are at the same
        //- position on all processors
        void unifyCoordinatesAtInterProcessorBoundaries();

        //- update coordinates of buffer layer points at processors
        //- that do not own the vertices
        void updateBufferLayerPoints();

        //- Disallow bitwise copy construct
        partTetMesh(const partTetMesh&);

        //- Disallow bitwise asignment
        void operator=(const partTetMesh&);

public:

    // Friend classes

        friend class partTetMeshSimplex;

    // Constructors
        //- construct from polyMeshGen and locked points
        partTetMesh(polyMeshGen& mesh, const labelLongList& lockedPoints);

        //- construct from polyMeshGen, locked points and the number of layers
        //- from the boundary
        partTetMesh
        (
            polyMeshGen& mesh,
            const labelLongList& lockedPoints,
            const direction nLayers
        );

        //- construct from polyMeshGen, bad faces and the number
        //- of additional layers
        partTetMesh
        (
            polyMeshGen& mesh,
            const labelLongList& lockedPoints,
            labelHashSet& badFaces,
            const direction additionalLayers = 0
        );

    // Enumerators

        enum vertexTypes
        {
            NONE = 0,
            SMOOTH = 1,
            FACECENTRE = 2,
            CELLCENTRE = 4,
            BOUNDARY = 8,
            INTERNALBOUNDARY = 16,
            PARALLELBOUNDARY = 32,
            LOCKED = 64
        };

    // Destructor
        ~partTetMesh();

    // Member functions
        //- return a const reference to polyMeshGen
        inline const polyMeshGen& origMesh() const;

        //- return mesh points
        inline const pointFieldPMG& points() const;

        //- labels of points inside the mesh that may be moved
        inline const labelLongList& movableInternalPoints() const;

        //- labels of points at the boundary that may be moved
        inline const labelLongList& movableBoundaryPoints() const;

        //- classification of vertices in the mesh
        inline const LongList<direction>& smoothVertex() const;

        //- find vertices with inverted tets attached and return the number
        //- of inverted vertices
        label findInvertedPoints
        (
            boolList& negativeNode,
            const boolList* activePointPtr = NULL
        ) const;

        //- write mesh to file
        void writeToVTK(const fileName& fName) const;

    // Modifiers
        //- move vertices to their new positions
        //- intended for SMP parallelisation
        void updateVerticesSMP(const List<LongList<labelledPoint> >&);

        //- add the LOCKED flag to mesh points which shall not move
        template<class labelListType>
        void lockPoints(const labelListType& l)
        {
            forAll(l, pI)
                smoothVertex_[l[pI]] |= LOCKED;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "partTetMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
