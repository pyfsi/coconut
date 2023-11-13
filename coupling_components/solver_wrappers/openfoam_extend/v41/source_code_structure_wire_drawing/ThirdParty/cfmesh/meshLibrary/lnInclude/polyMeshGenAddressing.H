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
    polyMeshGenAddressing

Description
    Cell-face mesh analysis engine

SourceFiles
    polyMeshGenAddressing.C
    polyMeshGenAddressingClear.C
    polyMeshGenAddressingCellCells.C
    polyMeshGenAddressingEdgeCells.C
    polyMeshGenAddressingPointCells.C
    polyMeshGenAddressingCells.C
    polyMeshGenAddressingEdgeFaces.C
    polyMeshGenAddressingPointFaces.C
    polyMeshGenAddressingCellEdges.C
    polyMeshGenAddressingPointEdges.C
    polyMeshGenAddressingPointPoints.C
    polyMeshGenAddressingEdges.C
    polyMeshGenAddressingCellCentresAndVols.C
    polyMeshGenAddressingFaceCentresAndAreas.C
    polyMeshGenAddressingUpdateGeometry.C
    polyMeshGenAddressingParallelAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenAddressing_H
#define polyMeshGenAddressing_H

#include "scalarLongList.H"
#include "vectorLongList.H"
#include "edgeLongList.H"
#include "polyMeshGenCells.H"
#include "boolList.H"
#include "HashSet.H"
#include "DynList.H"
#include "Map.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class polyMeshGenAddressing Declaration
\*---------------------------------------------------------------------------*/

class polyMeshGenAddressing
{
    // Mesh data

        //- reference to the mesh
        const polyMeshGenCells& mesh_;

    // Shapes

        //- Edges
        mutable edgeLongList* edgesPtr_;


    // Connectivity

        //- Cell-cells
        mutable VRWGraph* ccPtr_;

        //- Edge-cells
        mutable VRWGraph* ecPtr_;

        //- Point-cells
        mutable VRWGraph* pcPtr_;

        //- Edge-faces
        mutable VRWGraph* efPtr_;

        //- Point-faces
        mutable VRWGraph* pfPtr_;

        //- Cell-edges
        mutable VRWGraph* cePtr_;

        //- Face-edges
        mutable VRWGraph* fePtr_;

        //- Point-edges
        mutable VRWGraph* pePtr_;

        //- Point-points
        mutable VRWGraph* ppPtr_;

        //- Cell-points
        mutable VRWGraph* cpPtr_;


    // Geometric data

        //- Cell centres
        mutable vectorLongList* cellCentresPtr_;

        //- Face centres
        mutable vectorLongList* faceCentresPtr_;

        //- Cell volumes
        mutable scalarLongList* cellVolumesPtr_;

        //- Face areas
        mutable vectorLongList* faceAreasPtr_;

    // Parallel demand driven data

        //- global point labels
        mutable labelLongList* globalPointLabelPtr_;

        //- global face labels
        mutable labelLongList* globalFaceLabelPtr_;

        //- global cell labels
        mutable labelLongList* globalCellLabelPtr_;

        //- global edge labels
        mutable labelLongList* globalEdgeLabelPtr_;

    // Parallel addressing

        //- processors containing a vertex
        mutable VRWGraph* pProcsPtr_;

        //- global to local point adressing
        mutable Map<label>* globalToLocalPointAddressingPtr_;

        //- neighbour processors sharing a point with this processor
        mutable DynList<label>* pointNeiProcsPtr_;

        //- processors containing an edge
        mutable VRWGraph* eProcsPtr_;

        //- global to local edge addressing
        mutable Map<label>* globalToLocalEdgeAddressingPtr_;

        //- neighbour processors sharing an edge with the current processor
        mutable DynList<label>* edgeNeiProcsPtr_;


    // Private member functions

        //- Disallow construct as copy
        polyMeshGenAddressing(const polyMeshGenAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const polyMeshGenAddressing&);


    // Private member functions to calculate demand driven data

    // Topological calculations

        //- Calculate cell-cell addressing
        void calcCellCells() const;

        //- Calculate cell-edge addressing
        void calcEdgeCells() const;

        //- Calculate point-cell addressing
        void calcPointCells() const;

        //- Calculate edge-face addressing
        void calcEdgeFaces() const;

        //- Calculate face-edge addressing
        void calcFaceEdges() const;

        //- Calculate point-face addressing
        void calcPointFaces() const;

        //- Calculate edge list
        void calcCellEdges() const;

        //- Calculate point-edge addressing
        void calcPointEdges() const;

        //- Calculate point-point addressing
        void calcPointPoints() const;

        //- Calculate cell-point addressing
        void calcCellPoints() const;

        //- During edge calculation, a larger set of data is assembled.
        // Create and destroy as a set, using clearOutEdges()
        void calcEdges() const;
        void clearOutEdges() const;

    // Parallel calculations

        //- calculate global point labels
        void calcGlobalPointLabels() const;

        //- calculate global face labels
        void calcGlobalFaceLabels() const;

        //- calculate global cell labels
        void calcGlobalCellLabels() const;

        //- calculate global edge labels
        void calcGlobalEdgeLabels() const;

    // Geometrical calculations

        //- Calculate face centres and areas
        void calcFaceCentresAndAreas() const;
        void makeFaceCentresAndAreas
        (
            const pointFieldPMG& p,
            vectorLongList &fCtrs,
            vectorLongList &fAreas
        ) const;

        //- Calculate cell centres and volumes
        void calcCellCentresAndVols() const;
        void makeCellCentresAndVols
        (
            const vectorLongList &fCtrs,
            const vectorLongList &fAreas,
            vectorLongList &cellCtrs,
            scalarLongList &cellVols
        ) const;

        //- Calculate edge vectors
        void calcEdgeVectors() const;

    // Disallow bitwise construct
        //- Default constructor
        polyMeshGenAddressing();

public:

        // Static data

            ClassName("polyMeshGenAddressing");

    // Constructors

        //- Construct from components
        polyMeshGenAddressing(const polyMeshGenCells& mesh);


    // Destructor

        ~polyMeshGenAddressing();


    // Member Functions

        // Primitive mesh data

            //- Return mesh edges
            const edgeLongList& edges() const;

        // Return mesh connectivity

            const VRWGraph& cellCells() const;
            // faceCells given as faceOwner and faceNeighbour
            const VRWGraph& edgeCells() const;
            const VRWGraph& pointCells() const;

            // faceFaces considered unnecessary
            const VRWGraph& edgeFaces() const;
            const VRWGraph& pointFaces() const;

            const VRWGraph& cellEdges() const;
            const VRWGraph& faceEdges() const;
            // edgeEdges considered unnecessary
            const VRWGraph& pointEdges() const;
            const VRWGraph& pointPoints() const;
            const VRWGraph& cellPoints() const;

        // Parallel addressing

            const labelLongList& globalPointLabel() const;
            const labelLongList& globalFaceLabel() const;
            const labelLongList& globalCellLabel() const;
            const labelLongList& globalEdgeLabel() const;

            const VRWGraph& pointAtProcs() const;
            const Map<label>& globalToLocalPointAddressing() const;
            const DynList<label>& pointNeiProcs() const;
            const VRWGraph& edgeAtProcs() const;
            const Map<label>& globalToLocalEdgeAddressing() const;
            const DynList<label>& edgeNeiProcs() const;

        // Geometric data (raw!)

            const vectorLongList& cellCentres() const;
            const vectorLongList& faceCentres() const;
            const scalarLongList& cellVolumes() const;
            const vectorLongList& faceAreas() const;

    //  Storage management

        //- Print a list of all the currently allocated mesh data
        sizeType printAllocated() const;

        //- Clear geometry
        void clearGeom() const;

        //- Clear topological data
        void clearAddressing() const;

        //- Clear parallel addressing
        void clearParallelAddressing() const;

        //- Clear all geometry and addressing unnecessary for CFD
        void clearOut() const;

        //- Clear everything primitive, geometry and addressing
        void clearAll() const;

    // Update geometry data
        void updateGeometry(const boolList& changedFace) const;
        void updateGeometry() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
