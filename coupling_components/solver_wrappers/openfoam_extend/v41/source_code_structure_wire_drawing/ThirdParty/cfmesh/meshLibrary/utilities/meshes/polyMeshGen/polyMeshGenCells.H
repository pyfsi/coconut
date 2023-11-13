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
    polyMeshGenCells

Description
    Mesh with geometry and topology information

SourceFiles
    polyMeshGenCells.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenCells_H
#define polyMeshGenCells_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenFaces.H"
#include "cellListPMG.H"
#include "polyMeshGenCellZones.H"

namespace Foam
{

class polyMeshGenAddressing;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenCells
: public polyMeshGenFaces, public polyMeshGenCellZones
{
protected:

    // Protected data
        //- list of cells
        cellListPMG cells_;

        //- cell subsets
        std::map<label, meshSubset> cellSubsets_;

        //- primitive mesh which calculates addressing
        mutable polyMeshGenAddressing* addressingDataPtr_;

    // Private member functions
        //- calculate owner and neighbour
        void calculateOwnersAndNeighbours() const;

        //- calculate mesh addressing
        void calculateAddressingData() const;

        //- clear all pointer data
        void clearOut() const;

    // Disallow bitwise assignment
        void operator=(const polyMeshGenCells&);

        polyMeshGenCells(const polyMeshGenCells&);

public:

    friend class polyMeshGenModifier;

    // Constructors
        //- Construct from Time
        polyMeshGenCells
        (
            const Time&,
            const fileName instance = "constant",
            const fileName meshDir = "polyMesh"
        );

        //- Construct from components without the boundary
        polyMeshGenCells
        (
            const Time&,
            const pointField& points,
            const faceList& faces,
            const cellList& cells,
            const fileName instance = "constant",
            const fileName meshDir = "polyMesh"
        );

        //- Construct from components with the boundary
        polyMeshGenCells
        (
            const Time&,
            const pointField& points,
            const faceList& faces,
            const cellList& cells,
            const wordList& patchNames,
            const labelList& patchStart,
            const labelList& nFacesInPatch,
            const fileName instance = "constant",
            const fileName meshDir = "polyMesh"
        );

        //- Construct from components with the boundary and boundary types
        polyMeshGenCells
        (
            const Time&,
            const pointField& points,
            const faceList& faces,
            const cellList& cells,
            const wordList& patchNames,
            const wordList& patchTypes,
            const labelList& patchStart,
            const labelList& nFacesInPatch,
            const labelList& procPatchMyProcNo,
            const labelList& procPatchNeighbProcNo,
            const fileName instance = "constant",
            const fileName meshDir = "polyMesh"
        );

    // Destructor
        ~polyMeshGenCells();

    // Member functions
        //- access to cells
        inline const cellListPMG& cells() const;

        //- addressing which may be needed
        const polyMeshGenAddressing& addressingData() const;

        //- clear addressing data
        void clearAddressingData() const;

    // Subsets
        label addCellSubset(const word&);
        void removeCellSubset(const label);
        word cellSubsetName(const label) const;
        label cellSubsetIndex(const word&) const;
        inline void addCellToSubset(const label, const label);
        inline void removeCellFromSubset(const label, const label);
        inline void cellInSubsets(const label, DynList<label>&) const;
        inline void cellSubsetIndices(DynList<label>&) const;
        template<class ListType>
        inline void cellsInSubset(const label, ListType&) const;
        template<class ListType>
        inline void updateCellSubsets(const ListType&);
        inline void updateCellSubsets(const VRWGraph&);

    // Read mesh
        void read();

        void readFromLatestTime();

    // Write mesh
        void write() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenCellsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
