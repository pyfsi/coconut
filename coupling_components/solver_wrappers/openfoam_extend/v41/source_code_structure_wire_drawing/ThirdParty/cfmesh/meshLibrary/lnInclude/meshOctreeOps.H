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
    checkCellConnectionsOverFaces

Description
    Checks cells if all cells can be visited over face neighbours. Deletes

SourceFiles
    checkCellConnectionsOverFaces.C

\*---------------------------------------------------------------------------*/

#ifndef checkCellConnectionsOverFaces_H
#define checkCellConnectionsOverFaces_H

#include "meshOctree.H"
#include "meshOctreeAddressing.H"
#include "DynList.H"

#include <map>
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Namespace meshOctreeOps Declaration
\*---------------------------------------------------------------------------*/

namespace meshOctreeOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class meshOctreeNeighbourOperator
{
    // Private data
        //- const reference to meshOctree
        const meshOctree& octree_;

        //- leaf at inter-processor boundary
        mutable boolList interProcessorLeaves_;

public:

    // Constructors
        //- Construct from meshOctree
        meshOctreeNeighbourOperator(const meshOctree&);

    // Public member functions
        //- return the number of elements
        label size() const;

        //- find neighbours of a given element
        void operator()(const label, DynList<label>&) const;

        //- reduction operator needed for MPI runs
        template<class labelListType>
        void collectGroups
        (
            std::map<label, DynList<label> >&,
            const labelListType&,
            const DynList<label>&
        ) const;
};

class meshOctreeSelectorOperator
{
    // Private data
        //- const reference to octree box types
        const List<direction>& boxType_;

        //- type of selected boxes
        direction type_;

public:

    // Constructors
        //- construct from a list containing types of leaf boxes
        //- and a selected type
        meshOctreeSelectorOperator
        (
            const List<direction>& boxType,
            const direction type = meshOctreeAddressing::MESHCELL
        );

    // Public member functions
        //- selects all neighbours of a given type
        bool operator()(const label cellI) const;
};

class meshOctreeNonIntersectedSelectorOperator
{
    // Private data
        //- const reference to octree box types
        const boolList& intersectedLeaves_;

public:

    // Constructors
        //- construct from a list of marked intersected leaves
        meshOctreeNonIntersectedSelectorOperator
        (
            const boolList& intersectedLeaves
        );

    // Public member functions
        //- selects all neighbours of a given type
        bool operator()(const label cellI) const;
};

//- find the leaves that are in the groups with the connection to the outer
//- boundary of the octree's bounding box
template<class ListType>
void findOutsideGroups
(
    const meshOctree&,
    const ListType&,
    labelHashSet&
);

//- find the leaves that are in the groups inside the regions of marked boxes
template<class ListType>
void findInsideGroups
(
    const meshOctree&,
    const ListType&,
    labelHashSet&,
    const labelHashSet* outsideGroupsPtr = NULL
);

//- find the leaves that are enclosed within the inside regions of marked boxes
template<class ListType>
void findEnclosedGroups
(
    const meshOctree&,
    const ListType&,
    labelHashSet&,
    const labelHashSet* insideAndOutsideGroupsPtr = NULL
);

//- analyse the tree and find groups of boxes that are of the given type
template<class ListType>
void findGroupsByType
(
    const meshOctree&,
    const ListType&,
    const meshOctreeCubeBasic::typesOfCubes,
    labelHashSet&
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meshOctreeOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "meshOctreeOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
