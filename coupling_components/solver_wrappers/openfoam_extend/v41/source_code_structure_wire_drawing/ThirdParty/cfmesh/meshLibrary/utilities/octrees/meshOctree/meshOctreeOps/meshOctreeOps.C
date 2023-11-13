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

#include "meshOctreeOps.H"
#include "labelLongList.H"
#include "labelledPair.H"

# ifdef USE_OMP
#include <omp.h>
# endif


#include "helperFunctions.H"

//#define DEBUGCheck

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

namespace meshOctreeOps
{

meshOctreeNeighbourOperator::meshOctreeNeighbourOperator(const meshOctree& mo)
:
    octree_(mo),
    interProcessorLeaves_()
{
    if( Pstream::parRun() )
    {
        interProcessorLeaves_.setSize(octree_.numberOfLeaves());
        interProcessorLeaves_ = false;
    }
}

label meshOctreeNeighbourOperator::size() const
{
    return octree_.numberOfLeaves();
}

void meshOctreeNeighbourOperator::operator()
(
    const label leafI,
    DynList<label>& neighbourLeaves
) const
{
    neighbourLeaves.clear();

    for(label fI=0;fI<6;++fI)
    {
        DynList<label> neiLeaves;

        octree_.findNeighboursInDirection(leafI, fI, neiLeaves);

        forAll(neiLeaves, i)
        {
            const label nei = neiLeaves[i];

            if( nei < 0 )
            {
                if( nei == meshOctreeCube::OTHERPROC )
                    interProcessorLeaves_[leafI] = true;

                continue;
            }

            neighbourLeaves.append(nei);
        }
    }
}

//- selection operator for all leaves
meshOctreeSelectorOperator::meshOctreeSelectorOperator
(
    const List<direction>& boxType,
    const direction type)
:
    boxType_(boxType),
    type_(type)
{}

bool meshOctreeSelectorOperator::operator()(const label leafI) const
{
    return bool(boxType_[leafI] & type_);
}

//- selection operator for the intersected leaves
meshOctreeNonIntersectedSelectorOperator::
meshOctreeNonIntersectedSelectorOperator
(
    const boolList& intersectedLeaves
)
:
    intersectedLeaves_(intersectedLeaves)
{}

bool meshOctreeNonIntersectedSelectorOperator::operator()
(
    const label leafI
) const
{
    return !intersectedLeaves_[leafI];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meshConnectionsHelper

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
