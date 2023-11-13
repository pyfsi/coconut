/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "VRWGraphSMPModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcPointCells() const
{
    if( pcPtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcPointCells() const")
            << "pointCells already calculated"
            << abort(FatalError);
    }
    else
    {
        const VRWGraph& pointFaces = this->pointFaces();
        const labelLongList& owner = mesh_.owner();
        const labelLongList& neighbour = mesh_.neighbour();

        //- create the storage
        pcPtr_ = new VRWGraph();
        VRWGraph& pointCellAddr = *pcPtr_;
        pointCellAddr.setSize(mesh_.points().size());

        labelList npc(pointCellAddr.size());

        # ifdef USE_OMP
        # pragma omp parallel
        # endif
        {
            # ifdef USE_OMP
            # pragma omp for schedule(static, 1)
            # endif
            forAll(npc, i)
                npc[i] = 0;

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(pointFaces, pointI)
            {
                DynList<label> activeCells;

                forAllRow(pointFaces, pointI, pfI)
                {
                    const label faceI = pointFaces(pointI, pfI);

                    activeCells.appendIfNotIn(owner[faceI]);
                    if( neighbour[faceI] >= 0 )
                        activeCells.appendIfNotIn(neighbour[faceI]);
                }

                npc[pointI] = activeCells.size();
            }

            # ifdef USE_OMP
            # pragma omp single
            # endif
            {
                VRWGraphSMPModifier(pointCellAddr).setSizeAndRowSize(npc);
            }

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(pointFaces, pointI)
            {
                DynList<label> activeCells;

                forAllRow(pointFaces, pointI, pfI)
                {
                    const label faceI = pointFaces(pointI, pfI);

                    activeCells.appendIfNotIn(owner[faceI]);
                    if( neighbour[faceI] >= 0 )
                        activeCells.appendIfNotIn(neighbour[faceI]);
                }

                pointCellAddr.setRow(pointI, activeCells);
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const VRWGraph& polyMeshGenAddressing::pointCells() const
{
    if( !pcPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& polyMeshGenAddressing::pointCells() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcPointCells();
    }

    return *pcPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
