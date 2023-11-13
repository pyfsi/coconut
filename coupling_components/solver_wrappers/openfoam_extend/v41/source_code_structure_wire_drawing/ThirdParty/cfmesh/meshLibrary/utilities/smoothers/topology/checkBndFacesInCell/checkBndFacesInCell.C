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

#include "checkBndFacesInCell.H"
#include "labelLongList.H"
#include "labelledPair.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

void checkBndFacesInCell::findCellTypes()
{
    const cellListPMG& cells = mesh_.cells();
    cellType_.setSize(cells.size());

    //- find the range for boundary faces
    endBnd_ = -1;
    startBnd_ = mesh_.faces().size();
    forAll(mesh_.boundaries(), patchI)
    {
        const label start = mesh_.boundaries()[patchI].patchStart();
        const label end = start + mesh_.boundaries()[patchI].patchSize();

            endBnd_ = max(endBnd_, end);
            startBnd_ = min(startBnd_, start);
    }

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100)
    # endif
    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];

        direction cellType(INTERNAL);
        forAll(c, fI)
        {
            if( (c[fI] >= startBnd_) && (c[fI] < endBnd_) )
            {
                cellType = BOUNDARY;
                break;
            }
        }

        cellType_[cellI] = cellType;
    }
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //
// Constructors

checkBndFacesInCell::checkBndFacesInCell(polyMeshGen& mesh)
:
    mesh_(mesh),
    cellType_(mesh.cells().size()),
    startBnd_(),
    endBnd_()
{
    findCellTypes();
}

// * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * * * * //

checkBndFacesInCell::~checkBndFacesInCell()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool checkBndFacesInCell::checkCells()
{
    //- remove boundary cells whose are of boundary faces is more than twice
    //- the ares of internal faces
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();
    boolList removeCell(mesh_.cells().size());

    label nRemoved(0);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100) reduction(+ : nRemoved)
    # endif
    forAll(removeCell, cellI)
    {
        removeCell[cellI] = false;

        if( cellType_[cellI] & BOUNDARY )
        {
            scalar areaInternal(0.0);
            scalar areaBnd(0.0);

            const cell& c = cells[cellI];

            forAll(c, fI)
            {
                const vector v = faces[c[fI]].normal(points);

                if( c[fI] >= startBnd_ && c[fI] < endBnd_ )
                {
                    areaBnd += mag(v);
                }
                else
                {
                    areaInternal += mag(v);
                }
            }

            if( areaBnd > 2.01 * areaInternal )
            {
                removeCell[cellI] = true;
                ++nRemoved;
            }
        }
    }

    reduce(nRemoved, sumOp<label>());
    if( nRemoved == 0 )
        return false;

    Info << "Removing " << nRemoved
         << " cells with many boundary faces" << endl;

    polyMeshGenModifier(mesh_).removeCells(removeCell);

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
