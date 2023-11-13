/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "boundingBoxRemovalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "removeCells.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boundingBoxRemovalLaw, 0);
    addToRunTimeSelectionTable
    (
        cellRemovalLaw, boundingBoxRemovalLaw, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::boundingBoxRemovalLaw::boundingBoxRemovalLaw
(
    const word& name,
    fvMesh& mesh,
    const dictionary& dict
)
:
    cellRemovalLaw(name, mesh, dict),
    bb_
    (
        point(dict.lookup("boundingBoxMin")),
        point(dict.lookup("boundingBoxMax"))
    ),
    cellZoneID_
    (
        mesh.cellZones().findZoneID
        (
            dict.lookupOrDefault<word>("cellZone", "notSpecified")
        )
    )
{
    if (cellZoneID_ == -1)
    {
        InfoIn("boundingBoxRemovalLaw::boundingBoxRemovalLaw()")
            << "    No cellZone specified" << endl;
    }
    else
    {
        InfoIn("boundingBoxRemovalLaw::boundingBoxRemovalLaw()")
            << "    Using cellZone" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::boundingBoxRemovalLaw::~boundingBoxRemovalLaw()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::labelList Foam::boundingBoxRemovalLaw::cellsToRemove()
{
    // Set of cellIDs to be removed
    labelHashSet cellsToRemoveSet;

    // Find cells to remove

    const vectorField& C = mesh().C();

    if (cellZoneID_ != -1)
    {
        const labelList& zoneCells = mesh().cellZones()[cellZoneID_];

        forAll(zoneCells, cI)
        {
            const label cellID = zoneCells[cI];

            if (!bb_.contains(C[cellID]))
            {
                cellsToRemoveSet.insert(cellID);
            }
        }
    }
    else
    {
        forAll(C, cellID)
        {
            if (!bb_.contains(C[cellID]))
            {
                cellsToRemoveSet.insert(cellID);
            }
        }
    }

    return cellsToRemoveSet.toc();
}


// ************************************************************************* //
