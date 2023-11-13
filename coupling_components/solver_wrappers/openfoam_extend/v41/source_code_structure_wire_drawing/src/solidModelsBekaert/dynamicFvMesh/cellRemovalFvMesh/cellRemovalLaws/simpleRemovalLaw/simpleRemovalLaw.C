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

#include "simpleRemovalLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "removeCells.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleRemovalLaw, 0);
    addToRunTimeSelectionTable(cellRemovalLaw, simpleRemovalLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::simpleRemovalLaw::simpleRemovalLaw
(
    const word& name,
    fvMesh& mesh,
    const dictionary& dict
)
:
    cellRemovalLaw(name, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::simpleRemovalLaw::~simpleRemovalLaw()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::labelList Foam::simpleRemovalLaw::cellsToRemove()
{
    Info<< "simpleRemovalLaw::updateMesh()" << endl;

    // Remove the first ten cells

    labelList cellsToRemove(10, -1);

    forAll(cellsToRemove, cellI)
    {
        cellsToRemove[cellI] = cellI;
    }

    return cellsToRemove;
}


// ************************************************************************* //
