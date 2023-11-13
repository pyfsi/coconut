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

#include "sigmaEqLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "removeCells.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sigmaEqLaw, 0);
    addToRunTimeSelectionTable(cellRemovalLaw, sigmaEqLaw, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::sigmaEqLaw::sigmaEqLaw
(
    const word& name,
    fvMesh& mesh,
    const dictionary& dict
)
:
    cellRemovalLaw(name, mesh, dict),
    sigmaEqCrit_(readScalar(dict.lookup("sigmaEqCrit")))
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //

Foam::sigmaEqLaw::~sigmaEqLaw()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::labelList Foam::sigmaEqLaw::cellsToRemove()
{
    // Lookup gradient of displacement
    const volTensorField& gradU =
        mesh().lookupObject<volTensorField>("grad(U)");

    const volScalarField& mu = mesh().lookupObject<volScalarField>("mu");
    const volScalarField& lambda =
        mesh().lookupObject<volScalarField>("lambda");

    // Calculate current stress
    const volSymmTensorField sigma = mu*twoSymm(gradU) + lambda*I*tr(gradU);

    // Calculate current equivalent stress
    const volScalarField sigmaEq = sqrt((3.0/2.0)*magSqr(dev(sigma)));

    // Find cells with sigmaEq greater than the critical value
    const scalarField& sigmaEqI = sigmaEq.internalField();

    labelHashSet cellsToRemove;

    forAll(sigmaEqI, cellI)
    {
        if (sigmaEqI[cellI] > sigmaEqCrit_)
        {
            cellsToRemove.insert(cellI);
        }
    }

    return cellsToRemove.toc();
}


// ************************************************************************* //
