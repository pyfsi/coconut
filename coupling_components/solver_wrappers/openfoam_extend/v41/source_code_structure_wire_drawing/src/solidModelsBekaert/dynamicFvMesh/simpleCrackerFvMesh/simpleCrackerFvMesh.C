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

#include "simpleCrackerFvMesh.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "cohesiveZoneIncrementalFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleCrackerFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, simpleCrackerFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::simpleCrackerFvMesh::simpleCrackerFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io)
    // dict_
    // (
    //     IOdictionary
    //     (
    //         IOobject
    //         (
    //             "dynamicMeshDict",
    //             time().constant(),
    //             *this,
    //             IOobject::MUST_READ,
    //             IOobject::NO_WRITE
    //         )
    //     ).subDict(typeName + "Coeffs")
    // )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleCrackerFvMesh::~simpleCrackerFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::simpleCrackerFvMesh::update()
{
    // Lookup DU field
    const volVectorField& DU = this->lookupObject<volVectorField>("DU");

    // Find cohesive patch and update crack

    label nFacesToBreak = 0;

    forAll(DU.boundaryField(), patchI)
    {
        if
        (
            DU.boundaryField()[patchI].type()
         == cohesiveZoneIncrementalFvPatchVectorField::typeName
        )
        {
            // Cast patch and call update crack
            const cohesiveZoneIncrementalFvPatchVectorField& DUpatch =
                refCast<const cohesiveZoneIncrementalFvPatchVectorField>
                (
                    DU.boundaryField()[patchI]
                );

            // Cast away const-ness to call update crack
            cohesiveZoneIncrementalFvPatchVectorField& DUpatchC =
                const_cast<cohesiveZoneIncrementalFvPatchVectorField&>(DUpatch);

            nFacesToBreak = DUpatchC.updateCrack();

            break;
        }
    }

    Info<< nl << "Breaking " << nFacesToBreak << " faces" << nl << endl;

    return nFacesToBreak;
}


// ************************************************************************* //
