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

#include "damageLaw.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(damageLaw, 0);
defineRunTimeSelectionTable(damageLaw, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

damageLaw::damageLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const volScalarField& curMaterial
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    coupled_(dict.lookupOrDefault<Switch>("coupled", false)),
    damagePtr_(NULL),
    curMaterial_(curMaterial)
{
    // Create or lookup damage field
    if (mesh_.foundObject<volScalarField>("damage"))
    {
        // Lookup the field from the objectRegistry and set the pointer
        damagePtr_ =
            &const_cast<volScalarField&>
            (
                mesh_.lookupObject<volScalarField>("damage")
            );

        Info<< "Looking up damage field from the objectRegistry"
            << " and setting the field pointer" << endl;
    }
    else
    {
        // Create the field and read from the disk if present
        damagePtr_ =
            new volScalarField
            (
                IOobject
                (
                    "damage",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0.0)
            );

        Info<< "Creating damage field"<< endl;
    }

    Info<< "    coupled damage: " << coupled_ << endl;
}


// * * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //

damageLaw::~damageLaw()
{
    deleteDemandDrivenData(damagePtr_);
}


// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
