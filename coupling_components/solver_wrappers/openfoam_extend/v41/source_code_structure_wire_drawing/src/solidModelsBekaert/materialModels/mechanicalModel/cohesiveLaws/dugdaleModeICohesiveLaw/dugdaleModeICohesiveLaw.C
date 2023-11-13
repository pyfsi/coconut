/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "dugdaleModeICohesiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "cohesiveFvPatch.H"
#include "cohesiveLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dugdaleModeICohesiveLaw, 0);
    addToRunTimeSelectionTable
    (
        cohesiveLaw, dugdaleModeICohesiveLaw, dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::dugdaleModeICohesiveLaw::dugdaleModeICohesiveLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
  cohesiveLaw(name, mesh, dict),
  GIc_(dict.lookup("GIc")),
  sigmaMax_(dict.lookup("sigmaMax"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dugdaleModeICohesiveLaw::~dugdaleModeICohesiveLaw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::dugdaleModeICohesiveLaw::materials() const
{
  notImplemented(type() + "::materials()");

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "materials",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
}

Foam::tmp<Foam::surfaceScalarField>
Foam::dugdaleModeICohesiveLaw::sigmaMax() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "sigmaMax",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            sigmaMax_
        )
    );
}

Foam::tmp<Foam::surfaceScalarField>
Foam::dugdaleModeICohesiveLaw::tauMax() const
{
    notImplemented
    (
        "Foam::tmp<Foam::surfaceScalarField>"
        "Foam::dugdaleModeICohesiveLaw::tauMax() const"
    );

    // Keep the compiler happy; this is not used
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "tauMax",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            sigmaMax_
        )
    );
}

Foam::tmp<Foam::surfaceScalarField>
Foam::dugdaleModeICohesiveLaw::GIc() const
{
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "GIc",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            GIc_
        )
    );
}

Foam::tmp<Foam::surfaceScalarField>
Foam::dugdaleModeICohesiveLaw::GIIc() const
{
    notImplemented
    (
        "Foam::tmp<Foam::surfaceScalarField>"
        "Foam::dugdaleModeICohesiveLaw::GIIc() const"
    );

    // Keep the compiler happy; this is not used
    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "GIIc",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            GIc_
        )
    );
}


void Foam::dugdaleModeICohesiveLaw::damageTractions
(
    scalar& tN,
    scalar& tS,
    const scalar deltaN,
    const scalar deltaS,
    const scalar GI,
    const scalar GII,
    const label internalFaceID
) const
{
    // Set normal traction to sigmaMax
    tN = sigmaMax_.value();

    // Set shear traction to zero
    tS = 0.0;
}


void Foam::dugdaleModeICohesiveLaw::damageTractions
(
    scalar& tN,
    scalar& tS,
    const scalar deltaN,
    const scalar deltaS,
    const scalar GI,
    const scalar GII,
    const label internalFaceID,
    const scalarField& globalPatchMaterials
 ) const
{
    damageTractions(tN, tS, deltaN, deltaS, GI, GII, internalFaceID);
}


// ************************************************************************* //
