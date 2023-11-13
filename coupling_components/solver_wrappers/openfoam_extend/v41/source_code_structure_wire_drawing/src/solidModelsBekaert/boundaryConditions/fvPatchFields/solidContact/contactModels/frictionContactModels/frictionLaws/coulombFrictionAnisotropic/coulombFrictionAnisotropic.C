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

#include "coulombFrictionAnisotropic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coulombFrictionAnisotropic, 0);
    addToRunTimeSelectionTable
    (
        frictionLaw, coulombFrictionAnisotropic, dictionary
    );



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
coulombFrictionAnisotropic::coulombFrictionAnisotropic
(
    const word& name,
    const frictionContactModel& fricModel,
    const dictionary& dict
)
:
    frictionLaw(name, fricModel, dict),
    frictionLawDict_(dict.subDict("frictionLawDict")),
    primaryFrictionCoeff_
    (
        readScalar(frictionLawDict_.lookup("primaryFrictionCoeff"))
    ),
    secondaryFrictionCoeff_
    (
        readScalar(frictionLawDict_.lookup("secondaryFrictionCoeff"))
    ),
    primaryDirectionSqr_(tensor::zero)
{
    // Read primary direction
    vector primaryDirection =
        vector(frictionLawDict_.lookup("primaryDirection"));

    // Check direction has unit length
    const scalar magPrimDir = mag(primaryDirection);

    if (magPrimDir < SMALL)
    {
        FatalErrorIn("coulombFrictionAnisotropic::coulombFrictionAnisotropic")
            << "The primaryDirection should have a non-ezro length!"
            << abort(FatalError);
    }

    // Normalise
    primaryDirection /= magPrimDir;

    primaryDirectionSqr_ = sqr(primaryDirection);
}


// Construct as a copy
coulombFrictionAnisotropic::coulombFrictionAnisotropic
(
    const coulombFrictionAnisotropic& fricLaw
)
:
    frictionLaw(fricLaw),
    frictionLawDict_(fricLaw.frictionLawDict_),
    primaryFrictionCoeff_(fricLaw.primaryFrictionCoeff_),
    secondaryFrictionCoeff_(fricLaw.secondaryFrictionCoeff_),
    primaryDirectionSqr_(fricLaw.primaryDirectionSqr_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

coulombFrictionAnisotropic::~coulombFrictionAnisotropic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar coulombFrictionAnisotropic::slipTraction(const scalar pressure)
{
    FatalErrorIn
    (
        "scalar coulombFrictionAnisotropic::slipTraction(const scalar pressure)"
    )   << "slipTraction(const scalar pressure) should not be called for this "
        << "law; " << nl
        << "coulombFrictionAnisotropic::slipTraction"
        << "(const scalar pressure, const vector& slipDir)"
        << "should instead be called!"
        << abort(FatalError);

    return GREAT;
}


scalar coulombFrictionAnisotropic::slipTraction
(
    const scalar pressure, const vector& slipDir
)
{
    const scalar magSlipDir = mag(slipDir);

    if (magSlipDir < SMALL)
    {
        return GREAT;
    }

    // Normalise slipDir
    const vector s = slipDir/magSlipDir;

    // Calculate effective coefficient of friction
    const scalar effFricCoeff =
        mag
        (
            primaryFrictionCoeff_*(primaryDirectionSqr_ & s)
          + secondaryFrictionCoeff_*((I - primaryDirectionSqr_) & s)
        );

    return effFricCoeff*pressure;
}


void coulombFrictionAnisotropic::writeDict(Ostream& os) const
{
    word keyword("frictionLawDict");

    os.writeKeyword(keyword)
        << frictionLawDict_;
}

// ************************************************************************* //

} // end of namespace foam
