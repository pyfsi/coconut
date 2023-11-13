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

#include "stribeckCoulombFriction.H"
#include "addToRunTimeSelectionTable.H"
//#include "zeroGradientFvPatchFields.H"
//#include "frictionContactModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stribeckCoulombFriction, 0);
    addToRunTimeSelectionTable
    (
        frictionLaw, stribeckCoulombFriction, dictionary
    );



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar stribeckCoulombFriction::frictionCoeff
(
    const scalar pressure, const vector& velocity
) const
{
    const scalar viscosity =
        frictionLaw::frictionLawDict().lookupOrDefault
            (
                "viscosity",
                0.0009
            );
    const scalar surfaceRoughnessEQ =
        frictionLaw::frictionLawDict().lookupOrDefault
            (
                "surfaceRoughnessEQ",
                0.2e-6
            );
    const scalar minAxialFricCoeff =
        frictionLaw::frictionLawDict().lookupOrDefault
            (
                "minAxialFricCoeff",
                0.01
            );
    const scalar maxAxialFricCoeff =
        frictionLaw::frictionLawDict().lookupOrDefault
            (
                "maxAxialFricCoeff",
                0.15
            );
    const scalar minAxialFricCoeffLubNr =
        frictionLaw::frictionLawDict().lookupOrDefault
            (
                "minAxialFricCoeffLubNr",
                2.2e-5
            );
    const scalar axialHydDynSlope =
        frictionLaw::frictionLawDict().lookupOrDefault
            (
                "axialHydDynSlope",
                10.0
            );

    const scalar magVelocity = mag(velocity);

    const scalar lubNr = magVelocity*viscosity/(pressure*surfaceRoughnessEQ);

    scalar fricCoeff = 0.0;
    if(lubNr <= 2.0 * minAxialFricCoeffLubNr)
    {
        fricCoeff = minAxialFricCoeff + (maxAxialFricCoeff-minAxialFricCoeff)*
        (
            1 -
            tanh(1/(1e5*minAxialFricCoeffLubNr)*lubNr/minAxialFricCoeffLubNr)
        );

        if(fricCoeff > maxAxialFricCoeff || fricCoeff < minAxialFricCoeff)
        {
            FatalErrorIn("scalar stribeckCoulombFriction::frictionCoeff")
                << "Friction coefficient " << fricCoeff << "out of range"
                << abort(FatalError);
        }
    }
    else
    {
        fricCoeff = minAxialFricCoeff +
                    (lubNr - 2.0 * minAxialFricCoeffLubNr)*axialHydDynSlope;
    }
    return fricCoeff;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
stribeckCoulombFriction::stribeckCoulombFriction
(
    const word& name,
    const frictionContactModel& fricModel,
    const dictionary& dict
)
:
    frictionLaw(name, fricModel, dict),
    frictionLawDict_(dict.subDict("frictionLawDict"))
{}


// Construct as a copy
stribeckCoulombFriction::stribeckCoulombFriction
(
    const stribeckCoulombFriction& fricLaw
)
:
    frictionLaw(fricLaw),
    frictionLawDict_(fricLaw.frictionLawDict_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

stribeckCoulombFriction::~stribeckCoulombFriction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar stribeckCoulombFriction::slipTraction(const scalar pressure)
{
    vector tempVector = vector::zero;

    return frictionCoeff(pressure, tempVector)*pressure;
}


scalar stribeckCoulombFriction::slipTraction
(
    const scalar pressure,
    const vector&
)
{
    vector tempVector = vector::zero;

    return frictionCoeff(pressure, tempVector)*pressure;
}

scalar stribeckCoulombFriction::slipTraction
(
    const scalar pressure,
    const vector& slipDir,
    const vector& slaveFaceVelocity,
    const vector& masterFaceVelocity
)
{
    vector velocity = 0.5*(slaveFaceVelocity + masterFaceVelocity);

    return frictionCoeff(pressure, velocity)*pressure;
}

void stribeckCoulombFriction::writeDict(Ostream& os) const
{
    word keyword("frictionLawDict");

    os.writeKeyword(keyword)
        << frictionLawDict_;
}

// ************************************************************************* //

} // end of namespace foam
