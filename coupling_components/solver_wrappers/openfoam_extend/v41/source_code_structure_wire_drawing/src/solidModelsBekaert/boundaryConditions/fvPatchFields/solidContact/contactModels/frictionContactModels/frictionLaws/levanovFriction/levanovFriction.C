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

#include "levanovFriction.H"
#include "addToRunTimeSelectionTable.H"
//#include "zeroGradientFvPatchFields.H"
#include "frictionContactModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(levanovFriction, 0);
    addToRunTimeSelectionTable
    (
        frictionLaw, levanovFriction, dictionary
    );


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void levanovFriction::lookupSigmaCauchyField()
{
    if (sigmaCauchyPtr_)
    {
        FatalErrorIn("void levanovFriction::lookupSigmaCauchyField()")
            << "sigmaCauchy field pointer already set" << abort(FatalError);
    }

    const volSymmTensorField& sigmaCauchyRef =
        mesh_.objectRegistry::lookupObject<volSymmTensorField>("sigmaCauchy");

    sigmaCauchyPtr_ = &sigmaCauchyRef;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
levanovFriction::levanovFriction
(
    const word& name,
    const frictionContactModel& fricModel,
    const dictionary& dict
)
:
    frictionLaw(name, fricModel, dict),
    mesh_(fricModel.mesh()),
    frictionLawDict_(dict.subDict("frictionLawDict")),
    sigmaCauchyPtr_(NULL),
    k0_(readScalar(frictionLawDict_.lookup("k0"))),
    kf0_(readScalar(frictionLawDict_.lookup("kf0"))),
    C_(frictionLawDict_.lookupOrDefault<scalar>("C", 0.005)),
    a_(frictionLawDict_.lookupOrDefault<scalar>("a", 1.25))
{}


// Construct as a copy
levanovFriction::levanovFriction
(
    const levanovFriction& fricLaw
)
:
    frictionLaw(fricLaw),
    mesh_(fricLaw.mesh_),
    frictionLawDict_(fricLaw.frictionLawDict_),
    k0_(fricLaw.k0_),
    kf0_(fricLaw.kf0_),
    C_(fricLaw.C_),
    a_(fricLaw.a_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

levanovFriction::~levanovFriction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar levanovFriction::slipTraction
(
     const scalar contactPressure,         // Contact pressure
     const vector& faceSlip,               // Slip vector
     const vector& slaveFaceVelocity,      // Velocity of slave face
     const vector& masterFaceVelocity,     // Velocity of master face
     const label slavePatchIndex,          // Slave patch index
     const label faceIndex                 // Local slave face ID
)
{
   const symmTensor curSigma =
       sigmaCauchy().boundaryField()[slavePatchIndex][faceIndex];

    // Flow strength of the steel (i.e. axial yield stress measured from tensile
    // test)
    const scalar kf = max(kf0_, sqrt((3.0/2.0)*magSqr(dev(curSigma))));

    const scalar slipVelocity = mag(slaveFaceVelocity - masterFaceVelocity);

    // Slip factor
    const scalar fac_slip =
       2.0/Foam::mathematicalConstant::pi*Foam::atan(slipVelocity/C_);

    // Pressure factor
    const scalar fac_pres = 1.0 - Foam::exp(-a_*mag(contactPressure)/kf);

    // Return shear slip stress
    return k0_*fac_slip*fac_pres*kf/sqrt(3.0);
}

void levanovFriction::writeDict(Ostream& os) const
{
    word keyword("frictionLawDict");

    os.writeKeyword(keyword)
        << frictionLawDict_;
}

// ************************************************************************* //

} // end of namespace foam
