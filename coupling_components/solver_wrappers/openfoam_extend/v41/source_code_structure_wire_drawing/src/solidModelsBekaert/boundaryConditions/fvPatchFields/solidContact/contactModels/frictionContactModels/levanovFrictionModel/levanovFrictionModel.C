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

#include "levanovFrictionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(levanovFrictionModel, 0);
    addToRunTimeSelectionTable
    (
        frictionContactModel,
        levanovFrictionModel,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::levanovFrictionModel::levanovFrictionModel
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID
)
:
    frictionContactModel
    (
        name,
        patch,
        dict,
        masterPatchID,
        slavePatchID
    ),
    frictionContactModelDict_(dict.subDict(name + "FrictionModelDict")),
    mesh_(patch.boundaryMesh().mesh()),
    slaveTraction_(mesh().boundaryMesh()[slavePatchID].size(), vector::zero),
    prevSlaveTraction_(slaveTraction_),
    slip_(slaveTraction_),
    relaxFac_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>
        (
            "relaxationFactor", 0.1
        )
    ),
    k0_(readScalar(frictionContactModelDict_.lookup("k0"))),
    kf0_(readScalar(frictionContactModelDict_.lookup("kf0"))),
    C_(frictionContactModelDict_.lookupOrDefault<scalar>("C", 0.005)),
    a_(frictionContactModelDict_.lookupOrDefault<scalar>("a", 1.25))
{}


// Construct as a copy
Foam::levanovFrictionModel::levanovFrictionModel
(
    const levanovFrictionModel& fm
)
:
    frictionContactModel(fm),
    frictionContactModelDict_(fm.frictionContactModelDict_),
    mesh_(fm.mesh_),
    slaveTraction_(fm.slaveTraction_),
    prevSlaveTraction_(fm.prevSlaveTraction_),
    slip_(fm.slip_),
    relaxFac_(fm.relaxFac_),
    k0_(fm.k0_),
    kf0_(fm.kf0_),
    C_(fm.C_),
    a_(fm.a_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::levanovFrictionModel::correct
(
    const vectorField& slavePressure,
    const vectorField& slaveFaceNormals,
    const scalarField& areaInContact,
    const vectorField& slaveDD,
    const vectorField& masterDDInterpToSlave
)
{
    // Preliminaries
    const fvMesh& mesh = mesh_;

    // Calculate slave shear traction increments
    const scalarField magSlavePressure = mag(slavePressure);
    scalarField& stickSlip = stickSlipFaces();
    const scalarField oldStickSlip = stickSlip;
    scalarField slipTraction(magSlavePressure.size(), 0.0);

    // Pressure factor
    const scalarField fac_pres =
        1.0 - Foam::exp(-a_*mag(magSlavePressure)/kf0_);

    forAll(magSlavePressure, faceI)
    {
        if (areaInContact[faceI] > SMALL)
        {
            // Compute slip as the we need the difference of DD between the
            // master and slave
            slip_[faceI] = slaveDD[faceI] - masterDDInterpToSlave[faceI];

            // The shear traction direction is got by removing the normal
            // component of the DD
            //     (I - sqr(n)) removes the normal
            //    sqr(n) would remove the shear
            slip_[faceI] = (I - sqr(slaveFaceNormals[faceI])) & slip_[faceI];

            const scalar magSlip = mag(slip_[faceI]) + SMALL;

            const scalar deltaT = mesh.time().deltaTValue();
            const vector slaveVelocity = slaveDD[faceI]/deltaT;
            const vector masterVelocity = masterDDInterpToSlave[faceI]/deltaT;

            // Slip velocity
            const scalar slipVelocity = mag(slaveVelocity - masterVelocity);

            // Slip factor
            const scalar fac_slip =
               2.0/Foam::mathematicalConstant::pi*Foam::atan(slipVelocity/C_);

            // Return shear slip stress
            slaveTraction_[faceI] =
                k0_*fac_slip*fac_pres[faceI]*kf0_/sqrt(3.0)
               *(-slip_[faceI]/magSlip);
        }
        else
        {
            // No friction if pressure is negative or zero or face is not in
            // contact
            slaveTraction_[faceI] = vector::zero;
            slipTraction[faceI] = 0.0;
            stickSlip[faceI] = 0;
        }
    }

    // StickSlip field is just for visualisation but we will under-relax it to
    // allow us to see if a face is jumping between stick and slip
    stickSlip = relaxFac_*stickSlip + (1.0 - relaxFac_)*oldStickSlip;

    // Under-relax traction
    slaveTraction_ =
        relaxFac_*slaveTraction_ + (1.0 - relaxFac_)*prevSlaveTraction_;

    // Update the previous traction
    prevSlaveTraction_ = slaveTraction_;
}


void Foam::levanovFrictionModel::autoMap(const fvPatchFieldMapper& m)
{
    frictionContactModel::autoMap(m);

    if (debug)
    {
        InfoIn
        (
            "void levanovFrictionModel::autoMap(const fvPatchFieldMapper& m)"
        )   << "autoMap" << endl;
    }

    slaveTraction_.autoMap(m);
    prevSlaveTraction_.autoMap(m);
    slip_.autoMap(m);
}


void Foam::levanovFrictionModel::writeDict(Ostream& os) const
{
    word keyword(name()+"FrictionModelDict");
    os.writeKeyword(keyword)
        << frictionContactModelDict_;
}


// ************************************************************************* //
