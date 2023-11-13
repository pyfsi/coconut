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

#include "lubricatedPenaltyFriction.H"
#include "addToRunTimeSelectionTable.H"
#include "solidContactFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lubricatedPenaltyFriction, 0);
    addToRunTimeSelectionTable
    (
        frictionContactModel,
        lubricatedPenaltyFriction,
        dictionary
    );
}


// * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

void Foam::lubricatedPenaltyFriction::lookupFilmFields()
{
    // Lookup solidContact boundary field
    const volVectorField& dispField = mesh_.lookupObject<volVectorField>("DU");
    solidContactFvPatchVectorField& DUpatch =
        const_cast<solidContactFvPatchVectorField&>
        (
            refCast<const solidContactFvPatchVectorField>
            (
                dispField.boundaryField()[slavePatchID()]
            )
        );

    // Check if correct normal contact model is used
    if (DUpatch.normalModelForThisSlave().name() != "lubricatedContact")
    {
        FatalErrorIn
        (
            "Foam::lubricatedPenaltyFriction::lubricatedPenaltyFriction\n"
            "(\n"
            "    const word& name,\n"
            "    const fvPatch& patch,\n"
            "    const dictionary& dict,\n"
            "    const label masterPatchID,\n"
            "    const label slavePatchID\n"
            ") const"
        )   << "lubricatedContact normal contact model must be used"
            << exit(FatalError);
    }

    // Lookup film shear stress
    filmShearStressPtr_ =
        &mesh_.lookupObject<areaVectorField>("filmShearStress");

    filmHydrodynamicPressurePtr_ =
        &mesh_.lookupObject<areaScalarField>("filmHydrodynamicPressure");
}


void Foam::lubricatedPenaltyFriction::calcFrictionPenaltyFactor()
{
    // Set penalty factor using a similar method to the normal
    // contact where we approx penaltyFactor from mechanical properties
    // this can then be scaled using the penaltyScale

    const label masterPatchIndex =  masterPatchID();
    const label slavePatchIndex =  slavePatchID();

    // Lookup shear stiffness
    const volScalarField& mu = mesh_.lookupObject<volScalarField>("mu");

    // Avarage contact patch bulk modulus
    const scalar masterK = gAverage(mu.boundaryField()[masterPatchIndex]);
    const scalar slaveK = gAverage(mu.boundaryField()[slavePatchIndex]);

    // Average contact patch shear modulus
    const scalar modulus = 0.5*(masterK + slaveK);

    // average contact patch face area
    const scalar masterMagSf =
        gAverage(mesh_.magSf().boundaryField()[masterPatchIndex]);
    const scalar slaveMagSf =
        gAverage(mesh_.magSf().boundaryField()[slavePatchIndex]);
    const scalar faceArea = 0.5*(masterMagSf + slaveMagSf);

    // average contact patch cell volume
    scalarField masterV(mesh_.boundary()[masterPatchIndex].size(), 0.0);
    scalarField slaveV(mesh_.boundary()[slavePatchIndex].size(), 0.0);
    const volScalarField::DimensionedInternalField& V = mesh_.V();
    {
        const unallocLabelList& faceCells =
            mesh_.boundary()[masterPatchIndex].faceCells();
        forAll(mesh_.boundary()[masterPatchIndex], facei)
        {
            masterV[facei] = V[faceCells[facei]];
        }
    }
    {
        const unallocLabelList& faceCells =
            mesh_.boundary()[slavePatchIndex].faceCells();
        forAll(mesh_.boundary()[slavePatchIndex], facei)
        {
            slaveV[facei] = V[faceCells[facei]];
        }
    }
    const scalar cellVolume = 0.5*(gAverage(masterV) + gAverage(slaveV));

    // approximate penalty factor based on Hallquist et al.
    // we approximate penalty factor for traction instead of force
    frictionPenaltyFactorPtr_ =
        new scalar(frictionPenaltyScale_*modulus*faceArea/cellVolume);

    Info<< "    friction penalty factor: " << *frictionPenaltyFactorPtr_
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::lubricatedPenaltyFriction::lubricatedPenaltyFriction
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
    frictionLawPtr_(NULL),
    mesh_(patch.boundaryMesh().mesh()),
    writeDebugFile_
    (
        frictionContactModelDict_.lookupOrDefault<Switch>
        (
            "writeDebugFile",
            false
        )
    ),
    slaveTraction_(mesh().boundaryMesh()[slavePatchID].size(), vector::zero),
    prevSlaveTraction_(slaveTraction_),
    slip_(slaveTraction_),
    frictionPenaltyFactorPtr_(NULL),
    frictionPenaltyScale_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("penaltyScale", 1.0)
    ),
    relaxFac_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>
        (
            "relaxationFactor", 0.1
        )
    ),
    contactIterNum_(0),
    infoFreq_
    (
        frictionContactModelDict_.lookupOrDefault<int>("infoFrequency", 1)
    ),
    contactFilePtr_(NULL),
    filmShearStressPtr_(NULL),
    filmHydrodynamicPressurePtr_(NULL)
{
    // Create friction law
    frictionLawPtr_ =
        frictionLaw::New
        (
            frictionContactModelDict_.lookup("frictionLaw"),
            *this,
            frictionContactModelDict_
        ).ptr();

    // master proc open contact info file
    if (Pstream::master() && writeDebugFile_)
    {
        const word masterName = mesh_.boundary()[masterPatchID].name();
        const word slaveName = mesh_.boundary()[slavePatchID].name();

        const fileName contactFileDir = "contact";

        mkDir(contactFileDir);

        contactFilePtr_ =
            new OFstream
            (
                contactFileDir/
                "frictionContact_" + masterName + "_" + slaveName + ".txt"
            );

        OFstream& contactFile = *contactFilePtr_;

        int width = 20;
        contactFile
            << "time";
        contactFile.width(width);
        contactFile
            << "iterNum";
        contactFile.width(width);
        contactFile
            << "penaltyScale";
        contactFile.width(width);
        contactFile
            << "slipFaces";
        contactFile.width(width);
        contactFile
            << "stickFaces";
        contactFile.width(width);
        contactFile
            << "maxMagSlaveTraction";
        contactFile.width(width);
        contactFile
            << "maxMagSlip" << endl;
    }
}


// Construct as a copy
Foam::lubricatedPenaltyFriction::lubricatedPenaltyFriction
(
    const lubricatedPenaltyFriction& fm
)
:
    frictionContactModel(fm),
    frictionContactModelDict_(fm.frictionContactModelDict_),
    frictionLawPtr_(fm.frictionLawPtr_->clone().ptr()),
    mesh_(fm.mesh_),
    writeDebugFile_(fm.writeDebugFile_),
    slaveTraction_(fm.slaveTraction_),
    prevSlaveTraction_(fm.prevSlaveTraction_),
    slip_(fm.slip_),
    frictionPenaltyFactorPtr_(NULL),
    frictionPenaltyScale_(fm.frictionPenaltyScale_),
    relaxFac_(fm.relaxFac_),
    contactIterNum_(fm.contactIterNum_),
    infoFreq_(fm.infoFreq_),
    contactFilePtr_(NULL),
    filmShearStressPtr_(NULL),
    filmHydrodynamicPressurePtr_(NULL)
{
    if (fm.frictionPenaltyFactorPtr_)
    {
        frictionPenaltyFactorPtr_ = new scalar(*fm.frictionPenaltyFactorPtr_);
    }

    if (fm.contactFilePtr_)
    {
        contactFilePtr_ = new OFstream(*fm.contactFilePtr_);
    }

    if (fm.filmShearStressPtr_)
    {
        filmShearStressPtr_ = new areaVectorField(*fm.filmShearStressPtr_);
    }

    if (fm.filmHydrodynamicPressurePtr_)
    {
        filmHydrodynamicPressurePtr_ =
            new areaScalarField(*fm.filmHydrodynamicPressurePtr_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lubricatedPenaltyFriction::correct
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
    const label slavePatchIndex = slavePatchID();

    // Lookup film shear stress field from normal contact model
    if (!filmShearStressPtr_)
    {
        lookupFilmFields();
    }

    // Calculate slave shear traction increments
    const scalarField magSlavePressure = mag(slavePressure);
    label numSlipFaces = 0;
    label numStickFaces = 0;
    scalarField& stickSlip = stickSlipFaces();
    const scalarField oldStickSlip = stickSlip;
    scalar frictionPenaltyFac = frictionPenaltyFactor();
    scalar maxMagSlip = 0.0;
    scalarField slipTraction(magSlavePressure.size(), 0.0);

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

            slaveTraction_[faceI] = -frictionPenaltyFac*slip_[faceI];

            const scalar magSlip = mag(slip_[faceI]);
            maxMagSlip = max(maxMagSlip, magSlip);

            const scalar deltaT = mesh.time().deltaTValue();
            const vector slaveVelocity = slaveDD[faceI]/deltaT;
            const vector masterVelocity = masterDDInterpToSlave[faceI]/deltaT;

            // Traction to cause slipping i.e. the maximum shear traction the
            // face can hold for the given pressure, velocities, temperature,
            // etc.
            // Note: the actual friction law is implemented through the run-time
            // selectable frictionLaw

            // NOTE: Maybe the pressure going into frictionLaw should be
            // magSlavePressure - filmHydrodynamicPressure
            // VS, 2017-07-01
            scalar asperityContactPressureFaceI =
                max
                (
                    magSlavePressure[faceI]
                  - filmHydrodynamicPressurePtr_->internalField()[faceI]
                   *(1 - areaInContact[faceI]),
                    0.0
                );

            slipTraction[faceI] =
                frictionLawPtr_->slipTraction
                (
                    asperityContactPressureFaceI,    // Contact pressure
                    slip_[faceI],               // Slip vector
                    slaveVelocity,              // Velocity of slave face
                    masterVelocity,             // Velocity of master face
                    slavePatchIndex,            // Slave patch index
                    faceI                       // Local slave face ID
                );

            // Add lubricant shear stress
            slipTraction[faceI] +=
                mag(filmShearStressPtr_->internalField()[faceI]);

            if ((mag(slaveTraction_[faceI]) - slipTraction[faceI]) > SMALL)
            {
                // Analogous to plasticity
                // slip is a combination of elastic slip and plastic slip
                // elastic slip should be zero but is finite due to penalty
                // stiffness. plastic slip is the permanent deformation
                slaveTraction_[faceI] =
                    slipTraction[faceI]*(-slip_[faceI]/magSlip);

                numSlipFaces++;
                stickSlip[faceI] = 1;
            }
            else
            {
                numStickFaces++;
                stickSlip[faceI] = 2;
            }
        }
        else
        {
            // No asperity friction if pressure is negative or zero or face is
            // not in contact, only film shear stress
            slaveTraction_[faceI] = filmShearStressPtr_->internalField()[faceI];
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

    scalar maxMagSlaveTraction = 0.0;
    if (slaveTraction_.size() > 0)
    {
        maxMagSlaveTraction = max(mag(slaveTraction_));
    }
    reduce(maxMagSlaveTraction, maxOp<scalar>());
    reduce(numSlipFaces, sumOp<int>());
    reduce(numStickFaces, sumOp<int>());

    // Writes to contact info file
    if
    (
        Pstream::master()
     && (contactIterNum_++ %  infoFreq_ == 0)
     && writeDebugFile_
    )
    {
        OFstream& contactFile = *contactFilePtr_;
        int width = 20;
        contactFile
            << mesh.time().value();
        contactFile.width(width);
        contactFile
            << contactIterNum_;
        contactFile.width(width);
        contactFile
            << frictionPenaltyScale_;
        contactFile.width(width);
        contactFile
            << numSlipFaces;
        contactFile.width(width);
        contactFile
            << numStickFaces;
        contactFile.width(width);
        contactFile
            << maxMagSlaveTraction;
        contactFile.width(width);
        contactFile
            << maxMagSlip << endl;
    }
}


void Foam::lubricatedPenaltyFriction::writeDict(Ostream& os) const
{
    word keyword(name()+"FrictionModelDict");
    os.writeKeyword(keyword)
        << frictionContactModelDict_;
}


// ************************************************************************* //
