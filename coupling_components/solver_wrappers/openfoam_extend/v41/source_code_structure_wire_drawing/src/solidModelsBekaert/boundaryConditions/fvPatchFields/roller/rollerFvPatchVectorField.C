/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "rollerFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "mathematicalConstants.H"
#include "RodriguesRotation.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void rollerFvPatchVectorField::makeReferenceCentre() const
{
    if (origFaceCentresPtr_)
    {
        FatalErrorIn
        (
            "void rollerFvPatchVectorField::makeReferenceCentre() const"
        ) << "points already set" << abort(FatalError);
    }

    vectorField curFaceCentres = vectorField(patch().patch().faceCentres());

    // Remove current displacement from the current mesh position to find the
    // initial positions

    // Lookup total displacement field

    // U
    const fvPatchField<vector>& U =
        patch().lookupPatchField<volVectorField, vector>("U");

    // Set original positions

    origFaceCentresPtr_ = new vectorField(curFaceCentres - U);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rollerFvPatchVectorField::
rollerFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    displacement_(vector::zero),
    dispRampTime_(0.0),
    rpm_(0.0),
    currentRpm_(0.0),
    rotRampTime_(0.0),
    rotationAxis_(vector::zero),
    initialRotationOrigin_(vector::zero),
    currentAxisDisplacement_(vector::zero),
    origFaceCentresPtr_(NULL),
    curAngle_(0.0),
    curTimeIndex_(-1)
{}


rollerFvPatchVectorField::
rollerFvPatchVectorField
(
    const rollerFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    displacement_(ptf.displacement_),
    dispRampTime_(ptf.dispRampTime_),
    rpm_(ptf.rpm_),
    currentRpm_(ptf.currentRpm_),
    rotRampTime_(ptf.rotRampTime_),
    rotationAxis_(ptf.rotationAxis_),
    initialRotationOrigin_(ptf.initialRotationOrigin_),
    currentAxisDisplacement_(ptf.currentAxisDisplacement_),
    origFaceCentresPtr_(NULL),
    curAngle_(ptf.curAngle_),
    curTimeIndex_(ptf.curTimeIndex_)
{
    if (ptf.origFaceCentresPtr_)
    {
        origFaceCentresPtr_ = new vectorField(*ptf.origFaceCentresPtr_);
    }
}


rollerFvPatchVectorField::
rollerFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedDisplacementFvPatchVectorField(p, iF, dict),
    displacement_(dict.lookup("displacement")),
    dispRampTime_(readScalar(dict.lookup("dispRampTime"))),
    rpm_(readScalar(dict.lookup("rpm"))),
    currentRpm_(dict.lookupOrDefault<scalar>("currentRpm", 0.0)),
    rotRampTime_(readScalar(dict.lookup("rotationRampTime"))),
    rotationAxis_(dict.lookup("rotationAxis")),
    initialRotationOrigin_(vector::zero),
    currentAxisDisplacement_(vector::zero),
    origFaceCentresPtr_(NULL),
    curAngle_(dict.lookupOrDefault<scalar>("curAngle", 0.0)),
    curTimeIndex_(-1)
{
    Info<< patch().name() << ": " << type() << endl;

    // Set initial rotation origin
    if (dict.found("initialRotationOrigin"))
    {
        Info<< "    Reading initialRotationOrigin" << endl;
        initialRotationOrigin_ = dict.lookup("initialRotationOrigin");
    }
    else
    {
        Info<< "    Calculating initialRotationOrigin" << endl;
        const scalar gSumMagSf = gSum(patch().magSf());
        if (gSumMagSf > SMALL)
        {
            initialRotationOrigin_ =
                gSum(patch().Cf()*patch().magSf())/gSumMagSf;
        }
    }

    if (dict.found("currentAxisDisplacement"))
    {
        Info<< "Reading currentAxisDisplacement" << endl;
        currentAxisDisplacement_ = dict.lookup("currentAxisDisplacement");
    }

    if
    (
        dimensionedInternalField().name() != "DU"
     && dimensionedInternalField().name() != "DU_0"
     && dimensionedInternalField().name() != "DU_0_0"
     && dimensionedInternalField().name() != "DU_0_0_0"
    )
    {
        FatalError
            << "The roller boundary condition may only be used with a "
            << "DU based solver NOT with field "
            << dimensionedInternalField().name()
            << abort(FatalError);
    }

    Info<< "    initialRotationOrigin: " << initialRotationOrigin_ << endl;
}


rollerFvPatchVectorField::
rollerFvPatchVectorField
(
    const rollerFvPatchVectorField& ptf
)
:
    fixedDisplacementFvPatchVectorField(ptf),
    displacement_(ptf.displacement_),
    dispRampTime_(ptf.dispRampTime_),
    rpm_(ptf.rpm_),
    currentRpm_(ptf.currentRpm_),
    rotRampTime_(ptf.rotRampTime_),
    rotationAxis_(ptf.rotationAxis_),
    initialRotationOrigin_(ptf.initialRotationOrigin_),
    currentAxisDisplacement_(ptf.currentAxisDisplacement_),
    origFaceCentresPtr_(NULL),
    curAngle_(ptf.curAngle_),
    curTimeIndex_(ptf.curTimeIndex_)
{
    if (ptf.origFaceCentresPtr_)
    {
        origFaceCentresPtr_ = new vectorField(*ptf.origFaceCentresPtr_);
    }
}


rollerFvPatchVectorField::
rollerFvPatchVectorField
(
    const rollerFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(ptf, iF),
    displacement_(ptf.displacement_),
    dispRampTime_(ptf.dispRampTime_),
    rpm_(ptf.rpm_),
    currentRpm_(ptf.currentRpm_),
    rotRampTime_(ptf.rotRampTime_),
    rotationAxis_(ptf.rotationAxis_),
    initialRotationOrigin_(ptf.initialRotationOrigin_),
    currentAxisDisplacement_(ptf.currentAxisDisplacement_),
    origFaceCentresPtr_(NULL),
    curAngle_(ptf.curAngle_),
    curTimeIndex_(ptf.curTimeIndex_)
{
    if (ptf.origFaceCentresPtr_)
    {
        origFaceCentresPtr_ = new vectorField(*ptf.origFaceCentresPtr_);
    }
}


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::rollerFvPatchVectorField::~rollerFvPatchVectorField()
{
    deleteDemandDrivenData(origFaceCentresPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void Foam::rollerFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedDisplacementFvPatchVectorField::autoMap(m);

    // Force all data to be re-created when needed
    deleteDemandDrivenData(origFaceCentresPtr_);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::rollerFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedDisplacementFvPatchVectorField::rmap(ptf, addr);
}


void rollerFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // A rolling pass consists of three phases:
    //     1) displacement the roller to the required reduction;
    //     2) angular acceleration of the roller;
    //     3) rolling at a constant RPM;
    // Note: phase 1 and 2 overlap

    // Algorithm
    // if curTime < dispRampTime
    //     linearly ramp displacement
    // else if curTime < (dispRampTime + rotRampTime)
    //     linearly ramp angular velocity from zero to specified rpm
    // else
    //     hold angular velocity constant
    // end

    bool newStep = false;
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        curTimeIndex_ = this->db().time().timeIndex();
        newStep = true;
    }

    // Check if the current pass has been completed

    // curPassTime is the time spent on the current pass
    const scalar curTime = this->db().time().value();

    // Set patch displacements
    vectorField disp(this->size(), vector::zero);

    if (curTime <= (dispRampTime_ + SMALL))
    {
        if (newStep)
        {
            Info<< "Roller: wire compression phase" << endl;
        }

        // Linearly ramp displacement as roller is compressed into the wire
        disp = displacement_*curTime/dispRampTime_;
        currentAxisDisplacement_ = displacement_*curTime/dispRampTime_;
    }
    else
    {
        disp = displacement_;
        currentAxisDisplacement_ = displacement_;
    }

    // Rotation

    // Total time rotating
    // const scalar rotTime = curTime - dispRampTime_;
    const scalar rotTime = curTime;

    // Calculate the current total angle of rotation in degrees

    // Analytically calculate current angle
    if (rotTime <= (rotRampTime_ + SMALL))
    {
        if (newStep)
        {
            Info<< "Roller: rotational acceleration phase" << endl;
        }

        // Roller acceleration
        curAngle_ = Foam::pow(rotTime, 2)*rpm_*3.0/(rotRampTime_);

        // Current RPM
        currentRpm_ = (rotTime/rotRampTime_)*rpm_;
    }
    else
    {
        if (newStep)
        {
            Info<< "Roller: constant rotation phase" << endl;
        }

        // Constant angular velocity
        curAngle_ =
            rotRampTime_*rpm_*3.0 // acceleration
            + (rotTime - rotRampTime_)*rpm_*6.0; // constant velocity

        // Current RPM
        currentRpm_ = rpm_;
    }

    //Info<< "curAngle " << curAngle_ << endl;

    const tensor rotMat = RodriguesRotation(rotationAxis_, curAngle_);

    const vectorField newFaceCentres =
        (rotMat & (origFaceCentres() - initialRotationOrigin_))
        + initialRotationOrigin_;

    // Add rotation to total motion
    disp += newFaceCentres - origFaceCentres();

    // Old  U
    const fvPatchField<vector>& U =
        patch().lookupPatchField<volVectorField, vector>("U");

    // Subtract old displacement as we are setting the displacement increment
    disp -= U;

    fvPatchField<vector>::operator==(disp);

    fixedValueFvPatchVectorField::updateCoeffs();
}


const vectorField& rollerFvPatchVectorField::origFaceCentres() const
{
    if (!origFaceCentresPtr_)
    {
        makeReferenceCentre();
    }

    return *origFaceCentresPtr_;
}


void rollerFvPatchVectorField::write(Ostream& os) const
{
    fixedDisplacementFvPatchVectorField::write(os);

    os.writeKeyword("displacement")
        << displacement_ << token::END_STATEMENT << nl;
    os.writeKeyword("dispRampTime")
        << dispRampTime_ << token::END_STATEMENT << nl;
    os.writeKeyword("rpm")
        << rpm_ << token::END_STATEMENT << nl;
    os.writeKeyword("currentRpm")
        << currentRpm_ << token::END_STATEMENT << nl;
    os.writeKeyword("rotationRampTime")
        << rotRampTime_ << token::END_STATEMENT << nl;
    os.writeKeyword("rotationAxis")
        << rotationAxis_ << token::END_STATEMENT << nl;
    os.writeKeyword("initialRotationOrigin")
        << initialRotationOrigin_ << token::END_STATEMENT << nl;
    os.writeKeyword("currentAxisDisplacement")
        << currentAxisDisplacement_ << token::END_STATEMENT << nl;
    os.writeKeyword("curAngle")
        << curAngle_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    rollerFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
