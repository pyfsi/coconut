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

Class
    freeRollerNewFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "freeRollerNewFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "mathematicalConstants.H"
#include "RodriguesRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

label freeRollerNewFvPatchVectorField::zoneID() const
{
    const label zoneID =
        patch().boundaryMesh().mesh().cellZones().findZoneID(zoneName_);

    if (zoneID == -1)
    {
        FatalErrorIn
        (
            "const label freeRollerNewFvPatchVectorField::zoneID() const"
        )   << "zone: " << zoneName_ << " not found!"
            << abort(FatalError);
    }

    return zoneID;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freeRollerNewFvPatchVectorField::
freeRollerNewFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    //relaxFac_(1.0),
    displacement_(vector::zero),
    dispRampTime_(0.0),
    rotationAxis_(1, 0, 0),
    accelerationTime_(0.0),
    initialRotationOrigin_(vector::zero),
    currentAxisDisplacement_(vector::zero),
    totalAxisDisplacementReached_(false),
    zoneName_("undefined"),
    rpm_(0.0),
    curTimeIndex_(-1)
{}


freeRollerNewFvPatchVectorField::
freeRollerNewFvPatchVectorField
(
    const freeRollerNewFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    //relaxFac_(ptf.relaxFac_),
    displacement_(ptf.displacement_),
    dispRampTime_(ptf.dispRampTime_),
    rotationAxis_(ptf.rotationAxis_),
    accelerationTime_(ptf.accelerationTime_),
    initialRotationOrigin_(ptf.initialRotationOrigin_),
    currentAxisDisplacement_(ptf.currentAxisDisplacement_),
    totalAxisDisplacementReached_(ptf.totalAxisDisplacementReached_),
    zoneName_(ptf.zoneName_),
    rpm_(ptf.rpm_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


freeRollerNewFvPatchVectorField::
freeRollerNewFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    //relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 0.1)),
    displacement_(dict.lookup("displacement")),
    dispRampTime_(readScalar(dict.lookup("dispRampTime"))),
    rotationAxis_(dict.lookup("rotationAxis")),
    accelerationTime_(readScalar(dict.lookup("accelerationTime"))),
    initialRotationOrigin_(vector::zero),
    currentAxisDisplacement_(vector::zero),
    totalAxisDisplacementReached_
    (
        dict.lookupOrDefault<bool>("totalAxisDisplacementReached", false)
    ),
    zoneName_(dict.lookup("zoneName")),
    rpm_(readScalar(dict.lookup("rpm"))),
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

    // Reac current axis displacement
    if (dict.found("currentAxisDisplacement"))
    {
        Info<< "Reading currentAxisDisplacement" << endl;
        currentAxisDisplacement_ = dict.lookup("currentAxisDisplacement");
    }

    if
    (
        fieldName_ != "DU" && fieldName_ != "DU_0" && fieldName_ != "DU_0_0"
     && fieldName_ != "DU_0_0_0"
    )
    {
        FatalErrorIn("freeRollerNewFvPatchVectorField")
            << "The freeRoller boundary condition may only be used with a "
            << "DU based solver NOT with field " << fieldName_
            << abort(FatalError);
    }

    if (dispRampTime_ > accelerationTime_)
    {
        FatalErrorIn("freeRollerNewFvPatchVectorField")
            << "The accelerationTime should be greater than the dispRampTime!"
            << abort(FatalError);
    }

    if (dict.found("value"))
    {
        refValue() = vectorField("value", dict, p.size());
    }
    else
    {
        FatalErrorIn("freeRollerNewFvPatchVectorField")
            << "value entry not found for patch " << patch().name()
            << abort(FatalError);
    }

    // Normalise the rotationAxis
    if (mag(rotationAxis_) < SMALL)
    {
        FatalErrorIn("freeRollerNewFvPatchVectorField")
            << "The magnitude of the rotationAxis must not be zero"
            << abort(FatalError);
    }
    else
    {
        rotationAxis_ /= mag(rotationAxis_);
    }

    this->refGrad() = vector::zero;

    this->valueFraction() = sqr(patch().nf());

    Field<vector> normalValue = transform(valueFraction(), refValue());

    Field<vector> gradValue =
        this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

    Field<vector> transformGradValue =
        transform(I - valueFraction(), gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);
}


freeRollerNewFvPatchVectorField::
freeRollerNewFvPatchVectorField
(
    const freeRollerNewFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(ptf, iF),
    fieldName_(ptf.fieldName_),
    //relaxFac_(ptf.relaxFac_),
    displacement_(ptf.displacement_),
    dispRampTime_(ptf.dispRampTime_),
    rotationAxis_(ptf.rotationAxis_),
    accelerationTime_(ptf.accelerationTime_),
    initialRotationOrigin_(ptf.initialRotationOrigin_),
    currentAxisDisplacement_(ptf.currentAxisDisplacement_),
    totalAxisDisplacementReached_(ptf.totalAxisDisplacementReached_),
    zoneName_(ptf.zoneName_),
    rpm_(ptf.rpm_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::freeRollerNewFvPatchVectorField::~freeRollerNewFvPatchVectorField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void freeRollerNewFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void freeRollerNewFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(ptf, addr);

    //const freeRollerNewFvPatchVectorField& dmptf =
    //    refCast<const freeRollerNewFvPatchVectorField>(ptf);

    //totalDisp_.rmap(dmptf.totalDisp_, addr);
}


void freeRollerNewFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    bool newStep = false;

    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();
        newStep = true;

        if (currentAxisDisplacement_ == displacement_)
        {
            totalAxisDisplacementReached_ = true;
        }
    }

    // curPassTime is the time spent on the current pass
    const scalar curTime = db().time().value();

    // Set patch displacements
    if (curTime <= (dispRampTime_ - SMALL))
    {
        if (newStep)
        {
            Info<< type() << " : wire compression phase" << endl;

            // Old  U
            const fvPatchField<vector>& Uold =
                patch().lookupPatchField<volVectorField, vector>("U_0");

            // Linearly ramp displacement as freeRoller is compressed into the
            // wire
            currentAxisDisplacement_ = displacement_*curTime/dispRampTime_;
            refValue() = currentAxisDisplacement_ - Uold;

            // Set patch to fixedValue
            valueFraction() = I;
        }
    }
    else
    {
        if (newStep)
        {
            Info<< type() << " : free rolling phase" << endl;
        }

        nonLinearGeometry::nonLinearType nonLinear =
            nonLinearGeometry::nonLinearNames_.read
            (
                patch().boundaryMesh().mesh().solutionDict().subDict
                (
                    "solidMechanics"
                ).lookup("nonLinear")
            );

        const vectorField& n = patch().patch().faceNormals();

        // Set normal displacement increment to zero
        //refValue() = vector::zero;

        // Old  U
        const fvPatchField<vector>& Uold =
            patch().lookupPatchField<volVectorField, vector>("U_0");

        // Linearly ramp displacement as freeRoller is compressed into the
        // wire
        currentAxisDisplacement_ = displacement_;

        if (totalAxisDisplacementReached_)
        {
            refValue() = vector::zero;
        }
        else
        {
            refValue() = displacement_ - Uold;
        }

        // Lookup previous gradient
        const fvPatchTensorField& gradField =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(" + fieldName_ + ")"
            );

        // Lookup mechanical properties
        const fvPatchScalarField& mu =
            patch().lookupPatchField<volScalarField, scalar>("mu");

        const fvPatchScalarField& lambda =
            patch().lookupPatchField<volScalarField, scalar>("lambda");

        if (nonLinear == nonLinearGeometry::OFF)
        {
            // This option is required if we use a linear predictor step in the
            // solver

            // Fix old normal and rotationAxis directions
            // n ^ rotationAxis is the circumferential direction
            //valueFraction() = sqr(n);
            valueFraction() = I - sqr(n ^ rotationAxis_);

            // Set gradient to force zero shear traction
            refGrad() =
                (
                  - (n & (mu*gradField.T() - (mu + lambda)*gradField))
                  - n*lambda*tr(gradField)
                )/(2*mu + lambda);
        }
        else
        {
            // Calculate deformed unit normals

            const fvPatchField<tensor>& relFinv =
                patch().lookupPatchField<volTensorField, tensor>
                (
                    "relFinv"
                );

            const fvPatchField<scalar>& J =
                patch().lookupPatchField<volScalarField, scalar>("J");

            const vectorField nCurrent = J*relFinv.T() & n;

            // Average normal over the current time-step is the average of the
            // old the current: trapezoidal rule
            vectorField nAverage = 0.5*(n + nCurrent);
            nAverage /= mag(nAverage);

            // Set valueFraction to fix average normal direction
            //valueFraction() = sqr(nAverage);
            //valueFraction() =
            //    relaxFac_*sqr(nAverage) + (1.0 - relaxFac_)*valueFraction();
            //valueFraction() = sqr(nCurrent);
            // Fix normal and rotationAxis directions
            // n ^ rotationAxis is the circumferential direction
            //valueFraction() = I - sqr(nCurrent ^ rotationAxis_);
            valueFraction() = I - sqr(nAverage ^ rotationAxis_);

            // Set gradient to force zero shear traction in the nAverage
            // direction

            const fvPatchField<symmTensor>& tau =
                patch().lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "tauKirchhoff"
                );

            refGrad() =
                (
                    //-(nAverage & tau/J) + (2.0*mu + lambda)*(n & gradField)
                    -(nCurrent & tau/J) + (2.0*mu + lambda)*(n & gradField)
                )/(2*mu + lambda);
        }

        // Apply centrifugal force to the roller cells: in this way, we
        // artifically accelerate to the approximate expect RPM
        // If we didn't do this then it would take a long time for the free
        // roller to get up to speed
        if
        (
            newStep
         && (
                db().time().value()
             >= (accelerationTime_ - 0.5*db().time().deltaTValue() - SMALL)
             && db().time().value()
              < (accelerationTime_ + 0.5*db().time().deltaTValue() + SMALL)
            )
        )
        {
            InfoIn(type() + "::updateCoeffs()")
                << "Apply artificial centrifugal force to accelerate the free "
                << "roller: " << patch().name() << nl << endl;

            // We will calculate the current rotationOrigin
            const vector rotationOrigin =
                initialRotationOrigin_ + displacement_;

            // Take a reference to the mesh for convenience
            const fvMesh& mesh = patch().boundaryMesh().mesh();

            // Lookup displacement increment field
            volVectorField& DU =
                const_cast<volVectorField&>
                (
                    mesh.lookupObject<volVectorField>("DU")
                );
            //volVectorField& DU = DUnew.oldTime();

            // Lookup total displacement field
            volVectorField& U =
                const_cast<volVectorField&>
                (
                    mesh.lookupObject<volVectorField>("U")
                );
            //volVectorField& U = Unew.oldTime();

            // Const cast old fields as they will be overwritten
            volVectorField& DUo = const_cast<volVectorField&>(DU.oldTime());
            volVectorField& DUoo = const_cast<volVectorField&>(DUo.oldTime());
            volVectorField& Uo = const_cast<volVectorField&>(U.oldTime());
            volVectorField& Uoo = const_cast<volVectorField&>(Uo.oldTime());
            volVectorField& Uooo = const_cast<volVectorField&>(Uoo.oldTime());

            vectorField& DUI = DU.internalField();
            vectorField& DUoI = DUo.internalField();
            vectorField& DUooI = DUoo.internalField();
            //vectorField& UI = U.internalField();
            vectorField& UoI = Uo.internalField();
            vectorField& UooI = Uoo.internalField();
            vectorField& UoooI = Uooo.internalField();

            // Convert RPM to degrees per time-step
            const scalar rotationAngle =
                -rpm_*(360.0/60.0)*db().time().deltaTValue();

            // Rotation tensor
            const tensor rotMat =
                RodriguesRotation(rotationAxis_, rotationAngle);

            // Cell centres
            const volVectorField C = mesh.C();
            const vectorField& CI = C.internalField();

            // The free roller cell zone index
            const label zoneID = this->zoneID();
            const labelList& zoneCellIDs = mesh.cellZones()[zoneID];

            // Set displacement increment fields based on a const angular
            // velocity
            forAll(zoneCellIDs, cI)
            {
                const label cellID = zoneCellIDs[cI];

                const vector Cnew =
                    (rotMat.T() & (CI[cellID] - rotationOrigin))
                  + rotationOrigin;
                const vector Cold =
                    (rotMat & (CI[cellID] - rotationOrigin)) + rotationOrigin;
                const vector ColdOld =
                    (rotMat & (Cold - rotationOrigin)) + rotationOrigin;

                DUI[cellID] = Cnew - CI[cellID];
                DUoI[cellID] = CI[cellID] - Cold;
                DUooI[cellID] = Cold - ColdOld;

                UoI[cellID] = vector::zero;
                UooI[cellID] = -DUoI[cellID];
                UoooI[cellID] = UooI[cellID] - DUooI[cellID];
            }

            forAll(DU.boundaryField(), patchI)
            {
                vectorField& pDU = DU.boundaryField()[patchI];
                vectorField& pDUo = DUo.boundaryField()[patchI];
                vectorField& pDUoo = DUoo.boundaryField()[patchI];

                //vectorField& pU = U.boundaryField()[patchI];
                vectorField& pUo = Uo.boundaryField()[patchI];
                vectorField& pUoo = Uoo.boundaryField()[patchI];
                vectorField& pUooo = Uooo.boundaryField()[patchI];

                const vectorField& pC = C.boundaryField()[patchI];
                const labelList& faceCells =
                    mesh.boundaryMesh()[patchI].faceCells();

                forAll(pDU, faceI)
                {
                    const label cellID = faceCells[faceI];

                    if (mesh.cellZones().whichZone(cellID) == zoneID)
                    {
                        const vector Cnew =
                            (rotMat.T() & (pC[faceI] - rotationOrigin))
                            + rotationOrigin;
                        const vector Cold =
                            (rotMat & (pC[faceI] - rotationOrigin))
                            + rotationOrigin;
                        const vector ColdOld =
                            (rotMat & (Cold - rotationOrigin))
                            + rotationOrigin;

                        pDU[faceI] = Cnew - pC[faceI];
                        pDUo[faceI] = pC[faceI] - Cold;
                        pDUoo[faceI] = Cold - ColdOld;

                        pUo[faceI] = vector::zero;
                        pUoo[faceI] = -pDUo[faceI];
                        pUooo[faceI] = pUoo[faceI] - pDUoo[faceI];
                    }
                }
            }
        }

        // No need to call correctBoundaryConditions as we also set coupled
        // patches above
    }

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void freeRollerNewFvPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("displacement")
        << displacement_ << token::END_STATEMENT << nl;
    os.writeKeyword("dispRampTime")
        << dispRampTime_ << token::END_STATEMENT << nl;
    os.writeKeyword("rotationAxis")
        << rotationAxis_ << token::END_STATEMENT << nl;
    os.writeKeyword("accelerationTime")
        << accelerationTime_ << token::END_STATEMENT << nl;
    os.writeKeyword("initialRotationOrigin")
        << initialRotationOrigin_ << token::END_STATEMENT << nl;
    os.writeKeyword("currentAxisDisplacement")
        << currentAxisDisplacement_ << token::END_STATEMENT << nl;
    os.writeKeyword("totalAxisDisplacementReached")
        << totalAxisDisplacementReached_ << token::END_STATEMENT << nl;
    os.writeKeyword("zoneName")
        << zoneName_ << token::END_STATEMENT << nl;
    os.writeKeyword("rpm")
        << rpm_ << token::END_STATEMENT << nl;

    solidDirectionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    freeRollerNewFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
