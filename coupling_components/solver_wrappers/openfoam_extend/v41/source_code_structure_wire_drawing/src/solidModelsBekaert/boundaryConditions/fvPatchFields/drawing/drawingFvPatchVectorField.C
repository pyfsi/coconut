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

#include "foamTime.H"
#include "drawingFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

drawingFvPatchVectorField::
drawingFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    velocity_(0.0),
    velocityRampTime_(0.0)
{}


drawingFvPatchVectorField::
drawingFvPatchVectorField
(
    const drawingFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(ptf, p, iF, mapper),
    velocity_(ptf.velocity_),
    velocityRampTime_(ptf.velocityRampTime_)
{}


drawingFvPatchVectorField::
drawingFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    velocity_(readScalar(dict.lookup("velocity"))),
    velocityRampTime_(dict.lookupOrDefault<scalar>("velocityRampTime", 0.0))
{
    // Initialise refValue: this is updated by updateCoeffs()
    if (dict.found("refValue"))
    {
        refValue() = vectorField("refValue", dict, p.size());
    }
    else
    {
        refValue() = vectorField("value", dict, p.size());
    }

    // Initialise refGrad: this is updated by updateCoeffs()
    if (dict.found("refGradient"))
    {
        refGrad() = vectorField("refGradient", dict, p.size());
    }
    else
    {
        refGrad() = vector::zero;
    }

    // Do not set the value fraction field for old fields i.e. fields with a
    // name that ends in '_0'
    const word& fieldName = dimensionedInternalField().name();
    const int fieldNameSize = fieldName.size();
    bool oldField = false;
    if (fieldNameSize > 2)
    {
        if
        (
            fieldName[fieldNameSize - 2] == '_'
         && fieldName[fieldNameSize - 1] == '0'
        )
        {
            oldField = true;
        }
    }

    if (oldField)
    {
        // The reason we do this is to avoid the unnecesary hassle of updating
        // old field patches after topological changes
        valueFraction() = symmTensor::zero;
    }
    else
    {
        if (dict.found("valueFraction"))
        {
            valueFraction() = symmTensorField("valueFraction", dict, p.size());
        }
        else
        {
            valueFraction() = sqr(patch().nf());
        }
    }

    // Update the value on the patch
    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        // There is no non-orthogonal correction here, but updateCoeffs will update this
        Field<vector> normalValue = transform(valueFraction(), refValue());
        Field<vector> gradValue =
            patchInternalField() + refGrad()/patch().deltaCoeffs();
        Field<vector> transformGradValue =
            transform(I - valueFraction(), gradValue);

        Field<vector>::operator=(normalValue + transformGradValue);
    }
}


drawingFvPatchVectorField::
drawingFvPatchVectorField
(
    const drawingFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(ptf, iF),
    velocity_(ptf.velocity_),
    velocityRampTime_(ptf.velocityRampTime_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void drawingFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void drawingFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(ptf, addr);
}


void drawingFvPatchVectorField::setFaceDisplacementField
(
    vectorField& disp
) const
{
    // Use specified velocity to calculate the displacement increment
    scalar scaleFactor = 1.0;
    if (db().time().value() < velocityRampTime_)
    {
        scaleFactor = db().time().value()/velocityRampTime_;
    }

    disp = scaleFactor*velocity_*db().time().deltaTValue()*patch().nf();
}


void drawingFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Set displacement
    setFaceDisplacementField(refValue());

    // Set gradient to force zero shear traction
    refGrad() =
        tractionBoundaryGradient().snGrad
        (
            vectorField(patch().size(), vector::zero),
            scalarField(patch().size(), 0.0),
            dimensionedInternalField().name(),
            "U",
            patch(),
            false // incremental
        );

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


scalar drawingFvPatchVectorField::averageVelocity() const
{
    return velocity_;
}


// Write
void drawingFvPatchVectorField::write(Ostream& os) const
{
    os.writeKeyword("velocity")
        << velocity_ << token::END_STATEMENT << nl;
    os.writeKeyword("velocityRampTime")
        << velocityRampTime_ << token::END_STATEMENT << nl;

    solidDirectionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    drawingFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
