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
    fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField::
fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    fixedNormal_(vector::zero),
    dispTimeSeries_(),
    tracOrDispTimeSeries_()
{}


fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField::
fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
(
    const fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    traction_(ptf.traction_, mapper),
    pressure_(ptf.pressure_, mapper),
    fixedNormal_(ptf.fixedNormal_),
    dispTimeSeries_(ptf.dispTimeSeries_),
    tracOrDispTimeSeries_(ptf.tracOrDispTimeSeries_)
{}


fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField::
fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    traction_("traction", dict, p.size()),
    pressure_("pressure", dict, p.size()),
    fixedNormal_(dict.lookup("displacementFixedNormal")),
    dispTimeSeries_(dict.subDict("timeVaryingDisplacement")),
    tracOrDispTimeSeries_(dict.subDict("displacementOrTraction"))
{
    Info<< "fixedDisplacementZeroShearOrSolidTraction boundary condition"
        << endl;

    // Fix normal direction
    fixedNormal_ /= mag(fixedNormal_);

    this->refValue() = vector::zero;
    this->refGrad() = vector::zero;
    this->valueFraction() = symmTensor(1,0,0,1,0,1);

    Field<vector>::operator=(vector::zero);
}


fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField::
fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
(
    const fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    fieldName_(ptf.fieldName_),
    traction_(ptf.traction_),
    pressure_(ptf.pressure_),
    fixedNormal_(ptf.fixedNormal_),
    dispTimeSeries_(ptf.dispTimeSeries_),
    tracOrDispTimeSeries_(ptf.tracOrDispTimeSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);
}


void fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (mag(tracOrDispTimeSeries_(this->db().time().timeOutputValue())) < SMALL)
    {
        // traction boundary

        // set valueFraction to zero
        this->valueFraction() = symmTensor::zero;

        // set gradient to enfore specified traction
        refGrad() =
            tractionBoundaryGradient().snGrad
            (
                traction_,
                pressure_,
                fieldName_,
                "U",
                patch(),
                bool(fieldName_ == "DU")
            );
    }
    else
    {
        // fixed displacement zero shear

        // set valueFraction to fix normal
        this->valueFraction() = sqr(fixedNormal_);

        // force zero shear stresses
        refGrad() =
            tractionBoundaryGradient().snGrad
            (
                vectorField(traction_.size(), vector::zero),
                scalarField(traction_.size(), 0.0),
                fieldName_,
                "U",
                patch(),
                bool(fieldName_ == "DU")
            );

        // set displacement
        vectorField disp
        (
            patch().size(),
            dispTimeSeries_(this->db().time().timeOutputValue())
        );

        if (fieldName_ == "DU")
        {
            const fvPatchField<vector>& U =
                patch().lookupPatchField<volVectorField, vector>("U");
            disp -= U;
        }
        else if (fieldName_ != "U")
        {
            FatalError
                << "The displacement field should be U or DU"
                << abort(FatalError);
        }

        refValue() = disp;
    }

    directionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField::write
(
    Ostream& os
) const
{
    directionMixedFvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    pressure_.writeEntry("pressure", os);

    os.writeKeyword("timeVaryingDisplacement") << nl;
    os << token::BEGIN_BLOCK << nl;
    dispTimeSeries_.write(os);
    os << token::END_BLOCK << nl;

    os.writeKeyword("displacementOrTraction") << nl;
    os << token::BEGIN_BLOCK << nl;
    tracOrDispTimeSeries_.write(os);
    os << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementZeroShearOrSolidTractionFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
