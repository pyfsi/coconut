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
    fixedTangentialNormalPressureFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "fixedTangentialNormalPressureFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedTangentialNormalPressureFvPatchVectorField::
fixedTangentialNormalPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    pressure_(0.0),
    normal_(vector::zero),
    pressureTimeSeries_(),
    updateTangentialDisplacement_(false),
    updateTangentialDisplacementEndTime_(GREAT),
    updateTangentialDisplacementRelaxFactor_(1.0)
{}


fixedTangentialNormalPressureFvPatchVectorField::
fixedTangentialNormalPressureFvPatchVectorField
(
    const fixedTangentialNormalPressureFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(ptf, p, iF, mapper),
    pressure_(ptf.pressure_, mapper),
    normal_(ptf.normal_),
    pressureTimeSeries_(ptf.pressureTimeSeries_),
    updateTangentialDisplacement_(ptf.updateTangentialDisplacement_),
    updateTangentialDisplacementEndTime_
    (
        ptf.updateTangentialDisplacementEndTime_
    ),
    updateTangentialDisplacementRelaxFactor_
    (
        ptf.updateTangentialDisplacementRelaxFactor_
    )
{}


fixedTangentialNormalPressureFvPatchVectorField::
fixedTangentialNormalPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    pressure_
    (
        !dict.found("fileName")
      ? scalarField("pressure", dict, p.size())
      : scalarField(p.size(), 0.0)
    ),
    normal_(dict.lookup("normal")),
    pressureTimeSeries_
    (
        dict.found("fileName")
      ? interpolationTable<scalar>(dict)
      : interpolationTable<scalar>()
    ),
    updateTangentialDisplacement_
    (
        dict.lookupOrDefault<Switch>("updateTangentialDisplacement", false)
    ),
    updateTangentialDisplacementEndTime_
    (
        dict.lookupOrDefault<scalar>
        (
            "updateTangentialDisplacementEndTime",
            GREAT
        )
    ),
    updateTangentialDisplacementRelaxFactor_
    (
        dict.lookupOrDefault<scalar>("relaxationFactor", 0.01)
    )
{
    if (debug)
    {
        Info<< patch().name() << ": " << type() << endl;

        if (pressureTimeSeries_.size() > 0)
        {
            Info<<"    Pressure time series found" << endl;
        }

        if (updateTangentialDisplacement_)
        {
            Info<< "    updateTangentialDisplacement is active" << nl
                << "    updateTangentialDisplacementEndTime: "
                << updateTangentialDisplacementEndTime_ << nl
                << "    relaxationFactor: "
                << updateTangentialDisplacementRelaxFactor_ << endl;
        }
    }

    refGrad() = vector::zero;

    // Fixed tangential direction
    //vectorField n = patch().nf();
    //this->valueFraction() = I - sqr(n);
    valueFraction() = I - sqr(normal_);

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        Field<vector>::operator=(vectorField(p.size(), vector::zero));
    }

    refValue() = vector::zero;

    Field<vector> normalValue = transform(valueFraction(), refValue());

    Field<vector> gradValue =
        patchInternalField() + refGrad()/patch().deltaCoeffs();

    Field<vector> transformGradValue =
        transform(I - valueFraction(), gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);
}


fixedTangentialNormalPressureFvPatchVectorField::
fixedTangentialNormalPressureFvPatchVectorField
(
    const fixedTangentialNormalPressureFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(ptf, iF),
    pressure_(ptf.pressure_),
    normal_(ptf.normal_),
    pressureTimeSeries_(ptf.pressureTimeSeries_),
    updateTangentialDisplacement_(ptf.updateTangentialDisplacement_),
    updateTangentialDisplacementEndTime_
    (
        ptf.updateTangentialDisplacementEndTime_
    ),
    updateTangentialDisplacementRelaxFactor_
    (
        ptf.updateTangentialDisplacementRelaxFactor_
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedTangentialNormalPressureFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorField::autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedTangentialNormalPressureFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(ptf, addr);

    const fixedTangentialNormalPressureFvPatchVectorField& dmptf =
        refCast<const fixedTangentialNormalPressureFvPatchVectorField>(ptf);

    pressure_.rmap(dmptf.pressure_, addr);
}


void fixedTangentialNormalPressureFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Update pressure if there is a time series
    if (pressureTimeSeries_.size() > 0)
    {
        pressure_ = pressureTimeSeries_(db().time().timeOutputValue());
    }

    const word& fieldName = dimensionedInternalField().name();

    refGrad() =
        tractionBoundaryGradient().snGrad
        (
            vectorField(patch().size(), vector::zero), // surface traction
            pressure_,                  // surface pressure
            fieldName,                  // working field name
            "U",                        // total field name
            patch(),                    // polyPatch
            bool(fieldName == "DU")     // incremental
        );

    // Update tangential displacement if specified
    if (updateTangentialDisplacement_)
    {
        if (db().time().value() >= updateTangentialDisplacementEndTime_)
        {
            Info<< nl << type() << " : disabling updateTangentialDisplacement"
                << nl << endl;
            updateTangentialDisplacement_ = false;
        }
        else
        {
            // Calculate the average tangential displacement of the patch
            // internal field
            const vector UTangAvPif =
                (I - sqr(normal_)) & gAverage(patchInternalField());

            // Update the tangential displacement using under-relaxation
            // Note: the normal component is set to zero
            refValue() =
                updateTangentialDisplacementRelaxFactor_*UTangAvPif
              + (1.0 - updateTangentialDisplacementRelaxFactor_)
               *((I - sqr(normal_)) & refValue());
        }
    }

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void fixedTangentialNormalPressureFvPatchVectorField::write(Ostream& os) const
{
    if (pressureTimeSeries_.size() > 0)
    {
        pressureTimeSeries_.write(os);
    }
    else
    {
        pressure_.writeEntry("pressure", os);
    }
    os.writeKeyword("normal")
        << normal_ << token::END_STATEMENT << nl;
    os.writeKeyword("updateTangentialDisplacement")
        << updateTangentialDisplacement_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxationFactor")
        << updateTangentialDisplacementRelaxFactor_ << token::END_STATEMENT
        << nl;

    solidDirectionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedTangentialNormalPressureFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
