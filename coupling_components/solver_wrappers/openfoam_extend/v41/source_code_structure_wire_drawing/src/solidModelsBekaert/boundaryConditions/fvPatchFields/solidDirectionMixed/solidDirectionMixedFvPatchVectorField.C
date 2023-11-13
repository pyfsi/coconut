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

#include "solidDirectionMixedFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidDirectionMixedFvPatchVectorField::solidDirectionMixedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF)
{}


solidDirectionMixedFvPatchVectorField::solidDirectionMixedFvPatchVectorField
(
    const solidDirectionMixedFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper)
{}


solidDirectionMixedFvPatchVectorField::solidDirectionMixedFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
  directionMixedFvPatchVectorField(p, iF, dict)
{
    Field<vector> normalValue = transform(valueFraction(), refValue());

    Field<vector> gradValue =
        this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

    Field<vector> transformGradValue =
        transform(I - valueFraction(), gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);
}

solidDirectionMixedFvPatchVectorField::solidDirectionMixedFvPatchVectorField
(
    const solidDirectionMixedFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void solidDirectionMixedFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidDirectionMixedFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);
}


void solidDirectionMixedFvPatchVectorField::updateCoeffs()
{
    directionMixedFvPatchVectorField::updateCoeffs();
}

void solidDirectionMixedFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const Field<vector> pif = patchInternalField();
    const Field<vector> normalValue = transform(valueFraction(), refValue());

    if
    (
        db().foundObject<volTensorField>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        )
    )
    {
        // Lookup the gradient field
        const fvPatchField<tensor>& gradField =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(" + dimensionedInternalField().name() + ")"
            );

        // Patch unit normals
        const vectorField n = patch().nf();

        // Patch delta vectors
        const vectorField delta = patch().delta();

        // Non-orthogonal correction vectors
        const vectorField k = (I - sqr(n)) & delta;

        const Field<vector> gradValue =
            (refGrad()/patch().deltaCoeffs())
          + pif + (k & gradField.patchInternalField());

        const Field<vector> transformGradValue =
            transform(I - (valueFraction()), gradValue);

        Field<vector>::operator=(normalValue + transformGradValue);
    }
    else
    {
        const Field<vector> gradValue = (refGrad()/patch().deltaCoeffs()) + pif;

        const Field<vector> transformGradValue =
            transform(I - (valueFraction()), gradValue);

        Field<vector>::operator=(normalValue + transformGradValue);
    }

    fvPatchField<vector>::evaluate();
}


Foam::tmp<Foam::Field<vector> >
solidDirectionMixedFvPatchVectorField::snGrad() const
{
    const Field<vector> pif = patchInternalField();
    const Field<vector> normalValue = transform(valueFraction(), refValue());

    if
    (
        db().foundObject<volTensorField>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        )
    )
    {
        // Lookup the gradient field
        const fvPatchField<tensor>& gradField =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "grad(" + dimensionedInternalField().name() + ")"
            );

        // Patch unit normals
        const vectorField n = patch().nf();

        // Patch delta vectors
        const vectorField delta = patch().delta();

        // Non-orthogonal correction vectors
        const vectorField k = (I - sqr(n)) & delta;

        const Field<vector> gradValue =
            (refGrad()/patch().deltaCoeffs())
          + pif + (k & gradField.patchInternalField());

        const Field<vector> transformGradValue =
            transform(I - (valueFraction()), gradValue);

        const Field<vector> patchValue = normalValue + transformGradValue;

        return
        (
            patchValue
          - (patchInternalField() + (k & gradField.patchInternalField()))
        )*patch().deltaCoeffs();
    }
    else
    {
        const Field<vector> gradValue = (refGrad()/patch().deltaCoeffs()) + pif;

        const Field<vector> transformGradValue =
            transform(I - (valueFraction()), gradValue);

        const Field<vector> patchValue = normalValue + transformGradValue;

        // Return snGrad without non-orthogonal correction
        return
        (
            patchValue - patchInternalField()
        )*patch().deltaCoeffs();
    }
}

// Write
void solidDirectionMixedFvPatchVectorField::write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidDirectionMixedFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
