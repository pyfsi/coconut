/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fixedScalarFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedScalarFvPatchScalarField::fixedScalarFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    fieldName_(dimensionedInternalField().name())
{}


fixedScalarFvPatchScalarField::fixedScalarFvPatchScalarField
(
    const fixedScalarFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_)
{}


fixedScalarFvPatchScalarField::fixedScalarFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    fieldName_(dimensionedInternalField().name())
{}


fixedScalarFvPatchScalarField::fixedScalarFvPatchScalarField
(
    const fixedScalarFvPatchScalarField& pivpvf
)
:
    fixedValueFvPatchScalarField(pivpvf),
    fieldName_(pivpvf.fieldName_)
{}


fixedScalarFvPatchScalarField::fixedScalarFvPatchScalarField
(
    const fixedScalarFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pivpvf, iF),
    fieldName_(pivpvf.fieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::Field<scalar> > fixedScalarFvPatchScalarField::
snGrad() const
{
    //- fixedValue snGrad with no correction
    //  return (*this - patchInternalField())*this->patch().deltaCoeffs();

    // Lookup grad form solver
    const fvPatchField<vector>& gradField =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" +fieldName_ + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();

    // Correction vector
    vectorField k = delta - n*(n&delta);

    return
    (
        *this
        - (patchInternalField() + (k&gradField.patchInternalField()))
    )*this->patch().deltaCoeffs();
}

tmp<Field<scalar> > fixedScalarFvPatchScalarField::
gradientBoundaryCoeffs() const
{
    // Lookup grad from solver
    const fvPatchField<vector>& gradField =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" +fieldName_ + ")"
        );

    vectorField n = this->patch().nf();
    vectorField delta = this->patch().delta();

    // Correction vector
    vectorField k = delta - n*(n&delta);

    return this->patch().deltaCoeffs()
       *(*this - (k&gradField.patchInternalField()));
}

void fixedScalarFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchScalarField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fixedScalarFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
