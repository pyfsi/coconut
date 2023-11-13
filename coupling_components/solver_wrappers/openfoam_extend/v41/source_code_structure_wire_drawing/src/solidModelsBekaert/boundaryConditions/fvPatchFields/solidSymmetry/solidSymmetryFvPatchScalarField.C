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

#include "solidSymmetryFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidSymmetryFvPatchScalarField::solidSymmetryFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    symmetryFvPatchField<scalar>(p, iF)
{}


solidSymmetryFvPatchScalarField::solidSymmetryFvPatchScalarField
(
    const solidSymmetryFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    symmetryFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (!isType<symmetryFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "solidSymmetryFvPatchScalarField::"
            "solidSymmetryFvPatchScalarField\n"
            "(\n"
            "    const solidSymmetryFvPatchScalarField& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


solidSymmetryFvPatchScalarField::solidSymmetryFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    symmetryFvPatchField<scalar>(p, iF, dict)
{
    Info << "Symmetry boundary condition with non-orthogonal correction"
        << endl;

    if (!isType<symmetryFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "solidSymmetryFvPatchScalarField::"
            "solidSymmetryFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const Field<scalar>& field,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "\n    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


solidSymmetryFvPatchScalarField::solidSymmetryFvPatchScalarField
(
    const solidSymmetryFvPatchScalarField& ptf
)
:
    symmetryFvPatchField<scalar>(ptf)
{}


solidSymmetryFvPatchScalarField::solidSymmetryFvPatchScalarField
(
    const solidSymmetryFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    symmetryFvPatchField<scalar>(ptf, iF)
{}


// return gradient at boundary
tmp<Field<scalar> > solidSymmetryFvPatchScalarField::snGrad() const
{
    // Patch unit normals
    const vectorField n = patch().nf();

    if
    (
        db().foundObject<volVectorField>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        )
    )
    {
        // Lookup the gradient field
        const fvPatchField<vector>& gradField =
            patch().lookupPatchField<volVectorField, vector>
            (
                "grad(" + dimensionedInternalField().name() + ")"
            );

        // Patch delta vectors
        const vectorField delta = patch().delta();

        // Non-orthogonal correction vectors
        const vectorField k = (I - sqr(n)) & delta;

        // Patch internal field
        scalarField UP = patchInternalField();
        UP += (k & gradField.patchInternalField());

        // snGrad with non-orthogonal correction
        return
        (
            transform(I - 2.0*sqr(n), UP) - UP
        )*(patch().deltaCoeffs()/2.0);
    }
    else
    {
        // Patch internal field
        const scalarField UP = patchInternalField();

        // snGrad without non-orthogonal correction
        return
        (
            transform(I - 2.0*sqr(n), UP) - UP
        )*(patch().deltaCoeffs()/2.0);
    }
}


// Evaluate the field on the patch
void solidSymmetryFvPatchScalarField::
evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Patch unit normals
    const vectorField n = patch().nf();

    if
    (
        db().foundObject<volVectorField>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        )
    )
    {
        // Lookup the gradient field
        const fvPatchField<vector>& gradField =
            patch().lookupPatchField<volVectorField, vector>
            (
                "grad(" + dimensionedInternalField().name() + ")"
            );

        // Patch delta vectors
        const vectorField delta = patch().delta();

        // Non-orthogonal correction vectors
        const vectorField k = (I - sqr(n)) & delta;

        // Patch internal field
        scalarField UP = patchInternalField();
        UP += (k & gradField.patchInternalField());

        // Patch value with non-orthogonal correction
        Field<scalar>::operator=
        (
            (
                UP + transform(I - 2.0*sqr(n), UP)
            )/2.0
        );
    }
    else
    {
        // Patch internal field
        const scalarField UP = patchInternalField();

        // Patch value withuot non-orthogonal correction
        Field<scalar>::operator=
        (
            (
                UP + transform(I - 2.0*sqr(n), UP)
            )/2.0
        );
    }
}


// Write
void solidSymmetryFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, solidSymmetryFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
