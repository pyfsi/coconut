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

#include "fixedDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    constantTotalDisp_(vector::zero),
    dispSeries_()
{}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    constantTotalDisp_(ptf.constantTotalDisp_),
    dispSeries_(ptf.dispSeries_)
{}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    constantTotalDisp_(gAverage(*this)),
    dispSeries_()
{
    // Check displacement is not spatially varying
    // Hmnn: during full remeshing, it is not clear how to automap the patch
    // fields such as totalDisp, so instead we will use a scalar field for the
    // totalDisp, but this means a user cannot specify a spatially varying
    // initial field. There  may be a nicer way to do this.
    // if (max(mag(*this - constantTotalDisp_)) > SMALL)
    // {
    //     FatalErrorIn(type() + "::" + type() + "(...)")
    //         << "A spatially varying displacement is not implemented!"
    //         << abort(FatalError);
    // }

    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        if (debug)
        {
            Info<< "    displacement is time-varying" << endl;
        }

        dispSeries_ =
            interpolationTable<vector>(dict.subDict("displacementSeries"));

        fvPatchField<vector>::operator==
        (
            dispSeries_(db().time().timeOutputValue())
        );
    }
}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    constantTotalDisp_(pivpvf.constantTotalDisp_),
    dispSeries_(pivpvf.dispSeries_)
{}


fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    constantTotalDisp_(pivpvf.constantTotalDisp_),
    dispSeries_(pivpvf.dispSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedDisplacementFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}


void fixedDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField disp(patch().size(), constantTotalDisp_);

    if (dispSeries_.size())
    {
        disp = dispSeries_(this->db().time().timeOutputValue());
    }

    bool incremental = bool(dimensionedInternalField().name() == "DU");

    if (incremental)
    {
        if (patch().boundaryMesh().mesh().foundObject<volVectorField>("U_0"))
        {
            const fvPatchField<vector>& Uold =
              patch().lookupPatchField<volVectorField, vector>("U_0");

            disp -= Uold;
        }
        else if (patch().boundaryMesh().mesh().foundObject<volVectorField>("U"))
        {
            const fvPatchField<vector>& Uold =
              patch().lookupPatchField<volVectorField, vector>("U");

            disp -= Uold;
        }
    }

    fvPatchField<vector>::operator==
    (
        disp
    );

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<vector> >
fixedDisplacementFvPatchVectorField::snGrad() const
{
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

        // Return snGrad with non-orthogonal correction
        return
        (
            *this
          - (patchInternalField() + (k & gradField.patchInternalField()))
        )*patch().deltaCoeffs();
    }
    else
    {
        // No non-orthogonal correction
        return (*this - patchInternalField())*patch().deltaCoeffs();
    }
}

tmp<Field<vector> >
fixedDisplacementFvPatchVectorField::gradientBoundaryCoeffs() const
{
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

    return
    (
        patch().deltaCoeffs()
       *(*this - (k & gradField.patchInternalField()))
    );
}

void fixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    if (dispSeries_.size())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        dispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
