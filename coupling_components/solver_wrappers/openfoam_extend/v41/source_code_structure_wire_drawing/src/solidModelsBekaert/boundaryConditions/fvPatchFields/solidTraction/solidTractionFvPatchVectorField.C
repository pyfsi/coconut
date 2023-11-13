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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "solidTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"
#include "fvc.H"
#include "surfaceInterpolationScheme.H"
#include "fvBlockMatrix.H"
#include "IOReferencer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_()
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_()
{
    if (debug)
    {
        Info<< patch().name() << ": " << type() << endl;
    }

    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }

    // Check if traction is time-varying
    if (dict.found("tractionSeries"))
    {
        if (debug)
        {
            Info<< "traction is time-varying" << endl;
        }

        tractionSeries_ =
            interpolationTable<vector>(dict.subDict("tractionSeries"));
    }
    else
    {
        traction_ = vectorField("traction", dict, p.size());
    }

    // Check if pressure is time-varying
    if (dict.found("pressureSeries"))
    {
        if (debug)
        {
            Info<< "pressure is time-varying" << endl;
        }

        pressureSeries_ =
            interpolationTable<scalar>(dict.subDict("pressureSeries"));
    }
    else
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }
}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
    traction_(stpvf.traction_, mapper),
    pressure_(stpvf.pressure_, mapper),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_)
{}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf
)
:
    fixedGradientFvPatchVectorField(stpvf),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_)
{}


solidTractionFvPatchVectorField::
solidTractionFvPatchVectorField
(
    const solidTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(stpvf, iF),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const solidTractionFvPatchVectorField& dmptf =
        refCast<const solidTractionFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void solidTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (tractionSeries_.size())
    {
        traction_ = tractionSeries_(this->db().time().timeOutputValue());
    }

    if (pressureSeries_.size())
    {
        pressure_ = pressureSeries_(this->db().time().timeOutputValue());
    }

    const word fieldName = dimensionedInternalField().name();

    gradient() =
        tractionBoundaryGradient().snGrad
        (
            traction_,                  // surface traction
            pressure_,                  // surface pressure
            fieldName,                  // working field name
            "U",                        // total field name
            patch(),                    // polyPatch
            bool(fieldName == "DU")    // incremental
        );
	
	//Info << "Gradient() " << gradient() << nl << endl;

    // If a block coupled hybrid pressure-displacement approach is used then we
    // will include the pressure implicitly
    if
    (
        patch().boundaryMesh().mesh().foundObject
        <
            IOReferencer< fvBlockMatrix<vector4> >
        >
        (
            "DUpEqn"
        )
    )
    {
        // Lookup the block matrix
        fvBlockMatrix<vector4>& DUpEqn =
            const_cast<fvBlockMatrix<vector4>&>
            (
                patch().boundaryMesh().mesh().lookupObject
                <
                IOReferencer< fvBlockMatrix<vector4> >
                >
                (
                    "DUpEqn"
                )()
            );

        Field<tensor4>& diag = DUpEqn.diag().asSquare();
        Field<vector4>& source = DUpEqn.source();

        const unallocLabelList& faceCells = patch().faceCells();
        const vectorField& Sf = patch().Sf();

        // Lookup the pressure
        const volScalarField* pPtr = NULL;
        if (patch().boundaryMesh().mesh().foundObject<volScalarField>("Dp"))
        {
            pPtr =
                &patch().boundaryMesh().mesh().lookupObject<volScalarField>
                (
                    "Dp"
                );
        }
        else
        {
            pPtr =
                &patch().boundaryMesh().mesh().lookupObject<volScalarField>
                (
                    "p"
                );
        }
        const volScalarField& p = *pPtr;

        // Add coeffs to the matrix
        // To-do: we could try use the deformed normals here to further speed
        // up convergence
        forAll(faceCells, faceI)
        {
            const label cellID = faceCells[faceI];

            diag[cellID](0, 3) -= Sf[faceI][vector::X];
            diag[cellID](1, 3) -= Sf[faceI][vector::Y];
            diag[cellID](2, 3) -= Sf[faceI][vector::Z];

            source[cellID][0] -= p[cellID]*Sf[faceI][vector::X];
            source[cellID][1] -= p[cellID]*Sf[faceI][vector::Y];
            source[cellID][2] -= p[cellID]*Sf[faceI][vector::Z];
        }
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void solidTractionFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

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

        // Set value with non-orthogonal correction
        Field<vector>::operator=
        (
            patchInternalField()
          + (k & gradField.patchInternalField())
          + gradient()/patch().deltaCoeffs()
        );
    }
    else
    {
        // Set value without non-orthogonal correction
        Field<vector>::operator=
        (
            patchInternalField()
          + gradient()/patch().deltaCoeffs()
        );
    }

    fvPatchField<vector>::evaluate();
}

// Write
void solidTractionFvPatchVectorField::write(Ostream& os) const
{
    //fvPatchVectorField::write(os);
    fixedGradientFvPatchVectorField::write(os);

    if (tractionSeries_.size())
    {
        os.writeKeyword("tractionSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        tractionSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        traction_.writeEntry("traction", os);
    }

    if (pressureSeries_.size())
    {
        os.writeKeyword("pressureSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        pressureSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        pressure_.writeEntry("pressure", os);
    }

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidTractionFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
