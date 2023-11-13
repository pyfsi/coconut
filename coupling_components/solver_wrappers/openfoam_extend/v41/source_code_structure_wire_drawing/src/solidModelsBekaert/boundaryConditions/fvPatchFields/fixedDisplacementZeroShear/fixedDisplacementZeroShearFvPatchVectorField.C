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

#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    constantTotalDisp_(vector::zero),
    dispSeries_(),
    useCurrentNormals_(false),
    incremental_(iF.name() == "DU")
{}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidDirectionMixedFvPatchVectorField(ptf, p, iF, mapper),
    constantTotalDisp_(ptf.constantTotalDisp_),
    dispSeries_(ptf.dispSeries_),
    useCurrentNormals_(ptf.useCurrentNormals_),
    incremental_(ptf.incremental_)
{}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    constantTotalDisp_(gAverage(*this)),
    dispSeries_(),
    useCurrentNormals_(dict.lookupOrDefault<Switch>("useCurrentNormals", false)),
    incremental_(dimensionedInternalField().name() == "DU")
{
    // Check if displacement is time-varying
    if (dict.found("displacementSeries"))
    {
        if (debug)
        {
            Info<< "    displacement is time-varying" << endl;
        }

        dispSeries_ =
            interpolationTable<vector>(dict.subDict("displacementSeries"));

        // Not correct for incremental solver. Disable, as not needed here
        // refValue() = dispSeries_(this->db().time().timeOutputValue());
    }
    else
    {
        // We force this because otherwise the displacement field may be
        // spatially varying, which we do not deal with
        FatalErrorIn(type() + "::" + type() + "(...)")
            << "Currently the displacementSeries must be used!"
            << abort(FatalError);
    }

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

    // Check if useCurrentNormals is active
    if (useCurrentNormals_)
    {
        if (debug)
        {
            Info<< "Using deformed patch normals to set zero shear" << endl;
        }
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


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(ptf, iF),
    constantTotalDisp_(ptf.constantTotalDisp_),
    dispSeries_(ptf.dispSeries_),
    useCurrentNormals_(ptf.useCurrentNormals_),
    incremental_(ptf.incremental_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedDisplacementZeroShearFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedDisplacementZeroShearFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(ptf, addr);
}


void fixedDisplacementZeroShearFvPatchVectorField::setFaceDisplacementField
(
    vectorField& disp
) const
{
    if (dispSeries_.size())
    {
        disp = dispSeries_(db().time().timeOutputValue());

        if (incremental_)
        {
            const vector dispOld =
                dispSeries_
                (
                    db().time().timeOutputValue() - db().time().deltaTValue()
                );

            disp -= dispOld;
        }
    }
    else
    {
        FatalErrorIn(type() + "::setFaceDisplacementField(...)")
            << "Not implemented when dispSeries is not defined!"
            << abort(FatalError);
    }
}


void fixedDisplacementZeroShearFvPatchVectorField::updateCoeffs()
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
            incremental_
        );

    if (useCurrentNormals_)
    {
        // Calculate current normals
        const fvPatchField<tensor>& relFinv =
            patch().lookupPatchField<volTensorField, tensor>
            (
                "relFinv"
            );

        const fvPatchField<scalar>& J =
            patch().lookupPatchField<volScalarField, scalar>("J");

        const vectorField nCurrent = J*relFinv.T() & patch().nf();

        // Fix deformed normal directions
        valueFraction() = sqr(nCurrent);
    }

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void fixedDisplacementZeroShearFvPatchVectorField::write(Ostream& os) const
{
    if (dispSeries_.size())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        dispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    solidDirectionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementZeroShearFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
