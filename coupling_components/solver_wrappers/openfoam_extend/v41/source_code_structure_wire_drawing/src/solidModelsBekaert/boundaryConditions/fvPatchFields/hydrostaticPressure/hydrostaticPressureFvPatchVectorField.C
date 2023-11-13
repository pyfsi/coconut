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

#include "hydrostaticPressureFvPatchVectorField.H"
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

hydrostaticPressureFvPatchVectorField::
hydrostaticPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    rho_(1000.0),
    heightDirection_(0, 1, 0),
    referencePosition_(0.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


hydrostaticPressureFvPatchVectorField::
hydrostaticPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF),
    rho_(dict.lookupOrDefault<scalar>("rho", 1000)),
    heightDirection_
    (
        dict.lookupOrDefault<vector>("heightDirection", vector(0, 1, 0))
    ),
    referencePosition_(dict.lookupOrDefault<scalar>("referencePosition", 0))
{
    if (dict.found("pressure"))
    {
        pressure() = scalarField("pressure", dict, p.size());
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }
}


hydrostaticPressureFvPatchVectorField::
hydrostaticPressureFvPatchVectorField
(
    const hydrostaticPressureFvPatchVectorField& pvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(pvf, p, iF, mapper),
    rho_(pvf.rho_),
    heightDirection_(pvf.heightDirection_),
    referencePosition_(pvf.referencePosition_)
{}


hydrostaticPressureFvPatchVectorField::
hydrostaticPressureFvPatchVectorField
(
    const hydrostaticPressureFvPatchVectorField& pvf
)
:
    solidTractionFvPatchVectorField(pvf),
    rho_(pvf.rho_),
    heightDirection_(pvf.heightDirection_),
    referencePosition_(pvf.referencePosition_)
{}


hydrostaticPressureFvPatchVectorField::
hydrostaticPressureFvPatchVectorField
(
    const hydrostaticPressureFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(pvf, iF),
    rho_(pvf.rho_),
    heightDirection_(pvf.heightDirection_),
    referencePosition_(pvf.referencePosition_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hydrostaticPressureFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Set position field relative to the reference position
    const scalarField relativePosition =
        (heightDirection_ & patch().Cf()) - referencePosition_;

    // Set pressure field
    // Note: we multiply by -1 because all relative positions below the
    // reference position will be negative
    // And the neg function ensures that the pressure is zero for all positive
    // relative positions
    pressure() = neg(relativePosition)*(-rho_*9.81*relativePosition);

    solidTractionFvPatchVectorField::updateCoeffs();
}


void hydrostaticPressureFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);

    os.writeKeyword("rho")
        << rho_ << token::END_STATEMENT << endl;
    os.writeKeyword("heightDirection")
        << heightDirection_ << token::END_STATEMENT << endl;
    os.writeKeyword("referencePosition")
        << referencePosition_ << token::END_STATEMENT << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, hydrostaticPressureFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
