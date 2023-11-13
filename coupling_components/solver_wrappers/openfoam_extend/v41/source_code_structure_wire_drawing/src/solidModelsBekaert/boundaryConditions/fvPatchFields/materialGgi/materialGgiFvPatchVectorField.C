/*---------------------------------------------------------------------------* \
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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "materialGgiFvPatchVectorField.H"
#include "symmTransformField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"
#include "skewCorrectionVectors.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

materialGgiFvPatchVectorField::
materialGgiFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    ggiFvPatchField<vector>(p, iF)
{}


materialGgiFvPatchVectorField::
materialGgiFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    ggiFvPatchField<vector>(p, iF, dict)
{}


materialGgiFvPatchVectorField::
materialGgiFvPatchVectorField
(
    const materialGgiFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ggiFvPatchField<vector>(ptf, p, iF, mapper)
{}


materialGgiFvPatchVectorField::
materialGgiFvPatchVectorField
(
    const materialGgiFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    ggiFvPatchField<vector>(ptf, iF)
{}


materialGgiFvPatchVectorField::
materialGgiFvPatchVectorField
(
    const materialGgiFvPatchVectorField& ptf
)
:
    ggiFvPatchField<vector>(ptf)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<Field<vector> > materialGgiFvPatchVectorField::snGrad() const
{
    tmp<vectorField> tSnGrad
    (
        new vectorField
        (
            (this->patchNeighbourField() - this->patchInternalField())
           *this->patch().deltaCoeffs()
        )
    );

    return tSnGrad;
}


tmp<Field<vector> >
materialGgiFvPatchVectorField::patchNeighbourField() const
{
    tmp<Field<vector> > tpnf
    (
        new Field<vector>(ggiFvPatchField<vector>::patchNeighbourField())
    );

    return tpnf;
}


tmp<Field<vector> > materialGgiFvPatchVectorField::
modifiedTraction() const
{
    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    word DName = dimensionedInternalField().name();

    const volVectorField& D =
        mesh.lookupObject<volVectorField>(DName);
    const vectorField& DI = D.internalField();

    const volVectorField& prevD = D.prevIter();
    const vectorField& prevDI = prevD.internalField();

    vectorField ownD = this->patch().patchInternalField(DI);
    vectorField ownPrevD = this->patch().patchInternalField(prevDI);

    vectorField ownN = this->patch().nf();

    const vectorField& Cf = this->patch().Cf();
    const vectorField CP = this->patch().Cn();

    // scalarField ownDeltaN = mag(Cf-CP);
    scalarField ownDeltaN = mag((Cf-CP) & ownN);

    const surfaceVectorField& cauchyTraction =
        mesh.lookupObject<surfaceVectorField>("cauchyTraction");

    const volScalarField& mu =
        mesh.lookupObject<volScalarField>("mu");
    const volScalarField& lambda =
        mesh.lookupObject<volScalarField>("lambda");

    scalarField ownMu  =
        this->patch().patchInternalField
        (
            mu.internalField()
        );
    scalarField ownLambda =
        this->patch().patchInternalField
        (
            lambda.internalField()
        );
    scalarField ownK = (2*ownMu + ownLambda);

    // const symmTensorField& ownSigma =
    //   sigmaf.boundaryField()[this->patch().index()];
    vectorField ownCauchyTraction =
        cauchyTraction.boundaryField()[this->patch().index()];

    tmp<vectorField> tModifiedTraction
    (
        new vectorField
        (
            ownCauchyTraction
          + ownK*ownPrevD/ownDeltaN
          - ownK*ownD/ownDeltaN
        )
    );

    return tModifiedTraction;
}


tmp<scalarField> materialGgiFvPatchVectorField::iK() const
{
    const ggiFvPatch& ggiPatch =
        refCast<const ggiFvPatch>(this->patch());

    const fvMesh& mesh =
        this->patch().boundaryMesh().mesh();

    tmp<scalarField> tiK(new scalarField(this->patch().size(), 0.0));

    // if (ggiPatch.master())
    {
        const volScalarField& mu =
            mesh.lookupObject<volScalarField>("mu");
        const volScalarField& lambda =
            mesh.lookupObject<volScalarField>("lambda");

        const vectorField& Cf = this->patch().Cf();
        const vectorField CP = this->patch().Cn();
        const vectorField sCf = Cf;
        const vectorField CN = CP + ggiPatch.delta();

        vectorField ownN = patch().nf();
        vectorField ngbN = -ownN;

        // scalarField ngbDeltaN = mag(sCf-CN);
        // scalarField ownDeltaN = mag(Cf-CP);
        scalarField ngbDeltaN = mag((sCf-CN)&ngbN);
        scalarField ownDeltaN = mag((Cf-CP)&ownN);

        scalarField ownMu  =
            this->patch().patchInternalField
            (
                mu.internalField()
            );
        scalarField ngbMu  =
            this->patchNeighbourField
            (
                mu.internalField()
            );

        scalarField ownLambda =
            this->patch().patchInternalField
            (
                lambda.internalField()
            );
        scalarField ngbLambda =
	    this->patchNeighbourField
            (
                lambda.internalField()
            );

        scalarField ownK = (2*ownMu + ownLambda);
        scalarField ngbK = (2*ngbMu + ngbLambda);

        tiK() = (ownK/ownDeltaN + ngbK/ngbDeltaN);
    }
    // else
    // {
    //     tiK() =
    // }

    return tiK;
}


tmp<tensorField> materialGgiFvPatchVectorField::localGrad() const
{
    tmp<tensorField> tGrad
    (
        new tensorField(this->patch().size(), tensor::zero)
    );

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    word DName = dimensionedInternalField().name();

    const volTensorField& gradD =
        mesh.lookupObject<volTensorField>("grad(" + DName + ')');
    tGrad() =
      gradD.boundaryField()[this->patch().index()]
     .patchInternalField();

    // tGrad() = localExtrapolatedGrad();

    vectorField ownN = this->patch().nf();

    tGrad() -= ownN*(ownN&tGrad());
    tGrad() += ownN*localSnGrad();

    return tGrad;
}


tmp<Field<vector> > materialGgiFvPatchVectorField::localSnGrad() const
{
    const vectorField& interD = *this;

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    word DName = dimensionedInternalField().name();

    const volVectorField& D =
        mesh.lookupObject<volVectorField>(DName);
    const vectorField& DI = D.internalField();

    const volTensorField& gradD =
        mesh.lookupObject<volTensorField>("grad(" + DName + ')');
    const tensorField& gradDI = gradD.internalField();

    vectorField ownD = this->patch().patchInternalField(DI);
    tensorField ownGradD = this->patch().patchInternalField(gradDI);

    vectorField ownN = this->patch().nf();

    const vectorField& Cf = this->patch().Cf();
    const vectorField CP = this->patch().Cn();

    vectorField ownCorrVec = Cf - CP;
    ownCorrVec -= ownN*(ownN & ownCorrVec);

    ownD += (ownCorrVec & ownGradD);

    scalarField ownDeltaN = mag((Cf-CP)&ownN);

    tmp<vectorField> tLocalSnGrad
    (
        new vectorField((interD - ownD)/ownDeltaN)
    );

    return tLocalSnGrad;
}


void materialGgiFvPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const ggiFvPatch& ggiPatch =
        refCast<const ggiFvPatch>(this->patch());

    if (ggiPatch.master())
    {
        const fvMesh& mesh = this->patch().boundaryMesh().mesh();
        word DName = dimensionedInternalField().name();

        vectorField interD = *this;

        const volVectorField& D =
            mesh.lookupObject<volVectorField>(DName);

        // Get slave side
        materialGgiFvPatchVectorField& sD =
            const_cast<materialGgiFvPatchVectorField&>
            (
                refCast<const materialGgiFvPatchVectorField>
                (
                    D.boundaryField()[ggiPatch.shadowIndex()]
                )
            );

        vectorField ownModifiedTraction = this->modifiedTraction();

        vectorField ngbModifiedTraction =
           -ggiPatch.interpolate(sD.modifiedTraction());

        scalarField K = this->iK();

        interD += (ngbModifiedTraction - ownModifiedTraction)/K;

        Field<vector>::operator=(interD);

        sD.vectorField::operator=(ggiPatch.shadow().interpolate(interD));
    }
}


void materialGgiFvPatchVectorField::patchInterpolate
(
    GeometricField<vector, fvsPatchField, surfaceMesh>& fField,
    const scalarField& pL
) const
{
    fField.boundaryField()[this->patch().index()] = *this;

//     fField.boundaryField()[patchI] =
//         pL*this->patchInternalField()
//       + (1 - pL)*this->patchNeighbourField();
}


void materialGgiFvPatchVectorField::patchInterpolate
(
    GeometricField<vector, fvsPatchField, surfaceMesh>& fField,
    const scalarField& pL,
    const scalarField& pY
) const
{
    fField.boundaryField()[this->patch().index()] = *this;

//     fField.boundaryField()[patchI] =
//         pL*this->patchInternalField()
//       + pY*this->patchNeighbourField();
}

void materialGgiFvPatchVectorField::write(Ostream& os) const
{
    ggiFvPatchField<vector>::write(os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    materialGgiFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
