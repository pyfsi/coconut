/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "materialCouplingFvPatchVectorField.H"
#include "symmTransformField.H"
#include "harmonic.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::Field<vector>& materialCouplingFvPatchVectorField::originalPatchField() const
{
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        // Store original field for symmetric evaluation
        // Henrik Rusche, Aug/2011

        originalPatchField_ = *this;
        curTimeIndex_ = this->db().time().timeIndex();
    }

    return originalPatchField_;
}


tmp<Field<vector> > materialCouplingFvPatchVectorField::
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


tmp<scalarField> materialCouplingFvPatchVectorField::iK() const
{
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
        const vectorField CN = CP + regionCouplePatch_.delta();

        vectorField ownN = patch().nf();
        vectorField ngbN = -ownN;

        // scalarField ngbDeltaN = mag(sCf-CN);
        // scalarField ownDeltaN = mag(Cf-CP);
        scalarField ngbDeltaN = mag((sCf-CN)&ngbN);
        scalarField ownDeltaN = mag((Cf-CP)&ownN);

        scalarField ownMu =
            this->patch().patchInternalField
            (
                mu.internalField()
            );
        scalarField ngbMu =
            regionCouplePatch_.interpolate
            (
                lookupShadowPatchField<volScalarField, scalar>("mu")
               .patchInternalField()
            );
            // this->patchNeighbourField
            // (
            //     mu.internalField()
            // );

        scalarField ownLambda = 
            this->patch().patchInternalField
            (
                lambda.internalField()
            );
        scalarField ngbLambda = 
            regionCouplePatch_.interpolate
            (
                lookupShadowPatchField<volScalarField, scalar>("lambda")
               .patchInternalField()
            );
            // this->patchNeighbourField
            // (
            //     lambda.internalField()
            // );

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

materialCouplingFvPatchVectorField::materialCouplingFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    coupledFvPatchField<vector>(p, iF),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(iF.name()),
    matrixUpdateBuffer_(),
    originalPatchField_(),
    curTimeIndex_(-1)
{}


materialCouplingFvPatchVectorField::materialCouplingFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<vector>(p, iF, dict),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(dict.lookup("remoteField")),
    matrixUpdateBuffer_(),
    originalPatchField_(),
    curTimeIndex_(-1)
{
    if (!isType<regionCoupleFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "materialCouplingFvPatchVectorField::materialCouplingFvPatchVectorField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<vector, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not regionCouple type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }

    if (!dict.found("value"))
    {
        // Grab the internal value for initialisation. (?) HJ, 27/Feb/2009
        fvPatchField<vector>::operator=(this->patchInternalField()());
    }
}


materialCouplingFvPatchVectorField::materialCouplingFvPatchVectorField
(
    const materialCouplingFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<vector>(ptf, p, iF, mapper),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(p)),
    remoteFieldName_(ptf.remoteFieldName_),
    matrixUpdateBuffer_(),
    originalPatchField_(),
    curTimeIndex_(-1)
{
    if (!isType<regionCoupleFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "materialCouplingFvPatchVectorField::materialCouplingFvPatchVectorField\n"
            "(\n"
            "    const materialCouplingFvPatchVectorField& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<vector, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


materialCouplingFvPatchVectorField::materialCouplingFvPatchVectorField
(
    const materialCouplingFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    ggiLduInterfaceField(),
    coupledFvPatchField<vector>(ptf, iF),
    regionCouplePatch_(refCast<const regionCoupleFvPatch>(ptf.patch())),
    remoteFieldName_(ptf.remoteFieldName_),
    matrixUpdateBuffer_(),
    originalPatchField_(),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return shadow patch field
const materialCouplingFvPatchVectorField&
materialCouplingFvPatchVectorField::shadowPatchField() const
{
    // Lookup neighbour field

    return refCast<const materialCouplingFvPatchVectorField >
    (
        lookupShadowPatchField<volVectorField, vector>(remoteFieldName_)
    );
}


// Return neighbour field
tmp<Field<vector> > materialCouplingFvPatchVectorField::patchNeighbourField() const
{
    Field<vector> sField = shadowPatchField().patchInternalField();

    tmp<Field<vector> > tpnf
    (
         regionCouplePatch_.interpolate
         (
             shadowPatchField().patchInternalField()
         )
    );

    Field<vector>& pnf = tpnf();

    if (regionCouplePatch_.bridgeOverlap())
    {
        // Symmetry treatment used for overlap
        vectorField nHat = this->patch().nf();

        // Use mirrored neighbour field for interpolation
        // HJ, 21/Jan/2009
        Field<vector> bridgeField =
            transform(I - 2.0*sqr(nHat), this->patchInternalField());

#if FOAMEXTEND > 40
        regionCouplePatch_.setUncoveredFaces(bridgeField, pnf);
#else
        regionCouplePatch_.bridge(bridgeField, pnf);
#endif
    }

    return tpnf;
}


// Return neighbour field given internal cell data
tmp<Field<vector> > materialCouplingFvPatchVectorField::patchNeighbourField
(
    const word& name
) const
{
    // Lookup neighbour field

    return regionCouplePatch_.interpolate
    (
        lookupShadowPatchField<volVectorField, vector>(name)
       .patchInternalField()
    );

    // Note: this field is not bridged because local data does not exist
    // for named field.  HJ, 27/Sep/2011
}


void materialCouplingFvPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (debug)
    {
        Info << "In materialCouplingFvPatchVectorField::initEvaluate() on "
            << this->dimensionedInternalField().name()
            << " in " << this->patch().boundaryMesh().mesh().name()
            << " " << this->updated() << endl;
    }

    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Interpolation must happen at init

    // const ggiFvPatch& ggiPatch =
    //     refCast<const ggiFvPatch>(this->patch());

    if (regionCouplePatch_.master())
    {
        const fvMesh& mesh = this->patch().boundaryMesh().mesh();
        word DName = dimensionedInternalField().name();

        vectorField interD = *this;

        const volVectorField& D =
            mesh.lookupObject<volVectorField>(DName);

        // Get slave side
        materialCouplingFvPatchVectorField& sD =
            const_cast<materialCouplingFvPatchVectorField&>
            (
                refCast<const materialCouplingFvPatchVectorField>
                (
                    D.boundaryField()[regionCouplePatch_.shadowIndex()]
                )
            );

        vectorField ownModifiedTraction = this->modifiedTraction();

        vectorField ngbModifiedTraction =
           -regionCouplePatch_.interpolate(sD.modifiedTraction());

        scalarField K = this->iK();

        interD += (ngbModifiedTraction - ownModifiedTraction)/K;

        Field<vector>::operator=(interD);
        
        sD.vectorField::operator=
        (
            regionCouplePatch_.shadow().interpolate(interD)
        );
    }
}


void materialCouplingFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    // No interpolation allowed

    fvPatchField<vector>::evaluate();
}


tmp<Field<vector> > materialCouplingFvPatchVectorField::snGrad() const
{
    if (regionCouplePatch_.coupled())
    {
      // Info << "coupled snGrad" << endl;
        return coupledFvPatchField<vector>::snGrad();
    }
    else
    {
      // Info << "snGrad" << endl;
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
}

void materialCouplingFvPatchVectorField::patchInterpolate
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


void materialCouplingFvPatchVectorField::patchInterpolate
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

// Initialise neighbour processor internal cell data
void materialCouplingFvPatchVectorField::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{
    if (regionCouplePatch_.coupled())
    {
        // Prepare local matrix update buffer for the remote side.
        // Note that only remote side will have access to its psiInternal
        // as they are on different regions

        // Since interpolation needs to happen on the shadow, and within the
        // init, prepare interpolation for the other side.
        matrixUpdateBuffer_ =
            this->shadowPatchField().regionCouplePatch().interpolate
            (
                this->patch().patchInternalField(psiInternal)
            );
    }
    else
    {
        FatalErrorIn
        (
            "materialCouplingFvPatchVectorField::initInterfaceMatrixUpdate"
        )   << "init matrix update called in detached state"
            << abort(FatalError);
    }
}


// Return matrix product for coupled boundary
void materialCouplingFvPatchVectorField::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction ,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{
    if (regionCouplePatch_.coupled())
    {
        // Note: interpolation involves parallel communications and needs to
        // happen during init.  This changes the use of matrix update buffer
        // compared to earlier versions
        // HJ, 28/Sep/2011
        scalarField pnf = this->shadowPatchField().matrixUpdateBuffer();

        // Multiply the field by coefficients and add into the result
        const unallocLabelList& fc = regionCouplePatch_.faceCells();

        if (switchToLhs)
        {
            forAll(fc, elemI)
            {
                result[fc[elemI]] += coeffs[elemI]*pnf[elemI];
            }
        }
        else
        {
            forAll(fc, elemI)
            {
                result[fc[elemI]] -= coeffs[elemI]*pnf[elemI];
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "materialCouplingFvPatchVectorField::updateInterfaceMatrix"
        )   << "Matrix update called in detached state"
            << abort(FatalError);
    }
}


// Write
void materialCouplingFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeKeyword("remoteField")
        << remoteFieldName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        materialCouplingFvPatchVectorField
    );
}

// ************************************************************************* //
