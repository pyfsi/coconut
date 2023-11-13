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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "materialInterfaceFvPatchVectorField.H"
#include "symmTransformField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"

// #include "solidSolver.H"
// #include "unsULLSSolid.H"
// #include "pRveUnsULLSSolid.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

materialInterfaceFvPatchVectorField::
materialInterfaceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    ggiFvPatchField<vector>(p, iF),
    //cohesiveLawPtr_(NULL),
    // debonded_(p.size(), 0),
    // displacementJump_(p.size(), vector::zero),
    // traction_(p.size(), vector::zero),
    // normalTraction_(p.size(), 0),
    // initiationTraction_(p.size(), vector::zero),
    grad_(p.size(), vector::zero),
    curTimeIndex_(-1),
    iCorr_(0),
    // jumpConverged_(false),
    nIter_(p.size(), 0)
{}


materialInterfaceFvPatchVectorField::
materialInterfaceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    ggiFvPatchField<vector>(p, iF, dict),
    //cohesiveLawPtr_(NULL),
    // debonded_(p.size(), 0),
    // displacementJump_(p.size(), vector::zero),
    // traction_(p.size(), vector::zero),
    // normalTraction_(p.size(), 0),
    // initiationTraction_(p.size(), vector::zero),
    grad_(p.size(), vector::zero),
    curTimeIndex_(-1),
    iCorr_(0),
    // jumpConverged_(false),
    nIter_(p.size(), 0)
{
//     if (dict.found("debonded"))
//     {
//         debonded_ = Switch(dict.lookup("debonded"));
//     }

//     const ggiFvPatch& ggiPatch =
//         refCast<const ggiFvPatch>(this->patch());

//     if (ggiPatch.master())
    // {
    //     cohesiveLawPtr_ =
    //         simpleCohesiveLaw::New(dict.lookup("cohesiveLaw"), dict).ptr();
    // }

    // if (dict.found("displacementJump"))
    // {
    //     displacementJump_ = vectorField("displacementJump", dict, p.size());
    // }

    // if (dict.found("traction"))
    // {
    //     traction_ = vectorField("traction", dict, p.size());
    // }

    // if (dict.found("initiationTraction"))
    // {
    //     initiationTraction_ =
    //         vectorField("initiationTraction", dict, p.size());
    // }

//     Info << this->patch().name() << " " << cohesiveLawPtr_ << endl;
}


materialInterfaceFvPatchVectorField::
materialInterfaceFvPatchVectorField
(
    const materialInterfaceFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    ggiFvPatchField<vector>(ptf, p, iF, mapper),
    //cohesiveLawPtr_(ptf.cohesiveLawPtr_->clone().ptr()),
    // debonded_(ptf.debonded_, mapper),
    // displacementJump_(ptf.displacementJump_, mapper),
    // traction_(ptf.traction_, mapper),
    // normalTraction_(ptf.normalTraction_, mapper),
    // initiationTraction_(ptf.initiationTraction_, mapper),
    grad_(ptf.grad_, mapper),
    curTimeIndex_(-1),
    iCorr_(0),
    // jumpConverged_(false),
    nIter_(ptf.nIter_)
{}


materialInterfaceFvPatchVectorField::
materialInterfaceFvPatchVectorField
(
    const materialInterfaceFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    ggiFvPatchField<vector>(ptf, iF),
    //cohesiveLawPtr_(ptf.cohesiveLawPtr_->clone().ptr()),
    // debonded_(ptf.debonded_),
    // displacementJump_(ptf.displacementJump_),
    // traction_(ptf.traction_),
    // normalTraction_(ptf.normalTraction_),
    // initiationTraction_(ptf.initiationTraction_),
    grad_(ptf.grad_),
    curTimeIndex_(-1),
    iCorr_(0),
    // jumpConverged_(ptf.jumpConverged_),
    nIter_(ptf.nIter_)
{}


materialInterfaceFvPatchVectorField::
materialInterfaceFvPatchVectorField
(
    const materialInterfaceFvPatchVectorField& ptf
)
:
    ggiFvPatchField<vector>(ptf),
    //cohesiveLawPtr_(ptf.cohesiveLawPtr_->clone().ptr()),
    // debonded_(ptf.debonded_),
    // displacementJump_(ptf.displacementJump_),
    // traction_(ptf.traction_),
    // normalTraction_(ptf.normalTraction_),
    // initiationTraction_(ptf.initiationTraction_),
    grad_(ptf.grad_),
    curTimeIndex_(-1),
    iCorr_(0),
    // jumpConverged_(ptf.jumpConverged_),
    nIter_(ptf.nIter_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// const simpleCohesiveLaw& materialInterfaceFvPatchVectorField::law() const
// {
//     if (!cohesiveLawPtr_)
//     {
//         FatalErrorIn
//         (
//             "const cohesiveLaw& materialInterfaceFvPatchVectorField::law() "
//             "const"
//         )   << "Law pointer not set" << abort(FatalError);
//     }

//     return *cohesiveLawPtr_;
// }


tmp<Field<vector> > materialInterfaceFvPatchVectorField::snGrad() const
{
    const ggiFvPatch& ggiPatch =
        refCast<const ggiFvPatch>(this->patch());

    scalarField areaRatio(this->patch().size(), 1);

    if (!ggiPatch.master())
    {
        scalarField masterMagS =
            ggiPatch.interpolate
            (
                mag(ggiPatch.shadow().Sf())
            );
        scalarField slaveMagS = mag(ggiPatch.Sf());

        areaRatio = masterMagS/slaveMagS;
    }

    tmp<vectorField> tSnGrad
    (
        new vectorField
        (
            (this->patchNeighbourField() - this->patchInternalField())
           *this->patch().deltaCoeffs()*areaRatio
        )
    );

    // forAll(tSnGrad(), faceI)
    // {
    //     if (debonded_[faceI] > SMALL)
    //     {
    //         tSnGrad()[faceI] = grad_[faceI];
    //     }
    // }

    return tSnGrad;
}


tmp<Field<vector> >
materialInterfaceFvPatchVectorField::patchNeighbourField() const
{
    tmp<Field<vector> > tpnf
    (
        new Field<vector>(ggiFvPatchField<vector>::patchNeighbourField())
    );

//     Info << "materialInterfaceFvPatchVectorField::patchNeighbourField()"
//         << endl;

//     Field<vector>& pnf = tpnf();

//     const ggiFvPatch& ggiPatch =
//         refCast<const ggiFvPatch>(this->patch());

//     if (ggiPatch.master())
//     {
//         for (label faceI = 0; faceI < pnf.size(); faceI++)
//         {
//             pnf[faceI] -= displacementJump_[faceI];
//         }
//     }
//     else
//     {
//         for (label faceI = 0; faceI < pnf.size(); faceI++)
//         {
//             pnf[faceI] += displacementJump_[faceI];
//         }
//     }

    return tpnf;
}


tmp<Field<vector> > materialInterfaceFvPatchVectorField::localSnGrad() const
{
    const vectorField& interU = *this;

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();

    const volVectorField& U =
        mesh.lookupObject<volVectorField>("U");
    const vectorField& UI = U.internalField();

    const volTensorField& gradU =
        mesh.lookupObject<volTensorField>("grad(U)");
    const tensorField& gradUI = gradU.internalField();

    vectorField ownU = this->patch().patchInternalField(UI);
    tensorField ownGradU = this->patch().patchInternalField(gradUI);

    vectorField ownN = this->patch().nf();

    const vectorField& Cf = this->patch().Cf();
    const vectorField CP = this->patch().Cn();

    vectorField ownCorrVec = Cf - CP;
    ownCorrVec -= ownN*(ownN & ownCorrVec);

    ownU += (ownCorrVec & ownGradU);

    scalarField ownDeltaN = mag((Cf - CP) & ownN);

//     Info << "Using material interface snGrad: "
//         << average(ownDeltaN) << endl;

    tmp<vectorField> tLocalSnGrad
    (
        new vectorField((interU - ownU)/ownDeltaN)
    );

    // forAll(tLocalSnGrad(), faceI)
    // {
    //     if (debonded_[faceI] > SMALL)
    //     {
    //         tLocalSnGrad()[faceI] = grad_[faceI];
    //     }
    // }

    return tLocalSnGrad;
}


void materialInterfaceFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//     Info << "materialInterfaceFvPatchVectorField::updateCoeffs()" << endl;

    // const ggiFvPatch& ggiPatch =
    //   refCast<const ggiFvPatch>(this->patch());

    //const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    //word DDName = dimensionedInternalField().name();

//     const volVectorField& DD =
//         mesh.lookupObject<volVectorField>(DDName);

//     const volVectorField& D =
//         mesh.lookupObject<volVectorField>("D");

// //     const vectorField& DDI = DD.internalField();

//     const vectorField& interDD = *this;
//     vectorField shadowInterDD =
//         ggiPatch.interpolate
//         (
//             DD.boundaryField()[ggiPatch.shadowIndex()]
//         );

//     vectorField interD = D.boundaryField()[this->patch().index()];
//     vectorField shadowInterD =
//         ggiPatch.interpolate
//         (
//             D.boundaryField()[ggiPatch.shadowIndex()]
    //);

//     const surfaceScalarField& cTwoMuPlusLambda =
//         mesh.lookupObject<surfaceScalarField>("twoMuPlusLambda");

//     surfaceScalarField& twoMuPlusLambda =
//         const_cast<surfaceScalarField&>(cTwoMuPlusLambda);

//     nIter_ = 0;

    //if (curTimeIndex_ != this->db().time().timeIndex())
    //{
    //    curTimeIndex_ = this->db().time().timeIndex();
    //    iCorr_ = 0;
        //jumpConverged_ = false;

        // if (ggiPatch.master())
        // {
            //label nFaces = 0;

            //scalar avgNormalTraction = average(normalTraction_);
//             scalar maxDebonded = max(debonded_);

            //vectorField n = this->patch().nf();

//             if (maxDebonded < SMALL)
//             {
//                 if (avgNormalTraction > law().sigmaMax().value())
//                 {
//                     debonded_ = 1;
//                     initiationTraction_ = traction_;
//                     nFaces = traction_.size();

//                     displacementJump_ = 8.41638e-09*n;
//                 }
//             }

//             Info << "Displacement jump, max: " << max(n&displacementJump_)
//                 << ", avg: " << average(n&displacementJump_)
//                 << ", min: " << min(n&displacementJump_) << endl;

            // forAll(traction_, faceI)
            // {
            //     if (debonded_[faceI] < SMALL)
            //     {
            //         if (normalTraction_[faceI] > law().sigmaMax().value())
            //         {
            //             debonded_[faceI] = 1;
            //             initiationTraction_[faceI] = traction_[faceI];
            //             nFaces++;
            //         }
            //     }
            // }

            // Info << "Avg normal traction: " << avgNormalTraction << endl;
            // Info << "Number of debonded faces: " << nFaces << endl;

            // Check current relative separation distance

            // scalarField curSepDist =
            //     n
            //   & (
            //         shadowInterD
            //       - interD
            //     );

            // curSepDist /= law().deltaC().value() + SMALL;

            // Info << "Current relative separation distance ("
            //     << law().deltaC().value() << "), max: "
            //     << gMax(curSepDist) << ", avg: " << gAverage(curSepDist)
            //     << ", min: " << gMin(curSepDist) << endl;
        // }
        // else
        // {
        //     const materialInterfaceFvPatchVectorField& spDD =
        //         refCast<const materialInterfaceFvPatchVectorField>
        //         (
        //             DD.boundaryField()[ggiPatch.shadowIndex()]
        //         );

        //     debonded_ =
        //         ggiPatch.interpolate
        //         (
        //             spDD.debonded()
        //         );

        //     initiationTraction_ =
        //         ggiPatch.interpolate
        //         (
        //             -spDD.initiationTraction()
        //         );
        // }
        //}
    //else
    //{
        //if (ggiPatch.master())
        //{
            //label correctionFreq = 10;

            // Check for penetration
//             if (++iCorr_ % correctionFreq == 0)
//             {
// //                 Info << "Checking for penetration. "
// //                     << "Number of overlaped faces: ";

//                 vectorField n = this->patch().nf();

                // scalarField curSepDist =
                //     n
                //   & (
                //         (shadowInterD + shadowInterDD)
                //       - (interD + interDD)
                //     );

                // label nFaces = 0;

                // // Checking for penetration
                // forAll(curSepDist, faceI)
                // {
                //     if (debonded_[faceI] > SMALL)
                //     {
                //         if (curSepDist[faceI] < SMALL)
                //         {
                //             debonded_[faceI] = 0;
                //             nFaces++;
                //         }
                //     }
                // }

//                 Info << nFaces << endl;
            // }
            // }
        // else
        // {
            // const materialInterfaceFvPatchVectorField& spDD =
            //     refCast<const materialInterfaceFvPatchVectorField>
            //     (
            //         DD.boundaryField()[ggiPatch.shadowIndex()]
            //     );

            // debonded_ =
            //     ggiPatch.interpolate
            //     (
            //         spDD.debonded()
            //     );
        // }
        //}

    // Update normal gradient in cohezive zone
    // if (max(debonded_) > SMALL)
    // {
    //     word DDName = this->dimensionedInternalField().name();

    //     const fvsPatchField<scalar>& mu =
    //         patch().lookupPatchField<surfaceScalarField, scalar>
    //         (
    //             "muf"
    //         );

    //     const fvsPatchField<scalar>& lambda =
    //         patch().lookupPatchField<surfaceScalarField, scalar>
    //         (
    //             "lambdaf"
    //         );

    //     const fvsPatchField<tensor>& gradDD =
    //         patch().lookupPatchField<surfaceTensorField, tensor>
    //         (
    //             "grad" + DDName + "f"
    //         );

    //     const fvsPatchField<tensor>& F =
    //         patch().lookupPatchField<surfaceTensorField, tensor>
    //         (
    //             "Ff"
    //         );

    //     const symmTensorField& tau =
    //         patch().lookupPatchField<surfaceSymmTensorField, symmTensor>
    //         (
    //             "tauf"
    //         );

    //     vectorField n = this->patch().nf();

    //     scalarField J = det(F);

    //     tensorField relF = I + gradDD.T();
    //     tensorField invRelF = hinv(relF);
    //     scalarField relJ = det(relF);

    //     vectorField curN = (invRelF.T() & n);
    //     curN /= mag(curN);

    //     // Current separation distance
    //     scalarField curSepDist =
    //         n
    //       & (
    //             (shadowInterD + shadowInterDD)
    //           - (interD + interDD)
    //         );
    //     curSepDist = mag(curSepDist);

    //     // Current cohezive zone traction
    //     vectorField curCzmTraction = initiationTraction_;
    //     forAll(curCzmTraction, faceI)
    //     {
    //         curCzmTraction[faceI] *=
    //             law().traction(curSepDist[faceI])
    //            /law().sigmaMax().value();
    //     }

    //     // Normal gradient of displacement increment
    //     grad_ = curCzmTraction - (curN & tau/J)
    //       + (2.0*mu + lambda)*(n & gradDD);
    //     grad_ /= (2*mu + lambda);
    // }

    ggiFvPatchField<vector>::updateCoeffs();
}


void materialInterfaceFvPatchVectorField::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

//     Info << "materialInterfaceFvPatchVectorField::initEvaluate()" << endl;

    // Looking up solid solver
    // const solidSolver& solid =
    //     this->db().objectRegistry::lookupObject<solidSolver>
    //     (
    //         "solidProperties"
    //   );

    // tensor avgRelF = I;
    // if (isA<solidSolvers::pRveUnsULLSSolid>(solid))
    // {
    //     const solidSolvers::pRveUnsULLSSolid& pRveSolid =
    //         refCast<const solidSolvers::pRveUnsULLSSolid>(solid);

    //     avgRelF = pRveSolid.avgRelDeformationGradient();
    // }

    const ggiFvPatch& ggiPatch =
        refCast<const ggiFvPatch>(this->patch());

    const fvMesh& mesh = this->patch().boundaryMesh().mesh();
    //word DDName = dimensionedInternalField().name();

    //vectorField interDD = *this;
    vectorField interU = *this;

    if (ggiPatch.master())
    {
        // const volVectorField& DD =
        //     mesh.lookupObject<volVectorField>(DDName);
        const volVectorField& U =
            mesh.lookupObject<volVectorField>("U");
        const vectorField& UI = U.internalField();

        vectorField shadowInterU =
            ggiPatch.interpolate
            (
                U.boundaryField()[ggiPatch.shadowIndex()]
            );

        const volVectorField& prevU = U.prevIter();
        const vectorField& prevUI = prevU.internalField();

        const volTensorField& gradU =
            mesh.lookupObject<volTensorField>("grad(U)");
        const tensorField& gradUI = gradU.internalField();

        // const surfaceTensorField& gradDDf =
        //     mesh.lookupObject<surfaceTensorField>("grad" + DDName + 'f');

        const surfaceSymmTensorField& sigma =
           mesh.lookupObject<surfaceSymmTensorField>("sigmaf");
        // const volSymmTensorField& sigma =
        //     mesh.lookupObject<volSymmTensorField>("sigma");

        // const surfaceTensorField& Ff =
        //     mesh.lookupObject<surfaceTensorField>("Ff");

        const volScalarField& mu =
            mesh.lookupObject<volScalarField>("mu");
        const volScalarField& lambda =
            mesh.lookupObject<volScalarField>("lambda");

        const vectorField& Cf = this->patch().Cf();
        const vectorField CP = this->patch().Cn();
        const vectorField sCf =
            ggiPatch.interpolate
            (
                mesh.boundary()[ggiPatch.shadowIndex()].Cf()
            );
        const vectorField CN =
            ggiPatch.interpolate
            (
                mesh.boundary()[ggiPatch.shadowIndex()].Cn()
            );

        vectorField ownU = this->patch().patchInternalField(UI);
        vectorField ngbU = this->patchNeighbourField(UI);

        vectorField ownPrevU = this->patch().patchInternalField(prevUI);
        vectorField ngbPrevU = this->patchNeighbourField(prevUI);

        tensorField ownGradU = this->patch().patchInternalField(gradUI);
        tensorField ngbGradU = this->patchNeighbourField(gradUI);

        // tensorField ownGradUf =
        //     gradUf.boundaryField()[this->patch().index()];
        // tensorField ngbGradUf = this->patchNeighbourField(gradUf);

        // tensorField ownF =
        //     Ff.boundaryField()[this->patch().index()];
        // tensorField ngbF = this->patchNeighbourField(Ff);

        // scalarField ownJ = det(ownF);
        // scalarField ngbJ = det(ngbF);

        // tensorField ownRelF = avgRelF + ownGradDDf.T();
        // tensorField ngbRelF = avgRelF + ngbGradDDf.T();
//         tensorField ownRelF = I + ownGradDDf.T();
//         tensorField ngbRelF = I + ngbGradDDf.T();

        // scalarField ownRelJ = det(ownRelF);
        // scalarField ngbRelJ = det(ngbRelF);

        // tensorField invOwnRelF = hinv(ownRelF);
        // tensorField invNgbRelF = hinv(ngbRelF);

        vectorField ownN = patch().nf();
        vectorField ngbN =
            ggiPatch.interpolate
            (
                mesh.boundary()[ggiPatch.shadowIndex()].nf()
            );

        // vectorField ownCurN = (invOwnRelF.T() & ownN);
        // ownCurN /= mag(ownCurN);

        // vectorField ngbCurN = (invNgbRelF.T() & ngbN);
        // ngbCurN /= mag(ngbCurN);

        scalarField ngbDeltaN = mag((sCf - CN) & ngbN);
        scalarField ownDeltaN = mag((Cf - CP) & ownN);

        vectorField ownCorrVec = Cf - CP;
        ownCorrVec -= ownN*(ownN & ownCorrVec);
        vectorField ngbCorrVec = sCf - CN;
        ngbCorrVec -= ngbN*(ngbN & ngbCorrVec);

        ownU += (ownCorrVec & ownGradU);
        ngbU += (ngbCorrVec & ngbGradU);

        ownPrevU += (ownCorrVec & ownGradU);
        ngbPrevU += (ngbCorrVec & ngbGradU);

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

        scalarField ownLambda  =
            this->patch().patchInternalField
            (
                lambda.internalField()
            );
        scalarField ngbLambda  =
            this->patchNeighbourField
            (
                lambda.internalField()
            );

        // const symmTensorField& ownTau =
        //     tauf.boundaryField()[this->patch().index()];
        // symmTensorField ngbTau = this->patchNeighbourField(tauf);
        const symmTensorField& ownSigma =
            sigma.boundaryField()[this->patch().index()];
        symmTensorField ngbSigma = this->patchNeighbourField(sigma);

        // const vectorField& ownS = this->patch().Sf();
        // vectorField ngbS =
        //     ggiPatch.interpolate
        //     (
        //         mesh.boundary()[ggiPatch.shadowIndex()].Sf()
        //     );

        //vectorField ownCurS = (invOwnRelF.T() & ownS)*ownRelJ;
        //vectorField ngbCurS = (invNgbRelF.T() & ngbS)*ngbRelJ;

        //scalarField ownMagCurS = mag(ownCurS);
        ///scalarField ngbMagCurS = mag(ngbCurS);

        // This is Cauchy traction
        // vectorField ownCauchyTraction = (ownCurN & ownTau)/ownJ;
        // vectorField ngbCauchyTraction = -(ngbCurN & ngbTau)/ngbJ;
        vectorField ownCauchyTraction = (ownN & ownSigma);
        vectorField ngbCauchyTraction = -(ngbN & ngbSigma);

        forAll(interU, faceI)
        {
            // if (debonded_[faceI] < SMALL)
            // {
                //- Still attached face

                vector ngbModifiedTraction =
                    ngbCauchyTraction[faceI]
                  - (2*ngbMu[faceI] + ngbLambda[faceI])
                   *(ngbPrevU[faceI] - interU[faceI])
                   /ngbDeltaN[faceI];

                vector ownModifiedTraction =
                    ownCauchyTraction[faceI]
                  - (2*ownMu[faceI] + ownLambda[faceI])
                   *(interU[faceI] - ownPrevU[faceI])
                   /ownDeltaN[faceI];

                interU[faceI] =
                (
                    (2*ownMu[faceI] + ownLambda[faceI])
                   *ownU[faceI]*ngbDeltaN[faceI]
                  + (2*ngbMu[faceI] + ngbLambda[faceI])
                   *ngbU[faceI]*ownDeltaN[faceI]
                  + ownDeltaN[faceI]*ngbDeltaN[faceI]
                   *(ngbModifiedTraction - ownModifiedTraction)
                )
               /(
                    (2*ownMu[faceI] + ownLambda[faceI])
                   *ngbDeltaN[faceI]
                  + (2*ngbMu[faceI] + ngbLambda[faceI])
                   *ownDeltaN[faceI]
                );

                shadowInterU[faceI] = interU[faceI];
                //}
        }

        // displacementJump_ = shadowInterDD - interDD;

        //traction_ = ownCauchyTraction;

        //normalTraction_ = (ownCurN & ownCauchyTraction);
        //scalarField snt = -(ngbCurN & ngbCauchyTraction);

        //Field<vector>::operator=(interU);

//         Info << setprecision(12);

//         Info << "Master normal traction, max: " << max(normalTraction_)
//             << ", avg: " << average(normalTraction_)
//             << ", min: " << min(normalTraction_) << endl;

//         Info << "Slave normal traction, max: " << max(snt)
//             << ", avg: " << average(snt)
//             << ", min: " << min(snt) << endl;

//         Info << "Displacement jump, max: " << max(ownN&displacementJump_)
//             << ", avg: " << average(ownN&displacementJump_)
//             << ", min: " << min(ownN&displacementJump_) << endl;

        // Set slave side
        // materialInterfaceFvPatchVectorField& sDD =
        //     const_cast<materialInterfaceFvPatchVectorField&>
        //     (
        //         refCast<const materialInterfaceFvPatchVectorField>
        //         (
        //             DD.boundaryField()[ggiPatch.shadowIndex()]
        //         )
        //     );
        materialInterfaceFvPatchVectorField& sU =
            const_cast<materialInterfaceFvPatchVectorField&>
            (
                refCast<const materialInterfaceFvPatchVectorField>
                (
                    U.boundaryField()[ggiPatch.shadowIndex()]
                )
            );

        //sDD == (ggiPatch.shadow().interpolate(shadowInterDD));
        sU == (ggiPatch.shadow().interpolate(shadowInterU));
        //sDD.displacementJump() =
        //    ggiPatch.shadow().interpolate(displacementJump_);
    }

    // Calculate displacement for specified czm traction
    // if (max(debonded_) > SMALL)
    // {
    //     const fvPatchField<tensor>& gradDD =
    //         patch().lookupPatchField<volTensorField, tensor>
    //         (
    //             "grad(" + DDName + ")"
    //         );

    //     const vectorField n = this->patch().nf();
    //     const vectorField& Cf = this->patch().Cf();
    //     const vectorField CP = this->patch().Cn();

    //     vectorField delta = Cf - CP;
    //     vectorField k = delta - n*(n & delta);
    //     scalarField deltaN = (n & delta);

    //     vectorField newInterDD =
    //         this->patchInternalField()
    //       + (k & gradDD.patchInternalField())
    //       + grad_*deltaN;

    //     forAll(interDD, faceI)
    //     {
    //         if (debonded_[faceI] > SMALL)
    //         {
    //             interDD[faceI] = newInterDD[faceI];
    //         }
    //     }
    // }

    //Field<vector>::operator=(interDD);
    Field<vector>::operator=(interU);
}


// const Switch& materialInterfaceFvPatchVectorField::master() const
// {
// //     const fvMesh& mesh = this->patch().boundaryMesh().mesh();

//     return master_;
// }


void materialInterfaceFvPatchVectorField::manipulateMatrix
(
    fvMatrix<vector>& eqn
)
{
//     const fvMesh& mesh = this->patch().boundaryMesh().mesh();
//     const surfaceScalarField& twoMuPlusLambda =
//         mesh.lookupObject<surfaceScalarField>("twoMuPlusLambda");
//     const ggiFvPatch& ggiPatch =
//         refCast<const ggiFvPatch>(this->patch());

//     scalarField& diag = eqn.diag();



    // vectorField& source = eqn.source();

    // const scalarField& magSf = this->patch().magSf();

    // const fvsPatchField<scalar>& twoMuPlusLambda =
    //     patch().lookupPatchField<surfaceScalarField, scalar>
    //     (
    //         //"twoMuPlusLambda"
    //         "twoMuLambda"
    //     );

    //const unallocLabelList& faceCells = this->patch().faceCells();

    // forAll(faceCells, faceI)
    // {
    //     if (debonded_[faceI] > SMALL)
    //     {
    //         source[faceCells[faceI]] +=
    //             twoMuPlusLambda[faceI]*magSf[faceI]*grad_[faceI];
    //     }
    // }



    fvPatchField<vector>::manipulateMatrix(eqn);
}


void materialInterfaceFvPatchVectorField::patchInterpolate
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


void materialInterfaceFvPatchVectorField::patchInterpolate
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


void materialInterfaceFvPatchVectorField::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011

    const ggiFvPatch& ggiPatch = refCast<const ggiFvPatch>(this->patch());

    // Get shadow face-cells and assemble shadow field
    const unallocLabelList& sfc = ggiPatch.shadow().faceCells();

    scalarField sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = psiInternal[sfc[i]];
    }

    scalarField pnf = ggiPatch.interpolate(sField);

//     if
//     (
//         reinterpret_cast<const void*>(&psiInternal)
//      == reinterpret_cast<const void*>(&this->internalField())
//     )
//     {

//     nIter_[cmpt]++;

//     if (nIter_[cmpt]<3)
//     {
//         if (ggiPatch.master())
//         {
//             scalarField jump = displacementJump_.component(cmpt);

//             for (label faceI = 0; faceI < pnf.size(); faceI++)
//             {
//                 pnf[faceI] -= jump[faceI];
//             }
//         }
//         else
//         {
//             scalarField jump = displacementJump_.component(cmpt);

//             for (label faceI = 0; faceI < pnf.size(); faceI++)
//             {
//                 pnf[faceI] += jump[faceI];
//             }
//         }
//     }

//     }

    // Multiply the field by coefficients and add into the result
    const unallocLabelList& fc = ggiPatch.faceCells();

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


tmp<Field<vector> > materialInterfaceFvPatchVectorField::
gradientInternalCoeffs() const
{
    const ggiFvPatch& ggiPatch =
        refCast<const ggiFvPatch>(this->patch());

    scalarField areaRatio(this->patch().size(), 1);

    if (!ggiPatch.master())
    {
        scalarField masterMagS =
            ggiPatch.interpolate
            (
                mag(ggiPatch.shadow().Sf())
            );
        scalarField slaveMagS = mag(ggiPatch.Sf());

        areaRatio = masterMagS/slaveMagS;
    }

    tmp<Field<vector> > tgic
    (
        new vectorField
        (
            -pTraits<vector>::one*this->patch().deltaCoeffs()*areaRatio
        )
    );

    // if (max(debonded_) > SMALL)
    // {
    //     forAll(debonded_, faceI)
    //     {
    //         if (debonded_[faceI] > SMALL)
    //         {
    //             tgic()[faceI] = vector::zero;
    //         }
    //     }
    // }

    return tgic;

//     return -pTraits<vector>::one*this->patch().deltaCoeffs()*areaRatio;
}


tmp<Field<vector> > materialInterfaceFvPatchVectorField::
gradientBoundaryCoeffs() const
{
    tmp<Field<vector> > tgbc
    (
        new vectorField
        (
            -this->gradientInternalCoeffs()
        )
    );

//     if (max(debonded_) > SMALL)
//     {
//         forAll(debonded_, faceI)
//         {
//             if (debonded_[faceI] > SMALL)
//             {
//                 tgbc()[faceI] = grad_[faceI];
//             }
//         }
//     }

    return tgbc;
}


void materialInterfaceFvPatchVectorField::write(Ostream& os) const
{
    ggiFvPatchField<vector>::write(os);

    //traction_.writeEntry("traction", os);
    //initiationTraction_.writeEntry("initiationTraction", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    materialInterfaceFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
