/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |

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

#include "elastoPlasticBatheHill.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformGeometricField.H"
#include "fvc.H"
#include "mechanicalModel.H"
#include "thermalModel.H"
#include "logVolFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(elastoPlasticBatheHill, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, elastoPlasticBatheHill, dictionary
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    scalar elastoPlasticBatheHill::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

    scalar elastoPlasticBatheHill::sqrtThreeOverTwo_ = ::sqrt(3.0/2.0);

    Eigen::Matrix2d elastoPlasticBatheHill::I2_ = Eigen::Matrix2d::Identity();

} // End of namespace Foam


// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

Foam::volTensorField& Foam::elastoPlasticBatheHill::FpInv()
{
    if (!FpInvPtr_)
    {
        FatalErrorIn
        (
            "Foam::volTensorField& Foam::elastoPlasticBatheHill::FpInv()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *FpInvPtr_;
}


const Foam::volTensorField& Foam::elastoPlasticBatheHill::FpInv() const
{
    if (!FpInvPtr_)
    {
        FatalErrorIn
        (
            "Foam::volTensorField& Foam::elastoPlasticBatheHill::FpInv()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *FpInvPtr_;
}


Foam::volSymmTensorField& Foam::elastoPlasticBatheHill::Dp()
{
    if (!DpPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& Foam::elastoPlasticBatheHill::Dp()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *DpPtr_;
}


const Foam::volSymmTensorField& Foam::elastoPlasticBatheHill::Dp() const
{
    if (!DpPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& Foam::elastoPlasticBatheHill::Dp()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *DpPtr_;
}


Foam::volScalarField& Foam::elastoPlasticBatheHill::epsilonPEq()
{
    if (!epsilonPEqPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticBatheHill::epsilonPEq()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *epsilonPEqPtr_;
}


const Foam::volScalarField& Foam::elastoPlasticBatheHill::epsilonPEq() const
{
    if (!epsilonPEqPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticBatheHill::epsilonPEq()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *epsilonPEqPtr_;
}


Foam::volScalarField& Foam::elastoPlasticBatheHill::activeYield()
{
    if (!activeYieldPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticBatheHill::activeYield()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *activeYieldPtr_;
}


const Foam::volScalarField& Foam::elastoPlasticBatheHill::activeYield() const
{
    if (!activeYieldPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticBatheHill::activeYield()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *activeYieldPtr_;
}


Foam::symmTensor Foam::elastoPlasticBatheHill::calculateStress
(
    const symmTensor& Ee
)
{
    if (elasticIsotropy_)
    {
        // Hooke's law
        return
            symmTensor
            (
                Ee*2.0*mu_.value() + lambda_.value()*tr(Ee)*I
            );
    }
    else
    {
        // Convert to the tensor basis A1-5
        const tensor Ee1 = 0.5*(Ee && A1_)*A1_;
        const tensor Ee2 = (3.0/2.0)*(Ee && A2_)*A2_;
        const tensor Ee3 = 0.5*(Ee && A3_)*A3_;
        const tensor Ee4 = 0.5*(Ee && A4_)*A4_;
        const tensor Ee5 = 0.5*(Ee && A5_)*A5_;

        // Return stress
        return
            symmTensor
            (
                symm
                (
                    kappa_*tr(Ee)*I + 2.0*mu1_*Ee1 + 2.0*mu2_*Ee2
                    + 2.0*mu3_*Ee3  + 2.0*mu4_*Ee4 + 2.0*mu5_*Ee5
                    + beta1_*( (Ee && A1_)*A2_ + (Ee && A2_)*A1_ )
                    + beta5_*( (Ee && A1_)*I + tr(Ee)*A1_ )
                    + beta6_*( (Ee && A2_)*I + tr(Ee)*A2_ )
                )
            );
    }
}

Foam::symmTensor Foam::elastoPlasticBatheHill::logm(const symmTensor& T)
{
    // Finds the matrix log of the tensor T

    // Convert T to a MatrixXd
    Eigen::Matrix3d TM(3,3);
    TM  <<
        T.xx(), T.xy(), T.xz(),
        T.xy(), T.yy(), T.yz(),
        T.xz(), T.yz(), T.zz();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Eigen::Vector3d EVals;
    EVals = es.eigenvalues();

    Eigen::Vector3d LogEVals;
    LogEVals(0) = Foam::log(EVals(0));
    LogEVals(1) = Foam::log(EVals(1));
    LogEVals(2) = Foam::log(EVals(2));

    Eigen::Matrix3d D = LogEVals.asDiagonal();

    Eigen::Matrix3d EVecs;
    EVecs = es.eigenvectors();

    Eigen::Matrix3d resultM = EVecs*D*EVecs.inverse();

    return
        symmTensor
        (
            resultM(0,0), resultM(0,1), resultM(0,2),
                          resultM(1,1), resultM(1,2),
                                        resultM(2,2)
        );
}

Foam::symmTensor Foam::elastoPlasticBatheHill::expm(const symmTensor& T)
{
    // Calculate exponetial of a tensor

    // Convert T to a MatrixXd
    Eigen::Matrix3d TM(3,3);
    TM  <<
        T.xx(), T.xy(), T.xz(),
        T.xy(), T.yy(), T.yz(),
        T.xz(), T.yz(), T.zz();

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.compute(TM);

    Eigen::Vector3d EVals;
    EVals = es.eigenvalues();

    Eigen::Vector3d ExpEVals;
    ExpEVals(0) = Foam::exp(EVals(0));
    ExpEVals(1) = Foam::exp(EVals(1));
    ExpEVals(2) = Foam::exp(EVals(2));

    Eigen::Matrix3d D = ExpEVals.asDiagonal();

    Eigen::Matrix3d EVecs;
    EVecs = es.eigenvectors();

    Eigen::Matrix3d resultM = EVecs*D*EVecs.inverse();

    return
        symmTensor
        (
            resultM(0,0), resultM(0,1), resultM(0,2),
                          resultM(1,1), resultM(1,2),
                                        resultM(2,2)
        );
}


void Foam::elastoPlasticBatheHill::smallStrainReturnMap
(
    symmTensor& TTrial,
    tensor& FpInv,
    const tensor& FpInvOld,
    symmTensor& Dp,
    scalar& activeYield,
    scalar& epsilonPEq,
    const scalar& epsilonPEqOld
)
{
    symmTensor Sigma = dev(TTrial);

    // Convert to the tensor  basis A1-5
    const tensor Sigma1 = 0.5*(Sigma && A1_)*A1_;
    const tensor Sigma2 = (3.0/2.0)*(Sigma && A2_)*A2_;
    const tensor Sigma3 = 0.5*(Sigma && A3_)*A3_;
    const tensor Sigma4 = 0.5*(Sigma && A4_)*A4_;
    const tensor Sigma5 = 0.5*(Sigma && A5_)*A5_;

    scalar sigmaY =
        y0_ + hIso_*epsilonPEqOld
      + (yInf_ - y0_)*(1.0 - Foam::exp(-beta_*epsilonPEqOld));

    // Check the yield function
    const scalar phi =
        Foam::sqrt
        (
            alpha1_*Foam::sqr(mag(Sigma1)) + alpha2_*Foam::sqr(mag(Sigma2))
          + alpha3_*Foam::sqr(mag(Sigma3)) + alpha4_*Foam::sqr(mag(Sigma4))
          + alpha5_*Foam::sqr(mag(Sigma5))
          + gamma1_*(Sigma && A1_)*(Sigma && A2_)
        ) - sigmaY;

    if (phi > 0.0)
    {
        const scalar STR1 = Sigma && A1_;
        const scalar STR2 = Sigma && A2_;

        Eigen::Vector2d STR12;
        STR12
            << STR1, STR2;

        const scalar STR3 = 0.5*(Sigma && A3_);
        const scalar STR4 = 0.5*(Sigma && A4_);
        const scalar STR5 = 0.5*(Sigma && A5_);

        scalar deltaGamma = 0.0;
        scalar R = 1.0;
        scalar epsilonPEqNew = epsilonPEqOld;
        Eigen::Vector2d S12;
        scalar S3 = 0.0;
        scalar S4 = 0.0;
        scalar S5 = 0.0;

        for (int iteration = 0; iteration < 200; iteration++)
        {
            epsilonPEqNew = epsilonPEqOld + deltaGamma;

            sigmaY =
                y0_ + hIso_*epsilonPEqNew
              + (yInf_ - y0_)*(1.0 - Foam::exp(-beta_*epsilonPEqNew));

            const scalar DsigmaY =
                hIso_ + beta_*(yInf_ - y0_)*Foam::exp(-beta_*epsilonPEqNew);

            Eigen::Matrix2d F12;
            F12 = (sigmaY*I2_ + deltaGamma*Y12_).inverse();

            const scalar F3 = 1.0/(sigmaY + deltaGamma*Y3_);
            const scalar F4 = 1.0/(sigmaY + deltaGamma*Y4_);
            const scalar F5 = 1.0/(sigmaY + deltaGamma*Y5_);

            S12 = F12*STR12;

            S3 = F3*STR3;
            S4 = F4*STR4;
            S5 = F5*STR5;

            Eigen::Matrix2d X12;
            X12 = P12_*F12;

            const scalar X3 = 2.0*alpha3_*F3;
            const scalar X4 = 2.0*alpha4_*F4;
            const scalar X5 = 2.0*alpha5_*F5;

            Eigen::Vector2d xs12;
            xs12 = X12*S12;

            const scalar xs3 = X3*S3;
            const scalar xs4 = X4*S4;
            const scalar xs5 = X5*S5;

            R =
                S12.dot(P12_*S12)
              + 2.0*( alpha3_*Foam::sqr(S3)
              + alpha4_*Foam::sqr(S4)
              + alpha5_*Foam::sqr(S5) ) - 1.0;

            const scalar aHat =
                S12.dot((DsigmaY*I2_ + Y12_.transpose())*xs12)
              + S3*(DsigmaY + Y3_)*xs3
              + S4*(DsigmaY + Y4_)*xs4
              + S5*(DsigmaY + Y5_)*xs5;

            deltaGamma = deltaGamma + R/(2.0*aHat);

            if (mag(R) < yieldTol_)
            {
                break;
            }

            activeYield = R;
        }

        // update
        Eigen::Vector2d gs12;
        gs12 = deltaGamma*(P12_*S12);

        scalar gs3 = deltaGamma*(alpha3_*S3);
        scalar gs4 = deltaGamma*(alpha4_*S4);
        scalar gs5 = deltaGamma*(alpha5_*S5);

        Dp =
            symm
            (
                gs12(0)*A1_ + gs12(1)*A2_
              + gs3*A3_ + gs4*A4_ + gs5*A5_
            );

        //tensor n12_1 = D12_(0,0)*A1_ + D12_(0,1)*A2_ + beta5_*I;
        //tensor n12_2 = D12_(1,0)*A1_ + D12_(1,1)*A2_ + beta6_*I;

        Eigen::Vector2d Qgs12;
        Qgs12 = Q_*gs12;

        TTrial = TTrial - calculateStress(Dp);

        //scalar TEq = sqrtThreeOverTwo_*Foam::sqrt(dev(TTrial) && dev(TTrial));
        FpInv = FpInvOld & expm(-Dp);
        epsilonPEq = epsilonPEqNew;
        //activeYield = true;
    }
    else
    {
        // elastic
        Dp = symmTensor::zero;
        activeYield = false;
        FpInv = FpInvOld;
        epsilonPEq = epsilonPEqOld;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::elastoPlasticBatheHill::elastoPlasticBatheHill
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const label lawIndex
)
:
    mechanicalLaw(name, mesh, dict, lawIndex),
    elasticIsotropy_(dict.lookupOrDefault<Switch>("elasticIsotropy", true)),
    plasticIsotropy_(dict.lookupOrDefault<Switch>("plasticIsotropy", false)),
    rho_(dict.lookup("rho")),
    E_
    (
        elasticIsotropy_
      ? dimensionedScalar(dict.lookup("E"))
      : 1.0
    ),
    nu_
    (
        elasticIsotropy_
      ? dimensionedScalar(dict.lookup("nu"))
      : 1.0
    ),
    mu_(E_/(2.0*(1.0 + nu_))),
    lambda_
    (
        mesh.lookupObject<mechanicalModel>("mechanicalProperties").planeStress()
      ? nu_*E_/((1.0 + nu_)*(1.0 - nu_))
      : nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))
    ),
    K_
    (
        mesh.lookupObject<mechanicalModel>("mechanicalProperties").planeStress()
      ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
      : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
    ),
    E1_
    (
        elasticIsotropy_
      ? 1.0
      : readScalar(dict.lookup("E1"))
    ),
    E2_
    (
        elasticIsotropy_
      ? 1.0
      : readScalar(dict.lookup("E2"))
    ),
    E3_
    (
        elasticIsotropy_
      ? 1.0
      : readScalar(dict.lookup("E3"))
    ),
    G12_
    (
        elasticIsotropy_
      ? 1.0
      : readScalar(dict.lookup("G12"))
    ),
    G23_
    (
        elasticIsotropy_
      ? 1.0
      : readScalar(dict.lookup("G23"))
    ),
    G31_
    (
        elasticIsotropy_
      ? 1.0
      : readScalar(dict.lookup("G23"))
    ),
    nu21_
    (
        elasticIsotropy_
      ? 1.0
      : readScalar(dict.lookup("nu21"))
    ),
    nu23_
    (
        elasticIsotropy_
      ? 1.0
      : readScalar(dict.lookup("nu23"))
    ),
    nu31_
    (
        elasticIsotropy_
      ? 1.0
      : readScalar(dict.lookup("nu31"))
    ),
    nu12_(nu21_*E1_/E2_),
    nu32_(nu23_*E3_/E2_),
    nu13_(nu31_*E1_/E3_),
    kappa_(K_.value()),
    mu1_(mu_.value()),
    mu2_(mu_.value()),
    mu3_(mu_.value()),
    mu4_(mu_.value()),
    mu5_(mu_.value()),
    beta1_(0.0),
    beta5_(0.0),
    beta6_(0.0),
    N1_(dict.lookupOrDefault<vector>("N1", vector(1, 0, 0))),
    N2_(dict.lookupOrDefault<vector>("N2", vector(0, 1, 0))),
    N3_(dict.lookupOrDefault<vector>("N3", vector(0, 0, 1))),
    A0_(I),
    A1_(N2_*N2_ - N3_*N3_),
    A2_(N1_*N1_ - I/3.0),
    A3_(N1_*N2_ + N2_*N1_),
    A4_(N2_*N3_ + N3_*N2_),
    A5_(N3_*N1_ + N1_*N3_),
    Q_(2, 2),
    R11_
    (
        plasticIsotropy_
      ? 1.0
      : readScalar(dict.lookup("R11"))
    ),
    R22_
    (
        plasticIsotropy_
      ? 1.0
      : readScalar(dict.lookup("R22"))
    ),
    R33_
    (
        plasticIsotropy_
      ? 1.0
      : readScalar(dict.lookup("R33"))
    ),
    R12_
    (
        plasticIsotropy_
      ? 1.0/Foam::sqrt(3.0)
      : readScalar(dict.lookup("R12"))
    ),
    R23_
    (
        plasticIsotropy_
      ? 1.0/Foam::sqrt(3.0)
      : readScalar(dict.lookup("R23"))
    ),
    R31_
    (
        plasticIsotropy_
      ? 1.0/Foam::sqrt(3.0)
      : readScalar(dict.lookup("R31"))
    ),
    // P-Lu plasticity parameters
    alpha1_
    (
        (1.0/Foam::sqr(R22_)) + (1.0/Foam::sqr(R33_)) - (0.5/Foam::sqr(R11_))
    ),
    alpha2_((3.0/2.0)*(1.0/Foam::sqr(R11_))),
    alpha3_(0.5/Foam::sqr(R12_)),
    alpha4_(0.5/Foam::sqr(R23_)),
    alpha5_(0.5/Foam::sqr(R31_)),
    gamma1_((3.0/2.0)*((1.0/Foam::sqr(R33_)) - (1.0/Foam::sqr(R22_)))),
    DStar12_(2, 2),
    PStar12_(2, 2),
    D12_(2, 2),
    P12_(2, 2),
    Y12_(2, 2),
    Y3_(0.0),
    Y4_(0.0),
    Y5_(0.0),
    FpInvPtr_(NULL),
    DpPtr_(NULL),
    epsilonPEqPtr_(NULL),
    activeYieldPtr_(NULL),
    y0_(readScalar(dict.lookup("sigmaY0"))),
    beta_(readScalar(dict.lookup("beta"))),
    hIso_(readScalar(dict.lookup("hIso"))),
    yInf_(readScalar(dict.lookup("sigmaYInfinity"))),
    yieldTol_(1e-15),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    ),
    damageLawPtr_(NULL)
{
    // Create fields
    // We will check if the field already exists in the object registry (i.e. if
    // it has been created by another mechanical law); if it is not found then
    // the field is created; if it is found, then we will set a pointer to it
    // using const_cast. This may not be the most elegant or safe solution but
    // it is OK for now!

    SetFieldPtr<tensor>
    (
        FpInvPtr_,
        "FpInv",
        dimensionedTensor("one", dimless, tensor(I))
    );

    SetFieldPtr<symmTensor>
    (
        DpPtr_,
        "Dp",
        dimensionedSymmTensor("one", dimless, symmTensor::zero)
    );

    SetFieldPtr<scalar>
    (
        epsilonPEqPtr_,
        "epsilonPEq",
        dimensionedScalar("zero", dimless, 0.0)
    );

    SetFieldPtr<scalar>
    (
        activeYieldPtr_,
        "activeYield",
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Store old time for any field that we will need it for
    Dp().oldTime();
    FpInv().oldTime();
    epsilonPEq().oldTime();

    // On restart, some fields may be defaulted to zero when they should default
    // to I
    if (gMin(mag(FpInv().internalField())) < SMALL)
    {
        WarningIn("elastoPlasticBatheHill::elastoPlasticBatheHill()")
            << "Resetting zero in FpInv fields to I" << endl;

        tensorField& FpInvI = FpInv().internalField();

        forAll(FpInvI, cellI)
        {
            if (mag(FpInvI[cellI]) < SMALL)
            {
                FpInvI[cellI] = I;
            }
        }

        forAll(FpInv().boundaryField(), patchI)
        {
            forAll(FpInv().boundaryField()[patchI], faceI)
            {
                if (mag(FpInv().boundaryField()[patchI][faceI]) < SMALL)
                {
                    FpInv().boundaryField()[patchI][faceI] = I;
                }
            }
        }

        FpInv().correctBoundaryConditions();
    }

    if (!elasticIsotropy_)
    {
        Info<< "    Elasticity is orthotropic" << endl;

        Eigen::MatrixXd DInv = Eigen::MatrixXd::Zero(6,6);

        DInv(0,0) = 1.0/E1_;
        DInv(0,1) = -nu21_/E2_;
        DInv(0,2) = -nu31_/E3_;
        DInv(1,0) = -nu21_/E2_;
        DInv(1,1) = 1.0/E2_;
        DInv(1,2) = -nu23_/E2_;
        DInv(2,0) = -nu31_/E3_;
        DInv(2,1) = -nu23_/E2_;
        DInv(2,2) = 1.0/E3_;

        DInv(3,3) = 1.0/G12_;
        DInv(4,4) = 1.0/G23_;
        DInv(5,5) = 1.0/G31_;

        Eigen::MatrixXd DMat = DInv.inverse();

        scalar De11_ = DMat(0,0);
        scalar De22_ = DMat(1,1);
        scalar De33_ = DMat(2,2);
        scalar De44_ = DMat(3,3);
        scalar De55_ = DMat(4,4);
        scalar De66_ = DMat(5,5);

        scalar De23_ = DMat(1,2);
        scalar De12_ = DMat(0,1);
        scalar De31_ = DMat(2,0);


        // Converting to P-Lu parameters
        Eigen::VectorXd DeVec(6);

        DeVec
            << De11_, De12_, De31_, De22_, De23_, De33_;

        Eigen::MatrixXd EC(6,6);

        EC(0,0) = 4.0;
        EC(0,1) = 8.0;
        EC(0,2) = 8.0;
        EC(0,3) = 4.0;
        EC(0,4) = 8.0;
        EC(0,5) = 4.0;
        EC(1,0) = 0.0;
        EC(1,1) = 0.0;
        EC(1,2) = 0.0;
        EC(1,3) = 9.0;
        EC(1,4) = -18;
        EC(1,5) = 9.0;
        EC(2,0) = 12;
        EC(2,1) = -12;
        EC(2,2) = -12;
        EC(2,3) = 3.0;
        EC(2,4) = 6.0;
        EC(2,5) = 3.0;
        EC(3,0) = 0.0;
        EC(3,1) = 18;
        EC(3,2) = -18;
        EC(3,3) = -9;
        EC(3,4) = 0.0;
        EC(3,5) = 9.0;
        EC(4,0) = 0.0;
        EC(4,1) = 6.0;
        EC(4,2) = -6;
        EC(4,3) = 6.0;
        EC(4,4) = 0.0;
        EC(4,5) = -6;
        EC(5,0) = 12;
        EC(5,1) = 6.0;
        EC(5,2) = 6.0;
        EC(5,3) = -6.0;
        EC(5,4) = -12.0;
        EC(5,5) = -6.0;

        Eigen::VectorXd PLUE(6);

        PLUE = (EC*DeVec)/36.0;

        kappa_ = PLUE(0);
        mu1_ = PLUE(1);
        mu2_ = PLUE(2);
        beta1_ = PLUE(3);
        beta5_ = PLUE(4);
        beta6_ = PLUE(5);
        mu3_ = De44_;
        mu4_ = De55_;
        mu5_ = De66_;
    }

    // Const matrix
    Q_(0,0) = 2.0;
    Q_(0,1) = 0.0;
    Q_(0,0) = 0.0;
    Q_(1,1) = 2.0/3.0;

    // P-Lu matrices for algorithm given in ETH thesis
    DStar12_(0,0) = 2.0*mu1_;
    DStar12_(0,1) = 2.0*beta1_;
    DStar12_(1,0) = (2.0/3.0)*beta1_;
    DStar12_(1,1) = 2.0*mu2_;

    D12_(0,0) = mu1_;
    D12_(0,1) = beta1_;
    D12_(1,0) = beta1_;
    D12_(1,1) = 3.0*mu2_;

    PStar12_(0,0) = alpha1_;
    PStar12_(0,1) = gamma1_;
    PStar12_(1,0) = gamma1_/3.0;
    PStar12_(1,1) = alpha2_;

    P12_(0,0) = alpha1_/2.0;
    P12_(0,1) = gamma1_/2.0;
    P12_(1,0) = gamma1_/2.0;
    P12_(1,1) = (3.0/2.0)*alpha2_;

    Y12_ = DStar12_*PStar12_;
    Y3_ = 2.0*mu3_*alpha3_;
    Y4_ = 2.0*mu4_*alpha4_;
    Y5_ = 2.0*mu5_*alpha5_;

     // Check if a damage law is specified
    if (dict.found("damageLaw"))
    {
        Info<< "    Creating damage law" << endl;
        damageLawPtr_ =
            damageLaw::New
            (
                "law", mesh, dict.subDict("damageLaw"), curMaterial()
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::elastoPlasticBatheHill::~elastoPlasticBatheHill()
{
    CleanFieldPtr<tensor>(FpInvPtr_, "FpInv");
    CleanFieldPtr<symmTensor>(DpPtr_, "Dp");
    CleanFieldPtr<scalar>(epsilonPEqPtr_, "epsilonPEq");
    CleanFieldPtr<scalar>(activeYieldPtr_, "activeYield");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::elastoPlasticBatheHill::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            calculatedFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::elastoPlasticBatheHill::E() const
{
    if (elasticIsotropy_)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "E",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                E_,
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "E",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("E", dimPressure, (E1_ + E2_ + E3_)/3.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
}


Foam::tmp<Foam::volScalarField> Foam::elastoPlasticBatheHill::nu() const
{
    if (elasticIsotropy_)
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "nu",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                nu_,
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "nu",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("nu", dimless, (nu21_ + nu23_ + nu31_)/3.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::elastoPlasticBatheHill::RhieChowScaleFactor() const
{
    const scalar scaleFac =
        dict().lookupOrDefault<scalar>("RhieChowScaleFactor", 0.1);

    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
              "RhieChowScaleFactor",
              mesh().time().timeName(),
              mesh(),
              IOobject::NO_READ,
              IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, scaleFac)
        )
    );

    Info<< type() << ": RhieChowScaleFactor = " << scaleFac << endl;

    return tresult;
}


void Foam::elastoPlasticBatheHill::correct(volSymmTensorField& tau)
{
    const fvMesh& mesh = this->mesh();

    // Lookup relative deformation gradient from the solver
    const volTensorField& F =
        mesh.objectRegistry::lookupObject<volTensorField>("F");

    Dp().storePrevIter();

    // Take references to internal fields
    const tensorField& FI = F.internalField();
    symmTensorField& tauI = tau.internalField();
    symmTensorField& DpI = Dp().internalField();
    scalarField& epsilonPEqI = epsilonPEq().internalField();
    tensorField& FpInvI = FpInv().internalField();
    tensorField& FpInvIOld = FpInv().oldTime().internalField();
    const scalarField& epsilonPEqIOld = epsilonPEq().oldTime().internalField();
    scalarField& activeYieldI = activeYield().internalField();
    const scalarField& curMaterialI = curMaterial().internalField();

    // Storage for the subsequent loops
    tensor Fe = I;
    symmTensor Ce = I;
    symmTensor Ee = symmTensor::zero;
    symmTensor T = symmTensor::zero;
    tensor invCe = I;
    tensor mandel = tensor::zero;
    tensor SPK2 = tensor::zero;

    forAll(FI, cellI)
    {
        // Dont want to do the return map in the wrong material region
        if (curMaterialI[cellI] > SMALL)
        {
            // Pre-processor
            // Trial elastic def grad
            Fe = FI[cellI] & FpInvIOld[cellI];
            Ce = symm(Fe.T() & Fe);
            Ee = 0.5*logm(Ce);
            T = calculateStress(Ee);

            smallStrainReturnMap
            (
                T,
                FpInvI[cellI],
                FpInvIOld[cellI],
                DpI[cellI],
                activeYieldI[cellI],
                epsilonPEqI[cellI],
                epsilonPEqIOld[cellI]
            );

            // Post-processor
            Ee = Ee - DpI[cellI];
            invCe = inv(expm(2.0*Ee));
            mandel = T + ((Ee & T) - (T & Ee));
            SPK2 = 0.5*((invCe & mandel) + (mandel & invCe));
            Fe = FI[cellI] & FpInvI[cellI];
            tauI[cellI] = symm(Fe & (SPK2 & Fe.T()));
        }
    }

    forAll(F.boundaryField(), patchI)
    {
        // Note: we no longer treated processor boundaries differently, and we
        // should not called correctBoundaryConditions below
        // if
        // (
        //     !F.boundaryField()[patchI].coupled()
        //  && mesh.boundaryMesh()[patchI].type() != "empty"
        // )
        {
            const tensorField& FP = F.boundaryField()[patchI];
            symmTensorField& tauP = tau.boundaryField()[patchI];
            symmTensorField& DpP = Dp().boundaryField()[patchI];
            scalarField& epsilonPEqP = epsilonPEq().boundaryField()[patchI];
            tensorField& FpInvP = FpInv().boundaryField()[patchI];
            const tensorField& FpInvPOld =
                FpInv().oldTime().boundaryField()[patchI];
            const scalarField& epsilonPEqPOld =
                epsilonPEq().oldTime().boundaryField()[patchI];
            scalarField& activeYieldP = activeYield().boundaryField()[patchI];
            const scalarField& curMaterialP =
                curMaterial().boundaryField()[patchI];

            forAll(FP, faceI)
            {
                // Don't want to do the return map in the wrong material region
                if (curMaterialP[faceI] > SMALL)
                {
                    // Pre-processor
                    // Trial elastic def grad
                    Fe = FP[faceI] & FpInvPOld[faceI];
                    Ce = symm(Fe.T() & Fe);
                    Ee = 0.5*logm(Ce);
                    T = calculateStress(Ee);

                    smallStrainReturnMap
                    (
                        T,
                        FpInvP[faceI],
                        FpInvPOld[faceI],
                        DpP[faceI],
                        activeYieldP[faceI],
                        epsilonPEqP[faceI],
                        epsilonPEqPOld[faceI]
                    );

                    // Post-processor
                    Ee = Ee - DpP[faceI];
                    invCe = inv(expm(2.0*Ee));
                    mandel = T + ((Ee & T) - (T & Ee));
                    SPK2 = 0.5*((invCe & mandel) + (mandel & invCe));
                    Fe = FP[faceI] & FpInvP[faceI];
                    tauP[faceI] = symm(Fe & (SPK2 & Fe.T()));
                }
            }
        }
    }

    // Scale stress if damage law is coupled
    if (damageLawPtr_.valid())
    {
        if (damageLawPtr_->coupled())
        {
            // Update the damage field
            damageLawPtr_->updateDamage(sqrt((2.0/3.0)*magSqr(dev(Dp()))));

            // Scale the stress field
            tau = (1.0 - damageLawPtr_->damage())*tau;
        }
    }

    // We no longer treat rocessor boundaries differently
    // FpInv().correctBoundaryConditions();
    // Dp().correctBoundaryConditions();
    // epsilonPEq().correctBoundaryConditions();
    // activeYield().correctBoundaryConditions();
    // tau.correctBoundaryConditions();
}

Foam::scalar Foam::elastoPlasticBatheHill::residual()
{
    return
        gMax
        (
            mag
            (
                Dp().internalField()
              - Dp().prevIter().internalField()
            )
        )/gMax(SMALL + mag(Dp().prevIter().internalField()));
}



void Foam::elastoPlasticBatheHill::correct(surfaceSymmTensorField& tau)
{
    notImplemented
    (
        "elastoPlasticBatheHill::correct(surfaceSymmTensorField& tau)"
    );
}


Foam::tmp<Foam::volScalarField>
Foam::elastoPlasticBatheHill::plasticDissipationRate() const
{
    const volSymmTensorField& sigmaCauchy =
        mesh().lookupObject<volSymmTensorField>("sigmaCauchy");

    // We assume 90% conversion
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "plasticDissipationRate",
            max
            (
                dimensionedScalar("zero", dimForce/(dimArea*dimTime), 0.0),
                0.9*(sigmaCauchy && Dp())/mesh().time().deltaT()
            )
        )
    );
}


Foam::scalar Foam::elastoPlasticBatheHill::newDeltaT() const
{
    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time- step, where the return direction is given normalised deviatoric
    // strain. The error approximation is obtained using the difference between
    // the trapezoidal rule and the EUler backward method, as described in:

    // Nam-Sua Lee, Klaus-Jurgen Bathe, Error indicators and adaptive remeshing
    // in large deformation finite element analysis, Finite Elements in
    // Analysis and Design 16 (1994) 99-139.

    // Lookup the deformation gradient
    const volTensorField& F = mesh().lookupObject<volTensorField>("F");

    // Calculate the true (Hencky) strain
    volSymmTensorField epsilon = 0.5*log(symm(F.T() & F));

    // Calculate equivalent strain, for normalisation of the error
    volScalarField epsilonEq = sqrt((2.0/3.0)*magSqr(dev(epsilon)));

    // Calculate plastic return direction as normalised plastic strain increment
    // PC: soemthing strange is happening because Dp() should be equal to
    // Dp().oldTime as we are at the start of the time step.
    // const volSymmTensorField DEpsilonPOld = symm(Dp().oldTime());
    // const volSymmTensorField DEpsilonPOldOld =
    //    symm(Dp().oldTime().oldTime());
    const volSymmTensorField DEpsilonPOld = symm(Dp());
    const volSymmTensorField DEpsilonPOldOld = symm(Dp().oldTime());
    const volSymmTensorField plasticNOld =
        DEpsilonPOld/(mag(DEpsilonPOld) + SMALL);
    const volSymmTensorField plasticNOldOld =
        DEpsilonPOldOld/(mag(DEpsilonPOldOld) + SMALL);

    // Take reference to internal fields
    const symmTensorField& DEpsilonPI = DEpsilonPOld.internalField();
    const symmTensorField& plasticNIold = plasticNOld.internalField();
    const symmTensorField& plasticNIoldOld = plasticNOldOld.internalField();
    const scalarField& epsilonEqI = epsilonEq.internalField();

    // Calculate error field
    const symmTensorField DEpsilonPErrorI =
        Foam::sqrt(3.0/8.0)*DEpsilonPI*mag(plasticNIold - plasticNIoldOld)
       /(epsilonEqI + SMALL);

    // Max error
    const scalar maxMagDEpsilonPErr = gMax(mag(DEpsilonPErrorI));

    if (maxMagDEpsilonPErr > SMALL)
    {
        Info<< "    " << name() << ": max error = " << maxMagDEpsilonPErr
            << endl;

        // Calculate the time-step scaling factor, where maxDeltaErr_ is the
        // maximum allowed error
        const scalar scaleFac = maxDeltaErr_/maxMagDEpsilonPErr;

        // Return the new time-step size
        return scaleFac*mesh().time().deltaTValue();
    }

    return mesh().time().endTime().value();
}


void Foam::elastoPlasticBatheHill::updateYieldStress()
{
    // Force recalculation of curMaterial fields
    mechanicalLaw::updateYieldStress();

    // Count how many cells are actively yielding

    int numCellsYielding = 0;

    scalarField& activeYieldI = activeYield().internalField();

    forAll(activeYieldI, cellI)
    {
        if (magSqr(activeYieldI[cellI]) > 0.0)
        {
            numCellsYielding++;

            activeYieldI[cellI] = 1.0;
        }
        else
        {
            activeYieldI[cellI] = 0.0;
        }
    }


    forAll(activeYield().boundaryField(), patchI)
    {
        if
        (
            !activeYield().boundaryField()[patchI].coupled()
         && mesh().boundaryMesh()[patchI].type() != "empty"
        )
        {
            scalarField& activeYieldP = activeYield().boundaryField()[patchI];

            forAll(activeYieldP, faceI)
            {
                if (magSqr(activeYieldP[faceI]) > 0.0)
                {
                    activeYieldP[faceI] = 1.0;
                }
                else
                {
                    activeYieldP[faceI] = 0.0;
                }
            }
        }
    }

    activeYield().correctBoundaryConditions();

    reduce(numCellsYielding, sumOp<int>());

    Info<< nl << "There are " << numCellsYielding << " cells actively yielding"
        << nl << endl;

    if (damageLawPtr_.valid())
    {
        // Update the damage field
        damageLawPtr_->updateDamage(sqrt((2.0/3.0)*magSqr(dev(Dp()))));

        Info<< "Maximum damage: " << gMax(damageLawPtr_->damage()) << nl
            << "Minimum damage: " << gMin(damageLawPtr_->damage()) << nl
            << "Average damage: "<< gAverage(damageLawPtr_->damage()) << nl
            << endl;
    }

    if (mesh().time().outputTime())
    {
        volSymmTensorField yieldStress
        (
            IOobject
            (
                "yieldStress",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            (
                y0_ + hIso_*epsilonPEq()
              + (yInf_ - y0_)*(1.0 - Foam::exp(-beta_*epsilonPEq()))
            )*symmTensor(R11_, R12_, R31_, R22_, R23_, R33_)
        );

        yieldStress.write();
    }
}


// ************************************************************************* //
