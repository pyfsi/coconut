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

#include "elastoPlasticDistortional.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformGeometricField.H"
#include "fvc.H"
#include "mechanicalModel.H"
#include "logVolFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(elastoPlasticDistortional, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, elastoPlasticDistortional, dictionary
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    scalar elastoPlasticDistortional::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

    scalar elastoPlasticDistortional::sqrtOneOverThree_ = ::sqrt(1.0/3.0);

    scalar elastoPlasticDistortional::sqrtThreeOverTwo_ = ::sqrt(3.0/2.0);

} // End of namespace Foam


// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //


Foam::volTensorField& Foam::elastoPlasticDistortional::FpInv()
{
    if (!FpInvPtr_)
    {
        FatalErrorIn
        (
            "Foam::volTensorField& "
            "Foam::elastoPlasticDistortional::FpInv()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *FpInvPtr_;
}


const Foam::volTensorField& Foam::elastoPlasticDistortional::FpInv() const
{
    if (!FpInvPtr_)
    {
        FatalErrorIn
        (
            "Foam::volTensorField& "
            "Foam::elastoPlasticDistortional::FpInv()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *FpInvPtr_;
}


Foam::volSymmTensorField& Foam::elastoPlasticDistortional::henckyStrain()
{
    if (!henckyStrainPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& "
            "Foam::elastoPlasticDistortional::henckyStrain()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *henckyStrainPtr_;
}


const Foam::volSymmTensorField&
Foam::elastoPlasticDistortional::henckyStrain() const
{
    if (!henckyStrainPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& "
            "Foam::elastoPlasticDistortional::henckyStrain()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *henckyStrainPtr_;
}


Foam::volSymmTensorField& Foam::elastoPlasticDistortional::Dp()
{
    if (!DpPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& "
            "Foam::elastoPlasticDistortional::Dp()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *DpPtr_;
}


const Foam::volSymmTensorField& Foam::elastoPlasticDistortional::Dp() const
{
    if (!DpPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& "
            "Foam::elastoPlasticDistortional::Dp()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *DpPtr_;
}


Foam::volScalarField& Foam::elastoPlasticDistortional::epsilonPEq()
{
    if (!epsilonPEqPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& "
            "Foam::elastoPlasticDistortional::epsilonPEq()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *epsilonPEqPtr_;
}


const Foam::volScalarField&
Foam::elastoPlasticDistortional::epsilonPEq() const
{
    if (!epsilonPEqPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& "
            "Foam::elastoPlasticDistortional::epsilonPEq()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *epsilonPEqPtr_;
}


Foam::volSymmTensorField& Foam::elastoPlasticDistortional::backStrain()
{
    if (!backStrainPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& "
            "Foam::elastoPlasticDistortional::backStrain()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *backStrainPtr_;
}


const Foam::volSymmTensorField&
Foam::elastoPlasticDistortional::backStrain() const
{
    if (!backStrainPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& "
            "Foam::elastoPlasticDistortional::backStrain()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *backStrainPtr_;
}


Foam::volScalarField& Foam::elastoPlasticDistortional::activeYield()
{
    if (!activeYieldPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& "
            "Foam::elastoPlasticDistortional::activeYield()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *activeYieldPtr_;
}


const Foam::volScalarField&
Foam::elastoPlasticDistortional::activeYield() const
{
    if (!activeYieldPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& "
            "Foam::elastoPlasticDistortional::activeYield()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *activeYieldPtr_;
}


Foam::volTensorField&
Foam::elastoPlasticDistortional::aDist_1()
{
    if (!aDistPtr_1)
    {
        FatalErrorIn
        (
            "Foam::volDevSymmTensor4thOrderField& "
            "Foam::elastoPlasticDistortional::aDist_1()"
        )   << "pointer not set" << abort(FatalError);
    }
    return *aDistPtr_1;
}


Foam::volSymmTensorField&
Foam::elastoPlasticDistortional::aDist_2()
{
    if (!aDistPtr_2)
    {
        FatalErrorIn
        (
            "Foam::volDevSymmTensor4thOrderField& "
            "Foam::elastoPlasticDistortional::aDist_2()"
        )   << "pointer not set" << abort(FatalError);
    }
    return *aDistPtr_2;
}


Foam::volTensorField&
Foam::elastoPlasticDistortional::aLatt_1()
{
    if (!aLattPtr_1)
    {
        FatalErrorIn
        (
            "Foam::volDevSymmTensor4thOrderField& "
            "Foam::elastoPlasticDistortional::aLatt_1()"
        )   << "pointer not set" << abort(FatalError);
    }
    return *aLattPtr_1;
}


Foam::volSymmTensorField&
Foam::elastoPlasticDistortional::aLatt_2()
{
    if (!aLattPtr_2)
    {
        FatalErrorIn
        (
            "Foam::volDevSymmTensor4thOrderField& "
            "Foam::elastoPlasticDistortional::aLatt_2()"
        )   << "pointer not set" << abort(FatalError);
    }
    return *aLattPtr_2;
}


const Foam::volTensorField&
Foam::elastoPlasticDistortional::aDist_1() const
{
    if (!aDistPtr_1)
    {
        FatalErrorIn
        (
            "Foam::volDevSymmTensor4thOrderField& "
            "Foam::elastoPlasticDistortional::aDist_1()"
        )   << "pointer not set" << abort(FatalError);
    }
    return *aDistPtr_1;
}


const Foam::volSymmTensorField&
Foam::elastoPlasticDistortional::aDist_2() const
{
    if (!aDistPtr_2)
    {
        FatalErrorIn
        (
            "Foam::volDevSymmTensor4thOrderField& "
            "Foam::elastoPlasticDistortional::aDist_2()"
        )   << "pointer not set" << abort(FatalError);
    }
    return *aDistPtr_2;
}


const Foam::volTensorField&
Foam::elastoPlasticDistortional::aLatt_1() const
{
    if (!aLattPtr_1)
    {
        FatalErrorIn
        (
            "Foam::volDevSymmTensor4thOrderField& "
            "Foam::elastoPlasticDistortional::aLatt_1()"
        )   << "pointer not set" << abort(FatalError);
    }
    return *aLattPtr_1;
}


const Foam::volSymmTensorField&
Foam::elastoPlasticDistortional::aLatt_2() const
{
    if (!aLattPtr_2)
    {
        FatalErrorIn
        (
            "Foam::volDevSymmTensor4thOrderField& "
            "Foam::elastoPlasticDistortional::aLatt_2()"
        )   << "pointer not set" << abort(FatalError);
    }
    return *aLattPtr_2;
}


Foam::volSymmTensorField& Foam::elastoPlasticDistortional::yieldStress()
{
    if (!yieldStressPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& "
            "Foam::elastoPlasticDistortional::yieldStress()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *yieldStressPtr_;
}


const Foam::volSymmTensorField&
Foam::elastoPlasticDistortional::yieldStress() const
{
    if (!yieldStressPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& "
            "Foam::elastoPlasticDistortional::yieldStress()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *yieldStressPtr_;
}


Foam::volScalarField& Foam::elastoPlasticDistortional::plasticIncrement()
{
    if (!plasticIncrementPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& "
            "Foam::elastoPlasticDistortional::plasticIncrement()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *plasticIncrementPtr_;
}


const Foam::volScalarField&
Foam::elastoPlasticDistortional::plasticIncrement() const
{
    if (!plasticIncrementPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& "
            "Foam::elastoPlasticDistortional::plasticIncrement()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *plasticIncrementPtr_;
}


Foam::volScalarField& Foam::elastoPlasticDistortional::yieldResidual()
{
    if (!yieldResidualPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& "
            "Foam::elastoPlasticDistortional::yieldResidual()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *yieldResidualPtr_;
}


const Foam::volScalarField&
Foam::elastoPlasticDistortional::yieldResidual() const
{
    if (!yieldResidualPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& "
            "Foam::elastoPlasticDistortional::yieldResidual()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *yieldResidualPtr_;
}


Foam::devSymmTensor4thOrder
Foam::elastoPlasticDistortional::convertToTensors
(
    const tensor &T, const symmTensor &S
) const
{
    devSymmTensor4thOrder result = devSymmTensor4thOrder::zero;

    for (int i = 0; i < 9; i++)
    {
        result[i] = T[i];
    }

    for (int i = 0; i < 6; i++)
    {
        result[i + 9] = S[i];
    }

    return result;
}


Foam::tensor Foam::elastoPlasticDistortional::getTensor
(
    const devSymmTensor4thOrder& T
)
{
    return tensor(T[0],T[1],T[2],T[3],T[4],T[5],T[6],T[7],T[8]);
}


Foam::symmTensor Foam::elastoPlasticDistortional::getSymmTensor
(
    const devSymmTensor4thOrder& T
)
{
    return symmTensor(T[9],T[10],T[11],T[12],T[13],T[14]);
}


Foam::adSymmTensor7
Foam::elastoPlasticDistortional::approxExp(const adSymmTensor7& T)
{
    return symm(I + T + (T & T)/adScalar7(2.0) + (T & T & T)/adScalar7(6.0));
}


void Foam::elastoPlasticDistortional::functions
(
    List<adScalar7>& F,
    const List<adScalar7>& X,
    const parameters& P,
    updatedHistoryVariables& history,
    bool set_history
)
{
    if (debug)
    {
        Info<< "@@@@@@@@@@@@@@@@" << endl;
        Info<< "Functions call: " << endl;
        Info<< "@@@@@@@@@@@@@@@@" << endl;
    }


    adScalar7 deltaGamma = 0.0;
    adSymmTensor7 effectiveStress = adSymmTensor7::zero;

    unpack_X(X, deltaGamma, effectiveStress);

    if (debug)
    {
        Info<< "deltaGamma: " << deltaGamma << endl;
        Info<< "effectiveStress: " << effectiveStress << endl;
    }


    adScalar7 kappa = P.kappaOld + deltaGamma;
    adScalar7 Qi =
        y0_ + hIso_*kappa + (yInf_ - y0_)*(1.0 - expAD(-beta_*kappa));

    // new distortional hardening, the semi-explicit inclusion of old backstress
    // here really simplifies things
    adScalar7 magStress = sqrtAD(effectiveStress&&effectiveStress);


    adSymmTensor7 N = SphericalTensor<adScalar7>::I;
    if (magStress > SMALL)
    {
        N = (effectiveStress/magStress);
    }


    adDevSymmTensor4thOrder7 NN = N*N;

    adDevSymmTensor4thOrder7 A_dist =
        P.A_dist_old +
        Foam::operator*((expAD(-deltaGamma*cd_) - adScalar7(1.0)),
        Foam::operator*((P.A_dist_old && NN),NN));

    adDevSymmTensor4thOrder7 INN = adDevSymmTensor4thOrder7::one - NN;

    adDevSymmTensor4thOrder7 A_lat =
        P.A_lat_old
      + Foam::operator*
        (
            adScalar7((expAD(-deltaGamma*7.0*cl_)/7.0 - 1.0)
            *(INN && P.A_lat_old)),
            INN
        );

    if (debug)
    {
        Info<< "A_lat:" << A_lat << endl;
        Info<< "A_dist: " << A_dist << endl;
    }

    // need this for the return direction and yield function

    adScalar7 backPull = (N && (-ck_*(P.backStrainOld)));


    adScalar7 a = 1.0 - bd_ - bc_*backPull - bl_;
    adScalar7 b = bd_ + bc_*backPull;

    if (debug)
    {
        Info<< "a: " << a << endl;
        Info<< "b: " << b << endl;
    }

    adDevSymmTensor4thOrder7 H =
        Foam::operator*(a,adDevSymmTensor4thOrder7::one)
      + Foam::operator*(b,A_dist) + Foam::operator*(bl_,A_lat);

    // yield function and plastic flow
    // use this for handy indexing

    if (debug)
    {
        Info<< "H: " << H << endl;
    }

    adSymmTensor7 tempDirection = H&&effectiveStress;

    adScalar7 phiSqrt = sqrtAD(effectiveStress && (tempDirection) );

    if (debug)
    {
        Info<< "H: " << H << endl
            << "H && effectiveStress: "
            << adSymmTensor7(H&&effectiveStress) << endl;
    }


    adScalar7 phi = phiSqrt - Qi;

    if (debug)
    {
        Info<< "phi: " << phi << endl;
    }


    adSymmTensor7 Dp = symm(tempDirection/phiSqrt);


    // new backstrain
    adSymmTensor7 backStrain =
        (P.backStrainOld - deltaGamma*Dp)/(1 + deltaGamma*bk_*sqrtAD(Dp && Dp));

    //adSymmTensor7 backStrain = adSymmTensor7::zero;
    adSymmTensor7 backStress = -ck_*backStrain;

    // new stress
    adTensor7 FpNew = ((approxExp(deltaGamma*Dp)) & (inv(P.FpInvOld)));

    // do the function evaluations

    adSymmTensor7 Cp = symm(FpNew.T() & FpNew);

    adSymmTensor7 CpInv = inv(Cp);

    adSymmTensor7 C = symm(P.F.T() & P.F);

    adSymmTensor7 CInv = inv(C);
    adSymmTensor7 secondPiola = calculateStress(C,CInv,Cp,CpInv);
    adSymmTensor7 mandelStress =
        symm((FpNew & (CpInv & (C & secondPiola))) & FpNew.T());
    adSymmTensor7 devMandelStress =
        mandelStress - ((tr(mandelStress)/3.0)*adSymmTensor7(1,0,0,1,0,1));


    // the factor is to bring the stress, which is on the order of GPa,
    // back to around unity
    adSymmTensor7 stressEquation =
        (effectiveStress - devMandelStress + backStress)/1e6;

    if (debug)
    {
        Info<< "stressEquation: " << stressEquation << endl
            << "mandelStress: " << mandelStress << endl;
    }

    // set the functions
    F[0] = phi;
    F[1] = stressEquation.xx();
    F[2] = stressEquation.xy();
    F[3] = stressEquation.xz();
    F[4] = stressEquation.yy();
    F[5] = stressEquation.yz();
    F[6] = stressEquation.zz();


    // // also set the derived variables we want out
    if (set_history)
    {
        tensor FpNew_x = fadbadConvert7(FpNew);
        history.FpInv = inv(FpNew_x);
        history.backStrain = fadbadConvert7(backStrain);
        history.kappa = kappa.x();
        history.A_dist = fadbadConvert7(A_dist);
        history.A_lat = fadbadConvert7(A_lat);
        history.Dp = fadbadConvert7(deltaGamma*Dp);
    }

    history.residual = ((stressEquation&&stressEquation).x() + mag(phi.x()));

    if (debug)
    {
        Info<< "history.FpInv: " << history.FpInv << endl;
        Info<< "history.backStrain: " << history.backStrain << endl;
        Info<< "history.kappa: " << history.kappa << endl;
        Info<< "history.A_dist: " << history.A_dist << endl;
        Info<< "history.A_lat: " << history.A_lat << endl;
    }
}

void Foam::elastoPlasticDistortional::newtonLoop
(
    List<scalar>& X,
    const parameters& P,
    updatedHistoryVariables& history
)
{
    if (debug)
    {
        Info<< "################################################" << endl;
        Info<< "Newton loop: " << endl;
        Info<< "################################################" << endl;
    }

    List<adScalar7> F_ad;
    F_ad.setSize(7);
    List<adScalar7> X_ad;
    X_ad.setSize(7);

    Eigen::VectorXd Fv(7);
    Eigen::VectorXd Dx(7);
    Eigen::MatrixXd J(7,7);

    bool set_history = false;

    history.residual = GREAT;
    const int maxNewtonIterations = 100;

    List<scalar> X_init = X;

    for (int newtonIterI = 0; newtonIterI < maxNewtonIterations; newtonIterI++)
    {
        if (debug)
        {
            Info<< "Newton iteration: " << newtonIterI << endl
                << "X: " << X << endl;
        }

        for (int i = 0; i < 7; i++)
        {
            X_ad[i] = adScalar7(X[i]);
            X_ad[i].diff(i);
        }

        functions(F_ad, X_ad, P, history, set_history);

        if (debug)
        {
            Info<< "F_ad" << endl;
        }

        for (int i=0; i<7; i++)
        {
            Fv(i) = F_ad[i].x();
        }

        for (int i=0; i<7; i++)
        {
            for (int j=0; j<7; j++)
            {
                J(i,j) = F_ad[i].d(j);
            }
        }


        if (debug)
        {
            Info<< "Jacobian: " << endl;
            cout << J;
            Info<< endl << endl;
        }

        Dx = J.partialPivLu().solve(Fv);

        for (int j = 0; j<7; j++)
        {
            X[j] = X[j] - Dx(j);
        }

        if (set_history == true)
        {
            break;
        }

        if (history.residual < 1e-6)
        {
            set_history = true;
        }
        else if (newtonIterI == maxNewtonIterations-1)
        {
            FatalErrorIn
            (
                "void Foam::elastoPlasticDistortional::newtonLoop(...)"
            )   << "max newton iterations reached" << abort(FatalError);
        }
    }
}


Foam::scalar Foam::elastoPlasticDistortional::yieldFunction
(
    const symmTensor& devMandelStress,
    const scalar& kappa,
    const symmTensor& backStrain,
    const devSymmTensor4thOrder& A_dist,
    const devSymmTensor4thOrder& A_lat
) const
{
    if (debug)
    {
        Info<< "@@@@@@@@@@@@@@@@" << endl;
        Info<< "yieldFunction call: " << endl;
        Info<< "@@@@@@@@@@@@@@@@" << endl;
    }

    symmTensor backStress = -ck_*backStrain;
    symmTensor effectiveStress = devMandelStress - backStress;

    scalar magStress = sqrt(effectiveStress&&effectiveStress);

    symmTensor N = I;
    if (magStress > SMALL)
    {
        N = effectiveStress/magStress;
    }


    scalar backPull = N&&(backStress);


    scalar a = (1 - bd_ - bc_*(backPull) - bl_);
    scalar b = (bd_ + bc_*(backPull));

    if (debug)
    {
        Info<< "effectiveStress: " << effectiveStress << nl
            << "a: " << a << nl
            << "b: " << b << nl
            << "A_dist: " << A_dist << nl
            << "A_lat: " << A_lat << endl;
    }

    devSymmTensor4thOrder H =
        a*devSymmTensor4thOrder::one + b*A_dist + bl_*A_lat;


    if (debug)
    {
        Info<< "H: " << H << endl
            << "H && effectiveStress: " << symmTensor(H && effectiveStress)
            << endl;
    }

    scalar phiSqrd = (effectiveStress && (H && effectiveStress) );

    scalar Qi = y0_ + hIso_*kappa + (yInf_ - y0_)*(1.0 - exp(-beta_*kappa));


    scalar phi = sqrt(phiSqrd) - Qi;


    return phi;
}


void Foam::elastoPlasticDistortional::returnMap
(
    const tensor& F,
    const tensor& FpInvOld,
    const scalar& kappaOld,
    const symmTensor& backStrainOld,
    const devSymmTensor4thOrder& A_dist_old,
    const devSymmTensor4thOrder& A_lat_old,

    tensor& FpInv,
    scalar& kappa,
    symmTensor& backStrain,
    devSymmTensor4thOrder& A_dist,
    devSymmTensor4thOrder& A_lat,

    symmTensor& Dp,
    scalar& residual,
    scalar& plastic_increment,

    symmTensor& tau
)
{
    if (debug)
    {
        Info<< "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
            << endl;
        Info<< "Starting return map method" << endl;
        Info<< "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
            << endl;


        Info << "F:" << F << endl;
        Info << "FpInvOld:" << FpInvOld << endl;
        Info << "kappaOld:" << kappaOld << endl;
        Info << "backStrainOld:" << backStrainOld << endl;
        Info << "A_dist_old_1:" << getTensor(A_dist_old) << endl;
        Info << "A_dist_old_1:" << getSymmTensor(A_dist_old) << endl;
        Info << "A_lat_old_1:" << getTensor(A_lat_old) << endl;
        Info << "A_lat_old_1:" << getSymmTensor(A_lat_old) << endl;

    }

    // necessary kinematic items

    // Considering initial guess for fields
    // if (plastic_increment == 0)
    // {
    //     FpInv = FpInvOld;
    //     kappa = kappaOld;
    //     backStrain = backStrainOld;
    //     A_dist = A_dist_old;
    //     A_lat = A_lat_old;
    // }
    // else
    // {
    // }

    tensor Fp = inv(FpInvOld);
    symmTensor Cp = symm(Fp.T() & Fp);
    symmTensor CpInv = inv(Cp);
    symmTensor C = symm(F.T() & F);
    symmTensor CInv = inv(C);

    symmTensor secondPiola = calculateStress(C,CInv,Cp,CpInv);

    symmTensor mandelStress = symm((Fp & (CpInv & (C & secondPiola))) & Fp.T());
    symmTensor devMandelStress = dev(mandelStress);

    const symmTensor backStressOld = -ck_*backStrainOld;

    symmTensor effectiveStress = devMandelStress - backStressOld;


    scalar phi = yieldFunction
        (
            devMandelStress,
            kappaOld,
            backStrainOld,
            A_dist_old,
            A_lat_old
        );


    if (phi > 0.0)
    {
        if (debug)
        {
            Info<< "Yielding..." << endl;
            Info<< "yieldFunction: " << phi << endl;
            Info<< "A_dist_old: " << A_dist_old << endl;
            Info<< "A_lat_old: " << A_lat_old << endl;
        }


        // Uncomment this to allow an updated initial guess, might help
        // speed/convergence
        // but not guaranteed
        // if (plastic_increment > 0.0)
        // {
        //     Fp = inv(FpInv);
        //     Cp = symm(Fp.T() & Fp);
        //     CpInv = inv(Cp);

        //     secondPiola = calculateStress(C,CInv,Cp,CpInv);

        //     mandelStress = symm((Fp & (CpInv & (C & secondPiola))) & Fp.T());
        //     devMandelStress = dev(mandelStress);

        //     effectiveStress = devMandelStress + ck_*backStrain;
        // }


        parameters P;
        P.F = fadbadConvert7(F);
        P.FpInvOld = fadbadConvert7(FpInvOld);
        P.kappaOld = adScalar7(kappaOld);
        P.backStrainOld = fadbadConvert7(backStrainOld);
        P.A_dist_old = fadbadConvert7(A_dist_old);
        P.A_lat_old = fadbadConvert7(A_lat_old);


        updatedHistoryVariables history;
	history.FpInv = FpInvOld;
	history.kappa = kappaOld;
	history.backStrain = backStrainOld;
	history.A_dist = A_dist_old;
	history.A_lat = A_lat_old;

        List<scalar> X;
        pack_X(X,0.0,effectiveStress);
        //pack_X(X,plastic_increment,effectiveStress);


        newtonLoop(X,P,history);


        FpInv = history.FpInv;
        kappa = history.kappa;
        backStrain = history.backStrain;
        A_dist = history.A_dist;
        A_lat = history.A_lat;
        residual = history.residual;

        plastic_increment = X[0];
        if (debug)
        {
            Info<< "plastic_increment: " << plastic_increment << endl;

        }
    }
    else
    {
        // elastic
        Dp = symmTensor::zero;

        FpInv = FpInvOld;
        kappa = kappaOld;
        backStrain = backStrainOld;
        A_dist = A_dist_old;
        A_lat = A_lat_old;

        residual = 0.0;
        plastic_increment = 0.0;
    }

    // need to convert the mandel stress to tau now
    Fp = inv(FpInv);
    Cp = symm(Fp.T() & Fp);
    CpInv = inv(Cp);
    C = symm(F.T() & F);
    CInv = inv(C);

    secondPiola = calculateStress(C,CInv,Cp,CpInv);
    mandelStress = symm((Fp & (CpInv & (C & secondPiola))) & Fp.T());

    tau = symm(F & (secondPiola & F.T()));

    if (debug)
    {
        Info<< "doneso" << endl;
        Info<< "tau" << tau << endl;
        Info<< "secondPiola" << secondPiola << endl;
        Info<< "mandel: " << mandelStress << endl;
        Info<< "Fp" << Fp << endl;
        Info<< "det(Fp): " << det(Fp) << endl;

        Info<< "FpInv: " << inv(Fp) << endl;
        Info<< "A_dist: " << A_dist << endl;
        Info<< "A_lat: " << A_lat << endl;
        Info<< "backStrain: " << backStrain << endl;
        Info<< "kappa: " << kappa << endl;

        Info<< "Yield after: " <<
            yieldFunction
            (
                dev(mandelStress),
                kappa,
                backStrain,
                A_dist,
                A_lat
            ) << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::elastoPlasticDistortional::elastoPlasticDistortional
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const label lawIndex
)
:
    mechanicalLaw(name, mesh, dict, lawIndex),
    rho_(dict.lookup("rho")),
    E_
    (
      dimensionedScalar(dict.lookup("E"))
    ),
    nu_
    (
      dimensionedScalar(dict.lookup("nu"))
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
    FpInvPtr_(NULL),
    henckyStrainPtr_(NULL),
    DpPtr_(NULL),
    epsilonPEqPtr_(NULL),
    backStrainPtr_(NULL),
    activeYieldPtr_(NULL),
    yieldStressPtr_(NULL),
    aDistPtr_1(NULL),
    aDistPtr_2(NULL),
    aLattPtr_1(NULL),
    aLattPtr_2(NULL),
    y0_(readScalar(dict.lookup("sigmaY0"))),
    beta_(readScalar(dict.lookup("beta"))),
    hIso_(readScalar(dict.lookup("hIso"))),
    yInf_(readScalar(dict.lookup("sigmaYInfinity"))),
    yieldTol_(1e-6),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    ),

    // for paper comparison
    // ci_(readScalar(dict.lookup("ci"))),
    // bi_(readScalar(dict.lookup("bi"))),

    //kinematic
    ck_(readScalar(dict.lookup("ck"))),
    bk_(readScalar(dict.lookup("bk"))),

    //distortional
    cd_(readScalar(dict.lookup("cd"))),
    bd_(readScalar(dict.lookup("bd"))),
    cl_(readScalar(dict.lookup("cl"))),
    bl_(readScalar(dict.lookup("bl"))),
    bc_(readScalar(dict.lookup("bc"))),

    plasticIncrementPtr_(NULL),
    yieldResidualPtr_(NULL)
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
        dimensionedSymmTensor("one", dimless, symmTensor(I))
    );

    SetFieldPtr<symmTensor>
    (
        henckyStrainPtr_,
        "henckyStrain",
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    SetFieldPtr<scalar>
    (
        epsilonPEqPtr_,
        "epsilonPEq",
        dimensionedScalar("zero", dimless, 0.0)
    );

    SetFieldPtr<symmTensor>
    (
        backStrainPtr_,
        "backStrain",
        dimensionedSymmTensor("one", dimless, symmTensor::zero)
    );


    SetFieldPtr<scalar>
    (
        activeYieldPtr_,
        "activeYield",
        dimensionedScalar("zero", dimless, 0.0)
    );


    const tensor tensor_II = getTensor(devSymmTensor4thOrder::one);
    const symmTensor symmTensor_II = getSymmTensor(devSymmTensor4thOrder::one);

    SetFieldPtr<tensor>
    (
        aDistPtr_1,
        "aDist_1",
        dimensionedTensor("one",dimless,tensor_II)
    );

    SetFieldPtr<symmTensor>
    (
        aDistPtr_2,
        "aDist_2",
        dimensionedSymmTensor
        (
            "one",
            dimless,
            symmTensor_II
        )
    );

    SetFieldPtr<tensor>
    (
        aLattPtr_1,
        "aLatt_1",
        tensor_II
    );

    SetFieldPtr<symmTensor>
    (
        aLattPtr_2,
        "aLatt_2",
        dimensionedSymmTensor
        (
            "one",
            dimless,
            symmTensor_II
        )
    );

    SetFieldPtr<scalar>
    (
        plasticIncrementPtr_,
        "plasticIncrement",
        dimensionedScalar("zero", dimless, 0.0)
    );

    SetFieldPtr<scalar>
    (
        yieldResidualPtr_,
        "yieldResidual",
        dimensionedScalar("zero", dimless, 0.0)
    );

    SetFieldPtr<symmTensor>
    (
        yieldStressPtr_,
        "yieldStress",
        dimensionedSymmTensor("zero", dimless, symmTensor(y0_,y0_/sqrt(3.0),y0_/sqrt(3.0),y0_,y0_/sqrt(3.0),y0_))
    );


    if (debug)
    {
        Info << "Distortional hardening parameters: " << endl;
        Info<< "Boundary faces : " << endl;

        forAll(aDist_1().boundaryField(), patchI)
        {
            Info<< "patchI: " << patchI << endl;
            forAll(aDist_1().boundaryField()[patchI], faceI)
            {
                Info <<"faceI: " << faceI << endl;
                Info<< aDist_1().boundaryField()[patchI][faceI] << endl;
            }
        }

        // internal
        forAll(aDist_1().internalField(), cellI)
        {
            Info<< "cellI: " << cellI << endl;
            Info<< aDist_1().internalField()[cellI] << endl;
        }
    }

    if (debug)
    {
        Info << "After correct boundary conditions: " << endl;
        Info << "Distortional hardening parameters: " << endl;
        Info<< "Boundary faces : " << endl;

        forAll(aDist_1().boundaryField(), patchI)
        {
            Info<< "patchI: " << patchI << endl;
            forAll(aDist_1().boundaryField()[patchI], faceI)
            {
                Info<< aDist_1().boundaryField()[patchI][faceI] << endl;
            }
        }

        // internal
        forAll(aDist_1().internalField(), cellI)
        {
            Info<< "cellI: " << cellI << endl;
            Info<< aDist_1().internalField()[cellI] << endl;
        }
    }

    // // Tell fields to store their old time (this may not be necessary)
    FpInv().oldTime();
    epsilonPEq().oldTime();
    backStrain().oldTime();
    aDist_1().oldTime();
    aDist_2().oldTime();
    aLatt_1().oldTime();
    aLatt_2().oldTime();


    // The symmetry boundary conditions are being incorrectly applied to the tensor type
    // history variables after mapping fields from last pass, so in this step
    // the history variables are transferred from the internal cell to the symmetry patch

    // 2018-04-27: going to try it for all boundary conditions, just to see if its a processor
    // boundary condition problem that is causing the issue

    forAll(FpInv().boundaryField(), patchI)
      {
	// 2018-04-27: going to try it for all boundary conditions, just to see if its a processor
	// boundary condition problem that is causing the issue

	// 2018-04-27 seems to not make a difference if we do this to every patch or not

	if (mesh.boundaryMesh()[patchI].type() == "symmetryPlane")
	  {
	    forAll(FpInv().boundaryField()[patchI], faceI)
	      {
		label cellI = mesh.boundaryMesh()[patchI].faceCells()[faceI];
		FpInv().boundaryField()[patchI][faceI] = FpInv().internalField()[cellI];
		aDist_1().boundaryField()[patchI][faceI] = aDist_1().internalField()[cellI];
		aDist_2().boundaryField()[patchI][faceI] = aDist_2().internalField()[cellI];		
		aLatt_1().boundaryField()[patchI][faceI] = aLatt_1().internalField()[cellI];
		aLatt_2().boundaryField()[patchI][faceI] = aLatt_2().internalField()[cellI];
		backStrain().boundaryField()[patchI][faceI] = backStrain().internalField()[cellI];				
	      }
	  }
      }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::elastoPlasticDistortional::~elastoPlasticDistortional()
{
    CleanFieldPtr<tensor>(FpInvPtr_, "FpInv");
    CleanFieldPtr<symmTensor>(DpPtr_, "Dp");
    CleanFieldPtr<symmTensor>(henckyStrainPtr_, "henckyStrain");
    CleanFieldPtr<scalar>(epsilonPEqPtr_, "epsilonPEq");
    CleanFieldPtr<symmTensor>(backStrainPtr_, "backStrain");
    CleanFieldPtr<scalar>(activeYieldPtr_, "activeYield");
    CleanFieldPtr<tensor>(aDistPtr_1, "aDist_1");
    CleanFieldPtr<symmTensor>(aDistPtr_2, "aDist_2");
    CleanFieldPtr<tensor>(aLattPtr_1, "aLatt_1");
    CleanFieldPtr<symmTensor>(aLattPtr_2, "aLatt_2");
    CleanFieldPtr<scalar>(plasticIncrementPtr_, "plasticIncrement");
    CleanFieldPtr<scalar>(yieldResidualPtr_, "yieldResidual");
    CleanFieldPtr<symmTensor>(yieldStressPtr_, "yieldStress");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::elastoPlasticDistortional::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::elastoPlasticDistortional::E() const
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


Foam::tmp<Foam::volScalarField> Foam::elastoPlasticDistortional::nu() const
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


Foam::tmp<Foam::volScalarField>
Foam::elastoPlasticDistortional::RhieChowScaleFactor() const
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


void Foam::elastoPlasticDistortional::correct(volSymmTensorField& tau)
{
    const fvMesh& mesh = this->mesh();

    // Lookup relative deformation gradient from the solver
    const volTensorField& F = mesh.lookupObject<volTensorField>("F");

    // Store Dp as it is used to calculate the material residual
    Dp().storePrevIter();

    // Take references to internal fields
    const tensorField& FI = F.internalField();
    symmTensorField& tauI = tau.internalField();
    symmTensorField& DpI = Dp().internalField();
    scalarField& epsilonPEqI = epsilonPEq().internalField();
    tensorField& FpInvI = FpInv().internalField();
    const tensorField& FpInvOldI = FpInv().oldTime().internalField();
    const scalarField& epsilonPEqOldI = epsilonPEq().oldTime().internalField();
    //scalarField& activeYieldI = activeYield().internalField();
    scalarField& yieldResidualI = yieldResidual().internalField();
    const scalarField& curMaterialI = curMaterial().internalField();
    scalarField& plasticIncrementI = plasticIncrement().internalField();


    const symmTensorField& backStrainOldI =
        backStrain().oldTime().internalField();
    symmTensorField& backStrainI = backStrain().internalField();


    tensorField& aDist_1_I = aDist_1().internalField();
    symmTensorField& aDist_2_I = aDist_2().internalField();
    tensorField& aLatt_1_I = aLatt_1().internalField();
    symmTensorField& aLatt_2_I = aLatt_2().internalField();

    const tensorField& aDist_1_old_I = aDist_1().oldTime().internalField();
    const symmTensorField& aDist_2_old_I = aDist_2().oldTime().internalField();
    const tensorField& aLatt_1_old_I = aLatt_1().oldTime().internalField();
    const symmTensorField& aLatt_2_old_I = aLatt_2().oldTime().internalField();


    forAll(FI, cellI)
    {
        if (debug)
        {
            Info << "CellID: " << cellI << endl;
        }

        // Dont want to do the return map in the wrong material region
        if (curMaterialI[cellI] > SMALL)
        {
            devSymmTensor4thOrder aDist =
                convertToTensors(aDist_1_I[cellI],aDist_2_I[cellI]);
            const devSymmTensor4thOrder aDist_old =
                convertToTensors(aDist_1_old_I[cellI],aDist_2_old_I[cellI]);

            devSymmTensor4thOrder aLatt =
                convertToTensors(aLatt_1_I[cellI],aLatt_2_I[cellI]);
            const devSymmTensor4thOrder aLatt_old =
                convertToTensors(aLatt_1_old_I[cellI],aLatt_2_old_I[cellI]);

            returnMap
            (
                FI[cellI],
                FpInvOldI[cellI],
                epsilonPEqOldI[cellI],
                backStrainOldI[cellI],
                aDist_old,
                aLatt_old,

                FpInvI[cellI],
                epsilonPEqI[cellI],
                backStrainI[cellI],
                aDist,
                aLatt,

                DpI[cellI],
                yieldResidualI[cellI],
                plasticIncrementI[cellI],

                tauI[cellI]
            );

            aDist_1_I[cellI] = getTensor(aDist);
            aDist_2_I[cellI] = getSymmTensor(aDist);
            aLatt_1_I[cellI] = getTensor(aLatt);
            aLatt_2_I[cellI] = getSymmTensor(aLatt);
        }
    }


    forAll(F.boundaryField(), patchI)
    {
        // MC this could be the source of the fields being set to zero errors
        // PC, we update values on all patches including coupled
        // if (mesh.boundaryMesh()[patchI].type() != "empty")
        {
            const tensorField& FP = F.boundaryField()[patchI];
            symmTensorField& tauP = tau.boundaryField()[patchI];
            symmTensorField& DpP = Dp().boundaryField()[patchI];
            scalarField& epsilonPEqP = epsilonPEq().boundaryField()[patchI];
            tensorField& FpInvP = FpInv().boundaryField()[patchI];
            const tensorField& FpInvOldP =
                FpInv().oldTime().boundaryField()[patchI];
            const scalarField& epsilonPEqOldP =
                epsilonPEq().oldTime().boundaryField()[patchI];
            //scalarField& activeYieldP = activeYield().boundaryField()[patchI];
            scalarField& yieldResidualP =
                yieldResidual().boundaryField()[patchI];
            const scalarField& curMaterialP =
                curMaterial().boundaryField()[patchI];

            scalarField& plasticIncrementP =
                plasticIncrement().boundaryField()[patchI];

            const symmTensorField& backStrainOldP =
                backStrain().oldTime().boundaryField()[patchI];
            symmTensorField& backStrainP = backStrain().boundaryField()[patchI];


            tensorField& aDist_1_I = aDist_1().boundaryField()[patchI];
            symmTensorField& aDist_2_I = aDist_2().boundaryField()[patchI];
            tensorField& aLatt_1_I = aLatt_1().boundaryField()[patchI];
            symmTensorField& aLatt_2_I = aLatt_2().boundaryField()[patchI];

            const tensorField& aDist_1_old_I =
                aDist_1().oldTime().boundaryField()[patchI];
            const symmTensorField& aDist_2_old_I =
                aDist_2().oldTime().boundaryField()[patchI];
            const tensorField& aLatt_1_old_I =
                aLatt_1().oldTime().boundaryField()[patchI];
            const symmTensorField& aLatt_2_old_I =
                aLatt_2().oldTime().boundaryField()[patchI];

            forAll(FP, faceI)
            {
                if (debug)
                {
		  Info << "patchType: " << mesh.boundaryMesh()[patchI].type() << endl;
                    Info << "patchI: " << patchI << endl;
                    Info << "faceI: " << faceI << endl;
                }

                // Don't want to do the return map in the wrong material region
                if (curMaterialP[faceI] > SMALL)
                {
                    devSymmTensor4thOrder aDist =
                        convertToTensors(aDist_1_I[faceI],aDist_2_I[faceI]);
                    const devSymmTensor4thOrder aDist_old =
                        convertToTensors
                        (
                            aDist_1_old_I[faceI],aDist_2_old_I[faceI]
                        );

                    devSymmTensor4thOrder aLatt =
                        convertToTensors(aLatt_1_I[faceI],aLatt_2_I[faceI]);
                    const devSymmTensor4thOrder aLatt_old =
                        convertToTensors
                        (
                            aLatt_1_old_I[faceI],aLatt_2_old_I[faceI]
                        );

                    returnMap
                    (
                        FP[faceI],
                        FpInvOldP[faceI],
                        epsilonPEqOldP[faceI],
                        backStrainOldP[faceI],

                        aDist_old,
                        aLatt_old,

                        FpInvP[faceI],
                        epsilonPEqP[faceI],
                        backStrainP[faceI],

                        aDist,
                        aLatt,

                        DpP[faceI],
                        yieldResidualP[faceI],
                        plasticIncrementP[faceI],

                        tauP[faceI]
                    );

                    aDist_1_I[faceI] = getTensor(aDist);
                    aDist_2_I[faceI] = getSymmTensor(aDist);
                    aLatt_1_I[faceI] = getTensor(aLatt);
                    aLatt_2_I[faceI] = getSymmTensor(aLatt);
                }
            }
        }
    }

    // PC, do not correct boundary conditions, as instead we calculate values on
    // all patches above
    // aDist_1().correctBoundaryConditions();
    // aDist_2().correctBoundaryConditions();
    // aLatt_1().correctBoundaryConditions();
    // aLatt_2().correctBoundaryConditions();
    // epsilonPEq().correctBoundaryConditions();
    // backStrain().correctBoundaryConditions();
    // FpInv().correctBoundaryConditions();

    tau.correctBoundaryConditions();
}


Foam::scalar Foam::elastoPlasticDistortional::residual()
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


void Foam::elastoPlasticDistortional::correct(surfaceSymmTensorField& tau)
{
    notImplemented
    (
        "elastoPlasticDistortional::correct(surfaceSymmTensorField& tau)"
    );
}


Foam::tmp<Foam::volScalarField>
Foam::elastoPlasticDistortional::plasticDissipationRate() const
{
    const volSymmTensorField& sigmaCauchy =
        mesh().lookupObject<volSymmTensorField>("sigmaCauchy");

    // We assume 90% conversion
    // This might not be exactly correct but it is an OK estimate for now
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


Foam::scalar Foam::elastoPlasticDistortional::newDeltaT() const
{
    notImplemented
    (
        "Foam::scalar Foam::elastoPlasticDistortional::newDeltaT() const"
    );

    return mesh().time().endTime().value();
}


void Foam::elastoPlasticDistortional::updateYieldStress()
{
    // Force recalculation of curMaterial fields
    mechanicalLaw::updateYieldStress();

    // Check for FpInv determinant
    //scalar maxDet = gMax(det(FpInv().internalField()));
    //scalar minDet = gMin(det(FpInv().internalField()));
    //Info<< "max det: " << maxDet << endl;
    //Info<< "min det: " << minDet << endl;

    // Lookup relative deformation gradient from the solver
    const volTensorField& F = mesh().lookupObject<volTensorField>("F");

    // Update Hencky strain for visualisation
    henckyStrain() = 0.5*log(symm(F.T() & F));


    // Reset plasticIncrement so the trial stresses will be used next tie
    plasticIncrement() = 0.0;


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

    reduce(numCellsYielding, sumOp<int>());

    Info<< nl << "There are " << numCellsYielding << " cells actively yielding"
        << nl << endl;


    // now need to find the yield stress for all directions at each point
    // going to simply do this with a newton loop at each point
    scalar delta_stress = 1.0;

    forAll(yieldStress().internalField(), cellI)
    {

        if (curMaterial().internalField()[cellI] > SMALL)
        {
            for (int component=0; component<6; component++)
            {
                symmTensor stress_1 = symmTensor::zero;
                symmTensor stress_2 = symmTensor::zero;
                stress_1[component] = yieldStress().oldTime().internalField()[cellI][component];
                stress_2[component] = stress_1[component] + delta_stress;

                int max_its = 100;
                for (int newton_it=0; newton_it<max_its; newton_it++)
                {

                    scalar pt_1 = yieldFunction(stress_1,cellI);
                    scalar pt_2 = yieldFunction(stress_2,cellI);


                    if (pt_1 < 1e-2)
                        break;

                    scalar dPhiDStress = (pt_2 - pt_1)/delta_stress;

                    if (dPhiDStress < VSMALL)
                    {
                        break;
                    }

                    stress_1[component] -= pt_1/dPhiDStress;
                    stress_2[component] = stress_1[component] + delta_stress;

                    if(newton_it == max_its-1)
                    {
                        WarningIn
                        (
                            "elastoPlasticDistortional::updateYieldStress()"
                        )   << "Yield stress not calculated for cell: "
                            << cellI << endl;
                    }

                } // end of newton
		yieldStress().internalField()[cellI][component] =
		  stress_1[component];
            } // end of component
        }
    }
}


Foam::scalar Foam::elastoPlasticDistortional::yieldFunction
(
    const symmTensor& stress,
    const label cellI
) const
{
    // Take a reference to the mesh
    const fvMesh& mesh = this->mesh();

    // Lookup relative deformation gradient from the solver
    const tensor& F =
        mesh.lookupObject<volTensorField>("F").internalField()[cellI];
    const scalar J = det(F);
    const symmTensor C = symm(F.T() & F);

    const tensor& FpInv_local = (FpInv().internalField())[cellI];
    const tensor Fp = inv(FpInv_local);
    const symmTensor Cp = symm(Fp.T() & Fp);
    const symmTensor CpInv = inv(Cp);

    // Convert the input (cauchy) stress to the Mandel stress
    const symmTensor tau = stress*J;
    const symmTensor secondPiola = symm(inv(F) & (tau & inv(F.T())));
    const symmTensor mandelStress =
        symm((Fp & (CpInv & (C & secondPiola))) & Fp.T());

    devSymmTensor4thOrder a_dist =
        convertToTensors
        (aDist_1().internalField()[cellI],aDist_2().internalField()[cellI]);
    devSymmTensor4thOrder a_latt =
        convertToTensors
        (aLatt_1().internalField()[cellI],aLatt_2().internalField()[cellI]);

    return yieldFunction
    (
        dev(mandelStress),
        epsilonPEq().internalField()[cellI],
        backStrain().internalField()[cellI],
        a_dist,
        a_latt
    );
}


// ************************************************************************* //
