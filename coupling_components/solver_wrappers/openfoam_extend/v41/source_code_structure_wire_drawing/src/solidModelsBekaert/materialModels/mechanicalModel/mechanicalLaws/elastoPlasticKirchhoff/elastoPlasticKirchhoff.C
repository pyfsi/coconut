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

#include "elastoPlasticKirchhoff.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformGeometricField.H"
#include "fvc.H"
#include "mechanicalModel.H"
#include "thermalModel.H"
#include "logVolFields.H"
#include "ggiPolyPatch.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(elastoPlasticKirchhoff, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, elastoPlasticKirchhoff, dictionary
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

        // Tolerance for Newton loop
        scalar elastoPlasticKirchhoff::LoopTol_ = 1e-8;

        // Maximum number of iterations for Newton loop
        label elastoPlasticKirchhoff::MaxNewtonIter_ = 200;

        // finiteDiff is the delta for finite difference differentiation
        scalar elastoPlasticKirchhoff::finiteDiff_ = 0.25e-6;

        // Store sqrt(2/3) as we use it often
        scalar elastoPlasticKirchhoff::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::volScalarField& Foam::elastoPlasticKirchhoff::sigmaY()
{
    if (!sigmaYPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::sigmaY()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *sigmaYPtr_;
}


const Foam::volScalarField& Foam::elastoPlasticKirchhoff::sigmaY() const
{
    if (!sigmaYPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::sigmaY()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *sigmaYPtr_;
}


Foam::volScalarField& Foam::elastoPlasticKirchhoff::DSigmaY()
{
    if (!DSigmaYPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::DSigmaY()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *DSigmaYPtr_;
}


const Foam::volScalarField& Foam::elastoPlasticKirchhoff::DSigmaY() const
{
    if (!DSigmaYPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::DSigmaY()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *DSigmaYPtr_;
}


Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::epsilonP()
{
    if (!epsilonPPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::epsilonP()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *epsilonPPtr_;
}


const Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::epsilonP() const
{
    if (!epsilonPPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::epsilonP()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *epsilonPPtr_;
}


Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::DEpsilonP()
{
    if (!DEpsilonPPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& "
            "Foam::elastoPlasticKirchhoff::DEpsilonP()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *DEpsilonPPtr_;
}


const Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::DEpsilonP() const
{
    if (!DEpsilonPPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& "
            "Foam::elastoPlasticKirchhoff::DEpsilonP()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *DEpsilonPPtr_;
}


Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::bEbarTrial()
{
    if (!bEbarTrialPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& "
            "Foam::elastoPlasticKirchhoff::bEbarTrial()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *bEbarTrialPtr_;
}


const Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::bEbarTrial() const
{
    if (!bEbarTrialPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& "
            "Foam::elastoPlasticKirchhoff::bEbarTrial()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *bEbarTrialPtr_;
}


Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::bEbar()
{
    if (!bEbarPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::bEbar()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *bEbarPtr_;
}


const Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::bEbar() const
{
    if (!bEbarPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::bEbar()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *bEbarPtr_;
}


Foam::volScalarField& Foam::elastoPlasticKirchhoff::DEpsilonPEq()
{
    if (!DEpsilonPEqPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::DEpsilonPEq()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *DEpsilonPEqPtr_;
}


const Foam::volScalarField& Foam::elastoPlasticKirchhoff::DEpsilonPEq() const
{
    if (!DEpsilonPEqPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::DEpsilonPEq()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *DEpsilonPEqPtr_;
}


Foam::volScalarField& Foam::elastoPlasticKirchhoff::DLambda()
{
    if (!DLambdaPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::DLambda()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *DLambdaPtr_;
}


const Foam::volScalarField& Foam::elastoPlasticKirchhoff::DLambda() const
{
    if (!DLambdaPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::DLambda()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *DLambdaPtr_;
}


Foam::volScalarField& Foam::elastoPlasticKirchhoff::epsilonPEq()
{
    if (!epsilonPEqPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::epsilonPEq()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *epsilonPEqPtr_;
}


const Foam::volScalarField& Foam::elastoPlasticKirchhoff::epsilonPEq() const
{
    if (!epsilonPEqPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::epsilonPEq()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *epsilonPEqPtr_;
}


Foam::volScalarField& Foam::elastoPlasticKirchhoff::activeYield()
{
    if (!activeYieldPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::activeYield()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *activeYieldPtr_;
}


const Foam::volScalarField& Foam::elastoPlasticKirchhoff::activeYield() const
{
    if (!activeYieldPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::activeYield()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *activeYieldPtr_;
}


Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::plasticN()
{
    if (!plasticNPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::plasticN()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *plasticNPtr_;
}


const Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::plasticN() const
{
    if (!plasticNPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& Foam::elastoPlasticKirchhoff::plasticN()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *plasticNPtr_;
}


Foam::volScalarField& Foam::elastoPlasticKirchhoff::sigmaHyd()
{
    if (!sigmaHydPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::sigmaHyd()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *sigmaHydPtr_;
}


const Foam::volScalarField& Foam::elastoPlasticKirchhoff::sigmaHyd() const
{
    if (!sigmaHydPtr_)
    {
        FatalErrorIn
        (
            "Foam::volScalarField& Foam::elastoPlasticKirchhoff::sigmaHyd()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *sigmaHydPtr_;
}


Foam::scalar Foam::elastoPlasticKirchhoff::curYieldStress
(
    const scalar curEpsilonPEq,    // current plastic equivalent strain
    const scalar J,                // Current Jacobian
    const label cellID             // cellID is needed is case in multiMaterial
) const
{
    // We assume that the stress-strain curve was specifed as Cauchy stress vs
    // true strain, but we want the Kirchhoff (tau) yield stress,
    // so we multiply Cauchy stress by J as tauSigmaY = J*sigmaCauchySigmaY

    return J*stressPlasticStrainSeries_(max(curEpsilonPEq, SMALL));
}


Foam::scalar Foam::elastoPlasticKirchhoff::yieldFunction
(
    const scalar magSTrial,
    const scalar DLambda,
    const scalar muBar,
    const scalar J,
    const scalar epsilonPEq,
    const label cellID
) const
{
    // Evaluate current yield function
    // fy = mag(s) - sqrt(2/3)*curSigmaY
    // fy = mag(sTrial - 2*muBar*DLambda*plasticN) - ::sqrt(2.0/3.0)*curSigmaY;
    // fy = magSTrial - 2*muBar*DLambda - ::sqrt(2.0/3.0)*curSigmaY;
    // where
    // fy is the current value of the yield function - zero at convergence.
    // s is the current deviatoric component of tau
    // sTrial is trial version of s
    // plasticN is the return direction
    // DLambda is the current increment of plastic strain multiplier
    // curSigmaY is the current Kirchhoff yield stress which is typically a
    // function of total equivalent plastic strain (epsilonPEq + DEpsilonPEq)

    return
        magSTrial - 2*muBar*DLambda
      - sqrtTwoOverThree_
           *curYieldStress
            (
                epsilonPEq + sqrtTwoOverThree_*DLambda,
                J,
                cellID
            );
}


void Foam::elastoPlasticKirchhoff::newtonLoop
(
    scalar& DLambda,
    scalar& curSigmaY,
    const scalar magSTrial,
    const scalar muBar,
    const scalar J,
    const scalar epsilonPEq,
    const label cellID,
    const scalar maxMagDEpsilon
) const
{
    // Loop to determine DEpsilonPEq
    // using Newtion's method

    int i = 0;
    scalar fTrial =
        yieldFunction(magSTrial, DLambda, muBar, J, epsilonPEq, cellID);
    scalar residual = 1.0;
    do
    {
        // Numerically calculate derivative of yield function fTrial

        // First Order
        // Hauser 2009 suggested First Order is more efficient for Newton's
        // Method as only two function evaluaitons are required.

        // fTrial step is the the value of fTrial after a small finite
        // difference step
        const scalar fTrialStep  =
            yieldFunction
            (
                magSTrial, DLambda + finiteDiff_, muBar,  J, epsilonPEq, cellID
            );
        // Numerical derivative of fTrial
        const scalar fTrialDerivative = (fTrialStep - fTrial)/finiteDiff_;

        // Update DLambda
        residual = fTrial/fTrialDerivative;
        DLambda -= residual;

        residual /= maxMagDEpsilon; // Normalise wrt max strain increment

        // fTrial will go to zero at convergence
        fTrial =
            yieldFunction(magSTrial, DLambda, muBar,  J, epsilonPEq, cellID);

        if (i == MaxNewtonIter_)
        {
            Warning
                << "Plasticity Newton loop not converging, fx is "
                << mag(fTrial) << endl;
        }
    }
    while ((mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    // Update current yield stress
    // Bug-fix PC ZT 7-Nov-14
    // Note we divide by J to change the Kirchhoff yield stress to Cauchy yield
    // stress
    curSigmaY =
        curYieldStress
        (
            epsilonPEq + sqrtTwoOverThree_*DLambda,  J, cellID
        )/J;
}


Foam::tmp<Foam::volScalarField> Foam::elastoPlasticKirchhoff::Ibar
(
    const volSymmTensorField& devBEbar
)
{
    // From Simo & Hughes 1998:
    // but this incorrectly results in det(bEbar) =! 1
    //bEbar = (s/mu) + Ibar*I;

    // A method of calculating Ibar to enforce det(bEbar) == 1 is proposed
    // by solving a cubic equation.
    // Rubin and Attia, CALCULATION OF HYPERELASTIC RESPONSE OF FINITELY
    // DEFORMED ELASTIC-VISCOPLASTIC MATERIALS, INTERNATIONAL JOURNAL FOR
    // NUMERICAL METHODS IN ENGINEERING, VOL. 39,309-320(1996)
    // and
    // M. Hollenstein M. Jabareen M. B. Rubin, Modeling a smooth elastic-
    // inelastic transition with a strongly objective numerical integrator
    // needing no iteration, Comput Mech (2013) 52:649–667
    // DOI 10.1007/s00466-013-0838-7

    // Note: In Hollenstrain et al. (2013), they suggest that Eqn 49(a) in the
    // original Rubin and Attia paper should be used.

    // Method implemented below is translated from the SmoothMultiPhase fortran
    // subroutine of Rubin

    tmp<volScalarField> tIbar
    (
        new volScalarField
        (
            IOobject
            (
                "Ibar",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    volScalarField& Ibar = tIbar();

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar.internalField();
    const symmTensorField devBEbarI = devBEbar.internalField();

    // Calculate internal field
    forAll(IbarI, cellI)
    {
        const scalar detdevBepr = det(devBEbarI[cellI]);
        const scalar dotprod = devBEbarI[cellI] && devBEbarI[cellI];
        const scalar fac1 = 2.0*dotprod/3.0;

        scalar alpha1 = 0.0;

        if (mag(fac1) < SMALL)
        {
            alpha1 = 3.0;
        }
        else
        {
            const scalar fac2 = (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

            if (fac2 >= 1.0)
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cosh(Foam::acosh(fac2)/3.0);
            }
            else
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cos(Foam::acos(fac2)/3.0);
            }
        }

        IbarI[cellI] = alpha1/3.0;
    }

    // Calculate boundary field
    forAll(Ibar.boundaryField(), patchI)
    {
        if
        (
            !Ibar.boundaryField()[patchI].coupled()
         && Ibar.boundaryField()[patchI].type() != "empty"
        )
        {
            // Take reference to patch fields for efficiency
            scalarField& IbarP = Ibar.boundaryField()[patchI];
            const symmTensorField& devBEbarP =
                devBEbar.boundaryField()[patchI];

            forAll(IbarP, faceI)
            {
                const scalar detdevBepr = det(devBEbarP[faceI]);
                const scalar dotprod =
                    devBEbarP[faceI] && devBEbarP[faceI];
                const scalar fac1 = 2.0*dotprod/3.0;

                scalar alpha1 = 0.0;

                if (mag(fac1) < SMALL)
                {
                    alpha1 = 3.0;
                }
                else
                {
                    const scalar fac2 =
                        (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

                    if (fac2 >= 1.0)
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cosh(Foam::acosh(fac2)/3.0);
                    }
                    else
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cos(Foam::acos(fac2)/3.0);
                    }
                }

                IbarP[faceI] = alpha1/3.0;
            }
        }
    }

    Ibar.correctBoundaryConditions();

    return tIbar;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::elastoPlasticKirchhoff::elastoPlasticKirchhoff
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const label lawIndex
)
:
    mechanicalLaw(name, mesh, dict, lawIndex),
    restarted_
    (
        IOobject
        (
            "sigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        ).headerOk()
    ),
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    K_
    (
        mesh.lookupObject<mechanicalModel>("mechanicalProperties").planeStress()
      ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
      : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
    ),
    stressPlasticStrainSeries_(),
    sigmaYPtr_(NULL),
    DSigmaYPtr_(NULL),
    epsilonPPtr_(NULL),
    DEpsilonPPtr_(NULL),
    bEbarTrialPtr_(NULL),
    bEbarPtr_(NULL),
    DEpsilonPEqPtr_(NULL),
    DLambdaPtr_(NULL),
    epsilonPEqPtr_(NULL),
    activeYieldPtr_(NULL),
    plasticNPtr_(NULL),
    sigmaHydPtr_(NULL),
    smoothPressure_(dict.lookupOrDefault<Switch>("smoothPressure", false)),
    smoothFactor_(dict.lookupOrDefault<scalar>("smoothFactor", 10.0)),
    perfectPlasticity_(false),
    nonLinearPlasticity_(false),
    hardeningCoeff_("hardeningCoeff", dimPressure, 0.0),
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

    SetFieldPtr<scalar>
    (
        sigmaYPtr_, "sigmaY", dimensionedScalar("zero", dimPressure, 0.0)
    );

    SetFieldPtr<scalar>
    (
        DSigmaYPtr_, "DSigmaY", dimensionedScalar("zero", dimPressure, 0.0)
    );
    DSigmaY().writeOpt() = IOobject::NO_WRITE;

    SetFieldPtr<symmTensor>
    (
        epsilonPPtr_,
        "epsilonP",
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );
    epsilonP().writeOpt() = IOobject::NO_WRITE;

    SetFieldPtr<symmTensor>
    (
        DEpsilonPPtr_,
        "DEpsilonP",
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );
    DEpsilonP().writeOpt() = IOobject::NO_WRITE;

    SetFieldPtr<symmTensor>
    (
        bEbarTrialPtr_,
        "bEbarTrial",
        dimensionedSymmTensor("one", dimless, symmTensor(1,0,0,1,0,1))
    );
    bEbarTrial().writeOpt() = IOobject::NO_WRITE;

    SetFieldPtr<symmTensor>
    (
        bEbarPtr_,
        "bEbar",
        dimensionedSymmTensor("one", dimless, symmTensor(I))
    );

    SetFieldPtr<scalar>
    (
        DEpsilonPEqPtr_,
        "DEpsilonPEq",
        dimensionedScalar("zero", dimless, 0.0)
    );
    DEpsilonPEq().writeOpt() = IOobject::NO_WRITE;

    SetFieldPtr<scalar>
    (
        DLambdaPtr_, "DLambda", dimensionedScalar("zero", dimless, 0.0)
    );
    DLambda().writeOpt() = IOobject::NO_WRITE;

    SetFieldPtr<scalar>
    (
        epsilonPEqPtr_, "epsilonPEq", dimensionedScalar("zero", dimless, 0.0)
    );

    SetFieldPtr<scalar>
    (
        activeYieldPtr_, "activeYield", dimensionedScalar("zero", dimless, 0.0)
    );

    SetFieldPtr<symmTensor>
    (
        plasticNPtr_,
        "plasticN",
        dimensionedSymmTensor("zero", dimless, symmTensor(I))
    );
    plasticN().writeOpt() = IOobject::NO_WRITE;

    SetFieldPtr<scalar>
    (
        sigmaHydPtr_,
        "sigmaHyd",
        dimensionedScalar("zero", dimPressure, 0.0),
        "zeroGradient"
    );

    // Check if the yield stress vs. plastic strain is specified by a data
    // series or by a initialYieldStress and hardening coefficient

    // Calculate the sigmaY field for the current material (mechanical law) i.e.
    // will set temporarily set sigmaY to zero in all other materials
    // Note: curMaterial field is 1's and 0's (this is NOT the same as the
    // "materials" field

    const volScalarField curMatSigmaY = curMaterial()*sigmaY();

    if
    (
        dict.found("initialYieldStress")
     && dict.found("hardeningCoefficient")
     && !dict.found("fileName")
    )
    {
        // Reset yield stress if the case has not restarted
        if (gMax(curMatSigmaY) < SMALL)
        {
            // Set initial sigmaY just in the current material
            sigmaY() +=
                curMaterial()
               *dimensionedScalar(dict.lookup("initialYieldStress"));
        }

        hardeningCoeff_ =
            dimensionedScalar(dict.lookup("hardeningCoefficient"));

        if (mag(hardeningCoeff_.value()) < SMALL)
        {
            Info<< "    Perfect plasticity" << endl;

            // Set perfect plasticity flag
            perfectPlasticity_ = true;
        }
    }
    else if
    (
        !dict.found("initialYieldStress")
     && !dict.found("hardeningCoefficient")
     && dict.found("fileName")
    )
    {
        stressPlasticStrainSeries_ = interpolationTable<scalar>(dict);

        // Reset yield stress if the case has not restarted
        if (gMax(curMatSigmaY) < SMALL)
        {
            sigmaY() +=
                curMaterial()
               *dimensionedScalar
                (
                    "initialYieldStress",
                    dimPressure,
                    stressPlasticStrainSeries_[0].second()
                );
        }

        // Check how many points are specified
        if (stressPlasticStrainSeries_.size() == 1)
        {
            Info<< "    Perfect plasticity" << endl;

            // Set perfect plasticity flag
            perfectPlasticity_ = true;

            hardeningCoeff_ =
                dimensionedScalar("hardeningCoefficient", dimPressure, 0.0);
        }
        else if (stressPlasticStrainSeries_.size() == 2)
        {
            Info<< "    Linear hardening" << endl;

            hardeningCoeff_ =
                dimensionedScalar
                (
                    "hardeningCoefficient",
                    dimPressure,
                    (
                        stressPlasticStrainSeries_[1].second()
                      - stressPlasticStrainSeries_[0].second()
                    )
                   /(
                        stressPlasticStrainSeries_[1].first()
                      - stressPlasticStrainSeries_[0].first()
                    )
                );
        }
        else
        {
            Info<< "    Nonlinear hardening" << endl;

            // Set nonlinear plasticity flag
            nonLinearPlasticity_ = true;
        }
    }
    else
    {
        FatalErrorIn
        (
            "Foam::elastoPlasticKirchhoff::elastoPlasticKirchhoff\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "The is a problem with the material definition: "
            << "EITHER 'initialYieldStress' and 'hardeningCoefficient' should"
            << " be specified OR 'fileName' and 'outOfBounds'"
            << abort(FatalError);
    }

    // Note: on restart, sigmaY will be zero in elastic material cells, but this
    // is OK because we ignore these anyway

    // On restart, some fields may be defaulted to zero when they should default
    // to I
    if (gMin(mag(bEbar().internalField())) < SMALL)
    {
        WarningIn("elastoPlasticKirchhoff::elastoPlasticKirchhoff()")
            << "Resseting zero in bEbar fields to I" << endl;

        symmTensorField& bEbarI = bEbar().internalField();

        forAll(bEbarI, cellI)
        {
            if (mag(bEbarI[cellI]) < SMALL)
            {
                bEbarI[cellI] = 1.0*I;
            }
        }

        forAll(bEbar().boundaryField(), patchI)
        {
            symmTensorField& bEbarP = bEbar().boundaryField()[patchI];

            forAll(bEbarP, faceI)
            {
                if (mag(bEbarP[faceI]) < SMALL)
                {
                    bEbarP[faceI] = 1.0*I;
                }
            }
        }

        bEbar().correctBoundaryConditions();
    }

    // Force storage of old time for adjustable time-step calculations
    plasticN().oldTime();
    bEbar().oldTime();

    // Store previous iteration to allow under-relaxation and also calculation
    // of relative residuals
    DEpsilonP().storePrevIter();

    Info<< "    smoothPressure: " << smoothPressure_ << endl;
    Info<< "    smoothFactor: " << smoothFactor_ << endl;

    // Check if a damage law is specified
    if (dict.found("damageLaw"))
    {
        damageLawPtr_ =
            damageLaw::New
            (
                "law", mesh, dict.subDict("damageLaw"), curMaterial()
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::elastoPlasticKirchhoff::~elastoPlasticKirchhoff()
{
    CleanFieldPtr<scalar>(sigmaYPtr_, "sigmaY");
    CleanFieldPtr<scalar>(DSigmaYPtr_, "DSigmaY");
    CleanFieldPtr<symmTensor>(epsilonPPtr_, "epsilonP");
    CleanFieldPtr<symmTensor>(DEpsilonPPtr_, "DEpsilonP");
    CleanFieldPtr<symmTensor>(bEbarTrialPtr_, "bEbarTrial");
    CleanFieldPtr<symmTensor>(bEbarPtr_, "bEbar");
    CleanFieldPtr<scalar>(DEpsilonPEqPtr_, "DEpsilonPEq");
    CleanFieldPtr<scalar>(DLambdaPtr_, "DLambda");
    CleanFieldPtr<scalar>(epsilonPEqPtr_, "epsilonPEq");
    CleanFieldPtr<scalar>(activeYieldPtr_, "activeYield");
    CleanFieldPtr<symmTensor>(plasticNPtr_, "plasticN");
    CleanFieldPtr<scalar>(sigmaHydPtr_, "sigmaHyd");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::elastoPlasticKirchhoff::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::elastoPlasticKirchhoff::E() const
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


Foam::tmp<Foam::volScalarField> Foam::elastoPlasticKirchhoff::nu() const
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
Foam::elastoPlasticKirchhoff::RhieChowScaleFactor() const
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


void Foam::elastoPlasticKirchhoff::correct(volSymmTensorField& tau)
{
    const fvMesh& mesh = this->mesh();

    // Compute elastic predictor

    // Lookup relative deformation gradient  from the solver
    const volTensorField& relFbar =
        mesh.lookupObject<volTensorField>("relFbar");

    // Update bE trial
    bEbarTrial() = transform(relFbar, bEbar().oldTime());

    // Lookup the Jacobian of the deformation gradient
    const volScalarField& J = mesh.lookupObject<volScalarField>("J");

    // Calculate trial deviatoric stress
    volSymmTensorField sTrial = mu_*dev(bEbarTrial());

    const volScalarField Ibar = tr(bEbarTrial())/3.0;
    const volScalarField muBar = Ibar*mu_;

    // Check for plastic loading
    // and calculate increment of plastic equivalent strain
    // i.e. the plastic multiplier

    // Store previous iteration for under-relaxation and calculation of plastic
    // residual in the solver
    DEpsilonP().storePrevIter();

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(bEbarTrial().internalField())), SMALL);

    // Trial yield function
    // sigmaY is the Cauchy yield stress so we scale it by J
    const volScalarField fTrial = mag(sTrial) - sqrtTwoOverThree_*J*sigmaY();

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticN().internalField();
    scalarField& DSigmaYI = DSigmaY().internalField();
    scalarField& DLambdaI = DLambda().internalField();
    const scalarField& muBarI = muBar.internalField();
    const scalarField& JI = J.internalField();
    const scalarField& epsilonPEqI = epsilonPEq().internalField();
    const scalarField& sigmaYI = sigmaY().internalField();
    const scalarField& curMaterialI = curMaterial().internalField();

    // Calculate DLambda() and plasticN()
    forAll(fTrialI, cellI)
    {
        // Note: we must only calculate stress in the current mechanical laws
        // cells
        if (curMaterialI[cellI] > SMALL)
        {
            // Calculate return direction plasticN
            const scalar magS = mag(sTrialI[cellI]);
            if (magS > SMALL)
            {
                plasticNI[cellI] = sTrialI[cellI]/magS;
            }

            // Calculate DLambda/DEpsilonPEq
            if (fTrialI[cellI] < SMALL)
            {
                // elastic
                DSigmaYI[cellI] = 0.0;
                DLambdaI[cellI] = 0.0;
            }
            else
            {
                if (nonLinearPlasticity_)
                {
                    // Total equivalent plastic strain where t is start of
                    // time-step
                    // This is updated in loop below
                    scalar curSigmaY = 0.0;

                    // Calculates DEpsilonPEq using Newtons's method
                    newtonLoop
                    (
                        DLambdaI[cellI],
                        curSigmaY,
                        magS,
                        muBarI[cellI],
                        JI[cellI],
                        epsilonPEqI[cellI],
                        cellI,
                        maxMagBE
                    );

                    // Update increment of yield stress
                    DSigmaYI[cellI] = curSigmaY - sigmaYI[cellI];
                }
                else
                {
                    // Plastic modulus is linear
                    DLambdaI[cellI] = fTrialI[cellI]/(2*muBarI[cellI]);

                    if (!perfectPlasticity_)
                    {
                        DLambdaI[cellI] /=
                            1.0 + hardeningCoeff_.value()/(3*muBarI[cellI]);

                        // Update increment of yield stress
                        DSigmaYI[cellI] =
                            DLambdaI[cellI]*hardeningCoeff_.value();
			//Info<<"DSigmaYI[cellI] "<< DSigmaYI[cellI] << nl <<endl;
                    }
                }
            }
        }
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // No longer treat processor patches differently
        // if
        // (
        //     !fTrial.boundaryField()[patchI].coupled()
        //  || isA<ggiPolyPatch>(mesh.boundaryMesh()[patchI])
        // )
        {
            const labelList& faceCells = mesh.boundary()[patchI].faceCells();

            // Take references to the boundary patch fields for efficiency
            const scalarField& fTrialP = fTrial.boundaryField()[patchI];
            const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
            symmTensorField& plasticNP = plasticN().boundaryField()[patchI];
            scalarField& DSigmaYP = DSigmaY().boundaryField()[patchI];
            scalarField& DLambdaP = DLambda().boundaryField()[patchI];
            const scalarField& muBarP = muBar.boundaryField()[patchI];
            const scalarField& JP = J.boundaryField()[patchI];
            const scalarField& epsilonPEqP =
                epsilonPEq().boundaryField()[patchI];
            const scalarField& sigmaYP = sigmaY().boundaryField()[patchI];
            const scalarField& curMaterialP =
                curMaterial().boundaryField()[patchI];

            forAll(fTrialP, faceI)
            {
                // Note: we must only calculate stress in the current mechanical
                // laws cells
                if (curMaterialP[faceI] > SMALL)
                {
                    // Calculate direction plasticN
                    const scalar magS = mag(sTrialP[faceI]);
                    if (magS > SMALL)
                    {
                        plasticNP[faceI] = sTrialP[faceI]/magS;
                    }

                    // Calculate DEpsilonPEq
                    if (fTrialP[faceI] < SMALL)
                    {
                        // elasticity
                        DSigmaYP[faceI] = 0.0;
                        DLambdaP[faceI] = 0.0;
                    }
                    else
                    {
                        // yielding
                        if (nonLinearPlasticity_)
                        {
                            // Updated in loop below
                            scalar curSigmaY = 0.0;

                            // Calculate DEpsilonPEq and curSigmaY
                            newtonLoop
                            (
                                DLambdaP[faceI],
                                curSigmaY,
                                magS,
                                muBarP[faceI],
                                JP[faceI],
                                epsilonPEqP[faceI],
                                faceCells[faceI],
                                maxMagBE
                            );

                            // Update increment of yield stress
                            DSigmaYP[faceI] = curSigmaY - sigmaYP[faceI];
                        }
                        else
                        {
                            // Plastic modulus is linear
                            DLambdaP[faceI] =
                                fTrialP[faceI]/(2.0*muBarP[faceI]);

                            if (!perfectPlasticity_)
                            {
                                DLambdaP[faceI] /=
                                    1.0
                                  + hardeningCoeff_.value()/(3.0*muBarP[faceI]);

                                // Update increment of yield stress
                                DSigmaYP[faceI] =
                                    DLambdaP[faceI]*hardeningCoeff_.value();
				//Info<<"DsigmaYP[faceI] " <<DSigmaYP[faceI] << nl << endl;
                            }
                        }
                    }
                }
            }
        }
    }

    // No longer treat processor patches differently
    // DSigmaY().correctBoundaryConditions();
    // DLambda().correctBoundaryConditions();
    // plasticN().correctBoundaryConditions();

    // Update DEpsilonP and DEpsilonPEq
    DEpsilonPEq() = sqrtTwoOverThree_*DLambda();
    DEpsilonP() = Ibar*DLambda()*plasticN();
    DEpsilonP().relax();

    // Calculate new Kirchhoff stress

    // Calculate deviatoric stress term
    volSymmTensorField s = sTrial - 2*mu_*DEpsilonP();

    // Calculate hydrostatic pressure term

    // Take a copy of current sigmaHyd field
    // This is a quick fix for tha case when we have multiple materials
    sigmaHyd().storePrevIter();

    if (smoothPressure_)
    {
        // Calculate the hydrostatic pressure by solving a Laplace equation;
        // this ensures smoothness of the field and quells oscillations

        // Lookup the momentum equation inverse diagonal field
        const volScalarField& AD = mesh.lookupObject<volScalarField>("DUEqnA");

        // Pressure diffusivity field: multiple by
        // (2*mu + lambda) == (4.0/3.0)*mu + K
        const surfaceScalarField rDAf
        (
            "rDAf",
            smoothFactor_*fvc::interpolate
            (
                ((4.0/3.0)*mu_ + K_)/AD, "interpolate(grad(sigmaHyd))"
            )
        );

        const dimensionedScalar one("one", dimless, 1.0);

        // Construct the pressure equation
        fvScalarMatrix sigmaHydEqn
        (
            fvm::Sp(one, sigmaHyd())
          - fvm::laplacian(rDAf, sigmaHyd(), "laplacian(DDU,DU)")
         ==
            0.5*K_*(pow(J, 2.0) - 1.0)
          - fvc::div(rDAf*fvc::interpolate(fvc::grad(sigmaHyd())) & mesh.Sf())
        );

        // Under-relax the linear system
        sigmaHydEqn.relax(); //0.7);

        // Solve the pressure equation
        sigmaHydEqn.solve();

        // Under-relax the pressure field
        sigmaHyd().relax(); //0.2);

        // Update J to remove oscillations
        const_cast<volScalarField&>(J) = sqrt((2.0*sigmaHyd()/K_) + 1.0);
        // const_cast<volScalarField&>(J) = max(min(J, 1.1), 0.5);
        // if (max(J).value() > 1.10001)
        // {
        //     FatalError<< "stop" << abort(FatalError);
        // }
        const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("relJ")) =
            J/J.oldTime();
    }
    else
    {
        // Calculate pressure as a function of displacement
        sigmaHyd() = 0.5*K_*(pow(J, 2) - 1);
    }

    // Fix for multiple materials
    sigmaHyd() =
        curMaterial()*sigmaHyd() + (1.0 - curMaterial())*sigmaHyd().prevIter();

    // Add deviatoric stress term
    volSymmTensorField newTau = s + sigmaHyd()*I;

    // Add thermal component, if active
    if (mesh.foundObject<thermalModel>("thermalProperties"))
    {
        const volScalarField& threeKalpha =
            mesh.lookupObject<volScalarField>("(threeK*alpha)");

        const volScalarField& DT = mesh.lookupObject<volScalarField>("DT");

        newTau += threeKalpha*DT*symmTensor(I);
    }

    // Scale stress if damage law is coupled
    if (damageLawPtr_.valid())
    {
        if (damageLawPtr_->coupled())
        {
            // Update the damage field
            damageLawPtr_->updateDamage(DEpsilonPEq());

            // Scale the stress field
            newTau = (1.0 - damageLawPtr_->damage())*newTau;
        }
    }

    // Assign Kirchhoff stress
    // For now, to deal with multi-materials, we will multiply by curMaterial
    // index field so only cells in the current material are calculated:
    // we should be able to do this in a better way

    tau = curMaterial()*newTau + (1.0 - curMaterial())*tau;

    // Update bEbar consistently such that bEbar == 1
    const volSymmTensorField devBEbar = (s/mu_);
    const volSymmTensorField newBEbar = devBEbar + this->Ibar(devBEbar)*I;
    //const volSymmTensorField newBEbar = (s/mu_) + Ibar*I;

    bEbar() = curMaterial()*newBEbar + (1.0 - curMaterial())*bEbar();
}


void Foam::elastoPlasticKirchhoff::correct
(
    volSymmTensorField& tau, const int flag
)
{
    if (flag == 0)
    {
        correct(tau);
    }
    else
    {
        // Update stress but do not update plasticity
        // i.e. using trial deviatoric stress sTrial

        const fvMesh& mesh = this->mesh();

        // Compute elastic predictor

        // Lookup relative deformation gradient from the solver
        const volTensorField& relFbar =
            mesh.lookupObject<volTensorField>("relFbar");

        // Update bE trial
        bEbarTrial() = transform(relFbar, bEbar().oldTime());

        // Lookup the Jacobian of the deformation gradient
        const volScalarField& J = mesh.lookupObject<volScalarField>("J");

        // Calculate trial deviatoric stress
        volSymmTensorField sTrial = mu_*dev(bEbarTrial());

        const volScalarField Ibar = tr(bEbarTrial())/3.0;
        const volScalarField muBar = Ibar*mu_;

        // Calculate new Kirchhoff stress

        // Calculate deviatoric stress term
        volSymmTensorField s = sTrial - 2*mu_*DEpsilonP();

        // Add deviatoric stress term
        volSymmTensorField newTau("newTau", s);

        // Add hydrostatic pressure term
        if (mesh.foundObject<volScalarField>("p"))
        {
            // If the pressure field exists in the solver (e.g. from the hybrid
            // approach) then we will use it
            const volScalarField& p = mesh.lookupObject<volScalarField>("p");

            newTau -= p*symmTensor(I);
        }
        else
        {
            // Otherwise calculate pressure as a function of displacement
            newTau += 0.5*K_*(pow(J, 2) - 1)*symmTensor(I);
        }

        // Add thermal component, if active
        if (mesh.foundObject<thermalModel>("thermalProperties"))
        {
            const volScalarField& threeKalpha =
                mesh.lookupObject<volScalarField>("(threeK*alpha)");

            const volScalarField& DT = mesh.lookupObject<volScalarField>("DT");

            newTau += threeKalpha*DT*symmTensor(I);
        }

        // Assign Kirchhoff stress
        // For now, to deal with multi-materials, we will multiply by
        // curMaterial
        // index field so only cells in the current material are calculated:
        // we should be able to do this in a better way

        tau = curMaterial()*newTau + (1.0 - curMaterial())*tau;

        // Update bEbar
        const volSymmTensorField newBEbar = (s/mu_) + Ibar*I;
        bEbar() = curMaterial()*newBEbar + (1.0 - curMaterial())*bEbar();
    }
}


void Foam::elastoPlasticKirchhoff::correct(surfaceSymmTensorField& tau)
{
    notImplemented
    (
        "elastoPlasticKirchhoff::correct(surfaceSymmTensorField& tau)"
    );
}


Foam::scalar Foam::elastoPlasticKirchhoff::residual()
{
    // Calculate residual based on change in plastic strain increment
    return
        gMax
        (
            mag
            (
                DEpsilonP().internalField()
              - DEpsilonP().prevIter().internalField()
            )
        )/gMax(SMALL + mag(DEpsilonP().prevIter().internalField()));
}


void Foam::elastoPlasticKirchhoff::updateYieldStress()
{
    // Force recalculation of curMaterial fields
    mechanicalLaw::updateYieldStress();

   // Info<< nl << "Updating the yield stress" << endl;
    sigmaY() += DSigmaY();
   // sigmaY() = sigmaY().oldTime() + DSigmaY();//MATHIEU

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq()) << endl;
    epsilonPEq() += DEpsilonPEq();
    epsilonP() += DEpsilonP();

    // Count cells actively yielding
    int numCellsYielding = 0;

    scalarField& activeYieldI = activeYield().internalField();
    const scalarField& DEpsilonPEqI = DEpsilonPEq().internalField();

    forAll(activeYieldI, cellI)
    {
        if (DEpsilonPEqI[cellI] > SMALL)
        {
            activeYieldI[cellI] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYieldI[cellI] = 0.0;
        }
    }

    reduce(numCellsYielding, sumOp<int>());

    forAll(activeYield().boundaryField(), patchI)
    {
        // No longer treat processor patches differently
        //if (!activeYield().boundaryField()[patchI].coupled())
        {
            scalarField& activeYieldP = activeYield().boundaryField()[patchI];
            const scalarField& DEpsilonPEqP =
                DEpsilonPEq().boundaryField()[patchI];

            forAll(activeYieldP, faceI)
            {
                if (DEpsilonPEqP[faceI] > SMALL)
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

    // No longer treat processor patches differently
    //activeYield().correctBoundaryConditions();

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;

    if (damageLawPtr_.valid())
    {
        // Update the damage field
        damageLawPtr_->updateDamage(DEpsilonPEq());

        Info<< "Maximum damage: " << gMax(damageLawPtr_->damage()) << nl
            << "Minimum damage: " << gMin(damageLawPtr_->damage()) << nl
            << "Average damage: "<< gAverage(damageLawPtr_->damage()) << nl
            << endl;
    }
}

void Foam::elastoPlasticKirchhoff::resetYieldStress()//Mathieu
{
	mechanicalLaw::resetYieldStress();
      
        sigmaY() == sigmaY().oldTime();
}

Foam::tmp<Foam::volScalarField>
Foam::elastoPlasticKirchhoff::plasticDissipationRate() const
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
                0.9*(sigmaCauchy && DEpsilonP())/mesh().time().deltaT()
            )
        )
    );
}


Foam::scalar Foam::elastoPlasticKirchhoff::newDeltaT() const
{
    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time-step, where the return direction is given normalised deviatoric
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

    // Take reference to internal fields
    const symmTensorField& DEpsilonPI = DEpsilonP().internalField();
    const symmTensorField& plasticNIold = plasticN().oldTime().internalField();
    const symmTensorField& plasticNIoldOld =
        plasticN().oldTime().oldTime().internalField();
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

        if (maxMagDEpsilonPErr > 50*maxDeltaErr_)
        {
            WarningIn
            (
                "Foam::scalar Foam::elastoPlasticKirchhoff::newDeltaT() const"
            )   << "The error in the plastic strain is lover 50 times larger "
                << "than the specified value!\n    Consider starting the "
                << "simulation with a smaller initial time-step" << endl;
        }

        // Calculate the time-step scaling factor, where maxDeltaErr_ is the
        // maximum allowed error
        const scalar scaleFac = maxDeltaErr_/maxMagDEpsilonPErr;

        // Return the new time-step size
        return scaleFac*mesh().time().deltaTValue();
    }

    return mesh().time().endTime().value();
}


// ************************************************************************* //
