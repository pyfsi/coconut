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

#include "ceitMesoStressAtIntPointv2.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformGeometricField.H"
#include "fvc.H"
#include "mechanicalModel.H"
#include "thermalModel.H"
#include "logVolFields.H"
#include "ggiPolyPatch.H"
//#include "sqrtVolFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ceitMesoStressAtIntPointv2, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, ceitMesoStressAtIntPointv2, dictionary
    );

    // Declare fortran function prototypes
    extern "C"
    {
        // Note: all lowercase letters even if the fortran function has
        // uppercase letters
        //void mesostressatintpointv2_
        void mesostressatintpoint_2_22_
        (
            const double Fnew[9],
            const double Fold[9],
            const double props[11],
            const double statOld[25],
            double statNew[25],
            double sNew[6]
        );
    }
} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::ceitMesoStressAtIntPointv2::calculateStressAtPoint
(
    const tensor& Fnew,
    const tensor& Fold,
    const scalarList& props,
    const List<scalar>& statOld,
    List<scalar>& statNew,
    symmTensor& sNew
)
{

  // The implicit model expects a multidimensional array
  // need to reorder to column major indexing for fortran
    const double fortranFnew[9] =
    {
      Fnew.xx(), Fnew.yx(), Fnew.zx(),
      Fnew.xy(), Fnew.yy(), Fnew.zy(),
      Fnew.xz(), Fnew.yz(), Fnew.zz()
    };

    const double fortranFold[9] =
    {
      Fold.xx(), Fold.yx(), Fold.zx(),
      Fold.xy(), Fold.yy(), Fold.zy(),
      Fold.xz(), Fold.yz(), Fold.zz()
    };

    double fortranProps[props.size()];
    forAll(props, i)
    {
        fortranProps[i] = props[i];
    }

    double fortranStatOld[statOld.size()];
    double fortranStatNew[statNew.size()];
    forAll(statOld, i)
    {
        fortranStatOld[i] = statOld[i];
        fortranStatNew[i] = statNew[i];
    }

    double fortranSNew[6] =
    {
        sNew.xx(), sNew.yy(), sNew.zz(),
        sNew.xy(), sNew.xz(), sNew.yz()
    };

    mesostressatintpoint_2_22_
    (
        fortranFnew,
        fortranFold,
        fortranProps,
        fortranStatOld,
        fortranStatNew,
        fortranSNew
    );

    forAll(statNew, i)
    {
        statNew[i] = fortranStatNew[i];
    }


   sNew.xx() = fortranSNew[0];
    sNew.yy() = fortranSNew[1];
    sNew.zz() = fortranSNew[2];
    sNew.xy() = fortranSNew[3];
    sNew.xz() = fortranSNew[4];
    sNew.yz() = fortranSNew[5];

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::ceitMesoStressAtIntPointv2::ceitMesoStressAtIntPointv2
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
    props_(dict.lookup("props")),
    stats0_(dict.lookup("stats")),
    stats_(0),
    U_
    (
        IOobject
        (
            "Ustretch",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor(I))
    )
{
    // Initialise stats field
    stats_.setSize(25);
    forAll(stats_, fieldI)
    {
        stats_.set
        (
            fieldI,
            new volScalarField
            (
                IOobject
                (
                    "stateVariable" + Foam::name(fieldI + 1),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, stats0_[fieldI])
            )
        );

        stats_[fieldI].oldTime();
    }


    // Store U oldtime
    U_.oldTime();

    // Create fields
    // We will check if the field already exists in the object registry (i.e. if
    // it has been created by another mechanical law); if it is not found then
    // the field is created; if it is found, then we will set a pointer to it
    // using const_cast. This may not be the most elegant or safe solution but
    // it is OK for now!

    // SetFieldPtr<symmTensor>
    // (
    //     bEbarPtr_,
    //     "bEbar",
    //     dimensionedSymmTensor("one", dimless, symmTensor(I))
    // );

    // Check if the yield stress vs. plastic strain is specified by a data
    // series or by a initialYieldStress and hardening coefficient

    // On restart, some fields may be defaulted to zero when they should default
    // to I
    // if (gMin(mag(bEbar().internalField())) < SMALL)
    // {
    //     WarningIn("ceitMesoStressAtIntPointv2::ceitMesoStressAtIntPointv2()")
    //         << "Resseting zero in bEbar fields to I" << endl;

    //     symmTensorField& bEbarI = bEbar().internalField();

    //     forAll(bEbarI, cellI)
    //     {
    //         if (mag(bEbarI[cellI]) < SMALL)
    //         {
    //             bEbarI[cellI] = 1.0*I;
    //         }
    //     }

    //     forAll(bEbar().boundaryField(), patchI)
    //     {
    //         symmTensorField& bEbarP = bEbar().boundaryField()[patchI];

    //         forAll(bEbarP, faceI)
    //         {
    //             if (mag(bEbarP[faceI]) < SMALL)
    //             {
    //                 bEbarP[faceI] = 1.0*I;
    //             }
    //         }
    //     }

    //     bEbar().correctBoundaryConditions();
    // }

    // Force storage of old time for adjustable time-step calculations
    //bEbar().oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ceitMesoStressAtIntPointv2::~ceitMesoStressAtIntPointv2()
{
    //CleanFieldPtr<symmTensor>(bEbarPtr_, "bEbar");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::ceitMesoStressAtIntPointv2::rho() const
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


Foam::tmp<Foam::volScalarField> Foam::ceitMesoStressAtIntPointv2::E() const
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


Foam::tmp<Foam::volScalarField> Foam::ceitMesoStressAtIntPointv2::nu() const
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
Foam::ceitMesoStressAtIntPointv2::RhieChowScaleFactor() const
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


void Foam::ceitMesoStressAtIntPointv2::correct(volSymmTensorField& tau)
{
    const fvMesh& mesh = this->mesh();

    // Compute elastic predictor

    // Lookup deformation gradient from the solver
    const volTensorField& F = mesh.lookupObject<volTensorField>("F");

    // Required variables
    //const scalar deltaTValue = mesh.time().deltaTValue();
    //const scalar totalTime = mesh.time().value();

    // Update internal field

    // Take reference to variables
    const tensorField& FnewI = F.internalField();
    const tensorField& FoldI = F.oldTime().internalField();
    // WarningIn(type() + "::correct(...)")
    //     << "TODO: use sigmaCauchy instead of tau field" << endl;
    symmTensorField& sNewI = tau.internalField();

    forAll(sNewI, cellI)
    {
        // Copy stateOldI[cellI] to statOldCellI
        List<scalar> statOldCellI(stats_.size(), 0.0);
        List<scalar> statNewCellI(stats_.size(), 0.0);
        forAll(stats_, fieldI)
        {
	  statOldCellI[fieldI] = stats_[fieldI].oldTime()[cellI];
        }

        calculateStressAtPoint
        (
            FnewI[cellI],
            FoldI[cellI],
            props_,
            statOldCellI,
            statNewCellI,
            sNewI[cellI]
        );

        // Copy statNewCellI to stateNewI[cellI]
        forAll(stats_, fieldI)
        {
            stats_[fieldI][cellI] = statNewCellI[fieldI];
        }
    }

    // Update boundary field

    forAll(tau.boundaryField(), patchI)
    {
        // Take reference to variables
        const tensorField& FnewP = F.boundaryField()[patchI];
        const tensorField& FoldP = F.oldTime().boundaryField()[patchI];
        // TODO: use Cauchy stress instead of Kirchhoff stress
        symmTensorField& sNewP = tau.boundaryField()[patchI];

        forAll(sNewP, faceI)
        {
            // Copy stateOldP[faceI] to statOldFaceI
            List<scalar> statOldFaceI(stats_.size(), 0.0);
            List<scalar> statNewFaceI(stats_.size(), 0.0);
            forAll(stats_, fieldI)
            {
                statOldFaceI[fieldI] =
                    stats_[fieldI].oldTime().boundaryField()[patchI][faceI];
            }

            calculateStressAtPoint
            (
                FnewP[faceI],
                FoldP[faceI],
                props_,
                statOldFaceI,
                statNewFaceI,
                sNewP[faceI]
            );

            // Copy statNewFaceI to stateNewP[faceI]
            forAll(stats_, fieldI)
            {
                stats_[fieldI].boundaryField()[patchI][faceI] =
                    statNewFaceI[fieldI];
            }
        }
    }

    // WarningIn(type() + "::correct(...)")
    //     << "TODO: multiple by curMaterial field" << endl;
    //tau = curMaterial()*newTau + (1.0 - curMaterial())*tau;

    const volScalarField& J = mesh.lookupObject<volScalarField>("J");

    // const volSymmTensorField& sigma = mesh.lookupObject<volSymmTensorField>("sigmaCauchy");
    // sigma.storePrevIter();

    tau *= J;
}


void Foam::ceitMesoStressAtIntPointv2::correct
(
    volSymmTensorField& tau, const int flag
)
{
    notImplemented(type() + "::correct(...)");
}


void Foam::ceitMesoStressAtIntPointv2::correct(surfaceSymmTensorField& tau)
{
    notImplemented
    (
        "ceitMesoStressAtIntPointv2::correct(surfaceSymmTensorField& tau)"
    );
}


Foam::scalar Foam::ceitMesoStressAtIntPointv2::residual()
{
    //notImplemented(type() + "::residual()");

    return 0.0;
    // Calculate residual based on change in plastic strain increment
    // return
    //     gMax
    //     (
    //         mag
    //         (
    //             DEpsilonP().internalField()
    //           - DEpsilonP().prevIter().internalField()
    //         )
    //     )/gMax(SMALL + mag(DEpsilonP().prevIter().internalField()));

     // return
     //    gMax
     //    (
     //        mag
     //        (
     //            stats_[7].internalField()
     //          - stats_[7].prevIter().internalField()
     //        )
     //    )/gMax(SMALL + mag(stats_[].prevIter().internalField()));

}


void Foam::ceitMesoStressAtIntPointv2::updateYieldStress()
{
    //notImplemented(type() + "::updateYieldStress()");

    // Force recalculation of curMaterial fields
    mechanicalLaw::updateYieldStress();

    // Info<< nl << "Updating the yield stress" << endl;
    // sigmaY() += DSigmaY();

    // Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq()) << endl;
    // epsilonPEq() += DEpsilonPEq();
    // epsilonP() += DEpsilonP();

    // // Count cells actively yielding
    // int numCellsYielding = 0;

    // scalarField& activeYieldI = activeYield().internalField();
    // const scalarField& DEpsilonPEqI = DEpsilonPEq().internalField();

    // forAll(activeYieldI, cellI)
    // {
    //     if (DEpsilonPEqI[cellI] > SMALL)
    //     {
    //         activeYieldI[cellI] = 1.0;
    //         numCellsYielding++;
    //     }
    //     else
    //     {
    //         activeYieldI[cellI] = 0.0;
    //     }
    // }

    // reduce(numCellsYielding, sumOp<int>());

    // forAll(activeYield().boundaryField(), patchI)
    // {
    //     // No longer treat processor patches differently
    //     //if (!activeYield().boundaryField()[patchI].coupled())
    //     {
    //         scalarField& activeYieldP = activeYield().boundaryField()[patchI];
    //         const scalarField& DEpsilonPEqP =
    //             DEpsilonPEq().boundaryField()[patchI];

    //         forAll(activeYieldP, faceI)
    //         {
    //             if (DEpsilonPEqP[faceI] > SMALL)
    //             {
    //                 activeYieldP[faceI] = 1.0;
    //             }
    //             else
    //             {
    //                 activeYieldP[faceI] = 0.0;
    //             }
    //         }
    //     }
    // }

    // // No longer treat processor patches differently
    // //activeYield().correctBoundaryConditions();

    // Info<< "    " << numCellsYielding << " cells are actively yielding"
    //     << nl << endl;

    // if (damageLawPtr_.valid())
    // {
    //     // Update the damage field
    //     damageLawPtr_->updateDamage(DEpsilonPEq());

    //     Info<< "Maximum damage: " << gMax(damageLawPtr_->damage()) << nl
    //         << "Minimum damage: " << gMin(damageLawPtr_->damage()) << nl
    //         << "Average damage: "<< gAverage(damageLawPtr_->damage()) << nl
    //         << endl;
    // }

}


Foam::tmp<Foam::volScalarField>
Foam::ceitMesoStressAtIntPointv2::plasticDissipationRate() const
{
    notImplemented(type() + "::plasticDissipationRate()");

    // const volSymmTensorField& sigmaCauchy =
    //     mesh().lookupObject<volSymmTensorField>("sigmaCauchy");

    // // We assume 90% conversion
    // return tmp<volScalarField>
    // (
    //     new volScalarField
    //     (
    //         "plasticDissipationRate",
    //         max
    //         (
    //             dimensionedScalar("zero", dimForce/(dimArea*dimTime), 0.0),
    //             0.9*(sigmaCauchy && DEpsilonP())/mesh().time().deltaT()
    //         )
    //     )
    // );

    return rho();
}


Foam::scalar Foam::ceitMesoStressAtIntPointv2::newDeltaT() const
{
    notImplemented(type() + "::newDeltaT()");

    return 0.0;

    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time-step, where the return direction is given normalised deviatoric
    // strain. The error approximation is obtained using the difference between
    // the trapezoidal rule and the EUler backward method, as described in:

    // Nam-Sua Lee, Klaus-Jurgen Bathe, Error indicators and adaptive remeshing
    // in large deformation finite element analysis, Finite Elements in
    // Analysis and Design 16 (1994) 99-139.

    // // Lookup the deformation gradient
    // const volTensorField& F = mesh().lookupObject<volTensorField>("F");

    // // Calculate the true (Hencky) strain
    // volSymmTensorField epsilon = 0.5*log(symm(F.T() & F));

    // // Calculate equivalent strain, for normalisation of the error
    // volScalarField epsilonEq = sqrt((2.0/3.0)*magSqr(dev(epsilon)));

    // // Take reference to internal fields
    // const symmTensorField& DEpsilonPI = DEpsilonP().internalField();
    // const symmTensorField& plasticNIold = plasticN().oldTime().internalField();
    // const symmTensorField& plasticNIoldOld =
    //     plasticN().oldTime().oldTime().internalField();
    // const scalarField& epsilonEqI = epsilonEq.internalField();

    // // Calculate error field
    // const symmTensorField DEpsilonPErrorI =
    //     Foam::sqrt(3.0/8.0)*DEpsilonPI*mag(plasticNIold - plasticNIoldOld)
    //    /(epsilonEqI + SMALL);

    // // Max error
    // const scalar maxMagDEpsilonPErr = gMax(mag(DEpsilonPErrorI));

    // if (maxMagDEpsilonPErr > SMALL)
    // {
    //     Info<< "    " << name() << ": max error = " << maxMagDEpsilonPErr
    //         << endl;

    //     if (maxMagDEpsilonPErr > 50*maxDeltaErr_)
    //     {
    //         WarningIn
    //         (
    //             "Foam::scalar Foam::ceitMesoStressAtIntPointv2::newDeltaT() const"
    //         )   << "The error in the plastic strain is lover 50 times larger "
    //             << "than the specified value!\n    Consider starting the "
    //             << "simulation with a smaller initial time-step" << endl;
    //     }

    //     // Calculate the time-step scaling factor, where maxDeltaErr_ is the
    //     // maximum allowed error
    //     const scalar scaleFac = maxDeltaErr_/maxMagDEpsilonPErr;

    //     // Return the new time-step size
    //     return scaleFac*mesh().time().deltaTValue();
    // }

    // return mesh().time().endTime().value();
}


// ************************************************************************* //
