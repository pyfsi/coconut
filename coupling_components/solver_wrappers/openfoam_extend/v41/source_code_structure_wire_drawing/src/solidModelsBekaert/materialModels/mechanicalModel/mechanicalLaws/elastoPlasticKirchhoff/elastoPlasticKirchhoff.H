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

Class
    elastoPlasticKirchhoff

Description
    Hyper-elasto-plastic mechanical law as described by Simo & Hughes (1998)
    in Box 9.1.

    Neo-Hookean elasticity with J2 plasticity.

    The yield stress versus plastic strain can expressed as a piecewise linear
    series or by an initialYieldStress and hardeningCoeff.

    Reference:
    Simo & Hughes, Computational Inelasticity, 1998, Springer.

    The hydrostatic pressure can be calculated by solving a pressure equation
    (default beahviour) or directly from the displacement field. Solving the
    pressure equation can help alleviate oscillations in the pressure field due
    to the large incompressible plastic strains.

SourceFiles
    elastoPlasticKirchhoff.C

Author
    Philip Cardiff, UCD.
    Tian Tang, Bekaert.
    Peter De Jaeger, Bekaert.

\*---------------------------------------------------------------------------*/

#ifndef elastoPlasticKirchhoff_H
#define elastoPlasticKirchhoff_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"
#include "damageLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class elastoPlasticKirchhoff
:
    public mechanicalLaw
{
    // Private data

        //- If the case been restarted
        const bool restarted_;

        //- Initial density
        const dimensionedScalar rho_;

        // Elastic properties

        //- Young's modulus
        const dimensionedScalar E_;

        //- Poisson's ratio
        const dimensionedScalar nu_;

        //- Shear modulus
        const dimensionedScalar mu_;

        //- Bulk modulus
        const dimensionedScalar K_;

        //- Table of yield (Cauchy) stress versus true plastic strain
        interpolationTable<scalar> stressPlasticStrainSeries_;

        //- Cauchy yield stress
        mutable volScalarField* sigmaYPtr_;

        //- Incremental change in sigmaY
        mutable volScalarField* DSigmaYPtr_;

        //- Total plastic strain
        mutable volSymmTensorField* epsilonPPtr_;

        //- Incremental change of plastic strain
        mutable volSymmTensorField* DEpsilonPPtr_;

        //- Elastic left Cauchy-Green trial strain tensor with volumetric term
        //  removed. Trial means that plasticity is neglected.
        mutable volSymmTensorField* bEbarTrialPtr_;

        //- Elastic left Cauchy-Green strain tensor with volumetric term removed
        mutable volSymmTensorField* bEbarPtr_;

        //- Equivalent plastic strain increment
        mutable volScalarField* DEpsilonPEqPtr_;

        //- Plastic multiplier increment - plastric strain scaled by sqrt(2/3)
        mutable volScalarField* DLambdaPtr_;

        //- Equivalent plastic strain
        mutable volScalarField* epsilonPEqPtr_;

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise
        mutable volScalarField* activeYieldPtr_;

        //- plasticN is the return direction to the yield surface
        mutable volSymmTensorField* plasticNPtr_;

        //- Hydrostatic stress
        //  Note: sigmaHyd_ == tr(sigmaCauchy)/3.0 == -p,
        //  where p is the hydrostatic pressure
        mutable volScalarField* sigmaHydPtr_;

        //- Switch to enable solution of a pressure equation, which ensures a
        //  smooth sigmaHyd field
        Switch smoothPressure_;

        // Factor to scale the amount of smoothing of the sigmaHyd field
        const scalar smoothFactor_;

        // Check if the yield stress is independent of plastic strain
        //  i.e. perfect plasticity
        bool perfectPlasticity_;

        //- An iterative procedure is used when the yield stress is a nonlinear
        //  function of plastic strain
        bool nonLinearPlasticity_;

        //- Linear plastic modulus. It is only used if plasticity is linear,
        //  defined by two points on the stress-plastic strain curve
        dimensionedScalar hardeningCoeff_;

        //- Maximum allowed error in the plastic strain integration
        const scalar maxDeltaErr_;

        //- Optional damage law
        autoPtr<damageLaw> damageLawPtr_;

        //- Tolerance for Newton loop
        static scalar LoopTol_;

        //- Maximum number of iterations for Newton loop
        static label MaxNewtonIter_;

        //- finiteDiff is the delta for finite difference differentiation
        static scalar finiteDiff_;

        //- Store sqrt(2/3) as it is used often
        static scalar sqrtTwoOverThree_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        elastoPlasticKirchhoff(const elastoPlasticKirchhoff&);

        //- Disallow default bitwise assignment
        void operator=(const elastoPlasticKirchhoff&);

        //- Return a reference to sigmaY
        volScalarField& sigmaY();

        //- Return a reference to sigmaY
        const volScalarField& sigmaY() const;

        //- Return a reference to DSigmaY
        volScalarField& DSigmaY();

        //- Return a reference to DSigmaY
        const volScalarField& DSigmaY() const;

        //- Return a reference to epsilonP
        volSymmTensorField& epsilonP();

        //- Return a reference to epsilonP
        const volSymmTensorField& epsilonP() const;

        //- Return a reference to DEpsilonP
        volSymmTensorField& DEpsilonP();

        //- Return a reference to DEpsilonP
        const volSymmTensorField& DEpsilonP() const;

        //- Return a reference to bEbarTrial
        volSymmTensorField& bEbarTrial();

        //- Return a reference to bEbarTrial
        const volSymmTensorField& bEbarTrial() const;

        //- Return a reference to bEbar
        volSymmTensorField& bEbar();

        //- Return a reference to bEbar
        const volSymmTensorField& bEbar() const;

        //- Return a reference to DEpsilonPEq
        volScalarField& DEpsilonPEq();

        //- Return a reference to DEpsilonPEq
        const volScalarField& DEpsilonPEq() const;

        //- Return a reference to DLambda
        volScalarField& DLambda();

        //- Return a reference to DLambda
        const volScalarField& DLambda() const;

        //- Return a reference to epsilonPEq
        volScalarField& epsilonPEq();

        //- Return a reference to epsilonPEq
        const volScalarField& epsilonPEq() const;

        //- Return a reference to activeYield
        volScalarField& activeYield();

        //- Return a reference to activeYield
        const volScalarField& activeYield() const;

        //- Return a reference to plasticN
        volSymmTensorField& plasticN();

        //- Return a reference to plasticN
        const volSymmTensorField& plasticN() const;

        //- Return a reference to sigmaHyd
        volScalarField& sigmaHyd();

        //- Return a reference to sigmaHyd
        const volScalarField& sigmaHyd() const;

        // Return the current Kirchhoff yield stress
        scalar curYieldStress
        (
            const scalar curEpsilonPEq,
            const scalar J,
            const label cellID
        ) const;

        // Evaulate current value of the yield function
        scalar yieldFunction
        (
            const scalar magSTrial,
            const scalar DLambda,
            const scalar muBar,
            const scalar J,
            const scalar epsilonPEq,
            const label cellID
        ) const;

        // Iteratively calculate plastic multiplier increment (DLambda)
        // and current yield stress using Newton's method
        void newtonLoop
        (
            scalar& DLambda,
            scalar& curSigmaY,
            const scalar magSTrial,
            const scalar muBar,
            const scalar J,
            const scalar epsilonPEq,
            const label cellID,
            const scalar maxMagDEpsilon
        ) const;

        //- Calcualte Ibar such that det(bEbar) == 1
        tmp<volScalarField> Ibar
        (
            const volSymmTensorField& devBEbar
        );

public:

    //- Runtime type information
    TypeName("elastoPlasticKirchhoff");

    // Static data members


    // Constructors

        //- Construct from dictionary
        elastoPlasticKirchhoff
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const label lawIndex
        );


    // Destructor

        virtual ~elastoPlasticKirchhoff();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E() const;

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu() const;

        //- Return Rhie Chow scale factor
        virtual tmp<volScalarField> RhieChowScaleFactor() const;

        //- Update stress based on plasticity stress return vol field
        virtual void correct(volSymmTensorField& sigma);

        //- Update stress based on plasticity stress return
        virtual void correct(volSymmTensorField& sigma, const int flag);

        //- Update stress based on plasticity stress return surface field
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress: called at end of time-step
        virtual void updateYieldStress();

	//- MATHIEU: Reset the yield stress: called at the beginning of a partitioned FSI-calculation
	virtual void resetYieldStress();

        //- Rate of dissipated energy from plasticity
        virtual tmp<volScalarField> plasticDissipationRate() const;

        //- Return the desired new time-step
        virtual scalar newDeltaT() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
