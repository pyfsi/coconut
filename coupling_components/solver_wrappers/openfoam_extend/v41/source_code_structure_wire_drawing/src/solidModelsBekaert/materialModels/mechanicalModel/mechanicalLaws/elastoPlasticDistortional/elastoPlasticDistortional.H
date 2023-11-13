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
    Finite strain elasto-plastic model with distortional hardening based on the
    article by Bartels and Mosler (2017), International Journal of Plasticity.

SourceFiles
    elastoPlasticKirchhoff.C
    elastoPlasticKirchhoffTemplates.C

Author
    Michael Clancy, UCD.
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef elastoPlasticDistortional_H
#define elastoPlasticDistortional_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"
#include "Eigen/Dense"
#include "fadbadTensor7.H"
#include "devSymmTensor4thOrder.H"
#include "DevSymmTensor4thOrderTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Structs to hold model parameters in a compact manner
struct updatedHistoryVariables
{
    tensor FpInv;
    scalar kappa;
    symmTensor backStrain;
    devSymmTensor4thOrder A_dist;
    devSymmTensor4thOrder A_lat;
    scalar residual;
    symmTensor Dp;
};

struct parameters
{
    adTensor7 F;
    adTensor7 FpInvOld;
    adScalar7 kappaOld;
    adSymmTensor7 backStrainOld;
    adDevSymmTensor4thOrder7 A_dist_old;
    adDevSymmTensor4thOrder7 A_lat_old;
};

struct variables
{
    scalar deltaGamma;
    symmTensor effectiveStress;
};


/*---------------------------------------------------------------------------*\
                    Class Declaration
\*---------------------------------------------------------------------------*/

class elastoPlasticDistortional
:
public mechanicalLaw
{
    // Private data

        //- Initial density
        const dimensionedScalar rho_;

        // Isotropic elasticity properties

        //- Young's modulus
        const dimensionedScalar E_;

        //- Poisson's ratio
        const dimensionedScalar nu_;

        //- Shear modulus (Lame's second parameter)
        const dimensionedScalar mu_;

        //- Lame's first parameter
        const dimensionedScalar lambda_;

        //- Bulk modulus
        const dimensionedScalar K_;

        //- Inverse of the plastic deformation gradient
        mutable volTensorField* FpInvPtr_;

        //- Logarithmic strain
        mutable volSymmTensorField* henckyStrainPtr_;

        //- Increment of plastic strain
        mutable volSymmTensorField* DpPtr_;

        //- Equivalent plastic strain
        mutable volScalarField* epsilonPEqPtr_;

        //- Kinematic hardening strain
        mutable volSymmTensorField* backStrainPtr_;

        //- Indicator field for actively yielding regions (i.e. yielding in the
        //  current time-step
        //     1.0 indicates yielding
        //     0.0 indicates not yielding
        mutable volScalarField* activeYieldPtr_;

        //- Anistropic yield stress field
        //  This field is inly used for visualisation
        mutable volSymmTensorField* yieldStressPtr_;

        //- These are the distortional hardening internal variables.
        //  Since no devSymmTensor4thOrder type exists as a field,
        //  a symmTensor and tensor field are used together to contain
        //  the 15 independent components of each type of distortional hardening
        mutable volTensorField* aDistPtr_1;
        mutable volSymmTensorField* aDistPtr_2;
        mutable volTensorField* aLattPtr_1;
        mutable volSymmTensorField* aLattPtr_2;

        //- Reference yield stress
        const scalar y0_;

        // Isotropic hardening coefficient
        const scalar beta_;

        // Isotropic hardening slope
        const scalar hIso_;

        //- Asymptotic yield stress
        const scalar yInf_;

        //- Yieling tolerance
        const scalar yieldTol_;

        //- Maximum allowed error in the plastic strain integration
        const scalar maxDeltaErr_;

        //- Hardening laws

        //- isotropic hardening parameters used in the paper the
        // distortional hardening model is taken from, not currently used,
        // replaced with voce hardening law

        // const scalar ci_;
        // const scalar bi_;

        //- Kinematic hardening parameters
        const scalar ck_;
        const scalar bk_;

        //- Distortional hardening parameters for both latent and dynamic types
        const scalar cd_;
        const scalar bd_;

        const scalar cl_;
        const scalar bl_;

        //- This parameters affects the curvature (egg shaped)
        // of the yield-surface in the direction of the back stress
        const scalar bc_;


        //- History variables
        mutable volScalarField* plasticIncrementPtr_;
        mutable volScalarField* yieldResidualPtr_;

        //- Store sqrt(2/3) as it is used often
        static scalar sqrtTwoOverThree_;
        static scalar sqrtOneOverThree_;
        static scalar sqrtThreeOverTwo_;


    // Private Member Functions

        //- Return reference to epsilonPEq
        volTensorField& FpInv();

        //- Return reference to epsilonPEq
        const volTensorField& FpInv() const;

        //- Return reference to hencky strain
        volSymmTensorField& henckyStrain();
        const volSymmTensorField& henckyStrain() const;

        //- Return reference to Dp
        volSymmTensorField& Dp();

        //- Return reference to Dp
        const volSymmTensorField& Dp() const;

        //- Return reference to epsilonPEq
        volScalarField& epsilonPEq();

        //- Return reference to epsilonPEq
        const volScalarField& epsilonPEq() const;

        //- Return reference to backStrain
        volSymmTensorField& backStrain();
        const volSymmTensorField& backStrain() const;

        //- Return reference to plasticIncrement
        volScalarField& plasticIncrement();
        const volScalarField& plasticIncrement() const;

        //- Return reference to the yield function residual
        volScalarField& yieldResidual();
        const volScalarField& yieldResidual() const;

        //- Return reference to activeYield
        volScalarField& activeYield();
        const volScalarField& activeYield() const;

        //- Return reference to the distortional hardening parameters
        volTensorField& aDist_1();
        volSymmTensorField& aDist_2();
        volTensorField& aLatt_1();
        volSymmTensorField& aLatt_2();
        const volTensorField& aDist_1() const;
        const volSymmTensorField& aDist_2() const;
        const volTensorField& aLatt_1() const;
        const volSymmTensorField& aLatt_2() const;

        //- Return reference to yieldStress
        volSymmTensorField& yieldStress();
        const volSymmTensorField& yieldStress() const;


        //- These functions just fill and extract the tensor
        //  and symmTensor form devSymmtensor4thOrder objects
        devSymmTensor4thOrder convertToTensors
        (
            const tensor &T, const symmTensor &S
        ) const;
        tensor getTensor
        (
            const devSymmTensor4thOrder& T
        );
        symmTensor getSymmTensor
        (
            const devSymmTensor4thOrder& T
        );

        //- Calculate the mandel stress from the elastic deformation gradient
        //  will be isotropic, so will return a symmetric tensor
        template<class t>
        SymmTensor<t> calculateStress
        (
            const SymmTensor<t>& C,
            const SymmTensor<t>& CInv,
            const SymmTensor<t>& Cp,
            const SymmTensor<t>& CpInv
        );

        //- This exponential will be used to get the updated plastic
        //- configuration, because an exact solution would not be ammenable
        //- to automatic differentiation, a simple approximation will be used
        //- and the plastic incompressibility enforced afterwards
        adSymmTensor7 approxExp(const adSymmTensor7& T);

        //- Convenience to pack history variables into a vector
        template<class t>
        void pack_X
        (
            List<t>& X,
            const t& deltaGamma,
            const SymmTensor<t>& effectiveStress
        );

        //- Convenience to pack history variables into a vector
        template<class t>
        void unpack_X
        (
            const List<t>& X,
            t& deltaGamma,
            SymmTensor<t>& effectiveStress
        );

        //- Evaluates all the functions and puts it into a list, handy for
        //  black-box solution method (e.g. GSL solvers)
        inline void functions
        (
            List<adScalar7>& F,
            const List<adScalar7>& X,
            const parameters& P,
            updatedHistoryVariables& history,
            bool set_history
        );

        //- gsl solver not used by default, can cause issues when using the
        //  more advanced solution methods
        //  void gslSolver
        //  (
        //      List<scalar>& X_end,
        //      const List<scalar>& P_start,
        //      updatedHistoryVariables& history
        //  );

        //- Newton loop to solve return map equations
        void newtonLoop
        (
            List<scalar>& X,
            const parameters& P,
            updatedHistoryVariables& history
        );

        //- Yield-function
        scalar yieldFunction
        (
            const symmTensor& devMandelStress,
            const scalar& kappa,
            const symmTensor& backStrain,
            const devSymmTensor4thOrder& A_dist,
            const devSymmTensor4thOrder& A_lat
        ) const;

        //- Stress return map
        void returnMap
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
            scalar& plasticIncrement,

            symmTensor& tau
        );

        //- Disallow default bitwise copy construct
        elastoPlasticDistortional(const elastoPlasticDistortional&);

        //- Disallow default bitwise assignment
        void operator=(const elastoPlasticDistortional&);

    public:

        //- Runtime type information
        TypeName("elastoPlasticDistortional");

        // Constructors

        //- Construct from dictionary
        elastoPlasticDistortional
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const label lawIndex
        );

        // Destructor

        virtual ~elastoPlasticDistortional();


        // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E() const;

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu() const;

        //- Return Rhie Chow scale factor
        virtual tmp<volScalarField> RhieChowScaleFactor() const;

        //- Update stress based on mechanical law
        virtual void correct(volSymmTensorField& tau);

        //- Update stress based on mechanical law
        virtual void correct(surfaceSymmTensorField& tau);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress: called at end of time-step
        virtual void updateYieldStress();

        //- Rate of dissipated energy from plasticity
        virtual tmp<volScalarField> plasticDissipationRate() const;

        //- Return the desired new time-step
        virtual scalar newDeltaT() const;

        //- General yield function interface
        // TODO: take reference to stress
        virtual scalar yieldFunction
        (
            const symmTensor& stress,
            const label cellI
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "elastoPlasticDistortionalTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
