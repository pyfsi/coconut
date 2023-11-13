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
    elastoPlasticBatheHillEvolving

Description
    Large strain orthotropic elasto-plasticity.

    Large strain kinematics given by:

    Caminero, Montáns & Bathe, Modeling large strain anisotropic elasto-
    plasticity with logarithmic strain and stress measures, Computers and
    Structures, 2011, doi:10.1016/j.compstruc.2011.02.011.

    The above procedure allows the use of small-strain stress return maps in
    a large strain context; in the current implementation, the return map
    algorithm is given by:

    Mettler, Numerical Implementation of an Anisotropic Finite-Deformation
    Plasticity Model, 2012, Masters thesis, ETH, Zurich

    The work of Mettler is based on the work by:

    Papadopoulos and Lu, On the formulation and numerical solution of
    problems in anisotropic finite plasticity, Comput. Methods Appl. Mech.
    Engrg., 190  2001) 4889-4910.

SourceFiles
    elastoPlasticBatheHillEvolving.C

Author
    Michael Clancy, UCD. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef elastoPlasticBatheHillEvolving_H
#define elastoPlasticBatheHillEvolving_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"
#include "Eigen/Dense"
#include "damageLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class elastoPlasticBatheHillEvolving
:
    public mechanicalLaw
{
    // Private data

        //- Flag to switch between isotropic and orthotropic elasticity
        const Switch elasticIsotropy_;

        //- Flag to switch between isotropic and orthotropic plasticity
        const Switch plasticIsotropy_;

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

    // Orthotropic elasticity properties

        //- Young's modulii
        const scalar E1_;
        const scalar E2_;
        const scalar E3_;

        //- Shear modulii
        const scalar G12_;
        const scalar G23_;
        const scalar G31_;

        //- Poisson's ratio
        const scalar nu21_;
        const scalar nu23_;
        const scalar nu31_;
        const scalar nu12_;
        const scalar nu32_;
        const scalar nu13_;

        //- P-Lu elasticity parameters
        scalar kappa_;
        scalar mu1_;
        scalar mu2_;
        scalar mu3_;
        scalar mu4_;
        scalar mu5_;
        scalar beta1_;
        scalar beta5_;
        scalar beta6_;

        //- Tensor basis; material principal diretions
        const vector N1_;
        const vector N2_;
        const vector N3_;
        const tensor A0_;
        const tensor A1_;
        const tensor A2_;
        const tensor A3_;
        const tensor A4_;
        const tensor A5_;
        Eigen::Matrix2d Q_;

    // Hill plasticity parameters

        //- Hill yield surface parameters indicating ratio of yield stresses in
        //  the different directions
        const symmTensor RTensor0_;
        const symmTensor RTensorSat_;
        const symmTensor RTensorRate_;

        //- P-Lu plasticity parameters
        // non-const entries change due to evolving anisotropy
        scalar alpha1_;
        scalar alpha2_;
        scalar alpha3_;
        scalar alpha4_;
        scalar alpha5_;
        scalar gamma1_;

        //- P-Lu matrices for computation
        Eigen::Matrix2d DStar12_;
        Eigen::Matrix2d PStar12_;
        Eigen::Matrix2d D12_;
        Eigen::Matrix2d P12_;
        Eigen::Matrix2d Y12_;
        scalar Y3_;
        scalar Y4_;
        scalar Y5_;

        //- Hill 'R' coefficient tensor field
        mutable volSymmTensorField* RTensorPtr_;

        //- Inverse of the plastic deformation gradient
        mutable volTensorField* FpInvPtr_;

        //- Increment of plastic strain
        mutable volSymmTensorField* DpPtr_;

        //- Equivalent plastic strain
        mutable volScalarField* epsilonPEqPtr_;

        //- Indicator field for actively yielding regions (i.e. yielding in the
        //  current time-step
        //     1.0 indicates yielding
        //     0.0 indicates not yielding
        mutable volScalarField* activeYieldPtr_;

        //- Anistropic yield stress field
        //  This field is inly used for visualisation
        mutable volSymmTensorField* yieldStressPtr_;

        //- Reference yield stress
        const scalar y0_;

        // Isotropic hardening coefficient
        const scalar beta_;

        // Isotropic hardening slope
        const scalar hIso_;

        //- Asymptotic yield stress
        const scalar yInf_;

        // Yieling tolerance
        const scalar yieldTol_;

        //- Maximum allowed error in the plastic strain integration
        const scalar maxDeltaErr_;

         //- Optional damage law
        autoPtr<damageLaw> damageLawPtr_;

        // Store sqrt(2/3) as it is used often
        static scalar sqrtTwoOverThree_;

        static scalar sqrtOneOverThree_;

        // Store sqrt(3/2) as it is used often
        static scalar sqrtThreeOverTwo_;

        // Store Eigen identity as it is used often
        static Eigen::Matrix2d I2_;


    // Private Member Functions

        volSymmTensorField& RTensor();

        const volSymmTensorField& RTensor() const;


        // Return reference to epsilonPEq
        volTensorField& FpInv();

        // Return reference to epsilonPEq
        const volTensorField& FpInv() const;

        // Return reference to Dp
        volSymmTensorField& Dp();

        // Return reference to Dp
        const volSymmTensorField& Dp() const;

        // Return reference to epsilonPEq
        volScalarField& epsilonPEq();

        // Return reference to epsilonPEq
        const volScalarField& epsilonPEq() const;

        // Return reference to activeYield
        volScalarField& activeYield();

        // Return reference to yieldStress
        const volScalarField& activeYield() const;

        // Return reference to activeYield
        volSymmTensorField& yieldStress();

        // Return reference to yieldStress
        const volSymmTensorField& yieldStress() const;

        //- Calculate the stress from the strain
        symmTensor calculateStress(const symmTensor& Ee);

        //- Log of a symmTensor
        symmTensor logm(const symmTensor& T);

        //- Exponential of a symmTensor
        symmTensor expm(const symmTensor& T);

        //- Stress return map
        void smallStrainReturnMap
        (
            symmTensor& TTrial,
            tensor& FpInv,
            const tensor& FpInvOld,
            symmTensor& Dp,
            scalar& activeYield,
            scalar& epsilonPEq,
            const scalar& epsilonPEqOld,
            symmTensor& RTensor
        );

        //- Disallow default bitwise copy construct
        elastoPlasticBatheHillEvolving(const elastoPlasticBatheHillEvolving&);

        //- Disallow default bitwise assignment
        void operator=(const elastoPlasticBatheHillEvolving&);

public:

    //- Runtime type information
    TypeName("elastoPlasticBatheHillEvolving");


    // Constructors

        //- Construct from dictionary
        elastoPlasticBatheHillEvolving
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const label lawIndex
        );

    // Destructor

        virtual ~elastoPlasticBatheHillEvolving();


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
