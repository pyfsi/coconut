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
    elasticKirchhoff

Description
    Neo-Hookean elasticity, similar to model proposed in Simo & Hughes (1998)
    in Eqn 9.2.6:

        tau = J*p*I + mu*dev(bEbar)

    where
    J        Jacobian of the deformation gradient
    p        hydrostatic pressure
    I        indentity tensor
    mu       shear modulus
    bEbar    left Cauchy Green strain tensor with the volumetric term removed

    If thermal stresses are active then 3.0*K*alpha*deltaT*I is added.

    Reference
    Simo & Hughes, Computational Inelasticity, 1998, Springer.


SourceFiles
    elasticKirchhoff.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef elasticKirchhoff_H
#define elasticKirchhoff_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class elasticKirchhoff
:
    public mechanicalLaw
{
    // Private data

        // If the case been restarted
        const bool restarted_;

        // General properties

        // Initial density
        const dimensionedScalar rho_;

        // Elastic properties

        // Young's modulus
        const dimensionedScalar E_;

        // Poisson's ratio
        const dimensionedScalar nu_;

        // Shear modulus
        const dimensionedScalar mu_;

        // Bulk modulus
        const dimensionedScalar K_;

        // Elastic left Cauchy-Green strain tensor with volumetric term removed
        mutable volSymmTensorField* bEbarPtr_;

        // Store sqrt(2/3) as it is used often
        static scalar sqrtTwoOverThree_;


    // Private Member Functions

        //- Return a reference to bEbar
        volSymmTensorField& bEbar();

        //- Return a reference to bEbar
        const volSymmTensorField& bEbar() const;

        //- Disallow default bitwise copy construct
        elasticKirchhoff(const elasticKirchhoff&);

        //- Disallow default bitwise assignment
        void operator=(const elasticKirchhoff&);

public:

    //- Runtime type information
    TypeName("elasticKirchhoff");

    // Static data members


    // Constructors

        //- Construct from dictionary
        elasticKirchhoff
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const label lawIndex
        );


    // Destructor

        virtual ~elasticKirchhoff();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E() const;

        virtual tmp<volScalarField> E(const volScalarField& epsEq) const
        {
            return E();
        }

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu() const;

        //- Return density
        virtual tmp<volScalarField> rho(scalar t) const
        {
            return rho();
        }

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E(scalar t) const
        {
            return E();
        }

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu(scalar t) const
        {
            return nu();
        }

        //- Return creep compliance
        virtual tmp<volScalarField> J(scalar t) const
        {
            notImplemented(type() + "::J(scalar t)");

            return 1.0/E(t);
        }

        //- Return plastic modulus
        virtual tmp<volScalarField> Ep() const;

        virtual tmp<volScalarField> Ep(const volScalarField& sigmaEq) const
        {
            return Ep();
        }

        //- Correct the model
        virtual void correct()
        {
            notImplemented(type() + "::correct()");
        }

        //- Update stress based on plasticity stress return vol field
        virtual void correct(volSymmTensorField& sigma);

        //- Update stress based on plasticity stress return
        virtual void correct(volSymmTensorField& sigma, const int flag);

        //- Update stress based on plasticity stress return surface field
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Update the yield stress: called at end of time-step
        virtual void updateYieldStress()
        {}

        //- Rate of dissipated energy from plasticity
        virtual tmp<volScalarField> plasticDissipationRate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
