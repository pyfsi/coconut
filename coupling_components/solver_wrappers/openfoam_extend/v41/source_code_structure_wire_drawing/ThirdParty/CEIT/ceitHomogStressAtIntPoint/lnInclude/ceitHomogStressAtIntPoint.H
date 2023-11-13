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
    ceitHomogStressAtIntPoint

Description
    Wrapper class for CEIT homogStressAtIntPoint fortran sub-routine.
    The homogenised version of the mesoscopic model.

    Currently the procedure works by calling a specially created
    homogstressatintpoint_ fortran sub-routine.

SourceFiles
    ceitHomogStressAtIntPoint.C

Author
    Philip Cardiff, UCD.
    Peter De Jaeger, Bekaert.
    Michael Clancy, UCD.

\*---------------------------------------------------------------------------*/

#ifndef ceitHomogStressAtIntPoint_H
#define ceitHomogStressAtIntPoint_H

#include "mechanicalLaw.H"
#include "interpolationTable.H"
#include "damageLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class ceitHomogStressAtIntPoint
:
    public mechanicalLaw
{
    // Private data

        //- If the case been restarted
        const bool restarted_;

        //- Initial density
        const dimensionedScalar rho_;

        //- Material properties
        const List<scalar> props_;

        //- Initial global state variables
        // Only have lamellar spacing as a user defined
        // initial state variable for the Homogenised model
        //const List<scalar> stats0_;
        const scalar lamellarSpacing_;

        //- For the random initialisation of the vector field
        const scalar seed_;    

        //- List of internal state variable fields
        PtrList<volScalarField> stats_;

        //- Stretch tensor field
        volSymmTensorField U_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ceitHomogStressAtIntPoint(const ceitHomogStressAtIntPoint&);

        //- Disallow default bitwise assignment
        void operator=(const ceitHomogStressAtIntPoint&);

        // Wrapper function for fortran function to calculate stress at a point
        void calculateStressAtPoint
        (
            const tensor& Fnew,
            const tensor& Fold,
            const scalarList& props,
            const List<scalar>& statOld,
            List<scalar>& statNew,
            symmTensor& sNew
        );

public:

    //- Runtime type information
    TypeName("ceitHomogStressAtIntPoint");

    // Static data members


    // Constructors

        //- Construct from dictionary
        ceitHomogStressAtIntPoint
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const label lawIndex
        );


    // Destructor

        virtual ~ceitHomogStressAtIntPoint();


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
