/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    coulombFrictionAnisotropic

Description
    coulombFrictionAnisotropic friction law

    Shear traction is linearly related to the contact pressure via two friction
    coefficients.
    primaryFrictionCoeff is the friction coefficient in the primaryDirection,
    whereas secondaryFrictionCoeff is the friction coefficient in directions
    perpendicular to primaryDirection.

    The effective friction coefficient will be a mix of primaryFrictionCoeff and
    secondaryFrictionCoeff depending on the sliding direction.

SourceFiles
    coulombFrictionAnisotropic.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef coulombFrictionAnisotropic_H
#define coulombFrictionAnisotropic_H

#include "frictionLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coulombFrictionAnisotropic Declaration
\*---------------------------------------------------------------------------*/

class coulombFrictionAnisotropic
:
    public frictionLaw
{
    // Private data

         //- friction law coefficientdict
        dictionary frictionLawDict_;

        //- Coulomb friction coefficient in the primary direction
        const scalar primaryFrictionCoeff_;

        //- Coulomb friction coefficient in directions perpendicular to the
        //  primary direction
        const scalar secondaryFrictionCoeff_;

        //- Primary direction squared
        //  we store the square to save on repeat calculations of it
        tensor primaryDirectionSqr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        //coulombFrictionAnisotropic(const coulombFrictionAnisotropic&);

        //- Disallow default bitwise assignment
        void operator=(const coulombFrictionAnisotropic&);

public:

    //- Runtime type information
    TypeName("anisotropicCoulomb");

    // Static data members


    // Constructors

        //- Construct from dictionary
        coulombFrictionAnisotropic
        (
            const word& name,
            const frictionContactModel& fricModel,
            const dictionary& dict
        );

        //- Construct as a copy
        coulombFrictionAnisotropic(const coulombFrictionAnisotropic& fricLaw);

        //- Construct and return a clone
        virtual autoPtr<frictionLaw> clone() const
        {
            return
                autoPtr<frictionLaw>
                (
                    new coulombFrictionAnisotropic(*this)
                );
        }



    // Destructor

        virtual ~coulombFrictionAnisotropic();


    // Member Functions

        //- Return slip traction field
        // the traction that would cause slipping for the given pressure
        virtual scalar slipTraction(const scalar pressure);

        //- Return slip traction field given the pressure and slip direction
        //  This allows for anisotropic friction
        virtual scalar slipTraction
        (
            const scalar pressure, const vector& slipDir
        );

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
