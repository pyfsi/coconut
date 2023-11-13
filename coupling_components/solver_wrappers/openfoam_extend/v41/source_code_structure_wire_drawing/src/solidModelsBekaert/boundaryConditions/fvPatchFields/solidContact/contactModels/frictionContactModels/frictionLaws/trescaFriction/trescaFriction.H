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
    trescaFriction

Description
    trescaFriction friction law
    The slip traction is constant and equal to the shear yeild stress.

SourceFiles
    trescaFriction.C

Author
    Philip Cardiff, UCD.
    Tian Tang, Bekaert.

\*---------------------------------------------------------------------------*/

#ifndef trescaFriction_H
#define trescaFriction_H

#include "frictionLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class trescaFriction Declaration
\*---------------------------------------------------------------------------*/

class trescaFriction
:
    public frictionLaw
{
    // Private data

         //- friction law coefficientdict
        dictionary frictionLawDict_;

        //- shear stress limit
        scalar shearStressLimit_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        //trescaFriction(const trescaFriction&);

        //- Disallow default bitwise assignment
        void operator=(const trescaFriction&);

public:

    //- Runtime type information
    TypeName("tresca");

    // Static data members


    // Constructors

        //- Construct from dictionary
        trescaFriction
        (
            const word& name,
            const frictionContactModel& fricModel,
            const dictionary& dict
        );

        //- Construct as a copy
        trescaFriction(const trescaFriction& fricLaw);

        //- Construct and return a clone
        virtual autoPtr<frictionLaw> clone() const
        {
            return
                autoPtr<frictionLaw>
                (
                    new trescaFriction(*this)
                );
        }



    // Destructor

        virtual ~trescaFriction();


    // Member Functions

        //- Return slip traction field
        // the traction that would cause slipping
        // for the given pressure
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
