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
    levanovFriction

Description
    levanovFriction friction law.

    Shear traction is linearly related to the contact pressure via the friction
    coefficient, where the friction coefficient is specified as a time series vs
    contact pressure.

SourceFiles
    levanovFriction.C

Author
    Tian Tang Bekaert
    Peter De Jaeger Bekaert
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef levanovFriction_H
#define levanovFriction_H

#include "frictionLaw.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class levanovFriction Declaration
\*---------------------------------------------------------------------------*/

class levanovFriction
:
    public frictionLaw
{
    // Private data

         //- Reference to fvMesh
        const fvMesh& mesh_;

         //- friction law coefficientdict
        dictionary frictionLawDict_;

        //- stress field pointer
        const volSymmTensorField* sigmaCauchyPtr_;

        //- Lookup stress field
        void lookupSigmaCauchyField();

        //- Return temperature field
        const volSymmTensorField& sigmaCauchy()
        {
            if (!sigmaCauchyPtr_)
            {
                lookupSigmaCauchyField();
            }

            return *sigmaCauchyPtr_;
        }

    // Private Member Functions

        //- Surface constant taking into account the state of the roughness and
        //  the adhesion properties of the contacting surfaces. Have to be
        //  calibrated via measurements. k0 value has to be set between 0 and 1.
        const scalar k0_;

        //-Iinitial flow strength of the steel
        const scalar kf0_;

        //- Model constants with default values
        const scalar C_;
        const scalar a_;


        //- Disallow default bitwise assignment
        void operator=(const levanovFriction&);

public:

    //- Runtime type information
    TypeName("levanov");

    // Constructors

        //- Construct from dictionary
        levanovFriction
        (
            const word& name,
            const frictionContactModel& fricModel,
            const dictionary& dict
        );

        //- Construct as a copy
        levanovFriction(const levanovFriction& fricLaw);

        //- Construct and return a clone
        virtual autoPtr<frictionLaw> clone() const
        {
            return
                autoPtr<frictionLaw>
                (
                    new levanovFriction(*this)
                );
        }


    // Destructor

        virtual ~levanovFriction();


    // Member Functions



        //- Return slip traction field given the pressure and slip direction
        //  this allows for local variable friction stress
        virtual scalar slipTraction
        (
            const scalar contactPressure,         // Contact pressure
            const vector& faceSlip,               // Slip vector
            const vector& slaveFaceVelocity,      // Velocity of slave face
            const vector& masterFaceVelocity,     // Velocity of master face
            const label slavePatchIndex,          // Slave patch index
            const label faceIndex                 // Local slave face ID
        );

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
