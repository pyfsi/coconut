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
    LemaitreDamage

Description


SourceFiles
    LemaitreDamage.C

Author
    Philip Cardiff UCD
    Tian Tang Bekaert

\*---------------------------------------------------------------------------*/

#ifndef LemaitreDamage_H
#define LemaitreDamage_H

#include "damageLaw.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class LemaitreDamage
:
    public damageLaw
{
    // Private data

        // Critical damage parameter
        const dimensionedScalar h_;
        const dimensionedScalar S0_;
        const dimensionedScalar b_;
        const dimensionedScalar epsilonD_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        LemaitreDamage(const LemaitreDamage&);

        //- Disallow default bitwise assignment
        void operator=(const LemaitreDamage&);

public:

    //- Runtime type information
    TypeName("Lemaitre");

    // Static data members


    // Constructors

        //- Construct from dictionary
        LemaitreDamage
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const volScalarField& curMaterial
        );


    // Destructor

        virtual ~LemaitreDamage();


    // Member Functions

        //- Update damage field
        virtual void updateDamage();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
