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
    sigmaEqLaw

Description
    Cell removal law which selects cells with a equivalent stress greater than
    the specified critical equivalent stress value.

Author
    Philip Cardif UCD/UT

SourceFiles
    sigmaEqLaw.C

\*---------------------------------------------------------------------------*/

#ifndef sigmaEqLaw_H
#define sigmaEqLaw_H

#include "cellRemovalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sigmaEqLaw Declaration
\*---------------------------------------------------------------------------*/

class sigmaEqLaw
:
    public cellRemovalLaw
{
    // Private Data Member

        // Critical value for sigmaEq for cell removal
        const scalar sigmaEqCrit_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        sigmaEqLaw(const sigmaEqLaw&);

        //- Disallow default bitwise assignment
        void operator=(const sigmaEqLaw&);


public:

    //- Runtime type information
    TypeName("sigmaEq");

    // Static data members


    // Constructors

        //- Construct from dictionary
        sigmaEqLaw
        (
            const word& name,
            fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~sigmaEqLaw();


    // Member Functions

        // Update the mesh
        virtual labelList cellsToRemove();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
