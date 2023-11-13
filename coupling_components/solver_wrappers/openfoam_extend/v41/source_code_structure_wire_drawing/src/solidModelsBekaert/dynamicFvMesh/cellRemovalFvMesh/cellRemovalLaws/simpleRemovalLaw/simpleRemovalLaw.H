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
    simpleRemovalLaw

Description
    Simple cell removal law to test cell deletion; this law can also be used as
    as a template for new laws

Author
    Philip Cardif UCD/UT

SourceFiles
    simpleRemovalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef simpleRemovalLaw_H
#define simpleRemovalLaw_H

#include "cellRemovalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class simpleRemovalLaw Declaration
\*---------------------------------------------------------------------------*/

class simpleRemovalLaw
:
    public cellRemovalLaw
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        simpleRemovalLaw(const simpleRemovalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const simpleRemovalLaw&);


public:

    //- Runtime type information
    TypeName("simpleRemoval");

    // Static data members


    // Constructors

        //- Construct from dictionary
        simpleRemovalLaw
        (
            const word& name,
            fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~simpleRemovalLaw();


    // Member Functions

        // Update the mesh
        virtual labelList cellsToRemove();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
