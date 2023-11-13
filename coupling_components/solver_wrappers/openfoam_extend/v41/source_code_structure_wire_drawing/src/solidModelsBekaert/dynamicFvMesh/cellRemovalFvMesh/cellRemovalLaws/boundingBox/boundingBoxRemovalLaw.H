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
    boundingBoxRemovalLaw

Description
    Cell removal law that removes cells when they move outside a specified
    bounding box.

    Optional: specify cellZone to only consider the removal of cells within the
    cellZone

Author
    Philip Cardif UCD

SourceFiles
    boundingBoxRemovalLaw.C

\*---------------------------------------------------------------------------*/

#ifndef boundingBoxRemovalLaw_H
#define boundingBoxRemovalLaw_H

#include "cellRemovalLaw.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class boundingBoxRemovalLaw Declaration
\*---------------------------------------------------------------------------*/

class boundingBoxRemovalLaw
:
    public cellRemovalLaw
{
    // Private Data Members

        // Bounding box
        const boundBox bb_;

        // Optional: cell zone ID. -1 if not set
        const label cellZoneID_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        boundingBoxRemovalLaw(const boundingBoxRemovalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const boundingBoxRemovalLaw&);


public:

    //- Runtime type information
    TypeName("boundingBox");

    // Static data members


    // Constructors

        //- Construct from dictionary
        boundingBoxRemovalLaw
        (
            const word& name,
            fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~boundingBoxRemovalLaw();


    // Member Functions

        // Update the mesh
        virtual labelList cellsToRemove();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
