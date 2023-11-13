/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Class
    checkBndFacesInCell

Description
    Checks the number of boundary faces in each cell. Removes cells
    with dominant number of boundary faces.

SourceFiles
    checkBndFacesInCell.C

\*---------------------------------------------------------------------------*/

#ifndef checkBndFacesInCell_H
#define checkBndFacesInCell_H

#include "polyMeshGenModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class checkBndFacesInCell Declaration
\*---------------------------------------------------------------------------*/

class checkBndFacesInCell
{
    // Private data
        //- Reference to polyMeshGen
        polyMeshGen& mesh_;

        //- classify each cell whether it is at the boundary or internal
        labelList cellType_;

        //- start label of boundary faces
        label startBnd_;

        //- end label of boundary faces
        label endBnd_;

    // Private member functions
        //- classify cells
        void findCellTypes();

        //- Disallow default bitwise copy construct
        checkBndFacesInCell(const checkBndFacesInCell&);

        //- Disallow default bitwise assignment
        void operator=(const checkBndFacesInCell&);

    // Private enumerators
        enum cellTypes_
        {
            NONE = 0,
            INTERNAL = 1,
            BOUNDARY = 2
        };

public:

    // Constructors

        //- Construct from polyMeshGen
        checkBndFacesInCell(polyMeshGen& mesh);

    // Destructor

        ~checkBndFacesInCell();

    // Member Functions
        //- check if cells at the boundary have a dominant number of boundary
        //- faces and marks them for removal
        bool checkCells();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
