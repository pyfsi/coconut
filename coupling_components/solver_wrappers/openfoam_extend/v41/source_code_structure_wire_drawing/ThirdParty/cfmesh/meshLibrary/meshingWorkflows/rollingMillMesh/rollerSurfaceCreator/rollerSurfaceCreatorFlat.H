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
    rollerSurfaceCreator

Description


SourceFiles
    rollerSurfaceCreator.C

\*---------------------------------------------------------------------------*/

#ifndef rollerSurfaceCreatorFlat_H
#define rollerSurfaceCreatorFlat_H

#include "rollerSurfaceCreator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class rollerSurfaceCreatorFlat Declaration
\*---------------------------------------------------------------------------*/

class rollerSurfaceCreatorFlat
: public rollerSurfaceCreator
{
    // Private data

    // Private member functions
        //- create triangulated surface
        void createTriangulatedSurface();

        //- Disallow default bitwise assignment
        void operator=(const rollerSurfaceCreatorFlat&);

public:

    // Runtime type information
        TypeName("flatRoller");

    // Constructors

        //- Construct from components
        rollerSurfaceCreatorFlat
        (
            const word rollPosition,
            const direction symm,
            const rollingMillPatchNamesHandler& patchHandler,
            const dictionary& dict,
            const scalar tol = 1e-6
        );

        //- Copy construct
        rollerSurfaceCreatorFlat(const rollerSurfaceCreatorFlat&);

    // Return clone of the object
        virtual autoPtr<rollerSurfaceCreator> clone() const
        {
            return autoPtr<rollerSurfaceCreator>
            (
                new rollerSurfaceCreatorFlat(*this)
            );
        }


    // Member Functions


    // Generate surface mesh

        void generateGeometry();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
