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
    rollerSurfaceCreatorSurface

Description


SourceFiles
    rollerSurfaceCreatorSurface.C

\*---------------------------------------------------------------------------*/

#ifndef rollerSurfaceCreatorSurface_H
#define rollerSurfaceCreatorSurface_H

#include "rollerSurfaceCreator.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class rollerSurfaceCreatorSurface Declaration
\*---------------------------------------------------------------------------*/

class rollerSurfaceCreatorSurface
: public rollerSurfaceCreator
{
    // Private data

    // Private member functions
        //- read geometry from AutoCAD's DXF format
        //- and create the triangulated surface
        void readGeometry();

        //- Disallow default bitwise assignment
        void operator=(const rollerSurfaceCreatorSurface&);

public:

    // Runtime type information
        TypeName("surfaceFileRoller");

    // Constructors

        //- Construct from patch handler
        rollerSurfaceCreatorSurface
        (
            const word rollPosition,
            const direction symm,
            const rollingMillPatchNamesHandler& patchHandler,
            const dictionary& dict,
            const scalar tol = 1e-6
        );

        //- Copy construct
        rollerSurfaceCreatorSurface(const rollerSurfaceCreatorSurface&);


    // Member Functions

        void generateGeometry();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
