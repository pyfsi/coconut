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
    rollerSurfaceCreatorGeometry

Description


SourceFiles
    rollerSurfaceCreatorGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef rollerSurfaceCreatorGeometry_H
#define rollerSurfaceCreatorGeometry_H

#include "rollerSurfaceCreator.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMeshGen;
class rollingMillPatchNamesHandler;

/*---------------------------------------------------------------------------*\
                Class rollerSurfaceCreatorGeometry Declaration
\*---------------------------------------------------------------------------*/

class rollerSurfaceCreatorGeometry
: public rollerSurfaceCreator
{
    // Private data
        //- characteristic points in the geometry
        LongList<point> characteristicPoints_;
        word interpolationType_;

        //- tangent vectors
        vector dirBack_;
        vector dirFront_;

    // Private member functions
        //- read the geometry file and create the triangulated surface
        void parseGeometryFile();

        //- transform the profile into the x-y plane
        //- and create the profile from the set of points
        autoPtr<splineBase> createRollerProfile() const;

        //- add additional points to the profile to capture contactAreaLength
        void cutContactArea(LongList<point>&);

        //- create triangulated surface
        void createTriangulatedSurface();

        //- Disallow default bitwise assignment
        void operator=(const rollerSurfaceCreatorGeometry&);

public:

    // Runtime type information
        TypeName("geometryFileRoller");

    // Constructors

        //- Construct from patch handler
        rollerSurfaceCreatorGeometry
        (
            const word rollPosition,
            const direction symm,
            const rollingMillPatchNamesHandler& patchHandler,
            const dictionary& dict,
            const scalar tol = 1e-6
        );

        //- Copy construct
        rollerSurfaceCreatorGeometry(const rollerSurfaceCreatorGeometry&);


    // Member Functions

        void generateGeometry();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
