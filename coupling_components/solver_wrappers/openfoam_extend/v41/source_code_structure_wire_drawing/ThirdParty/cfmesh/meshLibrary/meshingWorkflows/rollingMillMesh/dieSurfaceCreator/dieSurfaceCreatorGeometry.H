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
    dieSurfaceCreatorGeometry

Description
    Generates inner profile of a die from a geometry file

SourceFiles
    dieSurfaceCreatorGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef dieSurfaceCreatorGeometry_H
#define dieSurfaceCreatorGeometry_H

#include "dieSurfaceCreator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class dieSurfaceCreatorGeometry Declaration
\*---------------------------------------------------------------------------*/

class dieSurfaceCreatorGeometry
: public dieSurfaceCreator
{
    // Private data
        //- fileName with the geometry data
        fileName fName_;

        //- characteristic points in the geometry
        LongList<point> characteristicPoints_;
        word interpolationType_;

        //- profile points
        pointField profilePoints_;

        //- tangent in the front director
        vector dirFront_;

        //- tangent in the back direction
        vector dirBack_;

        //- inlet diameter
        scalar inletDiameter_;

        //- outlet diameter
        scalar outletDiameter_;

        //- outer diameter
        scalar outerDiameter_;

    // Private member functions
        //- read the geometry file and create the triangulated surface
        void parseGeometryFile();

        //- transform the profile into the y-z plane
        //- and create the profile from the set of points
        autoPtr<splineBase> createProfile() const;

        //- create triangulated surface
        void createTriangulatedSurface();

        //- Disallow default bitwise assignment
        void operator=(const dieSurfaceCreatorGeometry&);

public:

    //- Runtime type information
        TypeName("geometryFile");

    // Constructors

        //- construct from patch names, dictionary and tolerance
        dieSurfaceCreatorGeometry
        (
            const rollingMillPatchNamesHandler& patchNames,
            const dictionary& dict,
            const scalar tol = 1e-6
        );

        //- copy construct
        dieSurfaceCreatorGeometry
        (
            const dieSurfaceCreatorGeometry&
        );

    // Return clone of the object
        virtual autoPtr<dieSurfaceCreator> clone() const
        {
            return autoPtr<dieSurfaceCreator>
            (
                new dieSurfaceCreatorGeometry(*this)
            );
        }

    // Member Functions

        //- average diameter at the upstream side
        scalar inletDiameter() const;

        //- average diameter at the downstream side
        scalar outletDiameter() const;

        //- outer diameter of a die
        scalar outerDiameter() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
