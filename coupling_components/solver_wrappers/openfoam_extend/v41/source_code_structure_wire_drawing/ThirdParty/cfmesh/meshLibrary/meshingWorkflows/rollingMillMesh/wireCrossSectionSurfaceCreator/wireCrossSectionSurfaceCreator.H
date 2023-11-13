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
    wireCrossSectionSurfaceCreator

Description
    Generates a ribbon from the points representing a cross section of a wire

SourceFiles
    wireCrossSectionSurfaceCreator.C

\*---------------------------------------------------------------------------*/

#ifndef wireCrossSectionSurfaceCreator_H
#define wireCrossSectionSurfaceCreator_H

#include "triSurf.H"
#include "tensor.H"
#include "splineBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMeshGen;
class rollingMillPatchNamesHandler;

/*---------------------------------------------------------------------------*\
                Class wireCrossSectionSurfaceCreator Declaration
\*---------------------------------------------------------------------------*/

class wireCrossSectionSurfaceCreator
{
    // Private data
        //- fileName with the geometry data
        fileName fName_;

        //- patch names
        const rollingMillPatchNamesHandler& patchHandler_;

        //- tolerance for geometry creation
        scalar geometryTol_;

        //- covariance matrix
        tensor transformationMatrix_;

        //- characteristic points in the geometry
        LongList<point> characteristicPoints_;
        word interpolationType_;

        //- generated triangulated surface
        triSurf* surfPtr_;

        //- check if the profile is symmetric. Store the information
        //- about the type of symmetry
        label typeOfSymmetry_;

        //- tangent in the front director
        vector dirFront_;

        //- tangent in the back direction
        vector dirBack_;

        //- shall the profile be centred at origin
        bool wireProfileCentering_;

    // Private member functions
        //- create transformation matrix into the x-y space
        void createTransformationMatrix();

        //- read the geometry file and create the triangulated surface
        void parseGeometryFile();

        //- transform the profile into the y-z plane
        //- and create the profile from the set of points
        autoPtr<splineBase> createProfile() const;

        //- create triangulated surface
        void createTriangulatedSurface();

        //- parse geometry from a DXF file
        void parseDXFFile();

        //- Disallow default bitwise copy construct
        wireCrossSectionSurfaceCreator
        (
            const wireCrossSectionSurfaceCreator&
        );

        //- Disallow default bitwise assignment
        void operator=(const wireCrossSectionSurfaceCreator&);

public:

    // Enumerators

        enum symmetryType
        {
            NONE=0,
            POSITIVE_Y = 1,
            NEGATIVE_Y = 2,
            POSITIVE_Z = 4,
            NEGATIVE_Z = 8,
            AXISYMMETRIC = 16
        };

    // Constructors

        //- Construct from fileName, tolerance and position
        wireCrossSectionSurfaceCreator
        (
            const dictionary& wireDict,
            const rollingMillPatchNamesHandler& patchHandler,
            const scalar tol = 1e-6
        );


    // Destructor

        ~wireCrossSectionSurfaceCreator();


    // Member Functions

        //- create a 2D surface in the y-z plane
        const triSurf& surface() const;

        //- create a 2D surface in the x-y plane used for meshing
        const triSurf* transformedSurface() const;

        //- snap the points onto the exact location at the parametric geometry
        void updatePointLocations(polyMeshGen& mesh) const;

        //- transform the mesh into requested position
        void transformToPosition(polyMeshGen& mesh) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
