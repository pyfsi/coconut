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

Namespace
    meshSurfaceDistanceFromGeometry

Description
    A class with a functions used to detect feature edges in the surface
    of the volume mesh and to detect the patches to which the boundary
    faces belong to

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceDistanceFromGeometry_H
#define meshSurfaceDistanceFromGeometry_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalarField.H"
#include "label.H"

namespace Foam
{

// Forward declarations
class polyMeshGen;
class meshSurfaceEngine;
class meshOctree;

/*---------------------------------------------------------------------------*\
            Class meshSurfaceDistanceFromGeometry Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceDistanceFromGeometry
{
    // Private data
        //- const reference to the meshSurfaceEngine
        const meshSurfaceEngine& surfaceEngine_;

        //- const reference to meshOctree
        const meshOctree& octree_;

        //- active points
        boolList activePoints_;

        //- active faces
        boolList activeFaces_;

        //- distance of surface vertices from the geometry
        scalarField boundaryPointDistance_;

        //- distance of face centres from the geometry
        scalarField faceCentreDistance_;

        //- point normal deviation
        scalarField pointNormalDeviationAngle_;

        //- face normal deviation angle
        scalarField faceNormalDeviationAngle_;

    // Private member functions
        //- mark active faces and point
        void findActiveFacesAndPoints(const wordList&);

        //- calculate point distance from the geometry
        void calculatePointDistances();

        //- calculate face distance from the geometry
        void calculateFaceCentreDistances();

        //- calculate point normal deviation
        void calculatePointNormalDeviationAngle();

        //- calculate face-normal deviation angle
        void calculateFaceNormalDeviationAngle();

        //- Disallow default bitwise copy construct
        meshSurfaceDistanceFromGeometry(const meshSurfaceDistanceFromGeometry&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceDistanceFromGeometry&);

public:

    // Constructors

        //- Construct from mesh surface and octree
        meshSurfaceDistanceFromGeometry
        (
            const meshSurfaceEngine& surfaceEngine,
            const meshOctree& octree,
            const wordList& avoidPatches = wordList()
        );


    // Destructor

        ~meshSurfaceDistanceFromGeometry();

    // Member Functions
        //- distance of boundary points from the geometry
        scalar boundaryPointDistance(const label bpI) const;

        //- distance of face centre from the geometry
        scalar boundaryFaceCentreDistance(const label bfI) const;

        //- deviation angle between the boundary point normal
        //- and the geometry normal (in radians)
        scalar angleDeviationAtBoundaryPoint(const label bpI) const;

        //- deviation angle between the boundary face normal
        //- and the geometry normal (in radians)
        scalar angleDeviationAtBoundaryFace(const label bfI) const;

        //- write data to VTK for display
        void writeToVTK(const fileName& fName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
