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
    profiledDieMeshGenerator

Description
    Creates 3D meshes for rolling mill components

SourceFiles
    profiledDieMeshGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef profiledDieMeshGenerator_H
#define profiledDieMeshGenerator_H

#include "polyMeshGen.H"
#include "vectorLongList.H"

#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;
class dieGeometryInfo;
class rollingMillPatchNamesHandler;
class profiledDieGeometryInterpolator;

/*---------------------------------------------------------------------------*\
                    Class profiledDieMeshGenerator Declaration
\*---------------------------------------------------------------------------*/

class profiledDieMeshGenerator
{
    // Private data
        //- reference to mesh
        polyMeshGen& mesh_;

        //- const reference to die cross sections
        const dieGeometryInfo& dieGeom_;

        //- handler of patch names
        const rollingMillPatchNamesHandler& patchNamesHandler_;

        //- geometry interpolator
        typedef std::shared_ptr<profiledDieGeometryInterpolator> geomPtr;
        mutable geomPtr interpolatorPtr_;

        //- surface engine
        mutable std::shared_ptr<meshSurfaceEngine> msePtr_;

        //- displacements of boundary points
        vectorLongList bndDisplacements_;

    // Private member functions
        //- return const reference to the surface engine
        const meshSurfaceEngine& surfaceEngine() const;

        //- return const reference to the interpolator
        const profiledDieGeometryInterpolator& interpolator() const;

        //- detect points at the surface of the mesh that shall not be modified
        void detectLockedBndPoints(labelLongList& lockedBndPoints) const;

        //- calculate displacements of points
        void calculateDisplacements();

        //- calculate displacements of remaining points
        void updateMeshPointsFVM(const std::map<label, vector>&);

        //- calculate displacement by using the FEM mesh motion
        void updateMeshPointsFEM(const std::map<label, vector>&);

        //- set up the mesh-motion equation and apply displacements
        //- to the mesh
        void updateMeshPoints();

    // Private member functions

        //- Disallow default bitwise copy construct
        profiledDieMeshGenerator(const profiledDieMeshGenerator&);

        //- Disallow default bitwise assignment
        void operator=(const profiledDieMeshGenerator&);

public:

    // Constructors

        //- Construct from time and the region name
        profiledDieMeshGenerator
        (
            polyMeshGen& mesh,
            const dieGeometryInfo&,
            const rollingMillPatchNamesHandler& patchHandler
        );

    // Destructor

        ~profiledDieMeshGenerator();


    // Member Functions
        //- generate profiled die
        void generateProfiledDie();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
