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

#ifndef rollerSurfaceCreator_H
#define rollerSurfaceCreator_H

#include "triSurf.H"
#include "tensor.H"
#include "splineBase.H"
#include "rollingMillMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMeshGen;

/*---------------------------------------------------------------------------*\
                    Class rollerSurfaceCreator Declaration
\*---------------------------------------------------------------------------*/

class rollerSurfaceCreator
{
protected:

    // Protected data
        //- roller patch names handler
        const rollingMillPatchNamesHandler& patchHandler_;

        //- file name (not used in all cases)
        fileName fName_;

        //- dictionary with roll data
        const dictionary& dict_;

        //- symmetry set for the the current pass
        const direction typeOfSymmetry_;

        //- tolerance for geometry creation
        scalar geometryTol_;

        //- covariance matrix
        tensor transformationMatrix_;

        //- inner diameter of a roll
        scalar innerDiameter_;

        //- outer diameter of a roll
        scalar outerDiameter_;

        //- roll width
        scalar rollWidth_;

        //- shift of the reference point from zero in the axial direction
        scalar rollAxialShift_;

        //- distance of the roll of the origin
        scalar rollRadialShift_;

        //- generated triangulated surface
        triSurf* surfPtr_;

        //- the name of the roll position (top, bottom, left, right)
        word rollPosition_;

        //- roller axis
        vector rollerAxis_;

        //- contact width centered at the origin
        scalar contactWidth_;

    // Private member functions
        //- create transformation matrix into the x-y space
        void createTransformationMatrix();

        //- merge identical vertices and triangles
        void mergeDuplicates();

        //- detect independent regions
        void detectIndependentRegions() const;

        //- cut the edges to have a hard point for contact area length
        void cutContactArea(triSurf& surf) const;

        //- position the roller such that the centre of the profile
        //- is at the origin
        void positionProfileInAxialDirection();

private:

    // Private member functions

        //- Disallow default bitwise assignment
        void operator=(const rollerSurfaceCreator&);

public:

    // Runtime type information
        TypeName("rollerProfile");

    // Constructors

        //- Construct from patch handler
        rollerSurfaceCreator
        (
            const word rollPosition,
            const direction symm,
            const rollingMillPatchNamesHandler& patchHandler,
            const dictionary& dict,
            const scalar tol = 1e-6
        );

        //- Copy construct
        rollerSurfaceCreator(const rollerSurfaceCreator&);

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            rollerSurfaceCreator,
            dictionary,
            (
                const word rollPosition,
                const rollingMillGeometryHandler::symmetryTypes_ symm,
                const rollingMillPatchNamesHandler& patchNames,
                const dictionary& dict,
                const scalar tol
            ),
            (rollPosition, symm, patchNames, dict, tol)
        );

    // Selectors

        //- Select constructed from dictionary
        static autoPtr<rollerSurfaceCreator> New
        (
            const word rollPosition,
            const rollingMillGeometryHandler::symmetryTypes_ symm,
            const rollingMillPatchNamesHandler& patchNames,
            const dictionary& dict,
            const scalar tol
        );

    // Destructor

        virtual ~rollerSurfaceCreator();


    // Member Functions
        //- get/set geometric tolerance
        const scalar& geometricTolerance() const;

        //- type of symmetry profile
        label symmetryType() const;

        //- get the position of the roll
        const word& rollPosition() const;

        //- rotation axis of the roller
        const vector& rotationAxis() const;

        //- origin of rotation
        point rotationOrigin() const;

        //- inner roll diameter
        const scalar& innerRollDiameter() const;

        //- outer diameter of the roll
        const scalar& outerRollDiameter() const;

        //- width of the roll (used for flat rolling only)
        const scalar& rollWidth() const;

        //- axial shift of the reference point from 0 in the axial direction
        const scalar& axialShift() const;

        //- radial shift of the roll from origin in the radial direction
        //- this value is read from the geometry
        const scalar& radialShift() const;

        //- width of the contact patch
        const scalar& contactWidth() const;

    // Generate surface mesh

        virtual void generateGeometry() = 0;

    // Access and transformation of geometry and meshes

        //- create a 2D surface in the y-z plane
        const triSurf& surface() const;

        //- create a 2D surface in the x-y plane used for meshing
        const triSurf* transformedSurface() const;

        //- create a 2D surface in the x-y plane used for contact evaluation
        //- it is positioned exactly where the roll shall be
        const triSurf* transformedSurfaceForContact() const;

        //- transform the mesh into requested position
        void transformToPosition(polyMeshGen& mesh) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
