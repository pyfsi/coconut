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
    casingAxialCrossSectionSurfaceCreator

Description
    Generates inner profile of a die

SourceFiles
    casingAxialCrossSectionSurfaceCreator.C

\*---------------------------------------------------------------------------*/

#ifndef casingAxialCrossSectionSurfaceCreator_H
#define casingAxialCrossSectionSurfaceCreator_H

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
            Class casingAxialCrossSectionSurfaceCreator Declaration
\*---------------------------------------------------------------------------*/

class casingAxialCrossSectionSurfaceCreator
{
protected:

    // Protected data
        //- dictionary
        const dictionary& dict_;

        //- reference to the patch names handler
        const rollingMillPatchNamesHandler& patchHandler_;

        //- tolerance for geometry creation
        scalar geometryTol_;

        //- pointer to the generated surface mesh
        triSurf* surfPtr_;

private:

    // Private member functions
        //- Disallow default bitwise copy construct
        casingAxialCrossSectionSurfaceCreator
        (
            const casingAxialCrossSectionSurfaceCreator&
        );

        //- Disallow default bitwise assignment
        void operator=(const casingAxialCrossSectionSurfaceCreator&);

public:

    // Runtime type information
        TypeName("casingAxialProfile");

    // Constructors

        //- construct from patch names, dictionary and tolerance
        casingAxialCrossSectionSurfaceCreator
        (
            const rollingMillPatchNamesHandler& patchNames,
            const dictionary& dict,
            const scalar tol = 1e-6
        );

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            casingAxialCrossSectionSurfaceCreator,
            dictionary,
            (
                const rollingMillPatchNamesHandler& patchNames,
                const dictionary& dict,
                const scalar tol
            ),
            (patchNames, dict, tol)
        );

    // Selectors

        //- Select constructed from dictionary
        static autoPtr<casingAxialCrossSectionSurfaceCreator> New
        (
            const rollingMillPatchNamesHandler& patchNames,
            const dictionary& dict,
            const scalar tol
        );

    // Destructor

        virtual ~casingAxialCrossSectionSurfaceCreator();


    // Member Functions

        //- create a 2D surface in the x-y plane
        const triSurf& surface() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
