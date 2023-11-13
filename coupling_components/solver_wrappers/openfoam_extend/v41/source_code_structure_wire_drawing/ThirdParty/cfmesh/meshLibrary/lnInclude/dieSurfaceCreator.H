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
    dieSurfaceCreator

Description
    Generates inner profile of a die

SourceFiles
    dieSurfaceCreator.C

\*---------------------------------------------------------------------------*/

#ifndef dieSurfaceCreator_H
#define dieSurfaceCreator_H

#include "triSurf.H"
#include "tensor.H"
#include "splineBase.H"

#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMeshGen;
class rollingMillPatchNamesHandler;

/*---------------------------------------------------------------------------*\
            Class dieSurfaceCreator Declaration
\*---------------------------------------------------------------------------*/

class dieSurfaceCreator
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

        //- cross-section profiles in case of profiled dies
        DynList<std::pair<scalar, std::shared_ptr<triSurf> > > crossSections_;

protected:

    // Protected member functions
        //- creates a profile with a given diameter
        std::shared_ptr<triSurf> createCircularProfile(const scalar d) const;

        //- position cross-sections of profiled dies at the origin
        void positionCrossSectionAtOrigin(triSurf&) const;

private:

    // Private member functions


        //- Disallow default bitwise copy construct
        dieSurfaceCreator
        (
            const dieSurfaceCreator&
        );

        //- Disallow default bitwise assignment
        void operator=(const dieSurfaceCreator&);

public:

    // Runtime type information
        TypeName("dieAxialProfile");

    // Constructors

        //- construct from patch names, dictionary and tolerance
        dieSurfaceCreator
        (
            const rollingMillPatchNamesHandler& patchNames,
            const dictionary& dict,
            const scalar tol = 1e-6
        );

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            dieSurfaceCreator,
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
        static autoPtr<dieSurfaceCreator> New
        (
            const rollingMillPatchNamesHandler& patchNames,
            const dictionary& dict,
            const scalar tol
        );

    // Destructor

        virtual ~dieSurfaceCreator();


    // Member Functions

        //- create a 2D surface in the x-y plane
        const triSurf& surface() const;

        //- profiles at given cross sections
        const DynList<std::pair<scalar, std::shared_ptr<triSurf> > >&
        crossSections() const;

        //- average diameter at the upstream side
        virtual scalar inletDiameter() const = 0;

        //- average diameter at the downstream side
        virtual scalar outletDiameter() const = 0;

        //- outer diameter of a die
        virtual scalar outerDiameter() const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
