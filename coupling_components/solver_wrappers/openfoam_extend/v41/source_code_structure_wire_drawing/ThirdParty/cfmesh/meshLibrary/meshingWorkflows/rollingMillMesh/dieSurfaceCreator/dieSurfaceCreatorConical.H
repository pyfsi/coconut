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
    dieSurfaceCreatorConical

Description
    Generates inner profile of a die

SourceFiles
    dieSurfaceCreatorConical.C

\*---------------------------------------------------------------------------*/

#ifndef dieSurfaceCreatorConical_H
#define dieSurfaceCreatorConical_H

#include "dieSurfaceCreator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class dieSurfaceCreatorConical Declaration
\*---------------------------------------------------------------------------*/

class dieSurfaceCreatorConical
: public dieSurfaceCreator
{
    // Private data
        //- inlet radius
        scalar inletRadius_;

        //- outlet radius
        scalar outletRadius_;

        //- axial length
        scalar axialLength_;

        //- outer diameter
        scalar outerRadius_;

    // Private member functions
        //- create surface triangulation
        void createTriangulatedSurface();

        //- Disallow default bitwise assignment
        void operator=(const dieSurfaceCreatorConical&);

public:

    //- Runtime type information
        TypeName("conical");

    // Constructors
        //- construct from patch names, dictionary and tolerance
        dieSurfaceCreatorConical
        (
            const rollingMillPatchNamesHandler& patchNames,
            const dictionary& dict,
            const scalar tol = 1e-6
        );

        //- copy construct
        dieSurfaceCreatorConical
        (
            const dieSurfaceCreatorConical&
        );

    // Return clone of the object
        virtual autoPtr<dieSurfaceCreator> clone() const
        {
            return autoPtr<dieSurfaceCreator>
            (
                new dieSurfaceCreatorConical(*this)
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
