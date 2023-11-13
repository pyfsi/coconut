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
    dieSurfaceCreatorBekaertSpecs

Description
    Generates inner profile of a die from Bekaert specs

SourceFiles
    dieSurfaceCreatorBekaertSpecs.C

\*---------------------------------------------------------------------------*/

#ifndef dieSurfaceCreatorBekaertSpecs_H
#define dieSurfaceCreatorBekaertSpecs_H

#include "dieSurfaceCreator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class dieSurfaceCreatorBekaertSpecs Declaration
\*---------------------------------------------------------------------------*/

class dieSurfaceCreatorBekaertSpecs
: public dieSurfaceCreator
{
    // Private data
        //- inlet diameter
        scalar inletDiameter_;

        //- scalar outlet diameter
        scalar outletDiameter_;

    // Private member functions
        //- create a line between two points
        void createLine
        (
            const point& firstPoint,
            const point& secondPoint,
            const label code,
            LongList<point>& points,
            edgeLongList& edges,
            labelLongList& edgePatchCode,
            labelLongList& featurePoints
        ) const;

        //- create edge segments across an arc
        void createArc
        (
            const point& firstPoint,
            const point& secondPoint,
            const point& centre,
            const label code,
            LongList<point>& points,
            edgeLongList& edges,
            labelLongList& edgePatchCode,
            labelLongList& featurePoints
        ) const;

        //- create geometry from parameters
        void createTriangulatedSurface
        (
            const scalar H,
            const scalar W,
            const scalar Con,
            const scalar TC,
            const scalar gamma2,
            const scalar BC,
            const scalar gamma1,
            const scalar CR,
            const scalar CRA,
            const scalar EBR,
            const scalar EBL,
            const scalar ECA,
            const scalar ECL,
            const scalar RCA,
            const scalar RCL,
            const scalar RCR,
            const scalar D,
            const scalar BL,
            const scalar BRA,
            const scalar BRL,
            const scalar BRR,
            const scalar XCA,
            const scalar XCL,
            const scalar XCR
        );

        //- Disallow default bitwise assignment
        void operator=(const dieSurfaceCreatorBekaertSpecs&);

public:

    //- Runtime type information
        TypeName("bekaertSpecs");

    // Constructors

        //- Construct from components defining the die
        dieSurfaceCreatorBekaertSpecs
        (
            const rollingMillPatchNamesHandler&,
            const dictionary& dict,
            const scalar tol=1e-6
        );

        //- copy construct
        dieSurfaceCreatorBekaertSpecs
        (
            const dieSurfaceCreatorBekaertSpecs&
        );

    // Return clone of the object
        virtual autoPtr<dieSurfaceCreator> clone() const
        {
            return autoPtr<dieSurfaceCreator>
            (
                new dieSurfaceCreatorBekaertSpecs(*this)
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
