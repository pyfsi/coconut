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
    casingAxialCrossSectionSurfaceCreatorBekaertSpecs

Description
    Generates inner profile of a die from Bekaert specs

SourceFiles
    casingAxialCrossSectionSurfaceCreatorBekaertSpecs.C

\*---------------------------------------------------------------------------*/

#ifndef casingAxialCrossSectionSurfaceCreatorBekaertSpecs_H
#define casingAxialCrossSectionSurfaceCreatorBekaertSpecs_H

#include "casingAxialCrossSectionSurfaceCreator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class casingAxialCrossSectionSurfaceCreatorBekaertSpecs Declaration
\*---------------------------------------------------------------------------*/

class casingAxialCrossSectionSurfaceCreatorBekaertSpecs
: public casingAxialCrossSectionSurfaceCreator
{
    // Private data

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
            const scalar diameter,
            const scalar height,
            const scalar BC,
            const scalar TC,
            const scalar EX_OD,
            const scalar N_P,
            const scalar EX_ID,
            const scalar N_X,
            const scalar EN_ID,
            const scalar EN_OD
        );

        //- Disallow default bitwise assignment
        void operator=
        (
            const casingAxialCrossSectionSurfaceCreatorBekaertSpecs&
        );

public:

    //- Runtime type information
        TypeName("bekaertSpecs");

    // Constructors

        //- Construct from components defining the die
        casingAxialCrossSectionSurfaceCreatorBekaertSpecs
        (
            const rollingMillPatchNamesHandler&,
            const dictionary& dict,
            const scalar tol=1e-6
        );

        //- copy construct
        casingAxialCrossSectionSurfaceCreatorBekaertSpecs
        (
            const casingAxialCrossSectionSurfaceCreatorBekaertSpecs&
        );

    // Return clone of the object
        virtual autoPtr<casingAxialCrossSectionSurfaceCreator> clone() const
        {
            return autoPtr<casingAxialCrossSectionSurfaceCreator>
            (
                new casingAxialCrossSectionSurfaceCreatorBekaertSpecs(*this)
            );
        }

    // Member Functions
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
