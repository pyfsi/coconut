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
    bezierSpline

Description
    Calculates Bezier spline through a set of points

SourceFiles
    bezierSpline.C

\*---------------------------------------------------------------------------*/

#ifndef bezierSpline_H
#define bezierSpline_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "splineBase.H"
#include "typeInfo.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class bezierSpline
:   public splineBase
{
    //- private data
    mutable LongList<point> copyPoints_;

public:

    //- Runtime type information
    TypeName("bezierSpline");

    // Constructors
        //- Null construct
        bezierSpline();

        //- Construct from components
        bezierSpline(const LongList<point>& points, const word& name);

    // Return clone of the object
        virtual autoPtr<splineBase> clone() const
        {
            return autoPtr<splineBase>
            (
                new bezierSpline
                (
                    this->points(),
                    this->type()
                )
            );
        }

    // Public member functions
        //- evaluate a point on the spline for the given parameter
        virtual point evaluate(const scalar t) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
