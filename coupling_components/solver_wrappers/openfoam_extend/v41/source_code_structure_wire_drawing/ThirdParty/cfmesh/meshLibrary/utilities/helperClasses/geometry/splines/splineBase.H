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
    splineBase

Description
    A base class for splines

SourceFiles
    splineBase.C

\*---------------------------------------------------------------------------*/

#ifndef splineBase_H
#define splineBase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "word.H"
#include "label.H"
#include "typeInfo.H"
#include "LongList.H"
#include "point.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class splineBase
{
protected:

    // Protected data
        //- control points of the spline
        mutable LongList<point> points_;

public:

    //- Runtime type information
    TypeName("splineBase");

    // Constructors
        //- Null construct
        splineBase();

        //- Construct from control points
        splineBase(const LongList<point>& points, const word& type);

    // Destructor

        virtual ~splineBase()
        {}

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            splineBase,
            typeOfSpline,
            (
                const LongList<point>& points,
                const word& name
            ),
            (points, name)
        );

    // Selectors

        static autoPtr<splineBase> New
        (
            const LongList<point>& points,
            const word& name
        );

    // Construct and return a clone

        autoPtr<splineBase> clone() const
        {
            return New
            (
                this->points(),
                this->type()
            );
        }

    // Member functions
        //- return the control points of the spline
        inline const LongList<point>& points() const
        {
            return points_;
        }

        //- evaluate a point on the spline for the given parameter
        virtual point evaluate(const scalar t) const = 0;

        //- find the nearest point on the spline to the given point
        virtual point nearestPointOnSpline
        (
            const point& p,
            const scalar tol=SMALL
        ) const;

        //- create a polyLine that differs from the spline less
        //- than the given tolerance
        virtual void createPolyLine
        (
            const scalar tol,
            LongList<point>& linePoints
        ) const;

    // Member operators
        friend Ostream& operator<<(Ostream&, const splineBase&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
