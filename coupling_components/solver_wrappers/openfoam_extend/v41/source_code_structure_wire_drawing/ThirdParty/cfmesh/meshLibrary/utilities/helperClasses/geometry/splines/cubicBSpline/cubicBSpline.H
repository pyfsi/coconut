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
    cubicBSpline

Description
    Implementation of a cubic b-spline using a cubic polynomial interpolation

SourceFiles
    cubicBSpline.C

\*---------------------------------------------------------------------------*/

#ifndef cubicBSpline_H
#define cubicBSpline_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "splineBase.H"
#include "typeInfo.H"
#include "vectorField.H"
#include "scalarLongList.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class cubicBSpline
:   public splineBase
{
    // Private data
        //- coefficients
        mutable List<FixedList<vector, 4> > coefficients_;

        //- tangents
        vector startTangent_;
        vector endTangent_;

        //- parameter on a curve for each given point
        mutable scalarLongList parameters_;

        //- has the user imposed requst on tangentionality
        bool startTangentRequest_;
        bool endTangentRequest_;

        //- is interpolation performed
        mutable bool done_;

    // Private member functions
        //- cleanup nearby points
        void cleanupNearbyPoints() const;

        //- calculate parameters
        void calculateParameters() const;

        //- evaluate the coefficients of the spline
        void calculateCoefficients() const;

public:

    //- Runtime type information
    TypeName("cubicBSpline");

    // Constructors
        //- Null construct
        cubicBSpline();

        //- Construct from points
        cubicBSpline(const LongList<point>& points, const word type);

        //- Construct from points and parameters assigned to each point
        cubicBSpline
        (
            const LongList<point>& points,
            const scalarLongList& params,
            const word type
        );

        //- copy constructor
        cubicBSpline(const cubicBSpline&);

    // Return clone of the object
        virtual autoPtr<splineBase> clone() const
        {
            return autoPtr<splineBase>
            (
                new cubicBSpline(*this)
            );
        }

    // Public member functions
        //- return the number of control points
        label numberOfControlPoints() const;

        //- return the coordinate a control point
        const point& controlPoint(const label) const;

        //- set the tangent at the first point
        void setTangentAtStartingPoint(const vector&);

        //- set the tangent at the second point
        void setTangentAtEndPoint(const vector&);

        //- get the parameter of a nearest spline point
        scalar evaluateParam(const point&) const;

        //- get the value of parameter assigned to a given control point
        scalar pointParam(const label) const;

        //- evaluate a point on the spline for the given parameter
        virtual point evaluate(const scalar t) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
