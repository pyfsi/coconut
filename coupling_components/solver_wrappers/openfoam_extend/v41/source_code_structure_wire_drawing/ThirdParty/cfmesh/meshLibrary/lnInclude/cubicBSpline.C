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

Description

\*---------------------------------------------------------------------------*/

#include "cubicBSpline.H"
#include "Ostream.H"
#include "Istream.H"
#include "token.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "tridiagonalMatrix.H"

//# define DEBUGSpline

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(cubicBSpline, 0);
addToRunTimeSelectionTable(splineBase, cubicBSpline, typeOfSpline);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void cubicBSpline::cleanupNearbyPoints() const
{
    //- check whether the two consecutive points are very close to one another
    //- remove points that are too close

    //- start by calculating the total length of the curve
    scalar sum(0.0);

    label nPoints = points_.size();
    for(label i=1;i<nPoints;++i)
        sum += mag(points_[i] - points_[i-1]);

    //- filter out vertices that are too close to the previous one
    label lastPoint(0);
    for(label i=1;i<points_.size();++i)
    {
        const point& lp = points_[lastPoint];

        if( (mag(points_[i] - lp) / sum) > SMALL )
        {
            ++lastPoint;

            if( i > lastPoint )
            {
                points_[lastPoint] = points_[i];
                parameters_[lastPoint] = parameters_[i];
            }
        }
    }

    points_.setSize(lastPoint+1);
    parameters_.setSize(lastPoint+1);

    //- set the parameter of the first point to zero
    parameters_[0] = 0.0;

    //- set the parameter of the last point to one
    parameters_[parameters_.size()-1] = 1.0;

    # ifdef DEBUGSpline
    Info << "Original number of spline points " << nPoints << endl;
    Info << "Number of spline points after filtering " << lastPoint << endl;
    Info << "Parameters " << parameters_ << endl;
    # endif
}

void cubicBSpline::calculateParameters() const
{
    scalar sum(0.0);

    parameters_.setSize(points_.size());
    parameters_[0] = 0.0;

    //- set the parameter value to every point on the spline
    const label nIntervals = points_.size() - 1;
    for(label intervalI=0;intervalI<nIntervals;++intervalI)
    {
        const scalar d = mag(points_[intervalI+1] - points_[intervalI]);

        parameters_[intervalI+1] = parameters_[intervalI] + d;
        sum += d;
    }

    //- normalize parameters
    forAll(parameters_, i)
        parameters_[i] /= sum;

    //- set the parameter of the first point to zero
    parameters_[0] = 0.0;

    //- set the parameter of the last point to one
    parameters_[parameters_.size()-1] = 1.0;
}

void cubicBSpline::calculateCoefficients() const
{
    if( parameters_.size() == 0 )
    {
        calculateParameters();
    }

    cleanupNearbyPoints();

    //- create the system matrix
    tridiagonalMatrix<vector> mat(points_.size());

    const label maxRow = points_.size() - 1;

    for(label j=1;j<maxRow;++j)
    {
        const scalar hj = parameters_[j+1] - parameters_[j] + ROOTVSMALL;
        const scalar hjj = parameters_[j] - parameters_[j-1] + ROOTVSMALL;

        const vector s =
            3.0 * (points_[j+1] - points_[j]) / hj -
            3.0 * (points_[j] - points_[j-1]) / hjj;

        mat.setSource(j, s);

        mat.setCoeff(j, j-1, hjj);
        mat.setCoeff(j, j, 2.0 * (hj + hjj));
        mat.setCoeff(j, j+1, hj);
    }

    //- boundary conditions
    if( startTangentRequest_ )
    {
        const scalar h0 = parameters_[1] - parameters_[0] + ROOTVSMALL;
        const vector s =
            3.0 * ((points_[1] - points_[0]) / h0) - (3.0 * startTangent_);

        mat.setCoeff(0, 0, 2.0 * h0);
        mat.setCoeff(0, 1, h0);
        mat.setSource(0, s);
    }
    else
    {
        mat.setCoeff(0, 0, 1.0);
        mat.setSource(0, vector::zero);
    }

    if( endTangentRequest_ )
    {
        const scalar hi = parameters_[maxRow]-parameters_[maxRow-1]+ROOTVSMALL;
        const vector s =
            3.0 * endTangent_ -
            3.0 * (points_[maxRow] - points_[maxRow-1]) / hi;

        mat.setCoeff(maxRow, maxRow-1, hi);
        mat.setCoeff(maxRow, maxRow, 2.0 * hi);
        mat.setSource(maxRow, s);
    }
    else
    {
        mat.setCoeff(maxRow, maxRow, 1.0);
        mat.setSource(maxRow, vector::zero);
    }

    //- solve the system
    const LongList<vector> c = mat.solve();

    //- calculate the coefficients
    coefficients_.setSize(maxRow);
    forAll(coefficients_, i)
    {
        const scalar h = parameters_[i+1] - parameters_[i] + ROOTVSMALL;

        FixedList<vector, 4>& coeffs = coefficients_[i];

        coeffs[0] = points_[i];
        coeffs[1] =
            ((points_[i+1] - points_[i]) / h) -
            ((c[i+1] + 2.0 * c[i]) * (h/3.));
        coeffs[2] = c[i];
        coeffs[3] = (c[i+1] - c[i]) / (3.0*h);
    }

    done_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

cubicBSpline::cubicBSpline()
:
    splineBase(),
    coefficients_(),
    startTangent_(),
    endTangent_(),
    parameters_(),
    startTangentRequest_(false),
    endTangentRequest_(false),
    done_(false)
{}

cubicBSpline::cubicBSpline(const LongList<point>& points, const word type)
:
    splineBase(points, type),
    coefficients_(points.size()-1, FixedList<vector, 4>()),
    startTangent_(),
    endTangent_(),
    parameters_(),
    startTangentRequest_(false),
    endTangentRequest_(false),
    done_(false)
{}

cubicBSpline::cubicBSpline
(
    const LongList<point>& points,
    const scalarLongList& params,
    const word type
)
:
    splineBase(points, type),
    coefficients_(points.size()-1, FixedList<vector, 4>()),
    startTangent_(),
    endTangent_(),
    parameters_(params),
    startTangentRequest_(false),
    endTangentRequest_(false),
    done_(false)
{}

cubicBSpline::cubicBSpline(const cubicBSpline& bs)
:
    splineBase(bs.points_, bs.type()),
    coefficients_(bs.coefficients_),
    startTangent_(bs.startTangent_),
    endTangent_(bs.endTangent_),
    parameters_(bs.parameters_),
    startTangentRequest_(bs.startTangentRequest_),
    endTangentRequest_(bs.endTangentRequest_),
    done_(bs.done_)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label cubicBSpline::numberOfControlPoints() const
{
    //- calculate coefficients for cubic spline interpolation
    if( !done_ )
        calculateCoefficients();

    return points_.size();
}

const point& cubicBSpline::controlPoint(const label pI) const
{
    //- calculate coefficients for cubic spline interpolation
    if( !done_ )
        calculateCoefficients();

    return points_[pI];
}

void cubicBSpline::setTangentAtStartingPoint(const vector& t0)
{
    endTangent_ = t0;
    endTangentRequest_ = true;
}

void cubicBSpline::setTangentAtEndPoint(const vector& t1)
{
    endTangent_ = t1;
    endTangentRequest_ = true;
}

scalar cubicBSpline::evaluateParam(const point& p) const
{
    //- calculate coefficients for cubic spline interpolation
    if( !done_ )
        calculateCoefficients();

    scalar distSq(VGREAT);
    scalar t = -1.0;

    //- detect the nearest point to the given point and return its parameter
    forAll(points_, pI)
    {
        const scalar dSq = magSqr(points_[pI] - p);

        if( dSq < distSq )
        {
            t = parameters_[pI];
            distSq = dSq;
        }
    }

    return t;
}

scalar cubicBSpline::pointParam(const label pI) const
{
    //- calculate coefficients for cubic spline interpolation
    if( !done_ )
        calculateCoefficients();

    if( pI < 0 || pI >= parameters_.size() )
    {
        FatalErrorIn
        (
            "scalar cubicBSpline::pointParam(const label) const"
        ) << "Point " << pI
          << " is not in a valid range" << abort(FatalError);
    }

    return parameters_[pI];
}

point cubicBSpline::evaluate(const scalar t) const
{
    //- calculate coefficients for cubic spline interpolation
    if( !done_ )
        calculateCoefficients();

    label minInterval(0), maxInterval(parameters_.size()-1);

    //- find the interval containing the requested value of the parameter
    while( (maxInterval-1) > minInterval )
    {
        const label i = (minInterval + maxInterval) / 2;

        if( parameters_[i] > t )
        {
            maxInterval = i;
        }
        else
        {
            minInterval = i;
        }
    }

    //- check if a valid interval was found
    if( parameters_[minInterval] > t || parameters_[minInterval+1] < t )
    {
        FatalErrorIn
        (
            "point cubicBSpline::evaluate(const scalar t) const"
        ) << "Paramemeter " << t << " is not inside the interval "
          << parameters_[minInterval] << ", " << parameters_[minInterval+1]
          << abort(FatalError);
    }

    //- calculate the point coordinates based on the interpolation polynomial
    point ret(vector::zero);

    const FixedList<vector, 4>& coeffs = coefficients_[minInterval];
    const scalar dt = t - parameters_[minInterval];

    ret += coeffs[0];
    scalar d = dt;
    ret += coeffs[1] * d;
    d *= dt;
    ret += coeffs[2] * d;
    d *= dt;
    ret += coeffs[3] * d;

    return ret;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
