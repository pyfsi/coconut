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

#include "splineBase.H"
#include "Ostream.H"
#include "Istream.H"
#include "token.H"
#include "IOPtrList.H"
#include "dictionary.H"
#include "helperFunctions.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(splineBase, 0);
defineRunTimeSelectionTable(splineBase, typeOfSpline);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<splineBase> splineBase::New
(
    const LongList<point>& points,
    const word& type
)
{
    typeOfSplineConstructorTable::iterator cstrIter =
        typeOfSplineConstructorTablePtr_->find(type);

    if( cstrIter == typeOfSplineConstructorTablePtr_->end() )
    {
        FatalErrorIn
        (
            "splineBase::New(const LongList<point>&, const word&)"
        )   << "Unknown splineBase type " << type << nl << nl
            << "Valid splineBase types are :" << nl
            << "[default: " << typeName_() << "]"
            << typeOfSplineConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<splineBase>(cstrIter()(points, type));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

splineBase::splineBase()
:
    points_()
{}

splineBase::splineBase(const LongList<point>& points, const word& type)
:
    points_(points)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

point splineBase::nearestPointOnSpline
(
    const point& p,
    const scalar tol
) const
{
    const scalar tolSqr = sqr(tol);

    std::map<scalar, scalar> distToParam;

    //- create an initial distribution of distances
    scalar dt = 0.25;
    scalar t = 0.0;
    while( t < (1.0+SMALL) )
    {
        t = min(1.0, t);
        distToParam[magSqr(evaluate(t)-p)] = t;
        t += dt;
    }

    Info << "1.Status of parameter " << distToParam.begin()->second << endl;
    Info << "1.Evaluations in map " << label(distToParam.size()) << endl;

    //- start the loop by trying to find the nearest point near
    //- the one currently classified as the nearest one
    label i(0);
    do
    {
        Info << "Starting iter " << ++i << endl;
        const scalar currt = distToParam.begin()->second;

        Info << "currt " << currt << endl;

        //- decrease the step
        dt *= 0.5;

        std::map<scalar, scalar>::const_iterator it = distToParam.begin();
        ++it;
        const scalar rSq = 9.0 * it->first;
        Info << "1. rSq " << rSq << endl;
        Info << "Num evaluations " << label(distToParam.size()) << endl;
        Info << "dt " << dt << endl;

        //- evaluate points in the vicinity of the nearest point
        DynList<std::pair<scalar, scalar> > newValues;
        for(it=distToParam.begin();it!=distToParam.end();++it)
        {
            if( it->first > rSq )
                break;

            Info << "curr rSq " << it->first << endl;
            const scalar tless = max(0.0, it->second - 0.5 * dt);
            const scalar tmore = min(1.0, it->second + 0.5 * dt);

            Info << "Orig t " << it->second << endl;
            Info << "tless " << tless << endl;
            Info << "tmore " << tmore << endl;

            //- calculate distances of adjacent points at the curve
            newValues.append(std::make_pair(magSqr(evaluate(tless)-p), tless));
            newValues.append(std::make_pair(magSqr(evaluate(tmore)-p), tmore));
        }

        //- insert new values into the map
        forAll(newValues, i)
            distToParam.insert(newValues[i]);

        Info << "min dist after " << distToParam.begin()->first << endl;

        if( magSqr(distToParam.begin()->second - currt) < tolSqr )
            break;

    } while( true );

    Info << "2.Status of parameter " << distToParam.begin()->second << endl;
    Info << "2.Evaluations in map " << label(distToParam.size()) << endl;
    Info << "points " << points_ << endl;

    return evaluate(distToParam.begin()->second);
}

void splineBase::createPolyLine
(
    const scalar tol,
    LongList<point>& pointsOnTheCurve
) const
{
    pointsOnTheCurve.clear();

    label nDiv = 10;

    point p, pnext;

    bool finished;
    do
    {
        Info << "Number of divisions " << nDiv << endl;
        finished = true;

        pointsOnTheCurve.setSize(nDiv+1);

        const scalar dt = 1.0 / nDiv;

        p = points_[0];
        pointsOnTheCurve[0] = p;

        //- set the last point to the last point on the spline
        pointsOnTheCurve[nDiv] = points_[points_.size()-1];

        for(label lI=1;lI<=nDiv;++lI)
        {
            //- evaluate the last point at the interval
            if( lI == (nDiv-1) )
            {
                pnext = points_[points_.size()-1];
            }
            else
            {
                pnext = evaluate(dt * lI);
            }

            pointsOnTheCurve[lI] = pnext;

            //- evaluate the point in the middle of the interval
            //- and calculate its distance from the line formed by
            //- the two consecutive points on the curve
            const point mp = evaluate((dt * lI) - 0.5 * dt);

            const point np = help::nearestPointOnTheEdgeExact(p, pnext, mp);

            if( magSqr(mp - np) > sqr(0.5 * tol) )
            {
                //- the error is larger than the tolerance
                //- increase the number of subdivisions
                nDiv *= 2;
                finished = false;
                break;
            }

            //- use the current point as the old point
            p = pnext;
        }
    } while( !finished );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const splineBase& sb)
{
    os.check("Ostream& operator<<(Ostream& f, const splineBase& sb");

    os << sb.points_;

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
