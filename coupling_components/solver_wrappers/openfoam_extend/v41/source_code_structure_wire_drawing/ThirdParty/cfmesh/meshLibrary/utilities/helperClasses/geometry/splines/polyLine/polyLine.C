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

#include "polyLine.H"
#include "Ostream.H"
#include "Istream.H"
#include "token.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(polyLine, 0);
addToRunTimeSelectionTable(splineBase, polyLine, typeOfSpline);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyLine::calculateCoefficients()
{
    startingCoeffs_.setSize(points_.size()-1);

    scalar s(0.0);
    forAll(startingCoeffs_, i)
    {
        const scalar l = mag(points_[i+1] - points_[i]);

        startingCoeffs_[i] = s;

        s += l;
    }

    startingCoeffs_ /= s;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyLine::polyLine()
:
    splineBase(),
    startingCoeffs_()
{
    calculateCoefficients();
}

polyLine::polyLine(const LongList<point>&points, const word& name)
:
    splineBase(points, name),
    startingCoeffs_()
{
    calculateCoefficients();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

point polyLine::evaluate(const scalar t) const
{
    scalar ct = max(0.0, t);
    ct = min(1.0, ct);

    forAllReverse(startingCoeffs_, i)
    {
        if( startingCoeffs_[i] <= ct )
        {
            if( i == (startingCoeffs_.size()-1) )
            {
                //- point is in the last interval
                const scalar lt =
                    (ct-startingCoeffs_[i]) / (1.0-startingCoeffs_[i]);

                const point np
                (
                    points_[i] * (1.0 - lt) + lt * points_[i+1]
                );

                return np;
            }
            else
            {
                const scalar lt =
                    (ct-startingCoeffs_[i]) /
                    (startingCoeffs_[i+1]-startingCoeffs_[i]);

                const point np
                (
                    points_[i] * (1.0 - lt) + lt * points_[i+1]
                );

                return np;
            }
        }
    }

    FatalErrorIn
    (
        "point polyLine::evaluate(const scalar) const"
    ) << "Invalid interpolation for parameter " << t << endl;
    return vector::zero;
}

void polyLine::createPolyLine
(
    const scalar tol,
    LongList<point>& linePoints
) const
{
    if( tol < VSMALL )
        Warning << "Invalid tolerance provided for the poly line" << endl;

    linePoints = points_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
