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

#include "bezierSpline.H"
#include "Ostream.H"
#include "Istream.H"
#include "token.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(bezierSpline, 0);
addToRunTimeSelectionTable(splineBase, bezierSpline, typeOfSpline);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bezierSpline::bezierSpline()
:
    splineBase(),
    copyPoints_()
{}

bezierSpline::bezierSpline(const LongList<point>& points, const word& type)
:
    splineBase(points, type),
    copyPoints_()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

point bezierSpline::evaluate(const scalar t) const
{
    const label nPoints = points_.size();

    copyPoints_ = points_;

    label i = nPoints - 1;
    while( i > 0 )
    {
        for(label k=0;k<i;++k)
            copyPoints_[k] += t * (copyPoints_[k+1] - copyPoints_[k] );

        --i;
    }

    return copyPoints_[0];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
