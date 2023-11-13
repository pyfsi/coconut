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
    Testing of quadric patch fitting 

\*---------------------------------------------------------------------------*/

#include "helperFunctions.H"
#include "quadricFitting.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    point p(0., 0., 0.);

    DynList<point> otherPoints(4);
    otherPoints[0] = point(1, 0, 1.0);
    otherPoints[1] = point(-1, 0, 1.0);
    otherPoints[2] = point(0, 1, 0);
    otherPoints[3] = point(0, -1, 0);

    quadricFitting qFit(p, vector(0, 0, 1), otherPoints);
    //quadricFitting qFit(p, otherPoints);
    Info << "Error " << qFit.cumulativeError() << endl;
    Info << "Normal " << qFit.normal() << endl;
    Info << "Max curv " << qFit.maxCurvature() << endl;
    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
