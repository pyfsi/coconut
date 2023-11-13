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
#include "triSurf.H"
#include "triSurfaceCurvatureEstimator.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    fileName surfName(argv[1]);
    Info << "Reading surface from file " << surfName << endl;

    triSurf surf(surfName);

    triSurfaceCurvatureEstimator curv(surf);

    forAll(surf, triI)
       Info << "triangle " << triI << " max curvature "
            << curv.maxCurvatureAtTriangle(triI) << endl;

    Info << "End\n" << endl;
    
    return 0;
}

// ************************************************************************* //
