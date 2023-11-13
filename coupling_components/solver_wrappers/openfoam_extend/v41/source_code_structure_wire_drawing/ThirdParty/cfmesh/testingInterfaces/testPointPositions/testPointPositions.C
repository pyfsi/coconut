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
    A testing interface for mesh quality optimisation

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "meshSurfaceEngine.H"
#include "polyMeshGenModifier.H"
#include "polyMeshGenChecks.H"
#include "meshOptimizer.H"
#include "meshSurfaceOptimizer.H"
#include "partTetMesh.H"
#include "tetMeshOptimisation.H"
#include "partTetMeshSimplex.H"
#include "volumeOptimizer.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"

#include "Random.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    pointFieldPMG pointSerial
    (
        IOobject
        (
            "points",
            runTime.constant(),
            "polyMesh",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    std::set<point> pointsInSet;
    forAll(pointSerial, i)
      pointsInSet.insert(pointSerial[i]);

    const label nProcs(8);
    
    for(label procI=0;procI<nProcs;++procI)
    {
        fileName fName = "/processor"+help::labelToText(procI)+"/constant/polyMesh";
        Info << "fName " << fName << endl;
        pointFieldPMG procPoints
        (
            IOobject
            (
                "points",
                "",
                fName,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        forAll(procPoints, i)
        {
            if( pointsInSet.find(procPoints[i]) == pointsInSet.end() )
	    {
	        Info << "Point " << i << " at proc " << procI << " is not found" << endl;
	    }
	}
    }
    
    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
