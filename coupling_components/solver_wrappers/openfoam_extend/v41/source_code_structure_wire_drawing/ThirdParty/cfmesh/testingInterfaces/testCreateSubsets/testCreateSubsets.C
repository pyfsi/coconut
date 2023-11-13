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

    polyMeshGen pmg(runTime);
    pmg.read();

    for(label i=0;i<500;++i)
    {
        const label cId = pmg.addCellSubset("cs_"+help::labelToText(i));
        const label fId = pmg.addFaceSubset("fs_"+help::labelToText(i));

        const label nObjects = 10 * pmg.cells().size();

        for(label j=0;j<nObjects;++j)
        {
            pmg.addCellToSubset(cId, j % pmg.cells().size());
            pmg.addFaceToSubset(fId, j % pmg.faces().size());
        }
    }

    pmg.write();

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
