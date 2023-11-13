/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2005-2007 Franjo Juretic
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    Test for smoothers

Description
    - reads the mesh and tries to untangle negative volume cells

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMeshGen.H"
#include "meshSurfaceEdgeExtractorFUN.H"
#include "triSurf.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "triSurfacePatchManipulator.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const fileName surfFile(meshDict.lookup("surfaceFile"));
    const triSurf* surfacePtr = new triSurf(surfFile);

    if( surfacePtr->featureEdges().size() != 0 )
    {
        //- create surface patches based on the feature edges
        //- and update the meshDict based on the given data
        triSurfacePatchManipulator manipulator(*surfacePtr);

        const triSurf* surfaceWithPatches =
            manipulator.surfaceWithPatches(&meshDict);

        //- delete the old surface and assign the new one
        deleteDemandDrivenData(surfacePtr);
        surfacePtr = surfaceWithPatches;
    }

    meshOctree mo(*surfacePtr);
    meshOctreeCreator(mo, meshDict).createOctreeWithRefinedBoundary(20, 30);

    polyMeshGen pmg(runTime);
    pmg.read();

    meshSurfaceEdgeExtractorFUN(pmg, mo, false);

    deleteDemandDrivenData(surfacePtr);

    pmg.write();

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
