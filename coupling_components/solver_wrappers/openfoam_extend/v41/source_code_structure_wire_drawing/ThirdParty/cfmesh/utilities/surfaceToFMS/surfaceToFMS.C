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
    Reads the specified surface and writes it in the fms format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fileName.H"
#include "triSurf.H"
#include "triSurfaceDetectFeatureEdges.H"
#include "triSurfacePatchManipulator.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("input surface file");

    argList::validOptions.insert("createPatches", "");

    argList args(argc, argv);

    const fileName inFileName(args.args()[1]);
    const fileName outFileName(inFileName.lessExt()+".fms");

    //- read the input surface mesh
    triSurf* surfacePtr = new triSurf(inFileName);

    //- generate feature edges
    triSurfaceDetectFeatureEdges edgeDetector(*surfacePtr);
    edgeDetector.detectFeatureEdges();

    //- create patches if the user asks for it
    if( args.options().found("createPatches") )
    {
        //- create surface patches based on the feature edges
        //- and update the meshDict based on the given data
        triSurfacePatchManipulator manipulator(*surfacePtr);

        const triSurf* surfaceWithPatchesPtr = manipulator.surfaceWithPatches();

        surfaceWithPatchesPtr->writeSurface(outFileName);

        deleteDemandDrivenData(surfaceWithPatchesPtr);
    }
    else
    {
        surfacePtr->writeSurface(outFileName);
    }

    deleteDemandDrivenData(surfacePtr);

    Info << "New surface generated at " << outFileName << "\n" << endl;
    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
