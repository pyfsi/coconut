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
    Finds feature edges and corners of a triangulated surface

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fileName.H"
#include "triSurf.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "demandDrivenData.H"
#include "triSurfaceImportSurfaceAsSubset.H"
#include <cstdlib>
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("master surface file");
    argList::validArgs.append("import surface file");

    argList::validOptions.insert
    (
        "maxNormalDeviationAngle",
        "Maximum angle deviation allowed between the surface facets"
    );
    argList::validOptions.insert
    (
        "maxDistanceTol",
        "Maximum distance between the surface facets"
    );

    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName importFileName(args.args()[2]);

    scalar maxAngleDeviation(5.0);
    if( args.options().found("maxNormalDeviationAngle") )
    {
        maxAngleDeviation = args.optionRead<scalar>("maxNormalDeviationAngle");
    }

    scalar maxDistanceDeviation(1e-6);
    if( args.options().found("maxDistanceTol") )
    {
        maxDistanceDeviation = args.optionRead<scalar>("maxDistanceTol");
    }

    Info << "Maximum deviation angle between surface normals is "
         << maxAngleDeviation << " deg" << endl;

    Info << "Reading master mesh" << endl;
    triSurf originalSurface(inFileName);

    Info << "Reading import surface" << endl;
    triSurf importedSurface(importFileName);

    Info << "Creating a subset in the master surface mesh" << endl;
    triSurfaceImportSurfaceAsSubset importSurf(originalSurface);

    const word importFileNameNoDir = importFileName.lessExt().name();
    Info << "Subset name " << importFileNameNoDir << endl;
    importSurf.addSurfaceAsSubset
    (
        importedSurface,
        importFileNameNoDir,
        maxAngleDeviation,
        maxDistanceDeviation
    );

    if( inFileName.ext() == "fms" )
    {
        Info << "Writting surface with the subset" << endl;
        originalSurface.writeSurface(inFileName);
    }
    else
    {
        fileName newName = inFileName.lessExt();
        newName.append(".fms");
        Warning << "Writting surface as " << newName
            << " to preserve the subset!!" << endl;

        originalSurface.writeSurface(newName);
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
