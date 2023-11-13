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
    Writes the mesh in fpma format readable by AVL's CfdWM

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurf.H"
#include "DXFReader.H"
#include "triSurfaceCopyParts.H"
#include "triSurfaceExtrude2DEdges.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("input surface file");
    argList::validOptions.insert("blockName", "scalar");
    argList args(argc, argv);

    const fileName inFileName(args.args()[1]);
    const fileName outName = inFileName.lessExt()+".fms";

    word blockName;
    if( args.options().found("blockName") )
    {
        blockName = word(IStringStream(args.options()["blockName"])());
    }

    triSurf surf;

    DXFReader(surf, inFileName, 0.1);

    if( !blockName.empty() )
    {
        triSurfaceCopyParts pCopy(surf);

        wordList w(1);
        w[0] = blockName;

        Info << "copying parts" << endl;
        triSurf* copySurfPtr = pCopy.copySurface(w);
        Info << "Extruding surface" << endl;
        const triSurf* eSurfPtr =
            triSurfaceExtrude2DEdges(*copySurfPtr).extrudeSurface();
        eSurfPtr->writeSurface(outName);
        delete eSurfPtr;
        delete copySurfPtr;
    }
    else
    {
        const triSurf* eSurfPtr =
            triSurfaceExtrude2DEdges(surf).extrudeSurface();
        eSurfPtr->writeSurface(outName);
        delete eSurfPtr;
    }

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
