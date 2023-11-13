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
    Reads the surface mesh, remove the selected facets
    and writes the modified mesh into a new file

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "edgeMesh.H"
#include "triSurfaceExtrude2DEdges.H"
#include "demandDrivenData.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input edge mesh file");
    argList::validArgs.append("output surface file");

    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);

    //- read the input surface
    triSurf origSurf;

    if( inFileName.ext() == "fms" )
    {
        Info << "Reading edges from " << inFileName << endl;
        origSurf.readSurface(inFileName);
    }
    else
    {
        Info << "Starting converting edge mesh" << endl;

        //- the the edgeMesh and convert it into the triSurf
        const edgeMesh eMesh(inFileName);

        const pointField& points = eMesh.points();
        const edgeList& edges = eMesh.edges();

        triSurfModifier sMod(origSurf);
        pointField& sPts = sMod.pointsAccess();
        edgeLongList& fEdges = sMod.featureEdgesAccess();

        sPts = points;

        fEdges.setSize(edges.size());

        forAll(edges, eI)
            fEdges[eI] = edges[eI];

        Info << "Points " << points << endl;
        Info << "Edges " << edges << endl;

        Info << "Finished converting edge mesh" << endl;
    }

    //- remove the selected facets
    triSurfaceExtrude2DEdges extruder(origSurf);

    const triSurf* newSurfacePtr = extruder.extrudeSurface();

    if( !newSurfacePtr )
        FatalError << "Extruding of the edge mesh failed!" << exit(FatalError);

    //- write the modified surface mesh
    newSurfacePtr->writeSurface(outFileName);

    deleteDemandDrivenData(newSurfacePtr);

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
