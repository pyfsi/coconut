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
#include "triSurfModifier.H"
#include "meshOctreeCreator.H"
#include "boundBox.H"
#include "OFstream.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");
    argList::validArgs.append("origin");
    argList::validArgs.append("symmetryNormal");

    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);

    if (outFileName == inFileName)
    {
        FatalErrorIn(args.executable())
            << "Output file " << outFileName
            << " would overwrite the input file."
            << exit(FatalError);
    }

    const vector origin(ITstream(args.args()[3], tokenList())());
    const vector normal(ITstream(args.args()[4], tokenList())());
    const plane pl(origin, normal);

    triSurf origSurface(inFileName);
    triSurfModifier sMod(origSurface);
    pointField& points = sMod.pointsAccess();
    List<geometricSurfacePatch>& patches = sMod.patchesAccess();

    const label currPatchId = patches.size();
    patches.setSize(currPatchId+1);
    patches[currPatchId] =
        geometricSurfacePatch("patch", "createdFacets", currPatchId);

    //- detect boundary edge and check their distance from the symmetry plane
    const edgeLongList& edges = origSurface.edges();
    const VRWGraph& edgeFaces = origSurface.edgeFacets();

    std::map<label, label> pointToNewLabel;

    forAll(edgeFaces, edgeI)
    {
        if( edgeFaces.sizeOfRow(edgeI) == 1 )
        {
            //- this edge is a boundary edge
            const edge& e = edges[edgeI];
            const point& sp = points[e.start()];
            const point& ep = points[e.end()];

            const point nsp = pl.nearestPoint(sp);
            const point nep = pl.nearestPoint(ep);

            if( pointToNewLabel.find(e.start()) == pointToNewLabel.end() )
            {
                pointToNewLabel[e.start()] = points.size();
                origSurface.appendVertex(nsp);
            }

            if( pointToNewLabel.find(e.end()) == pointToNewLabel.end() )
            {
                pointToNewLabel[e.end()] = points.size();
                origSurface.appendVertex(nep);
            }

            //- create a new face
            origSurface.appendTriangle
            (
                labelledTri
                (
                    pointToNewLabel[e.start()],
                    e.start(),
                    e.end(),
                    currPatchId
                )
            );

            origSurface.appendTriangle
            (
                labelledTri
                (
                    pointToNewLabel[e.start()],
                    e.end(),
                    pointToNewLabel[e.end()],
                    currPatchId
                )
            );
        }
    }

    //- write the surface
    origSurface.writeSurface(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
