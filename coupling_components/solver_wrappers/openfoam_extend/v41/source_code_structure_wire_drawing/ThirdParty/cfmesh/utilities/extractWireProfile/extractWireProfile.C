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
    Reads a mesh with many zones and translates the vertices in each zone
    by a given translation vector

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMeshGenModifier.H"
#include "rollingMillMesh.H"
#include "meshSurfaceEngine.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "helperFunctions.H"

#include "OFstream.H"

#include <map>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurf* extractWireProfile(polyMeshGen& mesh, const dictionary& meshDict)
{
    const PtrList<boundaryPatch>& boundaries = mesh.boundaries();

    //- extract the DOWNSTREAM and CONTACT patches
    const rollingMillPatchNamesHandler patchHandler(meshDict);

    const word downstreamName =
        patchHandler.patchNameForWire
        (
            rollingMillPatchNamesHandler::WIREDOWNSTREAM
        );

    label downstreamId(-1);
    forAll(boundaries, patchI)
    {
        if( boundaries[patchI].patchName().find(downstreamName) != word::npos )
        {
            downstreamId = patchI;
        }
    }

    if( downstreamId < 0 )
    {
        FatalError << "Cannot find downstream patch " << downstreamName
            << " in the mesh."
            << " Maybe you have tried to rename the patch in"
            << " rollingMillPatchNames dictionary."
            << exit(FatalError);
    }

    //- extract edges at the boundary between the patches
    meshSurfaceEngine mse(mesh);
    const edgeLongList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelLongList& facePatch = mse.boundaryFacePatches();

    edgeLongList profileEdges;
    labelLongList edgePatch;
    forAll(edgeFaces, beI)
    {
        if( edgeFaces.sizeOfRow(beI) != 2 )
            continue;

        const label patch0 = facePatch[edgeFaces(beI, 0)];
        const label patch1 = facePatch[edgeFaces(beI, 1)];

        if
        (
            ((patch0 != downstreamId) && (patch1 == downstreamId)) ||
            ((patch0 == downstreamId) && (patch1 != downstreamId))
        )
        {
            profileEdges.append(edges[beI]);
            edgePatch.append(patch0==downstreamId?patch1:patch0);
        }
    }

    //- create a surface mesh from profile edges
    LongList<point> sPts;
    std::map<label, label> newPointLabel;

    forAll(profileEdges, eI)
    {
        const edge& e = profileEdges[eI];

        forAll(e, pI)
        {
            if( newPointLabel.find(e[pI]) == newPointLabel.end() )
            {
                newPointLabel[e[pI]] = sPts.size();
                sPts.append(mesh.points()[e[pI]]);
            }

            profileEdges[eI][pI] = newPointLabel[e[pI]];
        }
    }

    //- create and extrude triangulated surface
    triSurf* surfPtr = new triSurf();
    triSurfModifier sMod(*surfPtr);

    pointField& points = sMod.pointsAccess();
    points.setSize(sPts.size());
    forAll(points, i)
    {
        points[i] = sPts[i];
        points[i].x() = points[i].z();
        points[i].z() = 0.0;
    }

    std::map<label, label> patchToSubset;
    forAll(profileEdges, i)
    {
        if( patchToSubset.find(edgePatch[i]) == patchToSubset.end() )
        {
            patchToSubset[edgePatch[i]] =
                surfPtr->addEdgeSubset(boundaries[edgePatch[i]].patchName());
        }

        surfPtr->addEdgeToSubset(patchToSubset[edgePatch[i]], i);
        surfPtr->appendFeatureEdge(profileEdges[i]);
    }

    triSurfaceExtrude2DEdges extruder(*surfPtr);
    triSurf* extrudedSurfPtr = new triSurf();
    extruder.extrudeSurface(*extrudedSurfPtr);
    deleteDemandDrivenData(surfPtr);

    pointField& extrudedPoints =
        triSurfModifier(*extrudedSurfPtr).pointsAccess();
    forAll(extrudedPoints, i)
    {
        point& p = extrudedPoints[i];

        const scalar z = p.z();
        p.z() = p.x();
        p.x() = z;
    }

    return extrudedSurfPtr;
}

int main(int argc, char *argv[])
{
    argList::validOptions.insert("region", "fileName");

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

    fileName instance("");
    word val;
    if( args.optionReadIfPresent("region", val) )
    {
        instance = val;
        Info<< "Reading wire from " << instance << nl << endl;
    }


    //- load the mesh from disk
    polyMeshGen pmg(runTime, instance);
    pmg.readFromLatestTime();

    Info << "The mesh has " << pmg.cells().size() << endl;

    //- extract the wire profile from the mesh
    triSurf* surfPtr = extractWireProfile(pmg, meshDict);

    Info << "Writing profile surface to disk" << endl;

    //- write file to disk
    surfPtr->writeSurface("wireProfile_"+runTime.caseName()+".fms");
    deleteDemandDrivenData(surfPtr);

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
