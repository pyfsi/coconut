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
    Performs point relocations in the mesh (smoothing) in order to
    improve quality measures. It does not make the mesh invalied.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMeshGenModifier.H"
#include "meshSurfaceEngine.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "triSurf.H"
#include "meshSurfaceDistanceFromGeometry.H"
#include "helperFunctions.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.clear();

    argList::validOptions.insert("distanceThreshold", "scalar");
    argList::validOptions.insert("normalDeviationAngleThreshold", "scalar");

#   include "setRootCase.H"
#   include "createTime.H"

    //- read the settings
    scalar distanceThreshold(VGREAT);
    scalar deviationAngleThreshold(5.0 * M_PI / 180.0);

    //- load the mesh from disk
    polyMeshGen pmg(runTime);
    pmg.read();

    //- create mesh surface engine
    const meshSurfaceEngine surfaceEngine(pmg);

    //- construct mesh dictionary
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

    //- load the surface mesh
    fileName surfaceFile = meshDict.lookup("surfaceFile");
    if( Pstream::parRun() )
        surfaceFile = ".."/surfaceFile;

    triSurf surf(runTime.path()/surfaceFile);

    //- create the octree
    meshOctree octree(surf);

    //- create the octree
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(20, 30);

    //- calculate the distance of the mesh from the input geometry
    meshSurfaceDistanceFromGeometry dist(surfaceEngine, octree);

    //- calculate the average at points
    const labelLongList& bPoints = surfaceEngine.boundaryPoints();
    const faceList::subList& bFaces = surfaceEngine.boundaryFaces();

    if( args.options().found("distanceThreshold") )
    {
        distanceThreshold =
            readScalar(IStringStream(args.options()["distanceThreshold"])());
    }
    else
    {
        //- take average deviation as threshold
        scalar avgPointDst(0.0), avgFaceDst(0.0);

        label nPoints = bPoints.size();
        if( Pstream::parRun() )
        {
            const labelLongList& gpl = surfaceEngine.globalBoundaryPointLabel();
            forAll(gpl, bpI)
                nPoints = max(gpl[bpI], nPoints);

            nPoints += 1;
        }

        reduce(nPoints, maxOp<label>());

        forAll(surfaceEngine.boundaryPoints(), bpI)
            avgPointDst += dist.boundaryPointDistance(bpI);

        avgPointDst = returnReduce(avgPointDst, sumOp<scalar>()) / nPoints;

        distanceThreshold = avgPointDst;

        //- calculate the average distance of face centres
        const label nFaces = returnReduce(bFaces.size(), sumOp<label>());

        forAll(bFaces, bfI)
            avgFaceDst += dist.boundaryFaceCentreDistance(bfI);

        avgFaceDst = returnReduce(avgFaceDst, sumOp<scalar>()) / nFaces;

        distanceThreshold = max(distanceThreshold, avgFaceDst);
    }

    if( args.options().found("normalDeviationAngleThreshold") )
    {
        deviationAngleThreshold =
            readScalar
            (
                IStringStream(args.options()["normalDeviationAngleThreshold"])()
            );
    }
    else
    {
        //- take average deviation as threshold
        scalar avgPointNormalDeviation(0.0), avgFaceNormalDeviation(0.0);

        label nPoints = bPoints.size();
        if( Pstream::parRun() )
        {
            const labelLongList& gpl = surfaceEngine.globalBoundaryPointLabel();
            forAll(gpl, bpI)
                nPoints = max(nPoints, gpl[bpI]);

            nPoints += 1;
        }
        reduce(nPoints, maxOp<label>());

        forAll(surfaceEngine.boundaryPoints(), bpI)
        {
            avgPointNormalDeviation += dist.angleDeviationAtBoundaryPoint(bpI);
        }

        avgPointNormalDeviation =
            returnReduce(avgPointNormalDeviation, sumOp<scalar>()) / nPoints;

        deviationAngleThreshold = avgPointNormalDeviation;

        //- calculate the average deviation of face normals
        const label nFaces = returnReduce(bFaces.size(), sumOp<label>());

        forAll(bFaces, bfI)
            avgFaceNormalDeviation += dist.angleDeviationAtBoundaryFace(bfI);

        avgFaceNormalDeviation =
            returnReduce(avgFaceNormalDeviation, sumOp<scalar>()) / nFaces;

        deviationAngleThreshold =
            max(deviationAngleThreshold, avgFaceNormalDeviation);
    }

    //- create subsets containing elements above the threshold
    const label pId =
        pmg.addPointSubset
        (
            "pointDistanceAbove_"+help::scalarToText(distanceThreshold)
        );
    const label pnId =
        pmg.addPointSubset
        (
            "pointNormalDeviationInRadians_" +
            help::scalarToText(deviationAngleThreshold)
        );

    //- find points above the threshold and store them into subsets
    forAll(bPoints, bpI)
    {
        if( dist.boundaryPointDistance(bpI) >= distanceThreshold )
            pmg.addPointToSubset(pId, bPoints[bpI]);
        if( dist.angleDeviationAtBoundaryPoint(bpI) >= deviationAngleThreshold )
            pmg.addPointToSubset(pnId, bPoints[bpI]);
    }

    //- calculate the average at face centres

    //- create subsets for faces above the threshold value
    const label fId =
        pmg.addFaceSubset
        (
            "faceDistanceAbove_"+help::scalarToText(distanceThreshold)
        );
    const label fnId =
        pmg.addFaceSubset
        (
            "faceNormalDeviationInRadians_" +
            help::scalarToText(deviationAngleThreshold)
        );

    //- find faces above the threshold
    forAll(bFaces, bfI)
    {
        if( dist.boundaryFaceCentreDistance(bfI) >= distanceThreshold )
            pmg.addFaceToSubset(fId, pmg.nInternalFaces() + bfI);
        if( dist.angleDeviationAtBoundaryFace(bfI) >= deviationAngleThreshold )
            pmg.addFaceToSubset(fnId, pmg.nInternalFaces() + bfI);
    }

    Info << "Writing mesh" << endl;
    pmg.write();

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
