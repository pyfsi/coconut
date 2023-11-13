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
#include "boolList.H"
#include "meshSurfaceEngine.H"
#include "polyMeshGenModifier.H"
#include "triSurfModifier.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "helperFunctions.H"

#include "Random.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("patch1");
    argList::validArgs.append("patch2");
    argList::validArgs.append("action");
    argList::validArgs.append("nIterations");

#   include "setRootCase.H"
#   include "createTime.H"

    //- read arguments
    const word patchName1(args.additionalArgs()[0]);
    const word patchName2(args.additionalArgs()[1]);
    const word action(args.additionalArgs()[2]);
    const int nIterations(readInt(IStringStream(args.additionalArgs()[3])()));

    polyMeshGen pmg(runTime);
    pmg.read();

    Info << "Detecting active patches" << endl;
    label patch1(-1), patch2(-1);
    forAll(pmg.boundaries(), patchI)
    {
        if( pmg.boundaries()[patchI].patchName() == patchName1 )
        {
            patch1 = patchI;
        }
        else if( pmg.boundaries()[patchI].patchName() == patchName2 )
        {
            patch2 = patchI;
        }
    }

    if( patch1 < 0 )
    {
        FatalError << "Patch " << patchName1 << " is not found"
            << exit(FatalError);
    }
    if( patch2 < 0 )
    {
        FatalError << "Patch " << patchName2 << " is not found"
            << exit(FatalError);
    }

    Info << "Creating mesh surface" << endl;
    meshSurfaceEngine mse(pmg);

    const pointFieldPMG& points = mse.points();
    const labelLongList& bp = mse.bp();
    const labelLongList& bPoints = mse.boundaryPoints();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelLongList& facePatch = mse.boundaryFacePatches();

    labelLongList activeFaceLabel(bFaces.size(), -1);
    labelLongList activePointLabel(bPoints.size(), -1);

    label nSurfPoints(0);
    label nSurfFaces(0);

    forAll(facePatch, bfI)
    {
        if( facePatch[bfI] == patch1 || facePatch[bfI] == patch2 )
        {
            activeFaceLabel[bfI] = nSurfFaces++;

            const face& bf = bFaces[bfI];
            forAll(bf, pI)
            {
                const label bpI = bp[bf[pI]];

                if( activePointLabel[bpI] < 0 )
                {
                    activePointLabel[bpI] = nSurfPoints++;
                }
            }
        }
    }

    //- create surface mesh
    Info << "Creating surface mesh" << endl;
    triSurf activeSurf;
    triSurfModifier sMod(activeSurf);
    pointField& sPts = sMod.pointsAccess();
    geometricSurfacePatchList& patches = sMod.patchesAccess();
    patches.setSize(2);
    patches[0].name() = pmg.boundaries()[patch1].patchName();
    patches[1].name() = pmg.boundaries()[patch2].patchName();

    sPts.setSize(nSurfFaces+nSurfPoints);
    forAll(activePointLabel, bpI)
    {
        if( activePointLabel[bpI] >= 0 )
        {
            sPts[activePointLabel[bpI]] = points[bPoints[bpI]];
        }
    }

    std::set<label> pointsInPatch1, pointsInPatch2;
    std::set<label> faceInPatch1, faceInPatch2;

    forAll(bFaces, bfI)
    {
        if( activeFaceLabel[bfI] >= 0 )
        {
            const label cLabel = nSurfPoints + activeFaceLabel[bfI];

            const face& bf = bFaces[bfI];

            label patchID(-1);

            if( facePatch[bfI] == patch1 )
            {
                patchID = 0;
                faceInPatch1.insert(bfI);

                forAll(bf, pI)
                    pointsInPatch1.insert(activePointLabel[bp[bf[pI]]]);
            }
            else if( facePatch[bfI] == patch2 )
            {
                patchID = 1;
                faceInPatch2.insert(bfI);

                forAll(bf, pI)
                    pointsInPatch2.insert(activePointLabel[bp[bf[pI]]]);
            }

            sPts[cLabel] = help::faceCentre(points, bFaces[bfI]);

            //- create triangles
            forAll(bf, eI)
            {
                activeSurf.appendTriangle
                (
                    labelledTri
                    (
                        activePointLabel[bp[bf[eI]]],
                        activePointLabel[bp[bf.nextLabel(eI)]],
                        cLabel,
                        patchID
                    )
                );
            }
        }
    }

    scalarList pointRange(sPts.size(), 0.0);
    scalarList faceRange(nSurfFaces, 0.0);

    forAll(activeSurf, tI)
    {
        const labelledTri& tri = activeSurf[tI];

        const scalar r1 = mag(sPts[tri[0]] - sPts[tri[2]]);
        const scalar r2 = mag(sPts[tri[1]] - sPts[tri[2]]);

        pointRange[tri[0]] = max(pointRange[tri[0]], r1);
        pointRange[tri[1]] = max(pointRange[tri[1]], r2);

        scalar& fr = faceRange[tri[2] - nSurfPoints];
        fr = max(fr, r1);
        fr = max(fr, r2);
    }

    Info << "Max point range " << max(pointRange) << endl;
    Info << "Max face range " << max(faceRange) << endl;

    //- construct octree for searches
    Info << "Creating octree" << endl;
    meshOctree octree(activeSurf);
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(20, 15);

    if( action == "pointDistances" )
    {
        Info << "Starting calculating point distances "
             << runTime.elapsedCpuTime() << endl;

        for(label iterI=0;iterI<nIterations;++iterI)
        {
            forAllConstIter(std::set<label>, pointsInPatch1, it)
            {
                const point& p = sPts[*it];

                point np;
                scalar dSq;
                label nt;
                octree.findNearestSurfacePointInRegion
                (
                    np,
                    dSq,
                    nt,
                    1,
                    p,
                    pointRange[*it]
                );
            }

            forAllConstIter(std::set<label>, pointsInPatch2, it)
            {
                const point& p = sPts[*it];

                point np;
                scalar dSq;
                label nt;
                octree.findNearestSurfacePointInRegion
                (
                    np,
                    dSq,
                    nt,
                    0,
                    p,
                    pointRange[*it]
                );
            }

            Info << "Iteration " << iterI << " calculated point distances "
             << runTime.elapsedCpuTime() << endl;
        }
    }
    else if( action == "fieldInterpolation" )
    {
        Info << "Starting calculating face distances "
             << runTime.elapsedCpuTime() << endl;

        for(label iterI=0;iterI<nIterations;++iterI)
        {
            forAllConstIter(std::set<label>, faceInPatch1, it)
            {
                const label fI = activeFaceLabel[*it];
                const point& p = sPts[fI+nSurfPoints];

                point np;
                scalar dSq;
                label nt;
                octree.findNearestSurfacePointInRegion
                (
                    np,
                    dSq,
                    nt,
                    1,
                    p,
                    faceRange[fI]
                );
            }

            forAllConstIter(std::set<label>, faceInPatch2, it)
            {
                const label fI = activeFaceLabel[*it];
                const point& p = sPts[fI+nSurfPoints];

                point np;
                scalar dSq;
                label nt;
                octree.findNearestSurfacePointInRegion
                (
                    np,
                    dSq,
                    nt,
                    0,
                    p,
                    faceRange[fI]
                );
            }

            Info << "Iteration " << iterI << " calculated face distances "
             << runTime.elapsedCpuTime() << endl;
        }
    }

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
