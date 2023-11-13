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

    polyMeshGen pmg(runTime);
    pmg.read();

    //- create a point subsets of points with backup coordinates
    const label backupId = pmg.addPointSubset("backupPoints");
    forAll(pmg.points(), pointI)
        if( pmg.hasBackup(pointI) )
            pmg.addPointToSubset(backupId, pointI);

    try
    {
        meshOptimizer mmOpt(pmg);
        mmOpt.optimizeMeshFVWithTopoChanges(10, 50, 5);
        //mmOpt.lockPoints(lockedPoints);
        //mmOpt.optimizeMeshFV(5, 10, 50, 5);
        //mmOpt.optimizeBoundaryLayer();
        //mmOpt.untangleMeshFV(10, 50, 5);
    }
    catch(...)
    {
        Pout << "Finished smoothing" << endl;
    }

    if( pmg.faceSubsetIndex("badFaces") >= 0 )
    {
        labelLongList facesInSubset;
        pmg.facesInSubset(pmg.faceSubsetIndex("badFaces"), facesInSubset);

        help::createSubsetOfCellsAttachedToFaces
        (
            pmg,
            facesInSubset,
            "badCells"
        );

        help::selectAdditionalLayersOfCells(pmg, 2, "badCells");
    }

    //polyMeshGenChecks::checkGeometry(pmg, true);

/*
    for(label i=0;i<10;++i)
    {
    labelHashSet badFaces;
    const label nBadFaces =
        polyMeshGenChecks::findBadFaces
        (
            pmg,
            badFaces,
            false
        );

    Info << "Number of bad-quality faces" << nBadFaces << endl;

    partTetMesh tm(pmg, labelLongList(), badFaces, 2);
    tetMeshOptimisation tmo(tm);

    tmo.optimiseUsingHeightOptimizer();

    tm.updateOrigMesh();
    tm.writeToVTK("tetMesh_"+help::labelToText(Pstream::myProcNo())+".vtk");
    }
*/
    pmg.write();
    return 0;

    if( pmg.hasPointsBackup() )
    {
        label nIdentical(0), noBackup(0);

        const label pId = pmg.addPointSubset("noBackup");

        forAll(pmg.points(), pointI)
        {
            if( !pmg.hasBackup(pointI) )
            {
                ++noBackup;
                pmg.addPointToSubset(pId, pointI);
                continue;
            }

            point pOrig;
            pmg.getOrigPoint(pointI, pOrig);

            if( mag(pOrig - pmg.points()[pointI]) < SMALL )
                ++nIdentical;
//            Info << "Point " << pointI << " coordinates "
//                 << pmg.points()[pointI] << " orig coordinates "
//                 << pOrig << endl;
        }

        Info << "Number of points with no backup " << noBackup << endl;
        Info << "Number of identical backup points " << nIdentical << endl;

        //- store boundary layer cells into a subset
        const label blId = pmg.addCellSubset("bndLayerCells");
        const labelLongList& owner = pmg.owner();
        forAll(pmg.boundaries(), patchI)
        {
            const label start = pmg.boundaries()[patchI].patchStart();
            const label end = start + pmg.boundaries()[patchI].patchSize();

            for(label faceI=start;faceI<end;++faceI)
                pmg.addCellToSubset(blId, owner[faceI]);
        }

        pmg.write();
        ::exit(0);
    }

//    {
//        //- create a tet mesh consisting of a single tet
//        polyMeshGenModifier meshModifier(pmg);

//        pointFieldPMG& points = meshModifier.pointsAccess();
//        faceListPMG& faces = meshModifier.facesAccess();
//        cellListPMG& cells = meshModifier.cellsAccess();
//        PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();

//        points.setSize(4);
//        points[0] = point(0., 0., 0.);
//        points[1] = point(1., 0., 0.);
//        points[2] = point(0.5, Foam::sqrt(3.0)/2.0, 0.);
//        points[3] = point(0., 0., -1.0);

//        faces.setSize(4);
//        face f(3);

//        f[0] = 0;
//        f[1] = 2;
//        f[2] = 1;
//        faces[0] = f;

//        f[0] = 0;
//        f[1] = 1;
//        f[2] = 3;
//        faces[1] = f;

//        f[0] = 1;
//        f[1] = 2;
//        f[2] = 3;
//        faces[2] = f;

//        f[0] = 2;
//        f[1] = 0;
//        f[2] = 3;
//        faces[3] = f;

//        cells.setSize(1);
//        cells[0].setSize(4);
//        for(label i=0;i<4;++i)
//            cells[0][i] = i;

//        boundaries.clear();
//        boundaries.setSize(1);
//        boundaries.set
//        (
//            0,
//            new boundaryPatch
//            (
//                "bnd",
//                "patch",
//                4,
//                0
//            )
//        );

//        meshModifier.clearAll();

//        polyMeshGenChecks::checkMesh(pmg, true);
//        Info << "Finished generating single tet" << endl;
//    }

/*
    {
        //- create a backup of point coordinates
        polyMeshGenModifier meshModifier(pmg);
        meshModifier.backupPoints();

        //- generate a boundary layer
        boundaryLayers(pmg).addLayerForPatch("layerPatch");

        //- add random noise to points
        pointFieldPMG& points = meshModifier.pointsAccess();

        point maxPoint(points[0]), minPoint(points[0]);
        for(label pI=1;pI<points.size();++pI)
        {
            maxPoint = max(maxPoint, points[pI]);
            minPoint = min(minPoint, points[pI]);
        }

        Random rand(points.size());
        forAll(points, pI)
        {
            points[pI] += mag(maxPoint - minPoint) * rand.vector01();
        }
    }
*/

//    labelLongList lockedPoints;
//    lockedPoints.setSize(3);
//    for(label i=0;i<3;++i)
//        lockedPoints[i] = i;
//    lockedPoints.append(0);
//    lockedPoints.append(210);
//    lockedPoints.append(30);
//    lockedPoints.append(80);

/*    Info << "Creating tetrahedral mesh " << endl;
    partTetMesh tetMesh(pmg, lockedPoints);
    Info << "Constructing tetMeshOptimisation" << endl;
    tetMeshOptimisation tmo(tetMesh);
    Info << "Starting smoothing the mesh" << endl;
    tmo.optimiseBoundaryVolumeOptimizer(100);
    //tmo.optimiseUsingHeightOptimizer(20000);
    //tmo.optimiseUsingVolumeOptimizer(200);

    tetMesh.updateOrigMesh();
*/
/*
    List<LongList<labelledPoint> > newPositions(1);
    label nIter(0);
    do
    {
        newPositions[0].clear();

        scalar maxDisp(0.0);

        forAll(tetMesh.nodeLabelInOrigMesh(), pI)
        {
            if( tetMesh.nodeLabelInOrigMesh()[pI] < 0 )
                continue;
            if( tetMesh.smoothVertex()[pI] & partTetMesh::LOCKED )
                continue;

            partTetMeshSimplex simplex(tetMesh, pI);

            const DynList<point, 128>& pts = simplex.pts();
            const DynList<partTet, 128>& tets = simplex.tets();

            const point op = simplex.centrePoint();

            symmTensor mat(symmTensor::zero), I(1.0, 0., 0., 1.0, 0., 1.0);
            vector source(vector::zero);

            scalar targetLength = 0.5;

            forAll(tets, tetI)
            {
                const partTet& t = tets[tetI];
                const FixedList<edge, 6> edges = t.edges();

                vector n = (pts[t[1]] - pts[t[0]]) ^ (pts[t[2]] - pts[t[0]]);
                n /= (mag(n) + VSMALL);

//                Info << "Tet " << t << endl;
//                Info << "n " << n << endl;

//                scalar maxEdge(0.0);
//                forAll(edges, edgeI)
//                {
//                    const edge& e = edges[edgeI];
//                    maxEdge = max(maxEdge, mag(pts[e[0]] - pts[e[1]]));
//                }

                const point fc = (pts[t[0]] + pts[t[1]] + pts[t[2]]) / 3.0;
                const scalar hIdeal = 0.81649 * targetLength;//maxEdge;
                //const scalar hIdeal = 1.0 * maxEdge;

                //disp += (fc + hIdeal * n) - op;
                mat += symm(n * n);
//                scalar l = (op - pts[t[0]]) & n;
//                mat += (targetLength - l) / (mag(l) + VSMALL) * pts[t[0]];
//                l = (op - pts[1])
                source +=
//                    (targetLength / lpts[t[0]] + pts[t[1]] + pts[t[2]] +
                    hIdeal * n + ((n * n) & fc);
            }

            Info << "Mat " << mat << endl;
            Info << "Source " << source << endl;
            const point nc = inv(mat) & source;

            //volumeOptimizer vOpt(simplex);
            //vOpt.optimizeNodePosition(1e-10);

            //const point nc = op + 0.5 * disp;//simplex.centrePoint();

            Info << "Displacement of point " << pI << " is " << (nc - op) << endl;

            maxDisp = max(maxDisp, mag(nc - op));

            newPositions[0].append(labelledPoint(pI, op + 0.5 *(nc - op)));
        }

        tetMesh.updateVerticesSMP(newPositions);

        if( maxDisp < 1e-4 )
            break;
    } while( true );
*/

    //tetMesh.createPolyMesh(pmg);
    polyMeshGenChecks::checkMesh(pmg);

//    tetMeshOptimisation tmo(tetMesh);

//    tmo.optimiseBoundaryVolumeOptimizer(100);
//    tetMesh.updateOrigMesh();

//    pmg.write();
//    return 0;

//    {
//        meshSurfaceEngine mse(pmg);

//        meshSurfaceOptimizer mso(mse);

//        mso.optimizeSurface();

//        pmg.write();
//        ::exit(0);
//    }

    meshOptimizer mOpt(pmg);
    //mOpt.optimizeMeshNearBoundaries();
    mOpt.untangleMeshFV(20, 100, 5);
    //mOpt.optimizeMeshFVBestQuality(100, 0.99);
    //mOpt.untangleMeshFV();
    //mOpt.optimizeMeshFV();

    //if( polyMeshGenChecks::checkGeometry(pmg, true) )
    //{
        //Info << "Reverting points" << endl;
        //polyMeshGenModifier(pmg).revertPoints();

        //Info << "Untangling again" << endl;
        //mOpt.untangleMeshFV();
    //}

   // polyMeshGenModifier(pmg).clearPointBackup();
    pmg.write();

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
