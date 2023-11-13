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
#include "meshOptimizer.H"
#include "polyMeshGen.H"
#include "triSurf.H"
#include "coordinateModifier.H"
#include "surfaceMeshGeometryModification.H"
#include "polyMeshGenGeometryModification.H"

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

    const dictionary& anisotropicDict = meshDict.subDict("anisotropicSources");

    coordinateModifier cm(anisotropicDict);

    const scalar cellSize = readScalar(meshDict.lookup("maxCellSize"));

    const wordList toc = anisotropicDict.toc();

    forAll(toc, i)
    {
        Info << "Constructing modifier " << toc[i] << endl;

        autoPtr<coordinateModification> cModPtr =
            coordinateModification::New
            (
                toc[i],
                anisotropicDict.subDict(toc[i])
            );

        PtrList<plane> bndPlanes;
        cModPtr->boundingPlanes(bndPlanes);

        Info << "Num planes " << bndPlanes.size() << endl;

        label pairI(0);
        while( (2*pairI+1) < bndPlanes.size() )
        {
            const plane& pl0 = bndPlanes[2*pairI];
            const plane& pl1 = bndPlanes[2*pairI+1];

            const scalar scalingDistance =
                ((pl1.refPoint() - pl0.refPoint()) & pl0.normal());

            Info << "Scaling distance " << scalingDistance << endl;

            const vector dv = pl0.normal() * cellSize;

            point p = pl0.refPoint();
            do
            {
                const scalar t =
                    ((p - pl0.refPoint()) & pl0.normal()) / scalingDistance;
                Info << "t " << t << endl;

                const point dm = cModPtr->displacement(p);
                const point dbm = cModPtr->backwardDisplacement(p);

                Info << "Displacement " << dm << endl;
                Info << "Backward displacement " << dbm << endl;

                p += dv;
            } while( ((p - pl1.refPoint()) & pl1.normal()) <= SMALL );

            ++pairI;
        }
    }

    //return 0;

    polyMeshGen pmg(runTime);

    Info << "Starting reading mesh" << endl;
    pmg.read();
    Info << "Finished reading mesh" << endl;

    const triSurf* surfPtr(NULL);
    surfPtr = new triSurf(fileName(meshDict.lookup("surfaceFile")));

    surfaceMeshGeometryModification surfaceModification(*surfPtr, meshDict);

    const triSurf* modSurfPtr = surfaceModification.modifyGeometry();
    modSurfPtr->writeSurface("modifiedSurface.stl");

    const triSurf* revertedSurfPtr =
        surfaceModification.revertGeometryModification();
    revertedSurfPtr->writeSurface("revertedSurface.stl");

//    deleteDemandDrivenData(modSurfPtr);
    deleteDemandDrivenData(revertedSurfPtr);

    surfaceMeshGeometryModification backSurfModification(*modSurfPtr, meshDict);

    const triSurf* backModSurfPtr =
        backSurfModification.revertGeometryModification();
    backModSurfPtr->writeSurface("backModSurf.stl");

    //- check if the starting surface and the modified surface are the same
    const pointField& pts = surfPtr->points();
    const pointField& bPts = backModSurfPtr->points();
    forAll(pts, pointI)
    {
        const scalar dst = mag(pts[pointI] - bPts[pointI]);

        if( dst > SMALL )
            Info << "Point " << pointI << " does not match by " << dst << endl;
    }

//    modSurfPtr->writeSurface("modifiedSurface.s");
    deleteDemandDrivenData(modSurfPtr);
    deleteDemandDrivenData(backModSurfPtr);

//    //- check the point field
//    const pointFieldPMG& points = pmg.points();
//    pointField ptsCopy(points.size());
//    forAll(points, pI)
//        ptsCopy[pI] = points[pI];

//    polyMeshGenGeometryModification meshModification(pmg, meshDict);
//    meshModification.modifyGeometry();
//    //meshModification.revertGeometryModification();

//    forAll(ptsCopy, pI)
//    {
//        const scalar dist = mag(ptsCopy[pI] - points[pI]);

//        if( dist > SMALL )
//            Info << "Point " << pI << " does not match by " << dist << endl;
//    }

//    pmg.write();

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
