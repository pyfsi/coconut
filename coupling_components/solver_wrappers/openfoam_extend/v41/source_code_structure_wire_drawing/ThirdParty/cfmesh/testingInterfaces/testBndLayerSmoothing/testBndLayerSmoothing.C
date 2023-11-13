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
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"
#include "polyMeshGen.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "meshOptimizer.H"
#include "boolList.H"

#include "boundaryLayerOptimisation.H"
#include "extrudeLayer.H"
#include "polyMeshGenChecks.H"
#include "HashSet.H"

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

    polyMeshGen pmg(runTime);

    Info << "Starting reading mesh" << endl;
    pmg.read();

    Info << "Finished reading mesh" << endl;

    Info << "Starting optimising mesh" << endl;
//    meshOptimizer mOpt(pmg);
//    mOpt.lockCellsInSubset("boundaryLayerCells");
//    mOpt.optimizeLowQualityFaces(15);
//    mOpt.untangleMeshFV(5, 50, 0);

    meshOptimizer mOpt(pmg);
    mOpt.optimizeBoundaryLayer();

    polyMeshGenChecks::checkTopology(pmg, true);

    forAll(pmg.points(), pointI)
    {
        const point& p = pmg.points()[pointI];

        if( help::isnan(p) )
            Info << "Vertex " << pointI << " is invalid " << p << endl;
    }

    pmg.write();

    return 0;

    forAll(pmg.points(), pointI)
    {
        const point& p = pmg.points()[pointI];

        if( help::isnan(p) )
            Info << "Vertex " << pointI << " is invalid " << p << endl;
    }

    //boundaryLayers(pmg).addLayerForAllPatches();
    //pmg.clearAddressingData();

    meshSurfaceEngine mse(pmg);
    boundaryLayerOptimisation blOpt(pmg, mse);

    const labelLongList& faceCell = mse.faceOwners();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const boolList& isBaseFace = blOpt.isBaseFace();

    const labelLongList& owner = pmg.owner();
    const labelLongList& neighbour = pmg.neighbour();
    const cellListPMG& cells = pmg.cells();
    const faceListPMG& faces = pmg.faces();

    Info << "Marking boundary layer cells" << endl;

    labelLongList layerCells;
    boolList layerCell(pmg.cells().size(), false);
    const label blCellsId = pmg.addCellSubset("boundaryLayerCells");
    forAll(isBaseFace, bfI)
    {
        pmg.addCellToSubset(blCellsId, faceCell[bfI]);
        layerCell[faceCell[bfI]] = true;
        layerCells.append(faceCell[bfI]);
    }

    //- marking faces at the inner boundary of the boundary layer
    LongList<labelPair> front;

    forAll(faceCell, bfI)
    {
        const cell& c = cells[faceCell[bfI]];

        const face& bf = bFaces[bfI];

        label faceOpposite(-1);

        forAll(c, fI)
            if( !help::shareAnEdge(faces[c[fI]], bf) )
                faceOpposite = c[fI];

        label cellI = owner[faceOpposite];
        if( cellI == faceCell[bfI] )
            cellI = neighbour[faceOpposite];

        if( layerCell[cellI] )
            continue;

        front.append(labelPair(faceOpposite, cellI));
    }

    //- find points in boundary layer
    boolList pointInBoundaryLayer(pmg.points().size(), false);

    forAll(cells, cellI)
    {
        if( layerCell[cellI] )
        {
            const cell& c = cells[cellI];
            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                forAll(f, pI)
                    pointInBoundaryLayer[f[pI]] = true;
            }
        }
    }

    //- check if there exist faces with all vertices in the boundary layer
    const label allPointsInLayerId = pmg.addCellSubset("allPointsInLayer");
    forAll(cells, cellI)
    {
        if( !layerCell[cellI] )
        {
            bool allInBndLayer(true);

            const cell& c = cells[cellI];
            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                forAll(f, pI)
                    allInBndLayer = false;
            }

            if( allInBndLayer )
            {
                Info << "Cell " << cellI
                     << " has all points in the layer" << endl;
                pmg.addCellToSubset(allPointsInLayerId, cellI);
            }
        }
    }

    //Info << "Optimising boundary layer" << endl;
    //blOpt.optimiseLayer();

    pmg.clearAddressingData();
/*
    //- refine boundary layers
    refineBoundaryLayers refLayers(pmg);

    refineBoundaryLayers::readSettings(meshDict, refLayers);

    refLayers.refineLayers();
*/
/*
    //- check bad quality cells in the layer
    boolList activeFace(pmg.faces().size(), false);
    if( blCellsId > -1 )
    {
        labelLongList cellsInLayer;
        pmg.cellsInSubset(blCellsId, cellsInLayer);

        forAll(cellsInLayer, i)
        {
            const cell& c = cells[cellsInLayer[i]];

            forAll(c, fI)
                activeFace[c[fI]] = true;
        }
    }

    labelHashSet badFaces;
    polyMeshGenChecks::findBadFaces(pmg, badFaces, true, &activeFace);

    if( returnReduce(badFaces.size(), sumOp<label>()) )
    {
        const labelLongList& own = pmg.owner();
        const labelLongList& nei = pmg.neighbour();

        label badCellId(-1);
        forAllConstIter(labelHashSet, badFaces, it)
        {
            if( badCellId < 0 )
                badCellId = pmg.addCellSubset("badBlCells");
            pmg.addCellToSubset(badCellId, own[it.key()]);

            if( nei[it.key()] >= 0 )
                pmg.addCellToSubset(badCellId, nei[it.key()]);
        }

        if( returnReduce(pmg.cellSubsetIndex("badBlCells")>=0, maxOp<bool>()) )
        {
            Info << "Found bad quality bl cells" << endl;
            pmg.write();
            returnReduce(1, sumOp<label>());
            ::exit(0);
        }
    }
*/
    Info << "Extruding layer of cells" << endl;
    extrudeLayer(pmg, front);
    pmg.clearAddressingData();


    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
