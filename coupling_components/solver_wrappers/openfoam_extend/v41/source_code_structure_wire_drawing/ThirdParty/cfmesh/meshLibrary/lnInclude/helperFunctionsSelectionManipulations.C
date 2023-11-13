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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "Pstream.H"

#include "helperFunctionsSelectionManipulations.H"
#include "helperFunctionsPar.H"
#include "polyMeshGen.H"
#include "polyMeshGenAddressing.H"

# ifdef USE_OMP
#include <omp.h>
# endif

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace help
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

label selectAdditionalLayersOfCells
(
    const polyMeshGen& mesh,
    const label nLayers,
    const labelLongList& activeCells,
    labelLongList& selectedCells
)
{
    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();
    const VRWGraph& pointCells = mesh.addressingData().pointCells();

    List<direction> useCell(cells.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        //- initialise useCell
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(useCell, cellI)
            useCell[cellI] = direction(0);

        //- select cells containing at least one vertex of the bad faces
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(activeCells, i)
            useCell[activeCells[i]] = 1;

        //- add additional layer of cells
        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            //- calculate the number of tasks
            # ifdef USE_OMP
            const label nTasks = 5 * omp_get_num_threads();
            # else
            const label nTasks = 1;
            # endif

            //- calculate chunk size for each task
            const label chSize = max(useCell.size() / nTasks, 1);

            for(direction layerI=1;layerI<(nLayers+1);++layerI)
            {
                for(label taskI=0;taskI<nTasks;++taskI)
                {
                    # ifdef USE_OMP
                    # pragma omp task default(shared) firstprivate(taskI)
                    # endif
                    {
                        const label sc = taskI * chSize;
                        label ec = min(sc + chSize, useCell.size());
                        if( taskI == (nTasks - 1) )
                            ec = useCell.size();

                        for(label cI=sc;cI<ec;++cI)
                        {
                            if( useCell[cI] == layerI )
                            {
                                const cell& c = cells[cI];

                                forAll(c, fI)
                                {
                                    const face& f = faces[c[fI]];

                                    forAll(f, pI)
                                    {
                                        forAllRow(pointCells, f[pI], pcI)
                                        {
                                            const label cLabel =
                                                pointCells(f[pI], pcI);
                                            if( !useCell[cLabel] )
                                            {
                                                useCell[cLabel] = layerI + 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                # ifdef USE_OMP
                # pragma omp taskwait
                # pragma omp flush(useCell)
                # endif

                if( Pstream::parRun() )
                {
                    const polyMeshGenAddressing& addr =
                        mesh.addressingData();
                    const VRWGraph& pProcs = addr.pointAtProcs();
                    const Map<label>& globalToLocal =
                        addr.globalToLocalPointAddressing();

                    std::map<label, labelLongList> eData;
                    forAllConstIter(Map<label>, globalToLocal, iter)
                    {
                        const label pointI = iter();

                        forAllRow(pProcs, pointI, procI)
                        {
                            const label neiProc = pProcs(pointI, procI);
                            if( neiProc == Pstream::myProcNo() )
                                continue;

                            if( eData.find(neiProc) == eData.end() )
                            {
                                eData.insert
                                (
                                    std::make_pair(neiProc, labelLongList())
                                );
                            }

                            forAllRow(pointCells, pointI, pcI)
                                if( useCell[pointCells(pointI, pcI)] == layerI )
                                {
                                    eData[neiProc].append(iter.key());
                                    break;
                                }
                       }
                    }

                    //- exchange data with other processors
                    labelLongList receivedData;
                    help::exchangeMap(eData, receivedData);

                    //- apply data locally
                    forAll(receivedData, i)
                    {
                        const label pointI = globalToLocal[receivedData[i]];

                        forAllRow(pointCells, pointI, pcI)
                        {
                            const label cLabel = pointCells(pointI, pcI);

                            if( !useCell[cLabel] )
                            {
                                useCell[cLabel] = layerI + 1;
                            }
                        }
                    }
                }
            }
        }
    }

    selectedCells.clear();
    forAll(useCell, i)
        if( useCell[i] )
            selectedCells.append(i);

    return selectedCells.size();
}

label selectAdditionalLayersOfCells
(
    polyMeshGen& mesh,
    const label nLayers,
    const label subsetId
)
{

    labelLongList selCells, aCells;

    mesh.cellsInSubset(subsetId, aCells);

    const label retVal =
        selectAdditionalLayersOfCells(mesh, nLayers, aCells, selCells);

    forAll(selCells, scI)
        mesh.addCellToSubset(subsetId, selCells[scI]);

    return retVal;
}

label selectAdditionalLayersOfCells
(
    polyMeshGen& mesh,
    const label nLayers,
    const word subsetName
)
{
    const label subsetId = mesh.cellSubsetIndex(subsetName);

    if( subsetId )
    {
        Warning << "Cell subset " << subsetName << " does not exist" << endl;
        return 0;
    }

    return selectAdditionalLayersOfCells(mesh, nLayers, subsetId);
}

label selectedFacesAttachedToPoints
(
    const polyMeshGen& mesh,
    const labelLongList& activePoints,
    labelLongList& selectedFaces
)
{
    selectedFaces.clear();

    const VRWGraph& pointFaces = mesh.addressingData().pointFaces();

    std::set<label> selFaces;

    forAll(activePoints, i)
    {
        const label pointI = activePoints[i];

        forAllRow(pointFaces, pointI, pfI)
        {
            selFaces.insert(pointFaces(pointI, pfI));
        }
    }

    forAllConstIter(std::set<label>, selFaces, it)
        selectedFaces.append(*it);

    return selectedFaces.size();
}

label createSubsetOfFacesAttachedToPoints
(
    polyMeshGen& mesh,
    const labelLongList& activePoints,
    const word subsetName
)
{
    labelLongList sf;
    const label retVal = selectedFacesAttachedToPoints(mesh, activePoints, sf);

    const label sId = mesh.addFaceSubset(subsetName);

    forAll(sf, i)
        mesh.addFaceToSubset(sId, sf[i]);

    return retVal;
}

label selectCellsAttachedToPoints
(
    const polyMeshGen& mesh,
    const labelLongList& activePoints,
    labelLongList& selectedCells
)
{
    const VRWGraph& pCells = mesh.addressingData().pointCells();

    std::set<label> selCells;

    forAll(activePoints, i)
    {
        const label pointI = activePoints[i];

        forAllRow(pCells, pointI, pcI)
            selCells.insert(pCells(pointI, pcI));
    }

    selectedCells.clear();

    forAllConstIter(std::set<label>, selCells, it)
        selectedCells.append(*it);

    return selectedCells.size();
}

label createSubsetOfCellsAttachedToPoints
(
    polyMeshGen& mesh,
    const labelLongList& activePoints,
    const word subsetName
)
{
    labelLongList sc;
    selectCellsAttachedToPoints(mesh, activePoints, sc);

    const label sId = mesh.addCellSubset(subsetName);

    forAll(sc, i)
        mesh.addCellToSubset(sId, sc[i]);

    return sc.size();
}

label selectCellsAttachedToFaces
(
    const polyMeshGen& mesh,
    const labelLongList& activeFaces,
    labelLongList& selectedCells
)
{
    const faceListPMG& faces = mesh.faces();
    const VRWGraph& pCells = mesh.addressingData().pointCells();

    std::set<label> selCells;

    forAll(activeFaces, i)
    {
        const label faceI = activeFaces[i];

        const face& f = faces[faceI];

        forAll(f, pI)
        {
            const label pointI = f[pI];

            forAllRow(pCells, pointI, pcI)
                selCells.insert(pCells(pointI, pcI));
        }
    }

    selectedCells.clear();
    forAllConstIter(std::set<label>, selCells, it)
        selectedCells.append(*it);

    return selectedCells.size();
}

label createSubsetOfCellsAttachedToFaces
(
    polyMeshGen& mesh,
    const labelLongList& activeFaces,
    const word subsetName
)
{
    labelLongList sc;
    selectCellsAttachedToFaces(mesh, activeFaces, sc);

    const label sId = mesh.addCellSubset(subsetName);

    forAll(sc, i)
        mesh.addCellToSubset(sId, sc[i]);

    return sc.size();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
