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

#include "checkNonMappableCellConnections.H"
#include "polyMeshGenModifier.H"
#include "helperFunctions.H"
#include "meshSurfaceEngine.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGCheck

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

void checkNonMappableCellConnections::findCellTypes()
{
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();

    meshSurfaceEngine mse(mesh_);
    const labelLongList& bp = mse.bp();
    const labelLongList& faceOwner = mse.faceOwners();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const edgeLongList& edges = mse.edges();

    cellType_.setSize(cells.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(cellType_, cellI)
            cellType_[cellI] = INTERNALCELL;

        //- find boundary cells
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(faceOwner, bfI)
            cellType_[faceOwner[bfI]] = BNDCELL;

        //- find boundary cells with all vertices at the boundary
        # ifdef USE_OMP
        # pragma omp for schedule(guided, 100)
        # endif
        for(label cellI=cells.size()-1;cellI>=0;--cellI)
        {
            if( cellType_[cellI] & INTERNALCELL )
                continue;

            const cell& c = cells[cellI];

            //- mark boundary cells with all vertices at the boundary
            bool allBoundary(true);

            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                //- check if all vertices are located at the surface of the mesh
                forAll(f, pI)
                    if( bp[f[pI]] == -1 )
                        allBoundary = false;
            }

            if( allBoundary )
            {
                cellType_[cellI] |= ALLBNDVERTEXCELL;
            }
            else if( !(cellType_[cellI] & BNDCELL) )
            {
                continue;
            }

            //- check if the internal faces are connected into a single group
            //- over their edges
            DynList<label> internalFaces;
            forAll(c, fI)
            {
                if( c[fI] < mesh_.nInternalFaces() )
                {
                    internalFaces.append(c[fI]);
                }
                else if( mesh_.faceIsInProcPatch(c[fI]) != -1 )
                {
                    internalFaces.append(c[fI]);
                }
            }

            bool hasInternalEdge(false);
            forAll(internalFaces, i)
            {
                const face& f = faces[internalFaces[i]];

                forAll(f, eI)
                {
                    const edge e = f.faceEdge(eI);

                    const label bps = bp[e.start()];
                    const label bpe = bp[e.end()];

                    bool isInternal(true);
                    if( bps != -1 )
                    {
                        forAllRow(bpEdges, bps, peI)
                        {
                            if( edges[bpEdges(bps, peI)] == e )
                            {
                                isInternal = false;
                                break;
                            }
                        }
                    }
                    else if( bpe != -1 )
                    {
                        forAllRow(bpEdges, bpe, peI)
                        {
                            if( edges[bpEdges(bpe, peI)] == e )
                            {
                                isInternal = false;
                                break;
                            }
                        }
                    }

                    if( isInternal )
                    {
                        hasInternalEdge = true;
                        break;
                    }
                }
            }

            if( hasInternalEdge )
            {
                cellType_[cellI] |= HASINTERNALEDGE;
            }
        }
    }
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //
// Constructors

checkNonMappableCellConnections::checkNonMappableCellConnections
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    cellType_()
{}

// * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * * * * //

checkNonMappableCellConnections::~checkNonMappableCellConnections()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void checkNonMappableCellConnections::findCells(labelHashSet& badCells)
{
    badCells.clear();

    //- classify cell types
    findCellTypes();

    //- select ALLBNDVERTEXCELL and INTERNALFACEGROUP cells
    //- with at least one INTERNALCELL neighbour
    //- these cells do not need to stay in the mesh
    const cellListPMG& cells = mesh_.cells();
    const labelLongList& owner = mesh_.owner();
    const labelLongList& neighbour = mesh_.neighbour();
    const label nInternalFaces = mesh_.nInternalFaces();
    const PtrList<processorBoundaryPatch>& procBoundaries = mesh_.procBoundaries();

    labelListList otherProcType;
    if( Pstream::parRun() )
    {
        //- exchange cell types at processor boundaries
        otherProcType.setSize(procBoundaries.size());

        //- send data to other processors
        forAll(procBoundaries, patchI)
        {
            label start = procBoundaries[patchI].patchStart();
            labelList patchCellType(procBoundaries[patchI].patchSize());

            forAll(patchCellType, faceI)
                patchCellType[faceI] = cellType_[owner[start++]];

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                patchCellType.byteSize()
            );

            toOtherProc << patchCellType;
        }

        //- receive data from other processors
        forAll(procBoundaries, patchI)
        {
            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            labelList& otherTypes = otherProcType[patchI];
            fromOtherProc >> otherTypes;
        }
    }

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100)
    # endif
    forAll(cellType_, cellI)
    {
        if( cellType_[cellI] & ALLBNDVERTEXCELL )
        {
            //- mark cells that do not have a neighbour with at least one
            //- internal edge
            const cell& c = cells[cellI];

            bool hasNeighbourWithInternalEdge(false);
            label nNeiCells(0);

            forAll(c, fI)
            {
                const label faceI = c[fI];

                if( faceI < nInternalFaces )
                {
                    ++nNeiCells;

                    label nei = neighbour[c[fI]];
                    if( nei == cellI )
                        nei = owner[c[fI]];

                    if( cellType_[nei] & HASINTERNALEDGE )
                    {
                        hasNeighbourWithInternalEdge = true;
                        break;
                    }
                }
                else if( mesh_.faceIsInProcPatch(faceI) != -1 )
                {
                    ++nNeiCells;

                    const label patchI = mesh_.faceIsInProcPatch(faceI);
                    const label j = faceI - procBoundaries[patchI].patchStart();

                    if( otherProcType[patchI][j] & HASINTERNALEDGE )
                    {
                        hasNeighbourWithInternalEdge = true;
                        break;
                    }
                }
            }

            if( !hasNeighbourWithInternalEdge )
            {
                # ifdef USE_OMP
                # pragma omp critical(noNeighbourWithInternalEdge)
                # endif
                badCells.insert(cellI);
            }
        }
    }
}

void checkNonMappableCellConnections::lockPoints()
{
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();
    polyMeshGenModifier meshModifier(mesh_);

    labelHashSet badCells;

    findCells(badCells);

    forAllConstIter(labelHashSet, badCells, it)
    {
        const cell& c = cells[it.key()];

        forAll(c, fI)
        {
            const face& f = faces[c[fI]];

            forAll(f, pI)
                meshModifier.lockPoint(f[pI]);
        }
    }

    const label pId = mesh_.addPointSubset("lockedPoints");
    forAll(mesh_.points(), pointI)
        if( mesh_.isLockedPoint(pointI) )
            mesh_.addPointToSubset(pId, pointI);
}

bool checkNonMappableCellConnections::removeCells()
{
    labelHashSet badCells;

    label nRemoved;
    bool changed(false);

    do
    {
        findCells(badCells);

        nRemoved = badCells.size();
        reduce(nRemoved, sumOp<label>());

        Info << "Found " << nRemoved << " non-mappable cells" << endl;

        if( nRemoved != 0 )
        {
            boolList removeCell(mesh_.cells().size(), false);
            forAllConstIter(labelHashSet, badCells, it)
                removeCell[it.key()] = true;

            polyMeshGenModifier(mesh_).removeCells(removeCell);

            changed = true;
        }
    } while( nRemoved );

    return changed;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
