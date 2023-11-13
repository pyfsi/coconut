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

#include "demandDrivenData.H"
#include "partTetMesh.H"
#include "partTetMeshSimplex.H"
#include "polyMeshGenModifier.H"
#include "VRWGraphList.H"
#include "polyMeshGenAddressing.H"
#include "helperFunctions.H"
#include "OFstream.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

partTetMesh::partTetMesh(polyMeshGen& mesh, const labelLongList& lockedPoints)
:
    origMesh_(mesh),
    points_(mesh.points()),
    faces_(mesh.faces()),
    cells_(mesh.cells()),
    owner_(mesh.owner()),
    neighbour_(mesh.neighbour()),
    faceCentres_(mesh.addressingData().faceCentres()),
    cellCentres_(mesh.addressingData().cellCentres()),
    pointFaces_(mesh.addressingData().pointFaces()),
    smoothVertex_(),
    movableInternal_(),
    movableBoundary_(),
    bndFaceToPatch_(),
    parPoints_(),
    tetsAtPoint_(),
    bndTrianglesAtPoint_(),
    pointToParPoint_(),
    faceCentreToParPoint_(),
    cellCentreToParPoint_(),
    globalParPointLabel_(),
    globalToLocalParPointLabel_(),
    parPointAtOtherProcs_()
{
    List<direction> useCell(mesh.cells().size(), direction(1));

    boolList lockedPoint(mesh.points().size(), false);
    forAll(lockedPoints, i)
        lockedPoint[lockedPoints[i]] = true;

    createPointsAndTets(useCell, lockedPoint);
}

partTetMesh::partTetMesh
(
    polyMeshGen& mesh,
    const labelLongList& lockedPoints,
    const direction nLayers
)
:
    origMesh_(mesh),
    points_(mesh.points()),
    faces_(mesh.faces()),
    cells_(mesh.cells()),
    owner_(mesh.owner()),
    neighbour_(mesh.neighbour()),
    faceCentres_(mesh.addressingData().faceCentres()),
    cellCentres_(mesh.addressingData().cellCentres()),
    pointFaces_(mesh.addressingData().pointFaces()),
    smoothVertex_(),
    movableInternal_(),
    movableBoundary_(),
    bndFaceToPatch_(),
    parPoints_(),
    tetsAtPoint_(),
    bndTrianglesAtPoint_(),
    pointToParPoint_(),
    faceCentreToParPoint_(),
    cellCentreToParPoint_(),
    globalParPointLabel_(),
    globalToLocalParPointLabel_(),
    parPointAtOtherProcs_()
{
    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();
    const labelLongList& owner = mesh.owner();
    const PtrList<boundaryPatch>& boundaries = mesh.boundaries();
    const VRWGraph& pointCells = mesh.addressingData().pointCells();

    List<direction> useCell(cells.size(), direction(0));

    //- select cells containing at least one vertex of the bad faces
    forAll(boundaries, patchI)
    {
        const label start = boundaries[patchI].patchStart();
        const label size = boundaries[patchI].patchSize();

        for(label fI=0;fI<size;++fI)
        {
            useCell[owner[start+fI]] = 1;
        }
    }

    //- add additional layer of cells
    for(direction layerI=1;layerI<(nLayers+1);++layerI)
    {
        forAll(useCell, cI)
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
                            const label cLabel = pointCells(f[pI], pcI);
                            if( !useCell[cLabel] )
                                useCell[cLabel] = layerI + 1;
                        }
                    }
                }
            }

        if( Pstream::parRun() )
        {
            const labelLongList& globalPointLabel =
                mesh.addressingData().globalPointLabel();
            const VRWGraph& pProcs = mesh.addressingData().pointAtProcs();
            const Map<label>& globalToLocal =
                mesh.addressingData().globalToLocalPointAddressing();

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
                            eData[neiProc].append(globalPointLabel[pointI]);
                            break;
                        }
                }
            }

            //- exchange data with other processors
            labelLongList receivedData;
            help::exchangeMap(eData, receivedData);

            forAll(receivedData, i)
            {
                const label pointI = globalToLocal[receivedData[i]];

                forAllRow(pointCells, pointI, pcI)
                {
                    const label cLabel = pointCells(pointI, pcI);
                    if( !useCell[cLabel] )
                        useCell[cLabel] = layerI + 1;
                }
            }
        }
    }

    boolList lockedPoint(mesh.points().size(), false);
    forAll(lockedPoints, i)
        lockedPoint[lockedPoints[i]] = true;

    createPointsAndTets(useCell, lockedPoint);
}

partTetMesh::partTetMesh
(
    polyMeshGen& mesh,
    const labelLongList& lockedPoints,
    labelHashSet& badFaces,
    const direction additionalLayers
)
:
    origMesh_(mesh),
    points_(mesh.points()),
    faces_(mesh.faces()),
    cells_(mesh.cells()),
    owner_(mesh.owner()),
    neighbour_(mesh.neighbour()),
    faceCentres_(mesh.addressingData().faceCentres()),
    cellCentres_(mesh.addressingData().cellCentres()),
    pointFaces_(mesh.addressingData().pointFaces()),
    smoothVertex_(),
    movableInternal_(),
    movableBoundary_(),
    bndFaceToPatch_(),
    parPoints_(),
    tetsAtPoint_(),
    bndTrianglesAtPoint_(),
    pointToParPoint_(),
    faceCentreToParPoint_(),
    cellCentreToParPoint_(),
    globalParPointLabel_(),
    globalToLocalParPointLabel_(),
    parPointAtOtherProcs_()
{
    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();
    const VRWGraph& pointCells = mesh.addressingData().pointCells();

    if( Pstream::parRun() )
    {
        //- calculate MPI addressing outside of the openmp parallel region
        const polyMeshGenAddressing& addr =
            origMesh_.addressingData();
        addr.pointAtProcs();
        addr.globalToLocalPointAddressing();
    }

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
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(faces, faceI)
        {
            if( badFaces.found(faceI) )
            {
                const face& f = faces[faceI];

                forAll(f, pI)
                {
                    forAllRow(pointCells, f[pI], pcI)
                        useCell[pointCells(f[pI], pcI)] = 1;
                }
            }
        }

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

            for(direction layerI=1;layerI<(additionalLayers+1);++layerI)
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
                        origMesh_.addressingData();
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

    # ifdef DEBUGSmooth
    Info << "Adding cells into subset" << endl;
    label i(1);
    while
    (
        origMesh_.cellSubsetIndex("smoothCells_"+help::labelToText(i)) >= 0
    )
    {
        ++i;
    }
    const label cId =
        origMesh_.addCellSubset("smoothCells_"+help::labelToText(i));
    forAll(useCell, cellI)
        if( useCell[cellI] )
            origMesh_.addCellToSubset(cId, cellI);
    # endif

    //- mark locked points
    boolList lockedPoint(mesh.points().size(), false);
    forAll(lockedPoints, i)
        lockedPoint[lockedPoints[i]] = true;

    createPointsAndTets(useCell, lockedPoint);
}

partTetMesh::~partTetMesh()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


label partTetMesh::findInvertedPoints
(
    boolList& negativeNode,
    const boolList *activePointPtr
) const
{
    negativeNode.setSize(origMesh_.points().size());

    //- find points attached to the inverted tetrahedra
    label nInverted(0);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100) reduction(+:nInverted)
    # endif
    forAll(smoothVertex_, pointI)
    {
        if( activePointPtr && !activePointPtr->operator[](pointI) )
        {
            negativeNode[pointI] = false;
            continue;
        }

        negativeNode[pointI] = false;

        if( smoothVertex_[pointI] & (SMOOTH|BOUNDARY) )
        {
            partTetMeshSimplex simplex(*this, pointI);

            forAll(simplex.tets(), tI)
            {
                const scalar vol = simplex.tets()[tI].mag(simplex.pts());

                if( vol < VSMALL )
                {
                    negativeNode[pointI] = true;
                    ++nInverted;

                    break;
                }
            }
        }
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            origMesh_.addressingData().globalToLocalPointAddressing();
        const VRWGraph& pointAtProcs =
            origMesh_.addressingData().pointAtProcs();
        const DynList<label>& pNeiProcs =
            origMesh_.addressingData().pointNeiProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(pNeiProcs, i)
            exchangeData[pNeiProcs[i]].clear();

        //- send the global labels of inverted vertices
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label pointI = it();

            if( negativeNode[pointI] )
            {
                forAllRow(pointAtProcs, pointI, i)
                {
                    const label neiProc = pointAtProcs(pointI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append(it.key());
                }
            }
        }

        //- exchange data
        labelLongList receiveData;
        help::exchangeMap(exchangeData, receiveData);

        //- apply information locally
        forAll(receiveData, i)
            negativeNode[globalToLocal[receiveData[i]]] = true;
    }

    return returnReduce(nInverted, sumOp<label>());
}

void partTetMesh::updateVerticesSMP(const List<LongList<labelledPoint> >& np)
{
    boolList updateFace(owner_.size());

    polyMeshGenModifier meshModifier(origMesh_);

    # ifdef USE_OMP
    # pragma omp parallel num_threads(np.size())
    # endif
    {
        # ifdef USE_OMP
        const LongList<labelledPoint>& newPoints = np[omp_get_thread_num()];
        # else
        const LongList<labelledPoint>& newPoints = np[0];
        # endif

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(updateFace, fI)
            updateFace[fI] = false;

        forAll(newPoints, i)
        {
            const labelledPoint& lp = newPoints[i];
            const label pointI = lp.pointLabel();

            meshModifier.movePoint(pointI, lp.coordinates());

            forAllRow(pointFaces_, pointI, pfI)
                updateFace[pointFaces_(pointI, pfI)] = true;
        }
    }

    //- update face centres and cell centres
    origMesh_.addressingData().updateGeometry(updateFace);

    if( Pstream::parRun() )
    {
        unifyCoordinatesAtInterProcessorBoundaries();

        updateBufferLayerPoints();
    }
}

void partTetMesh::writeToVTK(const fileName&) const
{
    notImplemented("void partTetMesh::writeToVTK(const fileName&) const");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
