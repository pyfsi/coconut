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
#include "checkCellConnectionsOverFaces.H"
#include "meshOptimizer.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGenChecks.H"
#include "labelLongList.H"
#include "helperFunctions.H"

//#define DEBUGSmoothing

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * Private member functions * * * * * * * * * * * * * * * * * //

const meshSurfaceEngine& meshOptimizer::meshSurface() const
{
    if( !msePtr_ )
        msePtr_ = new meshSurfaceEngine(mesh_);

    return *msePtr_;
}

void meshOptimizer::clearSurface()
{
    deleteDemandDrivenData(msePtr_);
}

void meshOptimizer::calculatePointLocations()
{
    vertexLocation_.setSize(mesh_.points().size());
    vertexLocation_ = INSIDE;

    const meshSurfaceEngine& mse = meshSurface();
    const labelLongList& bPoints = mse.boundaryPoints();

    //- mark boundary vertices
    forAll(bPoints, bpI)
        vertexLocation_[bPoints[bpI]] = BOUNDARY;

    //- mark edge vertices
    meshSurfacePartitioner mPart(mse);
    forAllConstIter(labelHashSet, mPart.edgePoints(), it)
        vertexLocation_[bPoints[it.key()]] = EDGE;

    //- mark corner vertices
    forAllConstIter(labelHashSet, mPart.corners(), it)
        vertexLocation_[bPoints[it.key()]] = CORNER;

    if( Pstream::parRun() )
    {
        const polyMeshGenAddressing& addresing = mesh_.addressingData();
        const VRWGraph& pointAtProcs = addresing.pointAtProcs();

        forAll(pointAtProcs, pointI)
            if( pointAtProcs.sizeOfRow(pointI) != 0 )
                vertexLocation_[pointI] |= PARALLELBOUNDARY;
    }
}

void meshOptimizer::revertPointLocations
(
    const labelHashSet& badFaces,
    boolList& changedFace,
    const scalar alpha,
    labelLongList& lockedPoints
)
{
    //- copy locked points into a set
    std::set<label> lockedVertices;
    forAll(lockedPoints, i)
        lockedVertices.insert(lockedPoints[i]);

    //- mark points affected by bad faces
    const faceListPMG& faces = mesh_.faces();

    std::set<label> affectedPoints;
    forAllConstIter(labelHashSet, badFaces, bfIt)
    {
        const face& f = faces[bfIt.key()];

        forAll(f, pI)
            affectedPoints.insert(f[pI]);
    }

    if( Pstream::parRun() )
    {
        //- distribute information accross all processors
        const polyMeshGenAddressing& addressing = mesh_.addressingData();
        const Map<label>& globalToLocal =
            addressing.globalToLocalPointAddressing();
        const VRWGraph& pointAtProcs = addressing.pointAtProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(addressing.pointNeiProcs(), i)
            exchangeData[addressing.pointNeiProcs()[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label pointI = it();

            if( affectedPoints.find(pointI) != affectedPoints.end() )
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

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
            affectedPoints.insert(globalToLocal[receivedData[i]]);
    }

    //- move vertices towards their original coordinates
    const meshSurfaceEngine& mse = meshSurface();
    const labelLongList& bp = mse.bp();

    # ifdef DEBUGSmoothing
    label pId = mesh_.pointSubsetIndex("surfPointsNoBackup");
    if( pId < 0 )
        pId = mesh_.addPointSubset("surfPointsNoBackup");
    forAll(bp, pI)
    {
        if( bp[pI] < 0 )
            continue;

        if( !mesh_.hasBackup(pI) )
            mesh_.addPointToSubset(pId, pI);
    }

    label lockedId = mesh_.pointSubsetIndex("lockedPoints");
    if( lockedId < 0 )
        lockedId = mesh_.addPointSubset("lockedPoints");

    forAll(lockedPoints, i)
        mesh_.addPointToSubset(lockedId, lockedPoints[i]);

    Info << "Alpha " << alpha << endl;
    # endif

    const pointFieldPMG& points = mesh_.points();
    const VRWGraph& pointFaces = mesh_.addressingData().pointFaces();
    forAllConstIter(std::set<label>, affectedPoints, pIt)
    {
        if( !mesh_.hasBackup(*pIt) || bp[*pIt] < 0 )
            continue;

        point pOrig;
        mesh_.getOrigPoint(*pIt, pOrig);

        polyMeshGenModifier(mesh_).movePoint
        (
            *pIt,
            pOrig + alpha * (points[*pIt] - pOrig)
        );

        # ifdef DEBUGSmoothing
        Info << "Point " << *pIt << " orig coordinates " << pOrig
             << " current coordinates " << mesh_.points()[*pIt] << endl;
        # endif

        if( lockedVertices.find(*pIt) == lockedVertices.end() )
        {
            lockedVertices.insert(*pIt);
            lockedPoints.append(*pIt);
        }

        forAllRow(pointFaces, *pIt, pfI)
            changedFace[pointFaces(*pIt, pfI)] = true;
    }

    mesh_.addressingData().updateGeometry(changedFace);
}

void meshOptimizer::updateActiveFaces
(
    const labelHashSet& badFaces,
    boolList& activeFaces,
    const direction nAdditionalLayers
) const
{
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();

    const VRWGraph& pointCells = mesh_.addressingData().pointCells();

    activeFaces.setSize(faces.size());

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

                activeFaces[faceI] = true;
            }
            else
            {
                activeFaces[faceI] = false;
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

            for(direction layerI=1;layerI<(nAdditionalLayers+1);++layerI)
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
                                    //- mark the face as active
                                    activeFaces[c[fI]] = true;

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

                                                const cell& cn = cells[cLabel];

                                                //- mark faces as active
                                                forAll(cn, fI)
                                                    activeFaces[cn[fI]] = true;
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
                    const labelLongList& globalPointLabel =
                        mesh_.addressingData().globalPointLabel();
                    const VRWGraph& pProcs =
                        mesh_.addressingData().pointAtProcs();
                    const Map<label>& globalToLocal =
                        mesh_.addressingData().globalToLocalPointAddressing();

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
                                    eData[neiProc].append
                                    (
                                        globalPointLabel[pointI]
                                    );
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

                                const cell& cn = cells[cLabel];

                                //- mark faces as active
                                forAll(cn, fI)
                                    activeFaces[cn[fI]] = true;
                            }
                        }
                    }
                }
            }
        }
    }
}

bool meshOptimizer::removeTangledCells(const labelHashSet& badFaces)
{
    const labelLongList& owner = mesh_.owner();
    const labelLongList& neighbour = mesh_.neighbour();

    boolList removeCell(mesh_.cells().size(), false);

    bool changed(false);

    forAllConstIter(labelHashSet, badFaces, it)
    {
        changed = true;

        const label faceI = it.key();

        removeCell[owner[faceI]] = true;

        if( neighbour[faceI] >= 0 )
            removeCell[neighbour[faceI]] = true;
    }

    reduce(changed, maxOp<bool>());

    if( !changed )
        return false;

    //- remove tangled cells from the mesh
    polyMeshGenModifier meshModifier(mesh_);
    meshModifier.removeCells(removeCell);

    //- check connections over faces
    if( checkCellConnectionsOverFaces(mesh_).checkCellGroups() )
        changed = true;

    meshModifier.removeUnusedVertices();

    clearSurface();
    calculatePointLocations();

    return changed;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
meshOptimizer::meshOptimizer(polyMeshGen& mesh)
:
    mesh_(mesh),
    vertexLocation_(),
    lockedFaces_(),
    msePtr_(NULL),
    enforceConstraints_(false),
    badPointsSubsetName_()
{
    calculatePointLocations();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshOptimizer::~meshOptimizer()
{
    clearSurface();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOptimizer::enforceConstraints(const word subsetName)
{
    enforceConstraints_ = true;

    badPointsSubsetName_ = subsetName;
}

void meshOptimizer::lockCellsInSubset(const word& subsetName)
{
    //- lock the points in the cell subset with the given name
    label subsetI = mesh_.cellSubsetIndex(subsetName);
    if( subsetI >= 0 )
    {
        labelLongList lc;
        mesh_.cellsInSubset(subsetI, lc);

        lockCells(lc);
    }
    else
    {
        Warning << "Subset " << subsetName << " is not a cell subset!"
            << " Cannot lock cells!" << endl;
    }
}

void meshOptimizer::lockFacesInSubset(const word& subsetName)
{
    //- lock the points in the face subset with the given name
    label subsetI = mesh_.faceSubsetIndex(subsetName);
    if( subsetI >= 0 )
    {
        labelLongList lf;
        mesh_.facesInSubset(subsetI, lf);

        lockFaces(lf);
    }
    else
    {
        Warning << "Subset " << subsetName << " is not a face subset!"
            << " Cannot lock faces!" << endl;
    }
}

void meshOptimizer::lockPointsInSubset(const word& subsetName)
{
    //- lock the points in the point subset with the given name
    label subsetI = mesh_.pointSubsetIndex(subsetName);
    if( subsetI >= 0 )
    {
        labelLongList lp;
        mesh_.pointsInSubset(subsetI, lp);

        lockPoints(lp);
    }
    else
    {
        Warning << "Subset " << subsetName << " is not a point subset!"
            << " Cannot lock points!" << endl;
    }
}

void meshOptimizer::removeUserConstraints()
{
    lockedFaces_.setSize(0);

    //- unlock points
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(vertexLocation_, i)
    {
        if( vertexLocation_[i] & LOCKED )
            vertexLocation_[i] ^= LOCKED;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
