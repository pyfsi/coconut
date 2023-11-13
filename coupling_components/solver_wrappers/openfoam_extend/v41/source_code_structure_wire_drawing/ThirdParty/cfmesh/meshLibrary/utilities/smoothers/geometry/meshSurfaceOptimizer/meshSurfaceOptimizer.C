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
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfacePartitioner.H"
#include "meshSurfaceMapper.H"
#include "polyMeshGenChecks.H"
#include "labelledScalar.H"
#include "labelledPointScalar.H"
#include "meshOctree.H"

#include <map>

//# define DEBUGSmoothing

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceOptimizer::classifySurfaceVertices()
{
    const labelHashSet& corners = partitionerPtr_->corners();
    const labelHashSet& edgePoints = partitionerPtr_->edgePoints();

    //- set all vertices to partition
    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(vertexType_, bpI)
        {
            vertexType_[bpI] = PARTITION;
        }

        # ifdef USE_OMP
        # pragma omp sections
        # endif
        {
            # ifdef USE_OMP
            # pragma omp section
            # endif
            {
                //- set corners
                forAllConstIter(labelHashSet, corners, it)
                    vertexType_[it.key()] = CORNER;
            }

            # ifdef USE_OMP
            # pragma omp section
            # endif
            {
                //- set edges
                forAllConstIter(labelHashSet, edgePoints, it)
                    vertexType_[it.key()] = EDGE;
            }
        }
    }

    if( surfaceEngine_.mesh().hasLockedPoints() )
    {
        const labelLongList& bp = surfaceEngine_.bp();
        const std::set<label>& lockedPts = surfaceEngine_.mesh().lockedPoints();

        forAllConstIter(std::set<label>, lockedPts, it)
            vertexType_[bp[*it]] |= LOCKED;
    }

    if( Pstream::parRun() )
    {
        //- mark nodes at parallel boundaries
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();

        forAllConstIter(Map<label>, globalToLocal, iter)
        {
            const label bpI = iter();

            vertexType_[bpI] |= PROCBND;
        }
    }
}

void meshSurfaceOptimizer::projectPointsBackOnTheSurface
(
    const labelLongList& nodesToSmooth,
    pointField& newCoordinates
) const
{
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const pointFieldPMG& points = surfaceEngine_.points();

    # ifdef DEBUGSmoothing
    const scalar start = omp_get_wtime();
    # endif

    const VRWGraph* bpAtProcsPtr(NULL);

    if( Pstream::parRun() )
    {
        bpAtProcsPtr = &surfaceEngine_.bpAtProcs();
    }

    std::map<label, std::pair<vector, scalar> > procPoints;

    # ifdef USE_OMP
    # pragma omp parallel if( !omp_in_parallel() )
    # endif
    {
        # ifdef USE_OMP
        const label nTasks = 5 * omp_get_num_threads();
        # else
        const label nTasks = 1;
        # endif

        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            for(label taskI=0;taskI<nTasks;++taskI)
            {
                # ifdef USE_OMP
                # pragma omp task default(shared) firstprivate(taskI)
                # endif
                {
                    for(label nI=taskI;nI<nodesToSmooth.size();nI+=nTasks)
                    {
                        const label bpI = nodesToSmooth[nI];
                        const point& p = newCoordinates[nI];

                        point nearest(p);
                        scalar distSq(VGREAT);

                        forAllRow(pFaces, bpI, pfI)
                        {
                            const face& bf = bFaces[pFaces(bpI, pfI)];
                            const point np =
                                help::nearestPointOnFace(bf, points, p);

                            const scalar dSq = magSqr(np - p);
                            if( dSq < distSq )
                            {
                                distSq = dSq;
                                nearest = np;
                            }
                        }

                        if( bpAtProcsPtr && bpAtProcsPtr->sizeOfRow(bpI) )
                        {
                            # ifdef USE_OMP
                            # pragma omp critical(collectProcPoints)
                            # endif
                            {
                                //- point is shared by two or more processors
                                procPoints[bpI] =
                                    std::make_pair(nearest, distSq);
                            }
                        }
                        else
                        {
                            newCoordinates[nI] = nearest;
                        }
                    }
                }
            }
        }
    }

    if( Pstream::parRun() )
    {
        const labelLongList& globalPointLabel =
            surfaceEngine_.globalBoundaryPointLabel();
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();

        //- allocate the map used to store messages to other processors
        std::map<label, LongList<labelledPointScalar> > exchangeData;
        forAll(surfaceEngine_.bpNeiProcs(), i)
            exchangeData[surfaceEngine_.bpNeiProcs()[i]].clear();

        //- generate messages
        std::map<label, std::pair<vector, scalar> >::const_iterator it;
        for(it=procPoints.begin();it!=procPoints.end();++it)
        {
            const label bpI = it->first;
            const label gpI = globalPointLabel[bpI];

            forAllRow(*bpAtProcsPtr, bpI, j)
            {
                const label neiProc = bpAtProcsPtr->operator()(bpI, j);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append
                (
                    labelledPointScalar
                    (
                        gpI,
                        it->second.first,
                        it->second.second
                    )
                );
            }
        }

        LongList<labelledPointScalar> receiveData;
        help::exchangeMap(exchangeData, receiveData);

        //- find the nearest point at the surface of the mesh
        forAll(receiveData, i)
        {
            const labelledPointScalar& lps = receiveData[i];

            const label bpI = globalToLocal[lps.pointLabel()];

            //- take the point with the smallest distance
            if( procPoints[bpI].second > lps.scalarValue() )
            {
                procPoints[bpI] =
                    std::make_pair(lps.coordinates(), lps.scalarValue());
            }
        }

        //- update the coordinates of the points
        forAll(newCoordinates, nI)
        {
            if( procPoints.find(nodesToSmooth[nI]) != procPoints.end() )
            {
                newCoordinates[nI] = procPoints[nodesToSmooth[nI]].first;
            }
        }
    }

    # ifdef DEBUGSmoothing
    Info << "Time for reprojecting " << (omp_get_wtime() - start) << endl;
    # endif
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceOptimizer::meshSurfaceOptimizer(const meshSurfaceEngine& surface)
:
    surfaceEngine_(surface),
    surfaceModifier_(surface),
    vertexType_(surface.boundaryPoints().size()),
    partitionerPtr_(new meshSurfacePartitioner(surface)),
    deletePartitioner_(true),
    octreePtr_(NULL),
    triMeshPtr_(NULL),
    enforceConstraints_(false),
    badPointsSubsetName_("invertedBoundaryPoints")
{
    classifySurfaceVertices();
}

meshSurfaceOptimizer::meshSurfaceOptimizer(const meshSurfacePartitioner& mPart)
:
    surfaceEngine_(mPart.surfaceEngine()),
    surfaceModifier_(mPart.surfaceEngine()),
    vertexType_(surfaceEngine_.boundaryPoints().size()),
    partitionerPtr_(&mPart),
    deletePartitioner_(false),
    octreePtr_(NULL),
    triMeshPtr_(NULL),
    enforceConstraints_(false),
    badPointsSubsetName_("invertedBoundaryPoints")
{
    classifySurfaceVertices();
}

meshSurfaceOptimizer::meshSurfaceOptimizer
(
    const meshSurfaceEngine& surface,
    const meshOctree& octree
)
:
    surfaceEngine_(surface),
    surfaceModifier_(surface),
    vertexType_(surface.boundaryPoints().size()),
    partitionerPtr_(new meshSurfacePartitioner(surface)),
    deletePartitioner_(true),
    octreePtr_(&octree),
    triMeshPtr_(NULL),
    enforceConstraints_(false),
    badPointsSubsetName_("invertedBoundaryPoints")
{
    classifySurfaceVertices();
}

meshSurfaceOptimizer::meshSurfaceOptimizer
(
    const meshSurfacePartitioner& partitioner,
    const meshOctree& octree
)
:
    surfaceEngine_(partitioner.surfaceEngine()),
    surfaceModifier_(partitioner.surfaceEngine()),
    vertexType_(surfaceEngine_.boundaryPoints().size()),
    partitionerPtr_(&partitioner),
    deletePartitioner_(false),
    octreePtr_(&octree),
    triMeshPtr_(NULL),
    enforceConstraints_(false),
    badPointsSubsetName_("invertedBoundaryPoints")
{
    classifySurfaceVertices();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceOptimizer::~meshSurfaceOptimizer()
{
    deleteDemandDrivenData(triMeshPtr_);

    if( deletePartitioner_ )
        deleteDemandDrivenData(partitionerPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceOptimizer::removeUserConstraints()
{
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(vertexType_, bpI)
        if( vertexType_[bpI] & LOCKED )
            vertexType_[bpI] ^= LOCKED;
}

void meshSurfaceOptimizer::enforceConstraints(const word subsetName)
{
    enforceConstraints_ = true;

    badPointsSubsetName_ = subsetName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
