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

#include "meshOctreeOps.H"
#include "labelLongList.H"
#include "labelledPair.H"
#include "helperFunctions.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGCheck

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace meshOctreeOps
{

template<class labelListType>
void meshOctreeNeighbourOperator::collectGroups
(
    std::map<label, DynList<label> >& neiGroups,
    const labelListType& elementInGroup,
    const DynList<label>& localGroupLabel
) const
{
    const labelList& neiProcs = octree_.neiProcs();
    const List<Pair<meshOctreeCubeCoordinates> >& neiRange = octree_.neiRange();

    std::map<label, labelLongList> exchangeData;
    forAll(neiProcs, i)
        exchangeData.insert(std::make_pair(neiProcs[i], labelLongList()));

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(interProcessorLeaves_, leafI)
    {
        if( !interProcessorLeaves_[leafI] )
        {
            const label groupI = elementInGroup[leafI];

            if( groupI < 0 )
                continue;

            for(label dir=0;dir<6;++dir)
            {
                DynList<label> neighbours;

                octree_.findNeighboursInDirection(leafI, dir, neighbours);

                forAll(neighbours, i)
                {
                    if( neighbours[i] == meshOctreeCube::OTHERPROC )
                    {
                        meshOctreeCubeCoordinates minCoord, maxCoord;

                        const meshOctreeCubeBasic& oc =
                            octree_.returnLeaf(leafI);

                        oc.neighbourRange(minCoord, maxCoord);

                        forAll(neiProcs, npI)
                        {
                            if( maxCoord >= neiRange[npI].first() )
                            {
                                if( minCoord <= neiRange[npI].second() )
                                {
                                    labelLongList& data = exchangeData[npI];

                                    data.append(oc.posX());
                                    data.append(oc.posY());
                                    data.append(oc.posZ());
                                    data.append(oc.level());

                                    data.append(localGroupLabel[groupI]);

                                    # ifdef OCTREE_DEBUG
                                    ++counter;
                                    # endif
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    labelLongList receiveData;
    help::exchangeMap(exchangeData, receiveData);

    label counter(0);

    while( counter < receiveData.size() )
    {
        const label posX = receiveData[counter++];
        const label posY = receiveData[counter++];
        const label posZ = receiveData[counter++];
        const label level = receiveData[counter++];

        meshOctreeCubeCoordinates cc(posX, posY, posZ, level);

        const label neiLocalGroupI = receiveData[counter++];

        for(label dir=0;dir<6;++dir)
        {
            DynList<label> neighbours;

            octree_.findNeighboursInDirection(cc, dir, neighbours);

            forAll(neighbours, i)
            {
                const label nei = neighbours[i];

                if( nei < 0 )
                    continue;

                const label groupI = elementInGroup[nei];

                if( groupI < 0 )
                    continue;

                DynList<label>& ng = neiGroups[localGroupLabel[groupI]];

                //- store the connection over the inter-processor boundary
                ng.appendIfNotIn(neiLocalGroupI);
            }
        }
    }
}

template<class T>
class insertOp
{
public:
    insertOp()
    {}

    T operator()(T& v1, const T& v2) const
    {
        v1 += v2;

        return v1;
    }
};

template<class ListType>
void findOutsideGroups
(
    const meshOctree& octree,
    const ListType& leafInGroup,
    labelHashSet& markedGroups
)
{
    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp single nowait
        # endif
        {
            # ifdef USE_OMP
            const label nTasks = 10 * omp_get_num_threads();
            # else
            const label nTasks = 1;
            # endif

            const label cSize = max(octree.numberOfLeaves() / nTasks, 1);

            for(label taskI=0;taskI<nTasks;++taskI)
            {
                # ifdef USE_OMP
                # pragma omp task shared(octree, markedGroups, leafInGroup) \
                firstprivate(taskI)
                # endif
                {
                    const label sl = taskI * cSize;
                    const label el =
                        taskI!=(nTasks-1)?sl+cSize:octree.numberOfLeaves();

                    //- find leaves that have no neighbours
                    labelHashSet localOutsideGroups;
                    for(label leafI=sl;leafI<el;++leafI)
                    {
                        for(label fI=0;fI<6;++fI)
                        {
                            DynList<label> neighbours;
                            octree.findNeighboursInDirection
                            (
                                leafI,
                                fI,
                                neighbours
                            );

                            forAll(neighbours, i)
                            {
                                if( neighbours[i] == -1 )
                                {
                                    localOutsideGroups.insert
                                    (
                                        leafInGroup[leafI]
                                    );
                                }
                            }
                        }
                    }

                    if( localOutsideGroups.size() != 0 )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(outsideGroups)
                        # endif
                        {
                            markedGroups += localOutsideGroups;
                        }
                    }
                }
            }
        }
    }

    //- transfer groups across all processors
    reduce(markedGroups, insertOp<labelHashSet>());
}

template<class ListType>
void findInsideGroups
(
    const meshOctree& octree,
    const ListType& leafInGroup,
    labelHashSet& markedGroups,
    const labelHashSet* outsideGroupsPtr
)
{
    markedGroups.clear();

    //- find the outside groups first or use the provided information
    labelHashSet outsideGroups;
    if( !outsideGroupsPtr )
    {
        findOutsideGroups(octree, leafInGroup, outsideGroups);
    }
    else
    {
        outsideGroups = *outsideGroupsPtr;
    }

    //- find octree leaves attached to the leaves marked as OUTSIDE
    LongList<meshOctreeCubeCoordinates> transferCoordinates;

    List<labelHashSet> hasOutsideAtTask;
    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp single nowait
        # endif
        {
            # ifdef USE_OMP
            const label nTasks = 10 * omp_get_num_threads();
            # else
            const label nTasks = 1;
            # endif

            hasOutsideAtTask.setSize(nTasks);
            const label cSize = max(octree.numberOfLeaves() / nTasks, 1);

            for(label taskI=0;taskI<nTasks;++taskI)
            {
                # ifdef USE_OMP
                # pragma omp task shared(octree, transferCoordinates) \
                shared(leafInGroup, hasOutsideAtTask, outsideGroups) \
                firstprivate(taskI)
                # endif
                {
                    const label sl = taskI * cSize;
                    const label el =
                        taskI!=(nTasks-1)?sl+cSize:octree.numberOfLeaves();

                    for(label leafI=sl;leafI<el;++leafI)
                    {
                        if
                        (
                            outsideGroups.find(leafInGroup[leafI]) ==
                            outsideGroups.end()
                        )
                            continue;

                        //- check if there exist a neighbour with a negative
                        //- group index. Set hasOutside to true then.
                        for(label fI=0;fI<6;++fI)
                        {
                            DynList<label> neighbours;

                            octree.findNeighboursInDirection
                            (
                                leafI,
                                fI,
                                neighbours
                            );

                            forAll(neighbours, i)
                            {
                                const label nei = neighbours[i];

                                if( nei == -1 )
                                    continue;
                                if( nei == meshOctreeCube::OTHERPROC )
                                {
                                    # ifdef USE_OMP
                                    # pragma omp critical(transferCoordinates)
                                    # endif
                                    {
                                        transferCoordinates.append
                                        (
                                            octree.returnLeaf(leafI)
                                        );
                                    }
                                    continue;
                                }

                                if( leafInGroup[nei] == -1 )
                                    hasOutsideAtTask[taskI].insert(nei);
                            }
                        }
                    }
                }
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- check the neighbours of the transferred leaves
        LongList<meshOctreeCubeCoordinates> receivedLeaves;
        octree.exchangeRequestsWithNeighbourProcessors
        (
            transferCoordinates,
            receivedLeaves
        );

        # ifdef USE_OMP
        # pragma omp parallel
        # endif
        {
            # ifdef USE_OMP
            # pragma omp single nowait
            # endif
            {
                const label nTasks = hasOutsideAtTask.size();
                const label cSize = max(label(receivedLeaves.size()/nTasks), 1);
                const label nReceived = receivedLeaves.size();

                for(label taskI=0;taskI<nTasks;++taskI)
                {
                    # ifdef USE_OMP
                    # pragma omp task shared(hasOutsideAtTask, receivedLeaves) \
                    shared(octree, leafInGroup) firstprivate(taskI)
                    # endif
                    {
                        const label sc = taskI * cSize;
                        const label ec = taskI!=(nTasks-1)?sc+cSize:nReceived;

                        for(label i=sc;i<ec;++i)
                        {
                            const meshOctreeCubeCoordinates& cc =
                                receivedLeaves[i];

                            for(label fI=0;fI<6;++fI)
                            {
                                DynList<label> neighbours;

                                octree.findNeighboursInDirection
                                (
                                    cc,
                                    fI,
                                    neighbours
                                );

                                forAll(neighbours, i)
                                {
                                    const label nei = neighbours[i];

                                    if( nei == -1 )
                                        continue;
                                    if( nei == meshOctreeCube::OTHERPROC )
                                        continue;

                                    if( leafInGroup[nei] == -1 )
                                        hasOutsideAtTask[taskI].insert(nei);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //- do the another loop through the leaves marked at each task
    //- an find out their neighbours that are not in any of the outside groups
    transferCoordinates.clear();
    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp single nowait
        # endif
        {
            forAll(hasOutsideAtTask, taskI)
            {
                # ifdef USE_OMP
                # pragma omp task shared(hasOutsideAtTask, leafInGroup) \
                shared(markedGroups, transferCoordinates) \
                shared(octree, outsideGroups) firstprivate(taskI)
                # endif
                {
                    labelHashSet localInsideGroups;

                    const labelHashSet& currLeaves = hasOutsideAtTask[taskI];
                    forAllConstIter(labelHashSet, currLeaves, it)
                    {
                        for(label fI=0;fI<6;++fI)
                        {
                            DynList<label> neighbours;

                            octree.findNeighboursInDirection
                            (
                                it.key(),
                                fI,
                                neighbours
                            );

                            forAll(neighbours, i)
                            {
                                const label nei = neighbours[i];

                                if( nei == -1 )
                                    continue;
                                if( nei == meshOctreeCube::OTHERPROC )
                                {
                                    # ifdef USE_OMP
                                    # pragma omp critical(transferCoordinates1)
                                    # endif
                                    {
                                        transferCoordinates.append
                                        (
                                            octree.returnLeaf(it.key())
                                        );
                                    }
                                    continue;
                                }

                                const label neiGroup = leafInGroup[nei];
                                if( neiGroup < 0 )
                                    continue;
                                if
                                (
                                    outsideGroups.find(neiGroup) !=
                                    outsideGroups.end()
                                )
                                    continue;

                                if( neiGroup != -1 )
                                    localInsideGroups.insert(neiGroup);
                            }
                        }
                    }

                    //- merge the local into th global set of inside groups
                    # ifdef USE_OMP
                    # pragma omp critical(insideGroups)
                    # endif
                    {
                        markedGroups += localInsideGroups;
                    }
                }
            }
        }
    }

    if( Pstream::parRun() )
    {
        LongList<meshOctreeCubeCoordinates> receivedCoordinates;
        octree.exchangeRequestsWithNeighbourProcessors
        (
            transferCoordinates,
            receivedCoordinates
        );

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 10)
        # endif
        forAll(receivedCoordinates, i)
        {
            const meshOctreeCubeCoordinates& cc = receivedCoordinates[i];

            for(label fI=0;fI<6;++fI)
            {
                DynList<label> neighbours;

                octree.findNeighboursInDirection
                (
                    cc,
                    fI,
                    neighbours
                );

                forAll(neighbours, i)
                {
                    const label nei = neighbours[i];

                    if( nei == -1 )
                        continue;
                    if( nei == meshOctreeCube::OTHERPROC )
                        continue;

                    const label neiGroup = leafInGroup[nei];
                    if( neiGroup < 0 )
                        continue;
                    if
                    (
                        outsideGroups.find(neiGroup) !=
                        outsideGroups.end()
                    )
                        continue;

                    if( neiGroup != -1 )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(insertGroups)
                        # endif
                        markedGroups.insert(neiGroup);
                    }
                }
            }
        }
    }

    //- transfer groups across all processors
    reduce(markedGroups, insertOp<labelHashSet>());
}

template<class ListType>
void findEnclosedGroups
(
    const meshOctree& octree,
    const ListType& leafInGroup,
    labelHashSet& markedGroups,
    const labelHashSet* insideAndOutsideGroupsPtr
)
{
    markedGroups.clear();

    //- find the groups that are inside and outside of the given leaves
    labelHashSet insideAndOutsideGroups;
    if( insideAndOutsideGroupsPtr )
    {
        insideAndOutsideGroups = *insideAndOutsideGroupsPtr;
    }
    else
    {
        //- find outside and inside groups
        labelHashSet outGroups;
        findOutsideGroups(octree, leafInGroup, outGroups);

        labelHashSet inGroups;
        findInsideGroups(octree, leafInGroup, inGroups, &outGroups);

        insideAndOutsideGroups += outGroups;
        insideAndOutsideGroups += inGroups;
    }

    //- pass through the leaves and select the ones that are not in the
    //- insideAndOutsideGroups
    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        labelHashSet localGroups;

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(leafInGroup, leafI)
        {
            const label groupI = leafInGroup[leafI];

            if( groupI < 0 )
                continue;

            if
            (
                insideAndOutsideGroups.find(groupI) ==
                insideAndOutsideGroups.end()
            )
                localGroups.insert(groupI);
        }

        # ifdef USE_OMP
        # pragma omp critical(localEnclosedGroups)
        # endif
        {
            markedGroups += localGroups;
        }
    }

    //- transfer groups across all processors
    reduce(markedGroups, insertOp<labelHashSet >());
}

template<class ListType>
void findGroupsByType
(
    const meshOctree& octree,
    const ListType& leafInGroup,
    const meshOctreeCubeBasic::typesOfCubes currType,
    labelHashSet& markedGroups
)
{
    labelHashSet outside, inside, enclosed;

    if
    (
        currType &
        (
            meshOctreeCubeBasic::OUTSIDE |
            meshOctreeCubeBasic::INSIDE |
            meshOctreeCubeBasic::UNKNOWN
        )
    )
    {
        findOutsideGroups(octree, leafInGroup, outside);
    }

    if
    (
        currType &
        (
            meshOctreeCubeBasic::INSIDE |
            meshOctreeCubeBasic::UNKNOWN
        )
    )
    {
        findInsideGroups(octree, leafInGroup, inside, &outside);
    }

    if( currType & meshOctreeCubeBasic::UNKNOWN )
    {
        labelHashSet inAndOut;
        inAndOut += outside;
        inAndOut += inside;
        findEnclosedGroups(octree, leafInGroup, enclosed, &inAndOut);
    }

    //- copy the selected groups
    markedGroups.clear();
    if( currType & meshOctreeCubeBasic::OUTSIDE )
    {
        markedGroups += outside;
    }
    if( currType & meshOctreeCubeBasic::INSIDE )
    {
        markedGroups += inside;
    }
    if( currType & meshOctreeCubeBasic::UNKNOWN )
    {
        markedGroups += enclosed;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meshConnectionsHelper

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
