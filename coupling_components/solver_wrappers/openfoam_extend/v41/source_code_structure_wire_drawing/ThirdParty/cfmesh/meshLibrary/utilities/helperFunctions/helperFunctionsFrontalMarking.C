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
#include "helperFunctionsFrontalMarking.H"
#include "DynList.H"
#include "labelPair.H"
#include "HashSet.H"

# ifdef USE_OMP
#include <omp.h>
# endif

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace help
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

template<class labelListType, class neiOp, class filterOp>
void frontalMarking
(
    labelListType& result,
    const label startingIndex,
    const neiOp& neighbourCalculator,
    const filterOp& selector
)
{
    //- add the starting element
    result.clear();
    result.append(startingIndex);

    //- add the starting element to the front
    DynList<label, 128> front;
    front.append(startingIndex);

    //- store information which element were already visited
    labelHashSet alreadySelected;

    //- start with frontal marking
    DynList<label> neighbours;

    while( front.size() )
    {
        const label eLabel = front.removeLastElement();

        //- find neighbours of the current element
        neighbours.clear();
        neighbourCalculator(eLabel, neighbours);

        forAll(neighbours, neiI)
        {
            const label nei = neighbours[neiI];

            if( nei < 0 )
                continue;
            if( alreadySelected.found(nei) )
                continue;

            if( selector(nei) )
            {
                alreadySelected.insert(nei);

                front.append(nei);
                result.append(nei);
            }
        }
    }
}

class graphNeiOp
{
    // Private data
        //- const reference to VRWGraph
        const VRWGraph& neiGroups_;

public:

    // Constructors
        //- Construct from VRWGraph
        inline graphNeiOp(const VRWGraph& neiGroups)
        :
            neiGroups_(neiGroups)
        {}

    // Public member functions
        //- return the size of the graph
        inline label size() const
        {
            return neiGroups_.size();
        }

        //- operator for finding neighbours
        inline void operator()(const label groupI, DynList<label>& ng) const
        {
            ng.setSize(neiGroups_.sizeOfRow(groupI));
            forAllRow(neiGroups_, groupI, i)
                ng[i] = neiGroups_(groupI, i);
        }
};

class graphSelectorOp
{
    // Private data
        //- const reference to VRWGraph
        const VRWGraph& neiGroups_;

public:

    // Constructors
        //- Construct from VRWGraph
        inline graphSelectorOp(const VRWGraph& neiGroups)
        :
            neiGroups_(neiGroups)
        {}

    // Public member functions
        //- operator for selecting elements
        inline bool operator()(const label groupI) const
        {
            if( (groupI < 0) || (groupI >= neiGroups_.size()) )
                return false;

            return true;
        }
};

template<class labelListType, class neiOp, class filterOp>
label groupMarking
(
    labelListType& elementInGroup,
    const neiOp& neighbourCalculator,
    const filterOp& selector
)
{
    label nGroups(0);

    elementInGroup.setSize(neighbourCalculator.size());

    DynList<label> nGroupsAtTask;
    List<LongList<std::pair<label, label> > > taskComms;

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        //- initialize data
        # pragma omp for schedule(static)
        forAll(elementInGroup, i)
            elementInGroup[i] = -1;

        const label nTasks = 5 * omp_get_num_threads();

        # pragma omp single
        {
            taskComms.setSize(nTasks);
            nGroupsAtTask.setSize(nTasks);
        }
        # else
        //- initialize data
        elementInGroup = -1;

        const label nTasks = 1;

        taskComms.setSize(nTasks);
        nGroupsAtTask.setSize(1);
        # endif

        const label chunkSize = Foam::max(1, int(elementInGroup.size()/nTasks));

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
                    labelLongList front;
                    DynList<label> neighbours;

                    LongList<std::pair<label, label> >& taskCommPairs =
                        taskComms[taskI];
                    taskCommPairs.clear();

                    const label minEl = taskI * chunkSize;

                    label maxEl =
                        min(minEl + chunkSize, int(elementInGroup.size()));
                    if( taskI == (nTasks - 1) )
                        maxEl = elementInGroup.size();

                    label& groupI = nGroupsAtTask[taskI];
                    groupI = 0;

                    for(label elI=minEl;elI<maxEl;++elI)
                    {
                        if( elementInGroup[elI] != -1 )
                            continue;
                        if( !selector(elI) )
                            continue;

                        elementInGroup[elI] = groupI;

                        front.clear();
                        front.append(elI);

                        while( front.size() )
                        {
                            const label eLabel = front.removeLastElement();

                            neighbours.clear();
                            neighbourCalculator(eLabel, neighbours);

                            forAll(neighbours, neiI)
                            {
                                const label nei = neighbours[neiI];

                                if( nei < 0 )
                                    continue;

                                if( (nei < minEl) || (nei >= maxEl) )
                                {
                                    //- this is a communication interface
                                    //- between two tasks
                                    taskCommPairs.append
                                    (
                                        std::make_pair(elI, nei)
                                    );
                                }
                                else if
                                (
                                    selector(nei) && (elementInGroup[nei] == -1)
                                )
                                {
                                    //- this element is part of the same group
                                    elementInGroup[nei] = groupI;
                                    front.append(nei);
                                }
                            }
                        }

                        ++groupI;
                    }
                }
            }
        }

        //- each task shall fill the correct data
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
                    label startGroup(0);
                    for(label i=0;i<taskI;++i)
                        startGroup += nGroupsAtTask[i];

                    const label minEl = taskI * chunkSize;
                    label maxEl =
                        min(minEl + chunkSize, int(elementInGroup.size()));
                    if( taskI == (nTasks - 1) )
                        maxEl = elementInGroup.size();

                    for(label elI=minEl;elI<maxEl;++elI)
                    {
                        if( elementInGroup[elI] < 0 )
                        {
                            continue;
                        }

                        elementInGroup[elI] += startGroup;
                    }
                }
            }
        }
    }

    //- calculate the total number of groups
    nGroups = 0;
    forAll(nGroupsAtTask, taskI)
        nGroups += nGroupsAtTask[taskI];

    //- find connections between groups
    List<labelHashSet> neiGroupsHelper;
    neiGroupsHelper.setSize(nGroups);

    forAll(taskComms, taskI)
    {
        const LongList<std::pair<label, label> >& taskCommPairs =
            taskComms[taskI];

        forAll(taskCommPairs, i)
        {
            const std::pair<label, label>& lp = taskCommPairs[i];
            const label groupI = elementInGroup[lp.first];
            const label neiGroup = elementInGroup[lp.second];

            if( (groupI < 0) || (neiGroup < 0) )
            {
                continue;
            }

            neiGroupsHelper[groupI].insert(neiGroup);
            neiGroupsHelper[neiGroup].insert(groupI);
        }
    }

    //- copy data into a VRWGraph for further processing
    VRWGraph neighbouringGroups;
    neighbouringGroups.setSize(nGroups);

    forAll(neiGroupsHelper, i)
    {
        labelList helper = neiGroupsHelper[i].toc();

        sort(helper);

        neighbouringGroups[i] = helper;
    }
    neiGroupsHelper.setSize(0);

    //- start processing connections between the group and merge the connected
    //- ones into a new group
    DynList<label> globalGroupLabel;
    globalGroupLabel.setSize(nGroups);
    globalGroupLabel = -1;

    //- reduce the information about the groups
    DynList<label> connectedGroups;
    label counter(0);

    forAll(neighbouringGroups, groupI)
    {
        if( globalGroupLabel[groupI] != -1 )
            continue;

        connectedGroups.clear();
        frontalMarking
        (
            connectedGroups,
            groupI,
            graphNeiOp(neighbouringGroups),
            graphSelectorOp(neighbouringGroups)
        );

        forAll(connectedGroups, gI)
            globalGroupLabel[connectedGroups[gI]] = counter;

        ++counter;
    }

    nGroups = counter;

    forAll(neighbouringGroups, groupI)
    {
        if( globalGroupLabel[groupI] != -1 )
            continue;

        forAllRow(neighbouringGroups, groupI, ngI)
            globalGroupLabel[neighbouringGroups(groupI, ngI)] = counter;

        ++counter;
    }

    if( Pstream::parRun() )
    {
        //- reduce the groups over processors of an MPI run
        //- count the total number of groups over all processors
        labelList nGroupsAtProc(Pstream::nProcs());
        nGroupsAtProc[Pstream::myProcNo()] = nGroups;

        Pstream::gatherList(nGroupsAtProc);
        Pstream::scatterList(nGroupsAtProc);

        label startGroup(0), totalNumGroups(0);
        for(label procI=0;procI<Pstream::nProcs();++procI)
        {
            totalNumGroups += nGroupsAtProc[procI];

            if( procI < Pstream::myProcNo() )
                startGroup += nGroupsAtProc[procI];
        }

        //- translate group labels
        forAll(globalGroupLabel, groupI)
        {
            if( globalGroupLabel[groupI] < 0 )
                continue;

            globalGroupLabel[groupI] += startGroup;
        }

        //- find the neighbouring groups
        //- collect groups on other processors
        //- this operator implements the algorithm for exchanging data
        //- over processors and collects information which groups
        //- are connected over inter-processor boundaries
        std::map<label, DynList<label> > neiGroups;

        neighbourCalculator.collectGroups
        (
            neiGroups,
            elementInGroup,
            globalGroupLabel
        );

        //- create a graph of connections
        List<List<labelPair> > globalNeiGroups(Pstream::nProcs());

        DynList<labelPair> connsAtProc;
        for
        (
            std::map<label, DynList<label> >::const_iterator it =
            neiGroups.begin();
            it!=neiGroups.end();
            ++it
        )
        {
            const DynList<label>& ng = it->second;

            forAll(ng, i)
                connsAtProc.append(labelPair(it->first, ng[i]));
        }

        //- copy the connections into the global neighbour list
        globalNeiGroups[Pstream::myProcNo()].setSize(connsAtProc.size());

        forAll(connsAtProc, i)
            globalNeiGroups[Pstream::myProcNo()][i] = connsAtProc[i];

        //- communicate partial graphs to the master processor
        Pstream::gatherList(globalNeiGroups);

        labelList allGroupsNewLabel;
        if( Pstream::master() )
        {
            //- collect the graph of connections for the whole system
            VRWGraph allGroups(totalNumGroups);
            forAll(allGroups, i)
                allGroups[i].append(i);

            forAll(globalNeiGroups, procI)
            {
                const List<labelPair>& connections = globalNeiGroups[procI];

                forAll(connections, i)
                {
                    const labelPair& lp = connections[i];

                    allGroups.appendIfNotIn(lp.first(), lp.second());
                    allGroups.appendIfNotIn(lp.second(), lp.first());
                }
            }

            //- assign a global label to each group
            allGroupsNewLabel.setSize(totalNumGroups);
            allGroupsNewLabel = -1;
            counter = 0;

            forAll(allGroups, groupI)
            {
                if( allGroupsNewLabel[groupI] != -1 )
                    continue;

                DynList<label> connectedGroups;
                frontalMarking
                (
                    connectedGroups,
                    groupI,
                    graphNeiOp(allGroups),
                    graphSelectorOp(allGroups)
                );

                forAll(connectedGroups, gI)
                    allGroupsNewLabel[connectedGroups[gI]] = counter;

                ++counter;
            }

            nGroups = counter;
        }

        //- broadcast group labels from the master to other processors
        Pstream::scatter(nGroups);
        Pstream::scatter(allGroupsNewLabel);

        //- assign correct group labels
        forAll(globalGroupLabel, groupI)
        {
            const label ngI = globalGroupLabel[groupI];
            globalGroupLabel[groupI] = allGroupsNewLabel[ngI];
        }
    }

    //- set the global group label
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(elementInGroup, elI)
    {
        if( elementInGroup[elI] < 0 )
            continue;

        elementInGroup[elI] = globalGroupLabel[elementInGroup[elI]];
    }

    return nGroups;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
