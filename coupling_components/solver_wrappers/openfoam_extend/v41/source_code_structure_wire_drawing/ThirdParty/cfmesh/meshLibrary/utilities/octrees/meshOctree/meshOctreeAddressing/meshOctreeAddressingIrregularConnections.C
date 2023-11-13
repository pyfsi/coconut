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

#include "meshOctreeAddressing.H"
#include "meshOctree.H"
#include "demandDrivenData.H"
#include "helperFunctions.H"

#include <map>

//#define DEBUGAutoRef

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOctreeAddressing::checkAndFixIrregularConnections()
{
    Info << "Checking the surface of selected boxes" << endl;

    //- create sets with invalid codes
    HashSet<direction> problematicVertexCodes(8);
    for(label i=0;i<4;++i)
    {
        direction available(0);
        available |= (1 << i);
        available |= (1 << (7-i));

        problematicVertexCodes.insert(available);

        direction missing(255);

        missing ^= (1 << i);
        missing ^= (1 <<(7-i));

        problematicVertexCodes.insert(missing);
    }

    FixedList<FixedList<FixedList<label, 4>, 3>, 8> edgePositionsAtVertex;
    forAll(edgePositionsAtVertex, i)
    {
        label pos = 0;

        for(label fI=0;fI<6;++fI)
        {
            //- the face must not contain the node
            bool containsNode(false);
            for(label pI=0;pI<4;++pI)
                if( meshOctreeCube::faceNodes_[fI][pI] == i )
                {
                    containsNode = true;
                    break;
                }

            if( !containsNode )
            {
                for(label pI=0;pI<4;++pI)
                    edgePositionsAtVertex[i][pos][pI] =
                        meshOctreeCube::faceNodes_[fI][pI];

                ++pos;
            }
        }
    }

    FixedList<HashSet<direction>, 8> problematicEdgeCodes;
    forAll(edgePositionsAtVertex, vI)
    {
        forAll(edgePositionsAtVertex[vI], dir)
        {
            const FixedList<label, 4>& activePos =
                edgePositionsAtVertex[vI][dir];

            for(label i=0;i<2;++i)
            {
                direction vertexCode(0);

                vertexCode |= (1 << activePos[i]);
                vertexCode |= (1 << activePos[i+2]);

                problematicEdgeCodes[vI].insert(vertexCode);
            }
        }
    }

    //- get information about octree leaves used as mesh cells
    if( !boxTypePtr_ )
        findUsedBoxes();

    List<direction>& boxType = *boxTypePtr_;

    //- find problematic connections and skip template cells with many boundary
    //- pairs
    bool foundProblematic;
    do
    {
        foundProblematic = false;

        //- find boundary faces of leaves marked as MESHCELL
        List<FixedList<direction, 6> > bndFaces
        (
            boxType.size(),
            FixedList<direction, 6>(direction(0))
        );

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(boxType, leafI)
        {
            if( boxType[leafI] & MESHCELL )
            {
                if( Pstream::parRun() )
                {
                    const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);

                    if( oc.procNo() != Pstream::myProcNo() )
                        continue;
                }

                for(label dir=0;dir<6;++dir)
                {
                    DynList<label> neighbours;
                    octree_.findNeighboursInDirection(leafI, dir, neighbours);

                    if( neighbours.size() == 1 )
                    {
                        //- neighbour is at the same level
                        const label nei = neighbours[0];

                        if( nei < 0 )
                        {
                            continue;
                        }
                        else if( boxType[nei] & BOUNDARY )
                        {
                            bndFaces[leafI][dir] = 1;
                        }
                    }
                    else
                    {
                        //- neighbour is at higher level
                        bool allBnd(true), hasBnd(false);
                        forAll(neighbours, i)
                        {
                            const label nei = neighbours[i];

                            if( nei < 0 )
                            {
                                hasBnd = true;
                            }
                            else if( boxType[nei] & BOUNDARY )
                            {
                                hasBnd = true;
                            }
                            else
                            {
                                allBnd = false;
                            }
                        }

                        if( allBnd )
                        {
                            bndFaces[leafI][dir] = 1;
                        }
                        else if( hasBnd )
                        {
                            bndFaces[leafI][dir] = 2;
                        }
                    }
                }
            }
        }

        //- find problematic connections over vertices and edges
        LongList<labelPair> problematicConnections;

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(boxType, leafI)
        {
            if( boxType[leafI] & MESHCELL )
            {
                if( Pstream::parRun()   )
                {
                    const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);

                    if( oc.procNo() != Pstream::myProcNo() )
                        continue;
                }

                for(label vI=0;vI<8;++vI)
                {
                    FixedList<label, 8> neighbours;
                    octree_.findLeavesForCubeVertex(leafI, vI, neighbours);

                    direction vertexCode(0);

                    forAll(neighbours, i)
                    {
                        const label boxI = neighbours[i];

                        if( boxType[boxI] & MESHCELL )
                            vertexCode |= (1 << i);
                    }

                    //- check if the code matches some invalid code
                    if( problematicVertexCodes.found(vertexCode) )
                    {
                        # ifdef USE_OMP
                        # pragma omp critical(problematicConnection)
                        # endif
                        problematicConnections.append(labelPair(leafI, vI));
                    }

                    //- check if the code matches some invalid code
                    forAll(edgePositionsAtVertex[vI], dir)
                    {
                        const FixedList<label, 4>& activeNodes =
                            edgePositionsAtVertex[vI][dir];

                        direction edgeCode(0);
                        forAll(activeNodes, i)
                            if( boxType[neighbours[activeNodes[i]]] & MESHCELL )
                                edgeCode |= (1 << activeNodes[i]);

                        if( problematicEdgeCodes[vI].found(edgeCode) )
                        {
                            # ifdef USE_OMP
                            # pragma omp critical(problematicConnection)
                            # endif
                            problematicConnections.append(labelPair(leafI, vI));

                            break;
                        }
                    }
                }
            }
        }

        Info << "Found " << problematicConnections.size()
                 << " problematic connections " << endl;

        labelLongList changedLeaves;
        forAll(problematicConnections, i)
        {
            const labelPair& lp = problematicConnections[i];

            const FixedList<direction, 6>& bbFaces = bndFaces[lp.first()];

            //- do not use octree boxes with pairs of faces classified
            //- as bnd faces
            label nBndPairs(0);

            for(label f=0;f<3;++f)
            {
                bool containsNode(false);

                const label fI = 2 * f;

                for(label pI=0;pI<4;++pI)
                    if( meshOctreeCube::faceNodes_[fI][pI] == lp.second() )
                    {
                        containsNode = true;
                        break;
                    }

                if( containsNode && (bbFaces[fI] & 1) && (bbFaces[fI+1] & 1) )
                {
                    ++nBndPairs;
                }
            }

            if( nBndPairs > 1 )
            {
                boxType[lp.first()] = BOUNDARY;
                changedLeaves.append(lp.first());
                foundProblematic = true;
            }
        }

        reduce(foundProblematic, maxOp<bool>());

        if( Pstream::parRun() )
        {
            LongList<meshOctreeCubeCoordinates> dts;
            forAll(changedLeaves, i)
                dts.append(octree_.returnLeaf(changedLeaves[i]));

            LongList<meshOctreeCubeCoordinates> rcvLeaves;
            octree_.exchangeRequestsWithNeighbourProcessors(dts, rcvLeaves);

            forAll(rcvLeaves, i)
            {
                const label leafI =
                    octree_.findLeafLabelForPosition(rcvLeaves[i]);

                if( leafI < 0 )
                    continue;

                boxType[leafI] = BOUNDARY;
            }
        }

    } while( foundProblematic );


    clearNodeAddressing();
    clearOctreeFaces();
    clearAddressing();

    Info << "Finished checking the surface of selected boxes" << endl;
}

void meshOctreeAddressing::reviseLeaves()
{
    if( meshDict_.found("nonManifoldMeshing") )
    {
        const bool nonManifoldMeshing(meshDict_.lookup("nonManifoldMeshing"));

        if( nonManifoldMeshing )
            return;
    }

    if( meshDict_.found("allowDisconnectedDomains") )
    {
        const bool disconnectedDomains
        (
            meshDict_.lookup("allowDisconnectedDomains")
        );

        if( disconnectedDomains )
            return;
    }

    //- this is not yet operational
    return;

    const triSurf& surf = octree_.surface();

    const boundBox& rootBox = octree_.rootBox();

    //- give labels to cubes which will be used as mesh cells
    const List<direction>& boxType = this->boxType();

    //- check the number of independent groups of MESHCELL boxes
    boolList changeType(boxType.size(), false);

    DynList<label> ct, allNeighbours, neiInGroup, otherNeighbours;
    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100) \
    private(ct, allNeighbours, neiInGroup, otherNeighbours)
    # endif
    forAll(boxType, leafI)
    {
        if( !(boxType[leafI] & (MESHCELL|BOUNDARY)) )
            continue;
        if( !octree_.hasContainedTriangles(leafI) )
            continue;

        ct.clear();
        octree_.containedTriangles(leafI, ct);

        const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);

        if( Pstream::parRun() && (oc.procNo() != Pstream::myProcNo()) )
            continue;

        const point c = oc.centre(rootBox);

        //- find the nearest point at the surface
        point nearest(c);
        scalar distSq(VGREAT);

        forAll(ct, i)
        {
            const point np = help::nearestPointOnTheTriangle(ct[i], surf, c);

            const scalar dSq = magSqr(np - c);

            if( dSq < distSq )
            {
                nearest = np;
                distSq = dSq;
            }
        }

        //- types of interfaces
        FixedList<direction, 6> neiTypes(NONE);
        for(label dirI=0;dirI<6;++dirI)
        {
            DynList<label> neighbours;
            octree_.findNeighboursInDirection(leafI, dirI, neighbours);

            forAll(neighbours, i)
            {
                neiTypes[dirI] |= octree_.returnLeaf(neighbours[i]).cubeType();
            }
        }

        bool useCell(true);
        direction allTypes(0);

        FixedList<point, 6> faceCentres;
        oc.faceCentres(rootBox, faceCentres);
        const vector vec = nearest - c;
        forAll(faceCentres, dirI)
        {
            allTypes |= neiTypes[dirI];

            if
            (
                (
                    (((faceCentres[dirI] - c) & vec) > 0.0) &&
                    (neiTypes[dirI] & meshOctreeCube::INSIDE)
                ) ||
                (
                    (
                        (neiTypes[dirI] & meshOctreeCube::OUTSIDE) ||
                        (neiTypes[dirI] & meshOctreeCube::UNKNOWN)
                    ) &&
                    (((faceCentres[dirI] - c) & vec) < 0.0)
                )
            )
                useCell = false;
        }

        //- check if the cell changes its type
        if( useCell ^ (boxType[leafI] & MESHCELL) )
        {
            bool changeCellType(true);
            Pout << "Leaf changes its type" << endl;
            if
            (
                (allTypes & meshOctreeCube::INSIDE) &&
                !(allTypes & meshOctreeCube::OUTSIDE) &&
                !(allTypes & meshOctreeCube::UNKNOWN)
            )
            {
                //- check if all internal boxes are reachable over faces
                allNeighbours.clear();
                octree_.findAllLeafNeighbours(leafI, allNeighbours);

                label groupI(0);
                neiInGroup.setSize(allNeighbours.size());
                neiInGroup = -1;

                forAll(allNeighbours, i)
                {
                    const meshOctreeCubeBasic& oc =
                        octree_.returnLeaf(allNeighbours[i]);
                    if( !(oc.cubeType() & meshOctreeCube::INSIDE) )
                        continue;
                    if( neiInGroup[i] != -1 )
                        continue;

                    DynList<label> front;
                    front.append(i);

                    while( front.size() )
                    {
                        const label posJ = front.removeLastElement();

                        neiInGroup[posJ] = groupI;

                        const label leafJ = allNeighbours[posJ];

                        otherNeighbours.clear();
                        octree_.findNeighboursForLeaf(leafJ, otherNeighbours);

                        forAll(otherNeighbours, j)
                        {
                            const label otherPos =
                                allNeighbours.containsAtPosition
                                (
                                    otherNeighbours[j]
                                );
                            if( (otherPos >= 0) && (neiInGroup[otherPos] < 0) )
                            {
                                front.append(otherPos);
                            }
                        }
                    }

                    ++groupI;
                }

                if( groupI > 1 )
                {
                    Pout << "1.Number of cell groups " << groupI << endl;
                    Pout << "1.Nei groups " << neiInGroup << endl;
                    changeCellType = false;
                }
            }
            else if
            (
                !(allTypes & meshOctreeCube::INSIDE) &&
                (allTypes & meshOctreeCube::OUTSIDE) &&
                !(allTypes & meshOctreeCube::UNKNOWN)
            )
            {
                //- check if all outside boxes are reachable over faces
                allNeighbours.clear();
                octree_.findAllLeafNeighbours(leafI, allNeighbours);

                label groupI(0);
                neiInGroup.setSize(allNeighbours.size());
                neiInGroup = -1;

                forAll(allNeighbours, i)
                {
                    const meshOctreeCubeBasic& oc =
                        octree_.returnLeaf(allNeighbours[i]);
                    if( !(oc.cubeType() & meshOctreeCube::OUTSIDE) )
                        continue;
                    if( neiInGroup[i] != -1 )
                        continue;

                    DynList<label> front;
                    front.append(i);

                    while( front.size() )
                    {
                        const label posJ = front.removeLastElement();

                        neiInGroup[posJ] = groupI;

                        const label leafJ = allNeighbours[posJ];

                        otherNeighbours.clear();
                        octree_.findNeighboursForLeaf(leafJ, otherNeighbours);

                        forAll(otherNeighbours, j)
                        {
                            const label otherPos =
                                allNeighbours.containsAtPosition
                                (
                                    otherNeighbours[j]
                                );
                            if( (otherPos >= 0) && (neiInGroup[otherPos] < 0) )
                            {
                                front.append(otherPos);
                            }
                        }
                    }

                    ++groupI;
                }

                if( groupI > 1 )
                {
                    Pout << "2.Number of cell groups " << groupI << endl;
                    Pout << "2.Nei groups " << neiInGroup << endl;
                    changeCellType = false;
                }
            }
            else if
            (
                !(allTypes & meshOctreeCube::INSIDE) &&
                !(allTypes & meshOctreeCube::OUTSIDE) &&
                (allTypes & meshOctreeCube::UNKNOWN)
            )
            {
                //- check if all unknown boxes are reachable over faces
                allNeighbours.clear();
                octree_.findAllLeafNeighbours(leafI, allNeighbours);

                label groupI(0);
                neiInGroup.setSize(allNeighbours.size());
                neiInGroup = -1;

                forAll(allNeighbours, i)
                {
                    const meshOctreeCubeBasic& oc =
                        octree_.returnLeaf(allNeighbours[i]);
                    if( !(oc.cubeType() & meshOctreeCube::UNKNOWN) )
                        continue;
                    if( neiInGroup[i] != -1 )
                        continue;

                    DynList<label> front;
                    front.append(i);

                    while( front.size() )
                    {
                        const label posJ = front.removeLastElement();

                        neiInGroup[posJ] = groupI;

                        const label leafJ = allNeighbours[posJ];

                        otherNeighbours.clear();
                        octree_.findNeighboursForLeaf(leafJ, otherNeighbours);

                        forAll(otherNeighbours, j)
                        {
                            const label otherPos =
                                allNeighbours.containsAtPosition
                                (
                                    otherNeighbours[j]
                                );
                            if( (otherPos >= 0) && (neiInGroup[otherPos] < 0) )
                            {
                                front.append(otherPos);
                            }
                        }
                    }

                    ++groupI;
                }

                if( groupI > 1 )
                {
                    Pout << "3.Number of cell groups " << groupI << endl;
                    Pout << "3.Nei groups " << neiInGroup << endl;
                    changeCellType = false;
                }
            }

            if( changeCellType )
                changeType[leafI] = true;
        }
    }

    forAll(changeType, leafI)
    {
        if( changeType[leafI] )
        {
            if( !(boxType[leafI] & MESHCELL) )
            {
                setBoxType(leafI, MESHCELL);
            }
            else
            {
                setBoxType(leafI, NONE);
            }
        }
    }

    //- set BOUNDARY flag to boxes which do not have a MESHCELL flag
    markBOUNDARYLeaves();

    clearNodeAddressing();
    clearOctreeFaces();
    clearAddressing();
    clearParallelAddressing();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
