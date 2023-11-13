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
#include "demandDrivenData.H"
#include "IOdictionary.H"
#include "helperFunctions.H"
#include "triSurf.H"
#include "meshOctreeModifier.H"
#include "faceListPMG.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void meshOctreeAddressing::clearOut()
{
    clearNodeAddressing();
    clearBoxTypes();
    clearOctreeFaces();
    clearAddressing();
    clearParallelAddressing();
}

void meshOctreeAddressing::clearNodeAddressing()
{
    nNodes_ = 0;
    deleteDemandDrivenData(octreePointsPtr_);
    deleteDemandDrivenData(nodeLabelsPtr_);
    deleteDemandDrivenData(nodeLeavesPtr_);

    deleteDemandDrivenData(nodeTypePtr_);
}

void meshOctreeAddressing::clearBoxTypes()
{
    deleteDemandDrivenData(boxTypePtr_);
    deleteDemandDrivenData(boxGroupPtr_);
}

void meshOctreeAddressing::clearOctreeFaces()
{
    deleteDemandDrivenData(octreeFacesPtr_);
    deleteDemandDrivenData(octreeFacesOwnersPtr_);
    deleteDemandDrivenData(octreeFacesNeighboursPtr_);
}

void meshOctreeAddressing::clearAddressing()
{
    deleteDemandDrivenData(leafFacesPtr_);
    deleteDemandDrivenData(nodeFacesPtr_);
    deleteDemandDrivenData(leafLeavesPtr_);
    deleteDemandDrivenData(octreeEdgesPtr_);
    deleteDemandDrivenData(edgeLeavesPtr_);
    deleteDemandDrivenData(leafEdgesPtr_);
    deleteDemandDrivenData(nodeEdgesPtr_);
    deleteDemandDrivenData(faceEdgesPtr_);
    deleteDemandDrivenData(edgeFacesPtr_);
}

void meshOctreeAddressing::clearParallelAddressing()
{
    deleteDemandDrivenData(globalPointLabelPtr_);
    deleteDemandDrivenData(globalPointToLocalPtr_);
    deleteDemandDrivenData(pointProcsPtr_);
    deleteDemandDrivenData(globalFaceLabelPtr_);
    deleteDemandDrivenData(globalFaceToLocalPtr_);
    deleteDemandDrivenData(faceProcsPtr_);
    deleteDemandDrivenData(globalLeafLabelPtr_);
    deleteDemandDrivenData(globalLeafToLocalPtr_);
    deleteDemandDrivenData(leafAtProcsPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from octree and IOdictionary
meshOctreeAddressing::meshOctreeAddressing
(
    const meshOctree& mo,
    const dictionary& dict,
    bool useDATABoxes
)
:
    octree_(mo),
    meshDict_(dict),
    useDATABoxes_(useDATABoxes),
    nNodes_(0),
    nGroups_(-1),
    octreePointsPtr_(NULL),
    nodeLabelsPtr_(NULL),
    nodeLeavesPtr_(NULL),
    boxTypePtr_(NULL),
    boxGroupPtr_(NULL),
    nodeTypePtr_(NULL),
    octreeFacesPtr_(NULL),
    octreeFacesOwnersPtr_(NULL),
    octreeFacesNeighboursPtr_(NULL),
    leafFacesPtr_(NULL),
    nodeFacesPtr_(NULL),
    leafLeavesPtr_(NULL),
    octreeEdgesPtr_(NULL),
    edgeLeavesPtr_(NULL),
    leafEdgesPtr_(NULL),
    nodeEdgesPtr_(NULL),
    faceEdgesPtr_(NULL),
    edgeFacesPtr_(NULL),
    globalPointLabelPtr_(NULL),
    globalPointToLocalPtr_(NULL),
    pointProcsPtr_(NULL),
    globalFaceLabelPtr_(NULL),
    globalFaceToLocalPtr_(NULL),
    faceProcsPtr_(NULL),
    globalLeafLabelPtr_(NULL),
    globalLeafToLocalPtr_(NULL),
    leafAtProcsPtr_(NULL)
{
    if( !useDATABoxes && dict.found("keepCellsIntersectingBoundary") )
    {
        useDATABoxes_ = readBool(dict.lookup("keepCellsIntersectingBoundary"));
    }

    if( dict.found("nonManifoldMeshing") )
    {
        const bool nonManifoldMesh
        (
            readBool(dict.lookup("nonManifoldMeshing"))
        );

        if( nonManifoldMesh )
            useDATABoxes_ = true;
    }

    //- check for glued regions
    checkGluedRegions();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshOctreeAddressing::~meshOctreeAddressing()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

sizeType meshOctreeAddressing::printAllocated() const
{
    sizeType totalMemUsage(sizeof(meshOctreeAddressing));

    const label toMb = 1048576;

    if( octreePointsPtr_ )
    {
        Pout << "octreePointsPtr_ " << octreePointsPtr_->byteSize()/toMb
             << endl;
        totalMemUsage += octreePointsPtr_->byteSize();
    }

    if( nodeLabelsPtr_ )
    {
        Pout << "nodeLabelsPtr_ "
             << nodeLabelsPtr_->printAllocated()/toMb << endl;
        totalMemUsage += nodeLabelsPtr_->printAllocated();
    }

    if( nodeLeavesPtr_ )
    {
        Pout << "nodeLeavesPtr_ "
             << nodeLeavesPtr_->printAllocated()/toMb << endl;
        totalMemUsage += nodeLeavesPtr_->printAllocated();
    }

    if( boxTypePtr_ )
    {
        Pout << "boxTypePtr_ " << boxTypePtr_->byteSize()/toMb << endl;
        totalMemUsage += boxTypePtr_->byteSize();
    }

    if( boxGroupPtr_ )
    {
        Pout << "boxGroupPtr_ " << boxGroupPtr_->byteSize()/toMb << endl;
        totalMemUsage += boxGroupPtr_->byteSize();
    }

    if( nodeTypePtr_ )
    {
        Pout << "nodeTypePtr_ " << nodeTypePtr_->byteSize()/toMb << endl;
        totalMemUsage += nodeTypePtr_->byteSize();
    }

    if( octreeFacesPtr_ )
    {
        Pout << "octreeFacesPtr_ "
             << octreeFacesPtr_->printAllocated()/toMb << endl;
        totalMemUsage += octreeFacesPtr_->printAllocated();
    }

    if( octreeFacesOwnersPtr_ )
    {
        Pout << "octreeFacesOwnersPtr_ "
             << octreeFacesOwnersPtr_->printAllocated()/toMb << endl;
        totalMemUsage += octreeFacesOwnersPtr_->printAllocated();
    }

    if( octreeFacesNeighboursPtr_ )
    {
        Pout << "octreeFacesNeighboursPtr_ "
             << octreeFacesNeighboursPtr_->printAllocated()/toMb << endl;
        totalMemUsage += octreeFacesNeighboursPtr_->printAllocated();
    }

    if( leafFacesPtr_ )
    {
        Pout << "leafFacesPtr_ "
             << leafFacesPtr_->printAllocated()/toMb << endl;
        totalMemUsage += leafFacesPtr_->printAllocated();
    }

    if( nodeFacesPtr_ )
    {
        Pout << "nodeFacesPtr_ "
             << nodeFacesPtr_->printAllocated()/toMb << endl;
        totalMemUsage += nodeFacesPtr_->printAllocated();
    }

    if( leafLeavesPtr_ )
    {
        Pout << "leafLeavesPtr_ "
             << leafLeavesPtr_->printAllocated()/toMb << endl;
        totalMemUsage += leafLeavesPtr_->printAllocated();
    }

    if( octreeEdgesPtr_ )
    {
        Pout << "octreeEdgesPtr_ "
             << octreeEdgesPtr_->printAllocated()/toMb << endl;
        totalMemUsage += octreeEdgesPtr_->printAllocated();
    }

    if( edgeLeavesPtr_ )
    {
        Pout << "edgeLeavesPtr_ "
             << edgeLeavesPtr_->printAllocated()/toMb << endl;
        totalMemUsage += edgeLeavesPtr_->printAllocated();
    }

    if( leafEdgesPtr_ )
    {
        Pout << "leafEdgesPtr_ "
             << leafEdgesPtr_->printAllocated()/toMb << endl;
        totalMemUsage += leafEdgesPtr_->printAllocated();
    }

    if( nodeEdgesPtr_ )
    {
        Pout << "nodeEdgesPtr_ "
             << nodeEdgesPtr_->printAllocated()/toMb << endl;
        totalMemUsage += nodeEdgesPtr_->printAllocated();
    }

    if( faceEdgesPtr_ )
    {
        Pout << "faceEdgesPtr_ "
             << faceEdgesPtr_->printAllocated()/toMb << endl;
        totalMemUsage += faceEdgesPtr_->printAllocated();
    }

    if( edgeFacesPtr_ )
    {
        Pout << "edgeFacesPtr_ "
             << edgeFacesPtr_->printAllocated()/toMb << endl;
        totalMemUsage += edgeFacesPtr_->printAllocated();
    }

    if( globalPointLabelPtr_ )
    {
        Pout << "globalPointLabelPtr_ "
             << globalPointLabelPtr_->printAllocated()/toMb << endl;
        totalMemUsage += globalPointLabelPtr_->printAllocated();
    }

    if( globalPointToLocalPtr_ )
    {
        Pout << "globalPointToLocalPtr_ " << endl;
    }

    if( pointProcsPtr_ )
    {
        Pout << "pointProcsPtr_ "
             << pointProcsPtr_->printAllocated()/toMb << endl;
        totalMemUsage += pointProcsPtr_->printAllocated();
    }

    if( globalFaceLabelPtr_ )
    {
        Pout << "globalFaceLabelPtr_ "
             << globalFaceLabelPtr_->printAllocated()/toMb
             << endl;
        totalMemUsage += globalFaceLabelPtr_->printAllocated();
    }

    if( globalFaceToLocalPtr_ )
    {
        Pout << "globalFaceToLocalPtr_" << endl;
    }

    if( faceProcsPtr_ )
    {
        Pout << "faceProcsPtr_ "
             << faceProcsPtr_->printAllocated()/toMb << endl;
        totalMemUsage += faceProcsPtr_->printAllocated();
    }

    if( globalLeafLabelPtr_ )
    {
        Pout << "globalLeafLabelPtr_ "
             << globalLeafLabelPtr_->printAllocated()/toMb << endl;
        totalMemUsage += globalLeafLabelPtr_->printAllocated();
    }

    if( globalLeafToLocalPtr_ )
    {
        Pout << "globalLeafToLocalPtr_" << endl;
    }

    if( leafAtProcsPtr_ )
    {
        Pout << "leafAtProcsPtr_ "
             << leafAtProcsPtr_->printAllocated()/toMb << endl;
        totalMemUsage += leafAtProcsPtr_->printAllocated();
    }

    Pout << "meshOctreeAddressing: memory usage " << totalMemUsage/toMb
         << " Mb." << endl;

    return totalMemUsage;
}

bool meshOctreeAddressing::isIntersectedFace(const label fI) const
{
    const labelLongList& owner = octreeFaceOwner();
    const labelLongList& neighbour = octreeFaceNeighbour();

    if( neighbour[fI] < 0 )
        return false;

    Map<label> nAppearances;
    DynList<label> triangles;
    octree_.containedTriangles(owner[fI], triangles);
    forAll(triangles, triI)
    {
        if( nAppearances.found(triangles[triI]) )
        {
            ++nAppearances[triangles[triI]];
        }
        else
        {
            nAppearances.insert(triangles[triI], 1);
        }
    }

    triangles.clear();
    octree_.containedTriangles(neighbour[fI], triangles);
    forAll(triangles, triI)
    {
        if( nAppearances.found(triangles[triI]) )
        {
            ++nAppearances[triangles[triI]];
        }
        else
        {
            nAppearances.insert(triangles[triI], 1);
        }
    }

    forAllConstIter(Map<label>, nAppearances, iter)
    {
        if( iter() == 2 )
        {
            if
            (
                octree_.returnLeaf(owner[fI]).level() ==
                octree_.returnLeaf(neighbour[fI]).level()
            )
                return true;

            //- check intersection by geometric testing
            const triSurf& surf = octree_.surface();
            const pointField& points = this->octreePoints();
            const VRWGraph& faces = this->octreeFaces();

            face f(faces.sizeOfRow(fI));
            forAll(f, pI)
                f[pI] = faces(fI, pI);

            if(
                help::doFaceAndTriangleIntersect
                (
                    surf,
                    iter.key(),
                    f,
                    points
                )
            )
                return true;
        }
    }

    return false;
}

bool meshOctreeAddressing::isIntersectedEdge(const label eI) const
{
    const VRWGraph& edgeCubes = this->edgeLeaves();

    Map<label> nAppearances;
    DynList<label> triangles;
    bool sameLevel(true);

    forAllRow(edgeCubes, eI, i)
    {
        const label leafI = edgeCubes(eI, i);
        if( !octree_.hasContainedTriangles(leafI) )
            return false;

        if
        (
            octree_.returnLeaf(leafI).level() !=
            octree_.returnLeaf(edgeCubes(eI, 0)).level()
        )
            sameLevel = false;

        triangles.clear();
        octree_.containedTriangles(leafI, triangles);
        forAll(triangles, triI)
        {
            if( nAppearances.found(triangles[triI]) )
            {
                ++nAppearances[triangles[triI]];
            }
            else
            {
                nAppearances.insert(triangles[triI], 1);
            }
        }
    }

    forAllConstIter(Map<label>, nAppearances, iter)
    {
        if( iter() == edgeCubes.sizeOfRow(eI) )
        {
            if( sameLevel )
                return true;

            //- check for geometric intersection
            const LongList<edge>& edges = this->octreeEdges();
            const pointField& points = this->octreePoints();
            point intersection(vector::zero);

            if(
                help::triLineIntersection
                (
                    octree_.surface(),
                    iter.key(),
                    points[edges[eI].start()],
                    points[edges[eI].end()],
                    intersection
                )
            )
                return true;
        }
    }

    return false;
}

void meshOctreeAddressing::edgeIntersections
(
    const label eI,
    DynList<point>& intersections
) const
{
    intersections.clear();

    const LongList<edge>& edges = this->octreeEdges();
    const pointField& points = this->octreePoints();
    const VRWGraph& edgeCubes = this->edgeLeaves();
    const scalar tol =
        SMALL * mag(points[edges[eI].start()] - points[edges[eI].end()]);

    Map<label> nAppearances;
    DynList<label> triangles;

    forAllRow(edgeCubes, eI, i)
    {
        const label leafI = edgeCubes(eI, i);
        if( !octree_.hasContainedTriangles(leafI) )
            return;

        triangles.clear();
        octree_.containedTriangles(leafI, triangles);
        forAll(triangles, triI)
        {
            if( nAppearances.found(triangles[triI]) )
            {
                ++nAppearances[triangles[triI]];
            }
            else
            {
                nAppearances.insert(triangles[triI], 1);
            }
        }
    }

    point intersection(vector::zero);

    forAllConstIter(Map<label>, nAppearances, iter)
    {
        if( iter() == edgeCubes.sizeOfRow(eI) )
        {
            //- check for geometric intersection
            const bool intersectionExists =
                help::triLineIntersection
                (
                    octree_.surface(),
                    iter.key(),
                    points[edges[eI].start()],
                    points[edges[eI].end()],
                    intersection
                );

            if( intersectionExists )
            {
                bool store(true);
                forAll(intersections, i)
                    if( mag(intersections[i] - intersection) <= tol )
                        store = false;

                if( store )
                    intersections.append(intersection);
            }
        }
    }
}

void meshOctreeAddressing::cubesAroundEdge
(
    const label leafI,
    const direction eI,
    FixedList<label, 4>& edgeCubes
) const
{
    const VRWGraph& nl = this->nodeLabels();
    const label nodeI = nl(leafI, meshOctreeCubeCoordinates::edgeNodes_[eI][0]);
    const FRWGraph<label, 8>& pLeaves = this->nodeLeaves();

    switch( eI )
    {
        case 0: case 1: case 2: case 3:
        {
            edgeCubes[0] = pLeaves(nodeI, 1);
            edgeCubes[1] = pLeaves(nodeI, 3);
            edgeCubes[2] = pLeaves(nodeI, 5);
            edgeCubes[3] = pLeaves(nodeI, 7);
        } break;
        case 4: case 5: case 6: case 7:
        {
            edgeCubes[0] = pLeaves(nodeI, 2);
            edgeCubes[1] = pLeaves(nodeI, 3);
            edgeCubes[2] = pLeaves(nodeI, 6);
            edgeCubes[3] = pLeaves(nodeI, 7);
        } break;
        case 8: case 9: case 10: case 11:
        {
            edgeCubes[0] = pLeaves(nodeI, 4);
            edgeCubes[1] = pLeaves(nodeI, 5);
            edgeCubes[2] = pLeaves(nodeI, 6);
            edgeCubes[3] = pLeaves(nodeI, 7);
        } break;
        default:
        {
            FatalErrorIn
            (
                "void tetMeshExtractorOctree::cubesAroundEdge(const label,"
                "const direction, FixedList<label, 4>&)"
            ) << "Invalid edge specified!!" << abort(FatalError);
        } break;
    };
}

label meshOctreeAddressing::findEdgeCentre
(
    const label leafI,
    const direction eI
) const
{
    if( octree_.isQuadtree() && eI >= 8 )
        return -1;

    const meshOctreeCubeBasic& oc = octree_.returnLeaf(leafI);
    const VRWGraph& nl = this->nodeLabels();
    const label nodeI = nl(leafI, meshOctreeCubeCoordinates::edgeNodes_[eI][0]);
    const FRWGraph<label, 8>& pLeaves = this->nodeLeaves();

    const direction level = oc.level();

    label fI(-1);
    switch( eI )
    {
        case 0: case 1: case 2: case 3:
        {
            fI = 1;
        } break;
        case 4: case 5: case 6: case 7:
        {
            fI = 3;
        } break;
        case 8: case 9: case 10: case 11:
        {
            fI = 5;
        } break;
        default:
        {
            FatalErrorIn
            (
                "label meshOctreeAddressing::findEdgeCentre"
                "(const label leafI, const direction eI) const"
            ) << "Invalid edge specified!!" << abort(FatalError);
        } break;
    };

    for(label i=0;i<4;++i)
    {
        const label fNode = meshOctreeCubeCoordinates::faceNodes_[fI][i];

        const label leafJ = pLeaves(nodeI, fNode);

        if( leafJ < 0 )
            continue;

        if( octree_.returnLeaf(leafJ).level() > level )
        {
            const label shift = (i+2)%4;
            return nl(leafJ, meshOctreeCubeCoordinates::faceNodes_[fI][shift]);
        }
    }

    return -1;
}

void meshOctreeAddressing::createFaces
(
    faceList& faces,
    labelList& owner,
    labelList& neighbour
) const
{
    const VRWGraph& nodeLabels = this->nodeLabels();
    const List<direction>& boxType = this->boxType();
    this->nodeLeaves();

    label nFaces(0);
    labelList numFacesAtTask;
    List<labelLongList> faceSizeAtTask;

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        const label nTasks = 5 * omp_get_num_threads();
        const label chunkSize = max(boxType.size() / nTasks, 1);
        # else
        const label nTasks(1);
        const label chunkSize = boxType.size();
        # endif

        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            numFacesAtTask.setSize(nTasks);
            numFacesAtTask = 0;

            faceSizeAtTask.setSize(nTasks);

            for(label taskI=0;taskI<nTasks;++taskI)
            {
                # ifdef USE_OMP
                # pragma omp task default(shared) firstprivate(taskI)
                # endif
                {
                    const label sl = taskI * chunkSize;
                    label el = Foam::min(sl+chunkSize, boxType.size());
                    if( taskI == (nTasks-1) )
                        el = boxType.size();

                    for(label leafI=sl;leafI<el;++leafI)
                    {
                        const meshOctreeCubeBasic& oc =
                            octree_.returnLeaf(leafI);

                        if( boxType[leafI] & MESHCELL )
                        {
                            FixedList<label, 12> edgeCentreLabel(-1);
                            for(label i=0;i<12;++i)
                                edgeCentreLabel[i] = findEdgeCentre(leafI, i);

                            for(label fI=0;fI<6;++fI)
                            {
                                DynList<label> neighbours;
                                octree_.findNeighboursInDirection
                                (
                                    leafI,
                                    fI,
                                    neighbours
                                );

                                if( neighbours.size() != 1 )
                                    continue;

                                const label nei = neighbours[0];

                                //- stop if the neighbour is on other processor
                                if( nei == meshOctreeCubeBasic::OTHERPROC )
                                    continue;

                                const label* fNodes =
                                    meshOctreeCubeCoordinates::faceNodes_[fI];
                                const label* fEdges =
                                    meshOctreeCubeCoordinates::faceEdges_[fI];

                                //- create face
                                DynList<label, 8> f;
                                for(label pI=0;pI<4;++pI)
                                {
                                    const label nI = fNodes[pI];
                                    const label feI = fEdges[pI];

                                    f.append(nodeLabels(leafI, nI));

                                    if( edgeCentreLabel[feI] != -1 )
                                        f.append(edgeCentreLabel[feI]);
                                }

                                if( nei < 0 )
                                {
                                    //- face is at the boundary of the octree
                                    faceSizeAtTask[taskI].append(f.size());
                                    ++numFacesAtTask[taskI];
                                }
                                else if( boxType[nei] & MESHCELL )
                                {
                                    //- face is an internal face
                                    if( nei > leafI )
                                    {
                                        faceSizeAtTask[taskI].append(f.size());
                                        ++numFacesAtTask[taskI];
                                    }
                                    else if
                                    (
                                        octree_.returnLeaf(nei).level() <
                                        oc.level()
                                    )
                                    {
                                        //- append a reversed face
                                        label i(1);
                                        for(label j=f.size()-1;j>i;--j)
                                        {
                                            const label add = f[j];
                                            f[j] = f[i];
                                            f[i] = add;
                                            ++i;
                                        }

                                        faceSizeAtTask[taskI].append(f.size());
                                        ++numFacesAtTask[taskI];
                                    }
                                }
                                else if( boxType[nei] & BOUNDARY )
                                {
                                    //- face is at the boundary of the mesh
                                    faceSizeAtTask[taskI].append(f.size());
                                    ++numFacesAtTask[taskI];
                                }
                            }
                        }
                        else if( boxType[leafI] & BOUNDARY )
                        {
                            for(label fI=0;fI<6;++fI)
                            {
                                DynList<label> neighbours;
                                octree_.findNeighboursInDirection
                                (
                                    leafI,
                                    fI,
                                    neighbours
                                );



                                if( neighbours.size() != 1 )
                                    continue;
                                const label nei = neighbours[0];
                                if( nei < 0 )
                                    continue;

                                const label* fNodes =
                                    meshOctreeCubeCoordinates::faceNodes_[fI];

                                if(
                                    (boxType[nei] & MESHCELL) &&
                                    (
                                        octree_.returnLeaf(nei).level() <
                                        oc.level()
                                    )
                                )
                                {
                                    //- add a boundary face
                                    FixedList<label, 4> cf;
                                    for(label i=0;i<4;++i)
                                    {
                                        cf[i] = nodeLabels(leafI, fNodes[i]);
                                    }

                                    faceSizeAtTask[taskI].append(cf.size());
                                    ++numFacesAtTask[taskI];
                                }
                            }
                        }
                    }

                }
            }
        }

        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            //- count the total number of faces
            nFaces = sum(numFacesAtTask);

            //- set the size of all lists
            owner.setSize(nFaces);
            neighbour.setSize(nFaces);
            faces.setSize(nFaces);

            //- allocated all faces
            for(label taskI=0;taskI<nTasks;++taskI)
            {
                # ifdef USE_OMP
                # pragma omp task default(shared) firstprivate(taskI)
                # endif
                {
                    label sf(0);
                    for(label i=0;i<taskI;++i)
                        sf += numFacesAtTask[i];

                    forAll(faceSizeAtTask[taskI], j)
                        faces[sf+j].setSize(faceSizeAtTask[taskI][j]);
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            faceSizeAtTask.clear();

            for(label taskI=0;taskI<nTasks;++taskI)
            {
                # ifdef USE_OMP
                # pragma omp task default(shared) firstprivate(taskI)
                # endif
                {
                    const label sl = taskI * chunkSize;
                    label el = Foam::min(sl+chunkSize, boxType.size());
                    if( taskI == (nTasks-1) )
                        el = boxType.size();

                    label faceI = 0;
                    for(label i=0;i<taskI;++i)
                        faceI += numFacesAtTask[i];

                    for(label leafI=sl;leafI<el;++leafI)
                    {
                        const meshOctreeCubeBasic& oc =
                            octree_.returnLeaf(leafI);

                        if( boxType[leafI] & MESHCELL )
                        {
                            FixedList<label, 12> edgeCentreLabel(-1);
                            for(label i=0;i<12;++i)
                                edgeCentreLabel[i] = findEdgeCentre(leafI, i);

                            for(label fI=0;fI<6;++fI)
                            {
                                DynList<label> neighbours;
                                octree_.findNeighboursInDirection
                                (
                                    leafI,
                                    fI,
                                    neighbours
                                );

                                if( neighbours.size() != 1 )
                                    continue;

                                const label nei = neighbours[0];

                                //- stop if the neighbour is on other processor
                                if( nei == meshOctreeCubeBasic::OTHERPROC )
                                    continue;

                                const label* fNodes =
                                    meshOctreeCubeCoordinates::faceNodes_[fI];
                                const label* fEdges =
                                    meshOctreeCubeCoordinates::faceEdges_[fI];

                                //- create face
                                DynList<label, 8> f;
                                for(label pI=0;pI<4;++pI)
                                {
                                    const label nI = fNodes[pI];
                                    const label feI = fEdges[pI];

                                    f.append(nodeLabels(leafI, nI));

                                    if( edgeCentreLabel[feI] != -1 )
                                        f.append(edgeCentreLabel[feI]);
                                }

                                if( nei < 0 )
                                {
                                    //- face is at the boundary of the octree
                                    forAll(faces[faceI], pI)
                                        faces[faceI][pI] = f[pI];
                                    owner[faceI] = leafI;
                                    neighbour[faceI] = -1;
                                    ++faceI;
                                }
                                else if( boxType[nei] & MESHCELL )
                                {
                                    //- face is an internal face
                                    if( nei > leafI )
                                    {
                                        forAll(faces[faceI], pI)
                                            faces[faceI][pI] = f[pI];
                                        owner[faceI] = leafI;
                                        neighbour[faceI] = nei;
                                        ++faceI;
                                    }
                                    else if
                                    (
                                        octree_.returnLeaf(nei).level() <
                                        oc.level()
                                    )
                                    {
                                        //- append a reversed face
                                        label i(1);
                                        for(label j=f.size()-1;j>i;--j)
                                        {
                                            const label add = f[j];
                                            f[j] = f[i];
                                            f[i] = add;
                                            ++i;
                                        }

                                        forAll(faces[faceI], pI)
                                            faces[faceI][pI] = f[pI];
                                        owner[faceI] = nei;
                                        neighbour[faceI] = leafI;
                                        ++faceI;
                                    }
                                }
                                else if( boxType[nei] & BOUNDARY )
                                {
                                    //- face is at the boundary of the mesh
                                    forAll(faces[faceI], pI)
                                        faces[faceI][pI] = f[pI];
                                    owner[faceI] = leafI;
                                    neighbour[faceI] = nei;
                                    ++faceI;
                                }
                            }
                        }
                        else if( boxType[leafI] & BOUNDARY )
                        {
                            for(label fI=0;fI<6;++fI)
                            {
                                DynList<label> neighbours;
                                octree_.findNeighboursInDirection
                                (
                                    leafI,
                                    fI,
                                    neighbours
                                );

                                if( neighbours.size() != 1 )
                                    continue;
                                const label nei = neighbours[0];
                                if( nei < 0 )
                                    continue;
                                if(
                                    (boxType[nei] & MESHCELL) &&
                                    (
                                        octree_.returnLeaf(nei).level() <
                                        oc.level()
                                    )
                                )
                                {
                                    //- add a boundary face
                                    const label* fNodes =
                                        meshOctreeCube::faceNodes_[fI];

                                    face& f = faces[faceI];
                                    f[0] = nodeLabels(leafI, fNodes[0]);
                                    f[1] = nodeLabels(leafI, fNodes[3]);
                                    f[2] = nodeLabels(leafI, fNodes[2]);
                                    f[3] = nodeLabels(leafI, fNodes[1]);

                                    owner[faceI] = nei;
                                    neighbour[faceI] = leafI;
                                    ++faceI;
                                }
                            }
                        }
                    }

                }
            }
        }
    }
}

void meshOctreeAddressing::clearPoints()
{
    deleteDemandDrivenData(octreePointsPtr_);
}

void meshOctreeAddressing::clearFacesAndNodeLabels()
{
    deleteDemandDrivenData(octreeFacesPtr_);
    deleteDemandDrivenData(octreeFacesOwnersPtr_);
    deleteDemandDrivenData(octreeFacesNeighboursPtr_);

    deleteDemandDrivenData(nodeLeavesPtr_);
    deleteDemandDrivenData(nodeLabelsPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
