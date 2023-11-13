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

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline label meshOctreeAddressing::numberOfNodes() const
{
    if( !nodeLabelsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline label meshOctreeAddressing::numberOfNodes() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        createNodeLabels();
    }

    return nNodes_;
}

inline const pointField& meshOctreeAddressing::octreePoints() const
{
    if( !octreePointsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const pointField&"
                " meshOctreeAddressing::octreePoints() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        createOctreePoints();
    }

    return *octreePointsPtr_;
}

inline const VRWGraph& meshOctreeAddressing::nodeLabels() const
{
    if( !nodeLabelsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph&"
                " meshOctreeAddressing::nodeLabels() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        createNodeLabels();
    }

    return *nodeLabelsPtr_;
}

inline const FRWGraph<label, 8>& meshOctreeAddressing::nodeLeaves() const
{
    if( !nodeLeavesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const FRWGraph<label, 8>&"
                " meshOctreeAddressing::nodeLeaves() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        createNodeLeaves();
    }

    return *nodeLeavesPtr_;
}

inline const List<direction>& meshOctreeAddressing::boxType() const
{
    if( !boxTypePtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const List<direction>&"
                " meshOctreeAddressing::boxType() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        findUsedBoxes();
    }

    return *boxTypePtr_;
}

inline void meshOctreeAddressing::setBoxType
(
    const label boxI,
    const direction type
)
{
    if( boxTypePtr_ )
    {
        List<direction>& boxType = *boxTypePtr_;
        boxType[boxI] |= type;
    }
}

inline labelList& meshOctreeAddressing::boxGroup() const
{
    if( !boxGroupPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline labelList& meshOctreeAddressing::boxGroup() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        findGroupsOfUsedBoxes();
    }

    return *boxGroupPtr_;
}

inline label meshOctreeAddressing::nGroups() const
{
    if( !boxGroupPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline label meshOctreeAddressing::nGroups() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        findGroupsOfUsedBoxes();
    }

    return nGroups_;
}

inline const meshOctree& meshOctreeAddressing::octree() const
{
    return octree_;
}

inline const List<direction>& meshOctreeAddressing::nodeType() const
{
    if( !nodeTypePtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const List<direction>&"
                " meshOctreeAddressing::nodeType() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calculateNodeType();
    }

    return *nodeTypePtr_;
}

inline const VRWGraph& meshOctreeAddressing::octreeFaces() const
{
    if( !octreeFacesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph&"
                " meshOctreeAddressing::octreeFaces() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        createOctreeFaces();
    }

    return *octreeFacesPtr_;
}

inline const labelLongList& meshOctreeAddressing::octreeFaceOwner() const
{
    if( !octreeFacesOwnersPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const labelLongList&"
                " meshOctreeAddressing::octreeFaceOwner() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        createOctreeFaces();
    }

    return *octreeFacesOwnersPtr_;
}

inline const labelLongList& meshOctreeAddressing::octreeFaceNeighbour() const
{
    if( !octreeFacesNeighboursPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const labelLongList&"
                " meshOctreeAddressing::octreeFaceNeighbour() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        createOctreeFaces();
    }

    return *octreeFacesNeighboursPtr_;
}

inline const VRWGraph& meshOctreeAddressing::leafFaces() const
{
    if( !leafFacesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph& meshOctreeAddressing::leafFaces() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calculateLeafFaces();
    }

    return *leafFacesPtr_;
}

inline const VRWGraph& meshOctreeAddressing::leafLeaves() const
{
    if( !leafLeavesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph&"
                " meshOctreeAddressing::leafLeaves() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calculateLeafLeaves();
    }

    return *leafLeavesPtr_;
}

inline const VRWGraph& meshOctreeAddressing::nodeFaces() const
{
    if( !nodeFacesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph& meshOctreeAddressing::nodeFaces() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calculateNodeFaces();
    }

    return *nodeFacesPtr_;
}

inline const LongList<edge>& meshOctreeAddressing::octreeEdges() const
{
    if( !octreeEdgesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const LongList<edge>&"
                " meshOctreeAddressing::octreeEdges() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        createOctreeEdges();
    }

    return *octreeEdgesPtr_;
}

inline const VRWGraph& meshOctreeAddressing::edgeLeaves() const
{
    if( !edgeLeavesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph&"
                " meshOctreeAddressing::edgeLeaves() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calculateEdgeLeaves();
    }

    return *edgeLeavesPtr_;
}

inline const VRWGraph& meshOctreeAddressing::leafEdges() const
{
    if( !leafEdgesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph& meshOctreeAddressing::leafEdges() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calculateLeafEdges();
    }

    return *leafEdgesPtr_;
}

inline const VRWGraph& meshOctreeAddressing::nodeEdges() const
{
    if( !nodeEdgesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph& meshOctreeAddressing::nodeEdges() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        createOctreeEdges();
    }

    return *nodeEdgesPtr_;
}

inline const VRWGraph& meshOctreeAddressing::faceEdges() const
{
    if( !faceEdgesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph& meshOctreeAddressing::faceEdges() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        createOctreeEdges();
    }

    return *faceEdgesPtr_;
}

inline const VRWGraph& meshOctreeAddressing::edgeFaces() const
{
    if( !edgeFacesPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph& meshOctreeAddressing::edgeFaces() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calculateEdgeFaces();
    }

    return *edgeFacesPtr_;
}

inline const labelLongList& meshOctreeAddressing::globalPointLabel() const
{
    if( !globalPointLabelPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const labelLongList&"
                " meshOctreeAddressing::globalPointLabel() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calcGlobalPointLabels();
    }

    return *globalPointLabelPtr_;
}

inline const Map<label>& meshOctreeAddressing::
globalToLocalPointAddressing() const
{
    if( !globalPointToLocalPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const Map<label>& meshOctreeAddressing::"
                "globalToLocalPointAddressing() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calcGlobalPointLabels();
    }

    return *globalPointToLocalPtr_;
}

inline const VRWGraph& meshOctreeAddressing::pointAtProcs() const
{
    if( !pointProcsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph&"
                " meshOctreeAddressing::pointAtProcs() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calcGlobalPointLabels();
    }

    return *pointProcsPtr_;
}

inline const labelLongList& meshOctreeAddressing::globalFaceLabel() const
{
    if( !globalFaceLabelPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const labelLongList&"
                " meshOctreeAddressing::globalFaceLabel() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calcGlobalFaceLabels();
    }

    return *globalFaceLabelPtr_;
}

inline const Map<label>& meshOctreeAddressing::
globalToLocalFaceAddressing() const
{
    if( !globalFaceToLocalPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const Map<label>& meshOctreeAddressing::"
                "globalToLocalFaceAddressing() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calcGlobalFaceLabels();
    }

    return *globalFaceToLocalPtr_;
}

inline const VRWGraph& meshOctreeAddressing::faceAtProcs() const
{
    if( !faceProcsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph&"
                " meshOctreeAddressing::faceAtProcs() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calcGlobalFaceLabels();
    }

    return *faceProcsPtr_;
}

inline const labelLongList& meshOctreeAddressing::globalLeafLabel() const
{
    if( !globalLeafLabelPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const labelLongList&"
                " meshOctreeAddressing::globalLeafLabel() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calcGlobalLeafLabels();
    }

    return *globalLeafLabelPtr_;
}

inline const VRWGraph& meshOctreeAddressing::leafAtProcs() const
{
    if( !leafAtProcsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const VRWGraph&"
                " meshOctreeAddressing::leafAtProcs() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calcGlobalLeafLabels();
    }

    return *leafAtProcsPtr_;
}

inline const Map<label>& meshOctreeAddressing::
globalToLocalLeafAddressing() const
{
    if( !globalLeafToLocalPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline const Map<label>& meshOctreeAddressing::"
                "globalToLocalLeafAddressing() const"
            ) << "Inside parallel region" << abort(FatalError);
        # endif

        calcGlobalLeafLabels();
    }

    return *globalLeafToLocalPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
