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

#include "polyMeshGen.H"
#include "demandDrivenData.H"
#include "polyMeshGenAddressing.H"
#include "OFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGen::polyMeshGen
(
    const Time& t,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenCells(t, instance, meshDir),
    metaDict_
    (
        IOobject
        (
            "meshMetaDict",
            instance,
            meshDir,
            runTime_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    )
{}

//- Construct from components without the boundary
polyMeshGen::polyMeshGen
(
    const Time& t,
    const pointField& points,
    const faceList& faces,
    const cellList& cells,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenCells(t, points, faces, cells, instance, meshDir),
    metaDict_
    (
        IOobject
        (
            "meshMetaDict",
            instance,
            meshDir,
            runTime_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    )
{}

//- Construct from components with the boundary
polyMeshGen::polyMeshGen
(
    const Time& t,
    const pointField& points,
    const faceList& faces,
    const cellList& cells,
    const wordList& patchNames,
    const labelList& patchStart,
    const labelList& nFacesInPatch,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenCells
    (
        t,
        points,
        faces,
        cells,
        patchNames,
        patchStart,
        nFacesInPatch,
        instance,
        meshDir
    ),
    metaDict_
    (
        IOobject
        (
            "meshMetaDict",
            instance,
            meshDir,
            runTime_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    )
{}

// PC, PDJ - 20-02-2019
//- Construct from components with the boundary and boundary patches
polyMeshGen::polyMeshGen
(
    const Time& t,
    const pointField& points,
    const faceList& faces,
    const cellList& cells,
    const wordList& patchNames,
    const wordList& patchTypes,
    const labelList& patchStart,
    const labelList& nFacesInPatch,
    const labelList& procPatchMyProcNo,
    const labelList& procPatchNeighbProcNo,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenCells
    (
        t,
        points,
        faces,
        cells,
        patchNames,
        patchTypes,
        patchStart,
        nFacesInPatch,
        procPatchMyProcNo,
        procPatchNeighbProcNo,
        instance,
        meshDir
    ),
    metaDict_
    (
        IOobject
        (
            "meshMetaDict",
            instance,
            meshDir,
            runTime_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    )
{}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Destructor
polyMeshGen::~polyMeshGen()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGen::read()
{
    polyMeshGenCells::read();
}

void polyMeshGen::readFromLatestTime()
{
    polyMeshGenCells::readFromLatestTime();
}

void polyMeshGen::write() const
{
    //- remove old mesh before writting
    const fileName meshDir = runTime_.path()/instance_/meshDir_;

    rm(meshDir/"points");
    rm(meshDir/"faces");
    rm(meshDir/"owner");
    rm(meshDir/"neighbour");
    rm(meshDir/"cells");
    rm(meshDir/"boundary");
    rm(meshDir/"pointZones");
    rm(meshDir/"faceZones");
    rm(meshDir/"cellZones");
    rm(meshDir/"meshModifiers");
    rm(meshDir/"parallelData");
    rm(meshDir/"meshMetaDict");
    rm(meshDir/"lockedPoints");
    rm(meshDir/"origPoints");

    // remove sets if they exist
    if (isDir(meshDir/"sets"))
    {
        rmDir(meshDir/"sets");
    }

    //- write the mesh
    polyMeshGenCells::write();

    //- write meta data
    OFstream fName(meshDir/"meshMetaDict");

    metaDict_.writeHeader(fName);
    metaDict_.writeData(fName);
}

void polyMeshGen::printAllocated() const
{
    sizeType totalMemUsage(0);

    totalMemUsage += sizeof(polyMeshGen);

    totalMemUsage += points_.byteSize();

    if( origPoints_.size() )
        totalMemUsage += (sizeof(label) + sizeof(point)) * origPoints_.size();

    if( lockedPoints_.size() )
        totalMemUsage += sizeof(label) * lockedPoints_.size();

    totalMemUsage += faces_.byteSize();

    totalMemUsage += sizeof(nIntFaces_);

    if( ownerPtr_ )
        totalMemUsage += ownerPtr_->byteSize();

    if( neighbourPtr_ )
        totalMemUsage += neighbourPtr_->byteSize();

    totalMemUsage += cells_.byteSize();

    forAll(boundaries_, patchI)
    {
        totalMemUsage += sizeof(boundaries_[patchI]);
        totalMemUsage += boundaries_[patchI].patchName().size();
        totalMemUsage += boundaries_[patchI].patchType().size();
    }

    forAll(procBoundaries_, patchI)
    {
        totalMemUsage += sizeof(procBoundaries_[patchI]);
        totalMemUsage += procBoundaries_[patchI].patchName().size();
        totalMemUsage += procBoundaries_[patchI].patchType().size();
    }

    //- point subsets
    for
    (
        std::map<label, meshSubset>::const_iterator it=pointSubsets_.begin();
        it!=pointSubsets_.end();
        ++it
    )
    {
        totalMemUsage += sizeof(label);
        totalMemUsage += it->second.byteSize();
    }

    //- face subsets
    for
    (
        std::map<label, meshSubset>::const_iterator it=faceSubsets_.begin();
        it!=faceSubsets_.end();
        ++it
    )
    {
        totalMemUsage += sizeof(label);
        totalMemUsage += it->second.byteSize();
    }

    for
    (
        std::map<label, meshSubset>::const_iterator it=cellSubsets_.begin();
        it!=cellSubsets_.end();
        ++it
    )
    {
        totalMemUsage += sizeof(label);
        totalMemUsage += it->second.byteSize();
    }

    if( addressingDataPtr_ )
        totalMemUsage += addressingDataPtr_->printAllocated();

    Pout << "Mesh uses " << totalMemUsage/1048576 << " Mb of memory" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
