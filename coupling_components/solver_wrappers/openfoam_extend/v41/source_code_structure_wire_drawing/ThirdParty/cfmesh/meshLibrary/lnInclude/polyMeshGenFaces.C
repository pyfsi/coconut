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

#include "polyMeshGenFaces.H"
#include "faceIOList.H"
#include "IOPtrList.H"
#include "IOobjectList.H"
#include "faceSet.H"
#include "demandDrivenData.H"
#include "stringListOps.H"

namespace Foam
{

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

void polyMeshGenFaces::clearOut() const
{
    deleteDemandDrivenData(ownerPtr_);
    deleteDemandDrivenData(neighbourPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGenFaces::polyMeshGenFaces
(
    const Time& runTime,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenPoints(runTime, instance, meshDir),
    polyMeshGenFaceZones(faces_),
    faces_
    (
        IOobject
        (
            "faces",
            instance,
            meshDir,
            runTime
        ),
        0
    ),
    procBoundaries_(),
    boundaries_(),
    faceSubsets_(),
    nIntFaces_(0),
    ownerPtr_(NULL),
    neighbourPtr_(NULL)
{}

//- Construct from components without the boundary
polyMeshGenFaces::polyMeshGenFaces
(
    const Time& runTime,
    const pointField& points,
    const faceList& faces,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenPoints(runTime, points, instance, meshDir),
    polyMeshGenFaceZones(faces_),
    faces_
    (
        IOobject
        (
            "faces",
            instance,
            meshDir,
            runTime
        ),
        faces
    ),
    procBoundaries_(),
    boundaries_(),
    faceSubsets_(),
    nIntFaces_(0),
    ownerPtr_(NULL),
    neighbourPtr_(NULL)
{}

//- Construct from components with the boundary
polyMeshGenFaces::polyMeshGenFaces
(
    const Time& runTime,
    const pointField& points,
    const faceList& faces,
    const wordList& patchNames,
    const labelList& patchStart,
    const labelList& nFacesInPatch,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenPoints(runTime, points, instance, meshDir),
    polyMeshGenFaceZones(faces_),
    faces_
    (
        IOobject
        (
            "faces",
            instance,
            meshDir,
            runTime
        ),
        faces
    ),
    procBoundaries_(),
    boundaries_(),
    faceSubsets_(),
    nIntFaces_(0),
    ownerPtr_(NULL),
    neighbourPtr_(NULL)
{
    if( Pstream::parRun() )
        FatalErrorIn
        (
            "polyMeshGenFaces::polyMeshGenFaces("
            "const Time& runTime,"
            "const pointField& points,"
            "const faceList& faces,"
            "const wordList& patchNames,"
            "const labelList& patchStart,"
            "const labelList& nFacesInPatch)"
        ) << "Cannot do this in parallel!" << exit(FatalError);

    boundaries_.setSize(patchNames.size());
    forAll(patchNames, patchI)
    {
        boundaries_.set
        (
            patchI,
            new boundaryPatch
            (
                patchNames[patchI],
                "patch",
                nFacesInPatch[patchI],
                patchStart[patchI]
            )
        );
    }
}


//- Construct from components with the boundary and boundary types
polyMeshGenFaces::polyMeshGenFaces
(
    const Time& runTime,
    const pointField& points,
    const faceList& faces,
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
    polyMeshGenPoints(runTime, points, instance, meshDir),
    polyMeshGenFaceZones(faces_),
    faces_
    (
        IOobject
        (
            "faces",
            instance,
            meshDir,
            runTime
        ),
        faces
    ),
    procBoundaries_(),
    boundaries_(),
    faceSubsets_(),
    nIntFaces_(0),
    ownerPtr_(NULL),
    neighbourPtr_(NULL)
{
    // PC, 13-Feb-19: adapted procedure from 'read' function to construct in
    // parallel

    // Count processor patches
    label i = 0;
    forAll(patchTypes, patchI)
    {
        if (patchTypes[patchI] == "processor")
        {
            ++i;
        }
    }

    // Set size of procBoundaries and boundaries
    procBoundaries_.setSize(i);
    boundaries_.setSize(patchTypes.size() - i);

    // Create non-processor patches
    i = 0;
    forAll(patchTypes, patchI)
    {
        if (patchTypes[patchI] != "processor")
        {
            boundaries_.set
            (
                i,
                new boundaryPatch
                (
                    patchNames[patchI],
                    patchTypes[patchI],
                    nFacesInPatch[patchI],
                    patchStart[patchI]
                )
            );
            ++i;
        }
    }

    // Create processor patches
    i = 0;
    forAll(patchTypes, patchI)
    {
        if (patchTypes[patchI] == "processor")
        {
            procBoundaries_.set
            (
                i,
                new processorBoundaryPatch
                (
                    patchNames[patchI],
                    patchTypes[patchI],
                    nFacesInPatch[patchI],
                    patchStart[patchI],
                    procPatchMyProcNo[patchI],
                    procPatchNeighbProcNo[patchI]
                )
            );
            ++i;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGenFaces::~polyMeshGenFaces()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label polyMeshGenFaces::faceIsInProcPatch(const label faceLabel) const
{
    const label i = procBoundaries_.size() - 1;
    if(
        (i < 0) ||
        (
            faceLabel >=
            (
                procBoundaries_[i].patchStart() +
                procBoundaries_[i].patchSize()
            )
        )
    )
        return -1;

    forAllReverse(procBoundaries_, patchI)
        if( faceLabel >= procBoundaries_[patchI].patchStart() )
            return patchI;

    return -1;
}

label polyMeshGenFaces::faceIsInPatch(const label faceLabel) const
{
    const label i = boundaries_.size() - 1;
    if( faceLabel >= (boundaries_[i].patchStart()+boundaries_[i].patchSize()) )
        return -1;

    forAllReverse(boundaries_, patchI)
        if( faceLabel >= boundaries_[patchI].patchStart() )
            return patchI;

    return -1;
}

wordList polyMeshGenFaces::patchNames() const
{
    wordList t(boundaries_.size());

    forAll(boundaries_, patchI)
    {
        t[patchI] = boundaries_[patchI].patchName();
    }

    return t;
}

label polyMeshGenFaces::getPatchID(const word& patchName) const
{
    forAll(boundaries_, patchI)
    {
        if(boundaries_.set(patchI))
        {
            if(boundaries_[patchI].patchName() == patchName)
            {
                return patchI;
            }
        }
    }

    // If the code gets here, it implies that the patch was not found.
    // return a -1 in this case
    return -1;
}

word polyMeshGenFaces::getPatchName(const label patchID) const
{
    if((patchID < 0) || (patchID >= boundaries_.size()))
    {
         FatalErrorIn
         (
             "polyMeshGenFaces::getPatchName(const label patchID) const"
         )   << "invalid patch ID supplied"
             << abort(FatalError);
    }

    return boundaries_[patchID].patchName();
}

labelList polyMeshGenFaces::findPatches(const word& patchName) const
{
    const wordList allPatches = patchNames();

    const labelList patchIDs = findStrings(patchName, allPatches);

    if(patchIDs.empty())
    {
        WarningIn("polyMeshGenFaces::findPatches(const word&)")
            << "Cannot find any patch names matching " << patchName << endl;
    }

    return patchIDs;
}

label polyMeshGenFaces::addFaceSubset(const word& setName)
{
    label id = faceSubsetIndex(setName);
    if( id >= 0 )
    {
        Warning << "Face subset " << setName << " already exists!" << endl;
        return id;
    }

    id = 0;
    for
    (
        std::map<label, meshSubset>::const_iterator it=faceSubsets_.begin();
        it!=faceSubsets_.end();
        ++it
    )
        id = Foam::max(id, it->first+1);

    faceSubsets_.insert
    (
        std::make_pair
        (
            id,
            meshSubset(setName, meshSubset::FACESUBSET)
        )
    );

    return id;
}

void polyMeshGenFaces::removeFaceSubset(const label setI)
{
    if( faceSubsets_.find(setI) == faceSubsets_.end() )
        return;

    faceSubsets_.erase(setI);
}

word polyMeshGenFaces::faceSubsetName(const label setI) const
{
    std::map<label, meshSubset>::const_iterator it =
        faceSubsets_.find(setI);
    if( it == faceSubsets_.end() )
    {
        Warning << "Subset " << setI << " is not a face subset" << endl;
        return word();
    }

    return it->second.name();
}

label polyMeshGenFaces::faceSubsetIndex(const word& setName) const
{
    std::map<label, meshSubset>::const_iterator it;
    for(it=faceSubsets_.begin();it!=faceSubsets_.end();++it)
    {
        if( it->second.name() == setName )
            return it->first;
    }

    return -1;
}

void polyMeshGenFaces::read()
{
    polyMeshGenPoints::read();

    faceIOList fcs
    (
        IOobject
        (
            "faces",
            instance_,
            meshDir_,
            runTime_,
            IOobject::MUST_READ
        )
    );
    faces_.transfer(fcs);

    deleteDemandDrivenData(ownerPtr_);
    deleteDemandDrivenData(neighbourPtr_);

    ownerPtr_ =
        new labelIOLongList
        (
            IOobject
            (
                "owner",
                instance_,
                meshDir_,
                runTime_,
                IOobject::MUST_READ
            )
        );

    neighbourPtr_ =
        new labelIOLongList
        (
            IOobject
            (
                "neighbour",
                instance_,
                meshDir_,
                runTime_,
                IOobject::MUST_READ
            )
        );

    if( neighbourPtr_->size() != ownerPtr_->size() )
    {
        const label origSize = neighbourPtr_->size();
        neighbourPtr_->setSize(ownerPtr_->size());
        for(label faceI=origSize;faceI<neighbourPtr_->size();++faceI)
            neighbourPtr_->operator[](faceI) = -1;
    }

    //- read boundary information
    IOPtrList<boundaryPatchBase> patches
    (
        IOobject
        (
            "boundary",
            instance_,
            meshDir_,
            runTime_,
            IOobject::MUST_READ
        )
    );

    label i(0);
    forAll(patches, patchI)
        if( patches[patchI].patchType() == "processor" )
            ++i;

    procBoundaries_.setSize(i);
    boundaries_.setSize(patches.size()-i);

    i=0;
    forAll(patches, patchI)
        if( patches[patchI].patchType() != "processor" )
        {
            boundaries_.set
            (
                i,
                new boundaryPatch
                (
                    patches[patchI].patchName(),
                    patches[patchI].dict()
                )
            );
            ++i;
        }

    i = 0;
    forAll(patches, patchI)
        if( patches[patchI].patchType() == "processor" )
        {
            procBoundaries_.set
            (
                i++,
                new processorBoundaryPatch
                (
                    patches[patchI].patchName(),
                    patches[patchI].dict()
                )
            );
        }

    nIntFaces_ = boundaries_[0].patchStart();

    //- read face subsets
    const fileName localDir = meshDir_ / "sets";
    fileName path = instance_.path();
    word name;
    if( path == "." )
    {
        path = instance_.name();
    }
    else
    {
        name = instance_.name();
    }

    IOobjectList allSets
    (
        runTime_,
        path,
        fileName(name/localDir)
    );

    wordList setNames = allSets.names("faceSet");
    forAll(setNames, setI)
    {
        IOobject* obj = allSets.lookup(setNames[setI]);

        faceSet fSet(*obj);
        const labelList content = fSet.toc();
        const label id = addFaceSubset(setNames[setI]);

        forAll(content, i)
            addFaceToSubset(id, content[i]);
    }

    readFaceZones(runTime_, instance_, meshDir_);
}

void polyMeshGenFaces::readFromLatestTime()
{
    polyMeshGenPoints::readFromLatestTime();

    const instantList instances = runTime_.times();

    fileName origInstance = instance_;

    forAllReverse(instances, instanceI)
    {
        instance_ = instances[instanceI].name() / origInstance;

        IOobject readFaces
        (
            "faces",
            instance_,
            meshDir_,
            runTime_,
            IOobject::MUST_READ
        );

        if( !readFaces.headerOk() )
            continue;

        polyMeshGenFaces::read();
        break;
    }

    instance_ = origInstance;
}

void polyMeshGenFaces::write() const
{
    polyMeshGenPoints::write();

    faces_.write();

    if( !ownerPtr_ || !neighbourPtr_ )
        calculateOwnersAndNeighbours();
    ownerPtr_->write();
    neighbourPtr_->write();

    //- write boundary data
    PtrList<boundaryPatchBase> ptchs
    (
        procBoundaries_.size() + boundaries_.size()
    );

    label i(0);

    //- ordinary patches come first
    forAll(boundaries_, patchI)
    {
        ptchs.set
        (
            i++,
            boundaryPatchBase::New
            (
                boundaries_[patchI].patchName(),
                boundaries_[patchI].dict()
            )
        );
    }

    //- processor patches are at the end
    forAll(procBoundaries_, patchI)
    {
        ptchs.set
        (
            i++,
            boundaryPatchBase::New
            (
                procBoundaries_[patchI].patchName(),
                procBoundaries_[patchI].dict()
            )
        );
    }

    IOPtrList<boundaryPatchBase> patches
    (
        IOobject
        (
            "boundary",
            instance_,
            meshDir_,
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ptchs
    );

    patches.write();

    //- write face subsets
    std::map<label, meshSubset>::const_iterator setIt;
    for(setIt=faceSubsets_.begin();setIt!=faceSubsets_.end();++setIt)
    {
        faceSet set
        (
            IOobject
            (
                setIt->second.name(),
                instance_,
                meshDir_/"sets",
                runTime_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        );

        labelLongList containedElements;
        setIt->second.containedElements(containedElements);

        forAll(containedElements, i)
            set.insert(containedElements[i]);

        set.write();
    }

    writeFaceZones(runTime_, instance_, meshDir_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
