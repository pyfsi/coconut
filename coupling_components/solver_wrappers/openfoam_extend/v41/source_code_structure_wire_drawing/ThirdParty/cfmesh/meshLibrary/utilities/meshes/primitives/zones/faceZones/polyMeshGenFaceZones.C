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

#include "polyMeshGenFaceZones.H"
#include "zoneIOList.H"
#include "polyMeshGen.H"

#include "VRWGraph.H"

#include "OFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenFaceZones::updateSize()
{
    if( faces_.size() != faceInZone_.size() )
    {
        const label prevSize = faceInZone_.size();
        faceInZone_.setSize(faces_.size());
        flipFace_.setSize(faces_.size());

        for(label fI=prevSize;fI<faceInZone_.size();++fI)
        {
            faceInZone_[fI] = -1;
            flipFace_[fI] = false;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGenFaceZones::polyMeshGenFaceZones(const faceListPMG& faces)
:
    faces_(faces),
    faceInZone_(),
    flipFace_(),
    nameToIndex_(),
    indexToName_(),
    index_(0)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGenFaceZones::~polyMeshGenFaceZones()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenFaceZones::updateFaceZones(const VRWGraph& newLabels)
{
    if( indexToName_.size() == 0 )
        return;

    labelLongList newFaceInZone(faces_.size(), -1);

    forAll(faceInZone_, faceI)
    {
        forAllRow(newLabels, faceI, i)
            newFaceInZone[newLabels(faceI, i)] = faceInZone_[faceI];
    }

    faceInZone_.transfer(newFaceInZone);
}

void polyMeshGenFaceZones::readFaceZones
(
    const Time& runTime,
    const fileName& instance,
    const fileName& meshDir
)
{
    zoneIOList readZones
    (
        IOobject
        (
            "faceZones",
            instance,
            meshDir,
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    index_ = 0;

    if( readZones.size() == 0 )
    {
        faceInZone_.clear();
        flipFace_.clear();
        nameToIndex_.clear();
        indexToName_.clear();
        flipFace_.clear();

        return;
    }

    updateSize();

    flipFace_.setSize(faces_.size());
    flipFace_ = false;

    forAll(readZones, zoneI)
    {
        const word zoneName = readZones[zoneI].name();

        const dictionary& dict = readZones[zoneI].dict();

        const label zoneId = addFaceZone(zoneName);

        const labelList indices(dict.lookup("faceLabels"));
        const boolList flipMap(dict.lookup("flipMap"));

        forAll(indices, i)
        {
            faceInZone_[indices[i]] = zoneId;
            flipFace_[indices[i]] = flipMap[i];
        }
    }
}

void polyMeshGenFaceZones::writeFaceZones
(
    const Time& runTime,
    const fileName& instance,
    const fileName& meshDir
) const
{
    if( indexToName_.size() == 0 )
        return;

    if( faceInZone_.size() != faces_.size() )
    {
        return;
        FatalError << "Wrong number of elements in the face zone list"
                   << exit(FatalError);
    }

    zoneIOList faceZones
    (
        IOobject
        (
            "faceZones",
            instance,
            meshDir,
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        indexToName_.size()
    );

    label zoneI(0);

    for
    (
        std::map<label, word>::const_iterator it=indexToName_.begin();
        it!=indexToName_.end();
        ++it
    )
    {
        const word zoneName = it->second;

        faceZones[zoneI].setName(zoneName);

        labelLongList facesInZone;
        LongList<bool> flipMap;
        forAll(faceInZone_, faceI)
        {
            if( faceInZone_[faceI] == it->first )
            {
                facesInZone.append(faceI);
                flipMap.append(flipFace_[faceI]);
            }
        }

        faceZones[zoneI].dict().add("type", "faceZone");
        faceZones[zoneI].dict().add("faceLabels", facesInZone, true);
        faceZones[zoneI].dict().add("flipMap", flipMap, true);

        ++zoneI;
    }

    faceZones.writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        runTime.writeCompression()
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
