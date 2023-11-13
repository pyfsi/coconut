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

#include "polyMeshGenPointZones.H"
#include "zoneIOList.H"
#include "polyMeshGen.H"
#include "OFstream.H"

#include "VRWGraph.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenPointZones::updateSize()
{
    if( points_.size() != pointInZone_.size() )
    {
        const label prevSize = pointInZone_.size();
        pointInZone_.setSize(points_.size());

        for(label pI=prevSize;pI<pointInZone_.size();++pI)
            pointInZone_[pI] = -1;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGenPointZones::polyMeshGenPointZones(const pointFieldPMG& points)
:
    points_(points),
    pointInZone_(),
    nameToIndex_(),
    indexToName_(),
    index_(0)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGenPointZones::~polyMeshGenPointZones()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenPointZones::updatePointZones(const VRWGraph& newLabels)
{
    if( indexToName_.size() == 0 )
        return;

    labelLongList newPointInZone(points_.size(), -1);

    forAll(pointInZone_, pointI)
    {
        forAllRow(newLabels, pointI, i)
            newPointInZone[newLabels(pointI, i)] = pointInZone_[pointI];
    }

    pointInZone_.transfer(newPointInZone);
}

void polyMeshGenPointZones::readPointZones
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
            "pointZones",
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
        pointInZone_.clear();
        nameToIndex_.clear();
        indexToName_.clear();

        return;
    }

    updateSize();

    forAll(readZones, zoneI)
    {
        const word zoneName = readZones[zoneI].name();

        const dictionary& dict = readZones[zoneI].dict();

        const label zoneId = addPointZone(zoneName);

        const labelList indices(dict.lookup("pointLabels"));

        forAll(indices, i)
            pointInZone_[indices[i]] = zoneId;
    }
}

void polyMeshGenPointZones::writePointZones
(
    const Time& runTime,
    const fileName& instance,
    const fileName& meshDir
) const
{
    if( indexToName_.size() == 0 )
        return;

    if( pointInZone_.size() != points_.size() )
        FatalError << "Wrong number of elements in the point zone list"
                   << exit(FatalError);

    zoneIOList pointZones
    (
        IOobject
        (
            "pointZones",
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

        pointZones[zoneI].setName(zoneName);

        labelLongList pointsInZone;
        forAll(pointInZone_, pointI)
            if( pointInZone_[pointI] == it->first )
                pointsInZone.append(pointI);

        pointZones[zoneI].dict().add("type", "pointZone");
        pointZones[zoneI].dict().add("pointLabels", pointsInZone, true);

        ++zoneI;
    }

    pointZones.writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        runTime.writeCompression()
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
