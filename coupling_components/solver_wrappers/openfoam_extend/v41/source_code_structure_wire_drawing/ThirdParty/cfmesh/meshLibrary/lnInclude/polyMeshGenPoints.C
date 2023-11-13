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

#include "polyMeshGenPoints.H"
#include "pointIOField.H"
#include "IOobjectList.H"
#include "pointSet.H"
#include "labelIOList.H"
#include "labelledPoint.H"

namespace Foam
{

defineCompoundTypeName(IOList<labelledPoint>, labelPointIOList);
defineTemplateTypeNameAndDebugWithName
(
    IOList<labelledPoint>,
    "labelledPointIOList",
    0
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGenPoints::polyMeshGenPoints
(
    const Time& runTime,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenPointZones(points_),
    runTime_(runTime),
    instance_(instance),
    meshDir_(meshDir),
    points_
    (
        IOobject
        (
            "points",
            instance_,
            meshDir_,
            runTime_
        ),
        0
    ),
    origPoints_(),
    lockedPoints_(),
    pointSubsets_()
{}

polyMeshGenPoints::polyMeshGenPoints
(
    const Time& runTime,
    const pointField& points,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenPointZones(points_),
    runTime_(runTime),
    instance_(instance),
    meshDir_(meshDir),
    points_
    (
        IOobject
        (
            "points",
            instance_,
            meshDir_,
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        points
    ),
    origPoints_(),
    lockedPoints_(),
    pointSubsets_()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGenPoints::~polyMeshGenPoints()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label polyMeshGenPoints::addPointSubset(const word& subsetName)
{
    label id = pointSubsetIndex(subsetName);
    if( id >= 0 )
    {
        Warning << "Point subset " << subsetName << " already exists!" << endl;
        return id;
    }

    id = 0;
    std::map<label, meshSubset>::const_iterator it;
    for(it=pointSubsets_.begin();it!=pointSubsets_.end();++it)
        id = Foam::max(id, it->first+1);

    pointSubsets_.insert
    (
        std::make_pair
        (
            id,
            meshSubset(subsetName, meshSubset::POINTSUBSET)
        )
    );

    return id;
}

void polyMeshGenPoints::removePointSubset(const label subsetID)
{
    if( pointSubsets_.find(subsetID) == pointSubsets_.end() )
        return;

    pointSubsets_.erase(subsetID);
}

word polyMeshGenPoints::pointSubsetName(const label subsetID) const
{
    std::map<label, meshSubset>::const_iterator it =
        pointSubsets_.find(subsetID);
    if( it == pointSubsets_.end() )
    {
        Warning << "Subset " << subsetID << " is not a point subset" << endl;
        return word();
    }

    return it->second.name();
}

label polyMeshGenPoints::pointSubsetIndex(const word& subsetName) const
{
    std::map<label, meshSubset>::const_iterator it;
    for(it=pointSubsets_.begin();it!=pointSubsets_.end();++it)
    {
        if( it->second.name() == subsetName )
            return it->first;
    }

    return -1;
}

void polyMeshGenPoints::read()
{
    pointIOField pts
    (
        IOobject
        (
            "points",
            instance_,
            meshDir_,
            runTime_,
            IOobject::MUST_READ
        )
    );
    points_.transfer(pts);

    //- find all fields in the polyMesh directory
    IOobjectList meshComponents
    (
        runTime_,
        instance_,
        meshDir_
    );

    if( meshComponents.lookup("origPoints") )
    {
        IOList<labelledPoint> origPoints
        (
            *meshComponents.lookup("origPoints")
        );

        forAll(origPoints, i)
        {
            const labelledPoint& lp = origPoints[i];
            origPoints_[lp.pointLabel()] = lp.coordinates();
        }
    }

    if( meshComponents.lookup("lockedPoints") )
    {
        labelIOList lockedPoints
        (
            *meshComponents.lookup("lockedPoints")
        );

        forAll(lockedPoints, i)
            lockedPoints_.insert(lockedPoints[i]);
    }

    //- read point subsets
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

    wordList setNames = allSets.names("pointSet");
    forAll(setNames, setI)
    {
        IOobject* obj = allSets.lookup(setNames[setI]);

        pointSet pSet(*obj);

        const labelList content = pSet.toc();
        const label id = addPointSubset(setNames[setI]);

        forAll(content, i)
            addPointToSubset(id, content[i]);
    }

    readPointZones(runTime_, instance_, meshDir_);
}

void polyMeshGenPoints::readFromLatestTime()
{
    const instantList instances = runTime_.times();

    fileName origInstance = instance_;

    forAllReverse(instances, instanceI)
    {
        instance_ = instances[instanceI].name() / origInstance;

        IOobject readPoints
        (
            "points",
            instance_,
            meshDir_,
            runTime_,
            IOobject::MUST_READ
        );

        if( !readPoints.headerOk() )
            continue;

        polyMeshGenPoints::read();
        break;
    }

    instance_ = origInstance;
}

void polyMeshGenPoints::write() const
{
    points_.write();

    if( hasPointsBackup() )
    {
        //- write original coordinates of points
        IOList<labelledPoint> origCoordinates
        (
            IOobject
            (
                "origPoints",
                instance_,
                meshDir_,
                runTime_
            ),
            origPoints_.size()
        );

        label i(0);
        for
        (
            std::map<label, point>::const_iterator it=origPoints_.begin();
            it!=origPoints_.end();
            ++it
        )
            origCoordinates[i++] = labelledPoint(it->first, it->second);

        origCoordinates.write();
    }

    if( hasLockedPoints() )
    {
        //- write the indices of locked points
        labelIOList lockedPoints
        (
            IOobject
            (
                "lockedPoints",
                instance_,
                meshDir_,
                runTime_
            ),
            lockedPoints_.size()
        );

        label i(0);
        forAllConstIter(std::set<label>, lockedPoints_, it)
            lockedPoints[i++] = *it;

        lockedPoints.write();
    }

    std::map<label, meshSubset>::const_iterator setIt;
    labelLongList containedElements;

    //- write point selections
    for(setIt=pointSubsets_.begin();setIt!=pointSubsets_.end();++setIt)
    {
        pointSet set
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

        setIt->second.containedElements(containedElements);

        forAll(containedElements, i)
            set.insert(containedElements[i]);

        set.write();
    }

    writePointZones(runTime_, instance_, meshDir_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
