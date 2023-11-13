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

#include "polyMeshGenCellZones.H"
#include "zoneIOList.H"
#include "dictionary.H"
#include "polyMeshGen.H"

#include "VRWGraph.H"

#include "OFstream.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenCellZones::updateSize()
{
    if( cells_.size() != cellInZone_.size() )
    {
        const label prevSize = cellInZone_.size();
        cellInZone_.setSize(cells_.size());

        for(label cI=prevSize;cI<cellInZone_.size();++cI)
            cellInZone_[cI] = -1;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGenCellZones::polyMeshGenCellZones(const cellListPMG& cells)
:
    cells_(cells),
    cellInZone_(),
    nameToIndex_(),
    indexToName_(),
    index_(0)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Destructor
polyMeshGenCellZones::~polyMeshGenCellZones()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenCellZones::updateCellZones(const VRWGraph& newLabels)
{
    if( indexToName_.size() == 0 )
        return;

    labelLongList newCellInZone(cells_.size(), -1);

    forAll(cellInZone_, cellI)
    {
        forAllRow(newLabels, cellI, i)
            newCellInZone[newLabels(cellI, i)] = cellInZone_[cellI];
    }

    cellInZone_.transfer(newCellInZone);
}

void polyMeshGenCellZones::readCellZones
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
            "cellZones",
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
        cellInZone_.clear();
        nameToIndex_.clear();
        indexToName_.clear();

        return;
    }

    cellInZone_.setSize(cells_.size());
    cellInZone_ = -1;

    forAll(readZones, zoneI)
    {
        const word zoneName = readZones[zoneI].name();

        const dictionary& dict = readZones[zoneI].dict();

        const label zoneId = addCellZone(zoneName);

        const labelList indices(dict.lookup("cellLabels"));

        forAll(indices, i)
            cellInZone_[indices[i]] = zoneId;
    }
}

void polyMeshGenCellZones::writeCellZones
(
    const Time& runTime,
    const fileName& instance,
    const fileName& meshDir
) const
{
    if( indexToName_.size() == 0 )
        return;

    if( cellInZone_.size() != cells_.size() )
    {
        return;
        FatalError << "Wrong number of elements in the cell zone list"
                   << exit(FatalError);
    }

    zoneIOList cellZones
    (
        IOobject
        (
            "cellZones",
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

        cellZones[zoneI].setName(zoneName);

        labelLongList cellsInZone;
        forAll(cellInZone_, cellI)
            if( cellInZone_[cellI] == it->first )
                cellsInZone.append(cellI);

        cellZones[zoneI].dict().add("type", "cellZone");
        cellZones[zoneI].dict().add("cellLabels", cellsInZone, true);

        ++zoneI;
    }

    cellZones.writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        runTime.writeCompression()
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
