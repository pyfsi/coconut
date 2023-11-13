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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenCellZones.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline label polyMeshGenCellZones::addCellZone(const word& name)
{
    updateSize();

    if( nameToIndex_.find(name) != nameToIndex_.end() )
        return nameToIndex_[name];

    nameToIndex_.insert(std::make_pair(name, index_));
    indexToName_.insert(std::make_pair(index_, name));
    ++index_;

    return nameToIndex_[name];
}

inline label polyMeshGenCellZones::removeCellZone(const word& name)
{
    if( nameToIndex_.find(name) == nameToIndex_.end() )
        return -1;

    updateSize();

    const label zoneId = nameToIndex_[name];

    indexToName_.erase(zoneId);
    nameToIndex_.erase(name);

    forAll(cellInZone_, cellI)
        if( cellInZone_[cellI] == zoneId )
            cellInZone_[cellI] = -1;

    return zoneId;
}

inline void polyMeshGenCellZones::cellZoneIndices(DynList<label>& indices) const
{
    indices.clear();

    for
    (
        std::map<label, word>::const_iterator it=indexToName_.begin();
        it!=indexToName_.end();
        ++it
    )
        indices.append(it->first);
}

inline void polyMeshGenCellZones::cellZoneNames(DynList<word>& names) const
{
    names.clear();

    for
    (
        std::map<label, word>::const_iterator it=indexToName_.begin();
        it!=indexToName_.end();
        ++it
    )
        names.append(it->second);
}

inline void polyMeshGenCellZones::addCellToZone
(
    const label zoneId,
    const label cellI
)
{
    if( indexToName_.find(zoneId) == indexToName_.end() )
        return;

    updateSize();

    cellInZone_[cellI] = zoneId;
}

inline void polyMeshGenCellZones::removeCellFromZone(const label cellI)
{
    if( (cellI < 0) || (cellI >= cellInZone_.size()) )
        FatalErrorIn
        (
            "inline void polyMeshGenCellZones::removeCellFromZone(const label)"
        ) << "Cell " << cellI << " is out of scope" << abort(FatalError);

    updateSize();

    cellInZone_[cellI] = -1;
}

inline label polyMeshGenCellZones::cellZoneIndex(const word& name) const
{
    std::map<word, label>::const_iterator it = nameToIndex_.find(name);

    if( it != nameToIndex_.end() )
        return it->second;

    return -1;
}

inline word polyMeshGenCellZones::cellZoneName(const label zoneId) const
{
    std::map<label, word>::const_iterator it = indexToName_.find(zoneId);

    if( it != indexToName_.end() )
        return it->second;

    WarningIn
    (
        "inline word polyMeshGenCellZones::cellZoneName(const label) const"
    ) << "Zone " << zoneId << " does not exist" << endl;

    return "";
}

inline bool polyMeshGenCellZones::isCellInZone
(
    const word& name,
    const label cellI
) const
{
    std::map<word, label>::const_iterator it = nameToIndex_.find(name);

    if( it == nameToIndex_.end() )
        return false;

    if( cellInZone_[cellI] == it->second )
        return true;

    return false;
}

inline bool polyMeshGenCellZones::isCellInZone
(
    const label zoneId,
    const label cellI
) const
{
    std::map<label, word>::const_iterator it = indexToName_.find(zoneId);

    if( it == indexToName_.end() )
        return false;

    if( cellInZone_[cellI] == zoneId )
        return true;

    return false;
}

inline label polyMeshGenCellZones::cellZone(const label cellI) const
{
    if( indexToName_.size() == 0 )
        return -1;

    return cellInZone_[cellI];
}

template<class ListType>
void polyMeshGenCellZones::updateCellZones(const ListType& newLabels)
{
    if( indexToName_.size() == 0 )
        return;

    updateSize();

    labelLongList newCellInZone(cells_.size(), -1);

    forAll(cellInZone_, cellI)
    {
        const label zoneId = cellInZone_[cellI];
        if( zoneId < 0 )
            continue;
        if( newLabels[cellI] < 0 )
            continue;

        newCellInZone[newLabels[cellI]] = zoneId;
    }

    cellInZone_.transfer(newCellInZone);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
