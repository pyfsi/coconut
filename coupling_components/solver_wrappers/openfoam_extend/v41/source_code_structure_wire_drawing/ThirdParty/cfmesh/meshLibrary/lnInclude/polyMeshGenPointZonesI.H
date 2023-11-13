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

#include "polyMeshGenPointZones.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline label polyMeshGenPointZones::addPointZone(const word& name)
{
    updateSize();

    if( nameToIndex_.find(name) != nameToIndex_.end() )
        return nameToIndex_[name];

    nameToIndex_.insert(std::make_pair(name, index_));
    indexToName_.insert(std::make_pair(index_, name));
    ++index_;

    return nameToIndex_[name];
}

inline label polyMeshGenPointZones::removePointZone(const word& name)
{
    if( nameToIndex_.find(name) == nameToIndex_.end() )
        return -1;

    updateSize();

    const label zoneId = nameToIndex_[name];

    indexToName_.erase(zoneId);
    nameToIndex_.erase(name);

    forAll(pointInZone_, pointI)
        if( pointInZone_[pointI] == zoneId )
            pointInZone_[pointI] = -1;

    return zoneId;
}

inline void polyMeshGenPointZones::pointZoneIndices
(
    DynList<label>& indices
) const
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

inline void polyMeshGenPointZones::pointZoneNames(DynList<word>& names) const
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

inline void polyMeshGenPointZones::addPointToZone
(
    const label zoneId,
    const label pointI
)
{
    if( indexToName_.find(zoneId) == indexToName_.end() )
        return;

    updateSize();

    pointInZone_[pointI] = zoneId;
}

inline void polyMeshGenPointZones::removePointFromZone(const label pointI)
{
    updateSize();

    if( (pointI < 0) || (pointI >= pointInZone_.size()) )
        FatalErrorIn
        (
            "inline void polyMeshGenPointZones::"
            "removePointFromZone(const label)"
        ) << "Point " << pointI << " is out of scope" << abort(FatalError);

    pointInZone_[pointI] = -1;
}

inline label polyMeshGenPointZones::pointZoneIndex(const word& name) const
{
    std::map<word, label>::const_iterator it = nameToIndex_.find(name);

    if( it != nameToIndex_.end() )
        return it->second;

    return -1;
}

inline word polyMeshGenPointZones::pointZoneName(const label zoneId) const
{
    std::map<label, word>::const_iterator it = indexToName_.find(zoneId);

    if( it != indexToName_.end() )
        return it->second;

    WarningIn
    (
        "inline word polyMeshGenPointZones::pointZoneName(const label) const"
    ) << "Zone " << zoneId << " does not exist" << endl;

    return "";
}

inline bool polyMeshGenPointZones::isPointInZone
(
    const word& name,
    const label pointI
) const
{
    std::map<word, label>::const_iterator it = nameToIndex_.find(name);

    if( it == nameToIndex_.end() )
        return false;

    if( pointInZone_[pointI] == it->second )
        return true;

    return false;
}

inline bool polyMeshGenPointZones::isPointInZone
(
    const label zoneId,
    const label pointI
) const
{
    std::map<label, word>::const_iterator it = indexToName_.find(zoneId);

    if( it == indexToName_.end() )
        return false;

    if( pointInZone_[pointI] == zoneId )
        return true;

    return false;
}

inline label polyMeshGenPointZones::pointZone(const label pointI) const
{
    if( indexToName_.size() == 0 )
        return -1;

    return pointInZone_[pointI];
}

template<class ListType>
void polyMeshGenPointZones::updatePointZones(const ListType& newLabels)
{
    if( indexToName_.size() == 0 )
        return;

    updateSize();

    labelLongList newPointInZone(points_.size(), -1);

    forAll(pointInZone_, pointI)
    {
        const label zoneId = pointInZone_[pointI];
        if( zoneId < 0 )
            continue;
        if( newLabels[pointI] < 0 )
            continue;

        newPointInZone[newLabels[pointI]] = zoneId;
    }

    pointInZone_.transfer(newPointInZone);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
