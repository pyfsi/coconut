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

#include "polyMeshGenFaceZones.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline label polyMeshGenFaceZones::addFaceZone(const word& name)
{
    updateSize();

    if( nameToIndex_.find(name) != nameToIndex_.end() )
        return nameToIndex_[name];

    nameToIndex_.insert(std::make_pair(name, index_));
    indexToName_.insert(std::make_pair(index_, name));
    ++index_;

    return nameToIndex_[name];
}

inline label polyMeshGenFaceZones::removeFaceZone(const word& name)
{
    if( nameToIndex_.find(name) == nameToIndex_.end() )
        return -1;

    updateSize();

    const label zoneId = nameToIndex_[name];

    indexToName_.erase(zoneId);
    nameToIndex_.erase(name);

    forAll(faceInZone_, faceI)
        if( faceInZone_[faceI] == zoneId )
            faceInZone_[faceI] = -1;

    return zoneId;
}

inline void polyMeshGenFaceZones::faceZoneIndices(DynList<label>& indices) const
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

inline void polyMeshGenFaceZones::faceZoneNames(DynList<word>& names) const
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

inline void polyMeshGenFaceZones::addFaceToZone
(
    const label zoneId,
    const label faceI
)
{
    if( indexToName_.find(zoneId) == indexToName_.end() )
        return;

    updateSize();

    faceInZone_[faceI] = zoneId;
}

inline void polyMeshGenFaceZones::removeFaceFromZone(const label faceI)
{
    updateSize();

    if( (faceI < 0) || (faceI >= faceInZone_.size()) )
        FatalErrorIn
        (
            "inline void polyMeshGenFaceZones::removeFaceFromZone(const label)"
        ) << "Face " << faceI << " is out of scope" << abort(FatalError);

    faceInZone_[faceI] = -1;
}

inline label polyMeshGenFaceZones::faceZoneIndex(const word& name) const
{
    std::map<word, label>::const_iterator it = nameToIndex_.find(name);

    if( it != nameToIndex_.end() )
        return it->second;

    return -1;
}

inline word polyMeshGenFaceZones::faceZoneName(const label zoneId) const
{
    std::map<label, word>::const_iterator it = indexToName_.find(zoneId);

    if( it != indexToName_.end() )
        return it->second;

    WarningIn
    (
        "inline word polyMeshGenFaceZones::faceZoneName(const label) const"
    ) << "Zone " << zoneId << " does not exist" << endl;

    return "";
}

inline bool polyMeshGenFaceZones::isFaceInZone
(
    const word& name,
    const label faceI
) const
{
    std::map<word, label>::const_iterator it = nameToIndex_.find(name);

    if( it == nameToIndex_.end() )
        return false;

    if( faceInZone_[faceI] == it->second )
        return true;

    return false;
}

inline bool polyMeshGenFaceZones::isFaceInZone
(
    const label zoneId,
    const label faceI
) const
{
    std::map<label, word>::const_iterator it = indexToName_.find(zoneId);

    if( it == indexToName_.end() )
        return false;

    if( faceInZone_[faceI] == zoneId )
        return true;

    return false;
}

inline label polyMeshGenFaceZones::faceZone(const label faceI) const
{
    if( indexToName_.size() == 0 )
        return -1;

    return faceInZone_[faceI];
}

inline bool polyMeshGenFaceZones::isFaceFlipped(const label faceI) const
{
    if( indexToName_.size() == 0 )
        return false;

    return flipFace_[faceI];
}

inline void polyMeshGenFaceZones::setFlip(const label faceI, const bool val)
{
    if( indexToName_.size() == 0 )
        return;

    updateSize();

    flipFace_[faceI] = val;
}

template<class ListType>
void polyMeshGenFaceZones::updateFaceZones(const ListType& newLabels)
{
    if( indexToName_.size() == 0 )
        return;

    updateSize();

    labelLongList newFaceInZone(faces_.size(), -1);

    forAll(faceInZone_, faceI)
    {
        const label zoneId = faceInZone_[faceI];
        if( zoneId < 0 )
            continue;
        if( newLabels[faceI] < 0 )
            continue;

        newFaceInZone[newLabels[faceI]] = zoneId;
    }

    faceInZone_.transfer(newFaceInZone);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
