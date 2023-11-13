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

#include "polyMeshGenPoints.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Time& polyMeshGenPoints::returnTime() const
{
    return runTime_;
}

inline const fileName& polyMeshGenPoints::instance() const
{
    return instance_;
}

inline fileName& polyMeshGenPoints::instance()
{
    return instance_;
}

inline const fileName& polyMeshGenPoints::meshDir() const
{
    return meshDir_;
}

inline fileName& polyMeshGenPoints::meshDir()
{
    return meshDir_;
}

inline const pointFieldPMG& polyMeshGenPoints::points() const
{
    return points_;
}

inline void polyMeshGenPoints::appendVertex(const point& p)
{
    points_.append(p);
}

inline bool polyMeshGenPoints::hasPointsBackup() const
{
    return (origPoints_.size() != 0);
}

inline bool polyMeshGenPoints::hasBackup(const label pointI) const
{
    const std::map<label, point>::const_iterator it = origPoints_.find(pointI);

    return it != origPoints_.end();
}

inline bool polyMeshGenPoints::getOrigPoint(const label pointI, point& p) const
{
    std::map<label, point>::const_iterator it = origPoints_.find(pointI);
    if( it != origPoints_.end() )
    {
        p = it->second;
        return true;
    }

    p = points_[pointI];
    return false;
}

inline const std::set<label>& polyMeshGenPoints::lockedPoints() const
{
    return lockedPoints_;
}

inline bool polyMeshGenPoints::hasLockedPoints() const
{
    return (lockedPoints_.size() != 0);
}

inline bool polyMeshGenPoints::isLockedPoint(const label pointI) const
{
    return (lockedPoints_.find(pointI) != lockedPoints_.end());
}

inline void polyMeshGenPoints::addPointToSubset
(
    const label setI,
    const label pointI
)
{
    std::map<label, meshSubset>::iterator it = pointSubsets_.find(setI);
    if( it == pointSubsets_.end() )
        return;

    it->second.addElement(pointI);
}

inline void polyMeshGenPoints::removePointFromSubset
(
    const label setI,
    const label pointI
)
{
    std::map<label, meshSubset>::iterator it = pointSubsets_.find(setI);
    if( it == pointSubsets_.end() )
        return;

    it->second.removeElement(pointI);
}

inline void polyMeshGenPoints::pointInSubsets
(
    const label pointI,
    DynList<label>& pointSubsets
) const
{
    pointSubsets.clear();

    std::map<label, meshSubset>::const_iterator it;
    for
    (
        it=pointSubsets_.begin();
        it!=pointSubsets_.end();
        ++it
    )
    {
        if( it->second.contains(pointI) )
            pointSubsets.append(it->first);
    }
}

inline void polyMeshGenPoints::pointSubsetIndices(DynList<label>& indices) const
{
    indices.clear();

    std::map<label, meshSubset>::const_iterator it;
    for
    (
        it=pointSubsets_.begin();
        it!=pointSubsets_.end();
        ++it
    )
        indices.append(it->first);
}

template<class ListType>
inline void polyMeshGenPoints::pointsInSubset
(
    const label setI,
    ListType& pointLabels
) const
{
    pointLabels.clear();

    std::map<label, meshSubset>::const_iterator it =
        pointSubsets_.find(setI);
    if( it == pointSubsets_.end() )
        return;

    it->second.containedElements(pointLabels);
}

template<class ListType>
inline void polyMeshGenPoints::updatePointSubsets(const ListType& newNodeLabels)
{
    for
    (
        std::map<label, meshSubset>::iterator it=pointSubsets_.begin();
        it!=pointSubsets_.end();
        ++it
    )
        it->second.updateSubset(newNodeLabels);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
