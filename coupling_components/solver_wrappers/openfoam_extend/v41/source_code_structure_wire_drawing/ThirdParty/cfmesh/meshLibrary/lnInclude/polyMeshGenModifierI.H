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

#include "polyMeshGenModifier.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline pointFieldPMG& polyMeshGenModifier::pointsAccess()
{
    return mesh_.points_;
}

inline faceListPMG& polyMeshGenModifier::facesAccess()
{
    return mesh_.faces_;
}

inline cellListPMG& polyMeshGenModifier::cellsAccess()
{
    return mesh_.cells_;
}

inline PtrList<processorBoundaryPatch>&
polyMeshGenModifier::procBoundariesAccess()
{
    return mesh_.procBoundaries_;
}

inline PtrList<boundaryPatch>& polyMeshGenModifier::boundariesAccess()
{
    return mesh_.boundaries_;
}

inline const point& polyMeshGenModifier::movePoint
(
    const label pointI,
    const point& p
)
{
    if( mesh_.lockedPoints_.find(pointI) == mesh_.lockedPoints_.end() )
    {
        mesh_.points_[pointI] = p;
    }
    else
    {
        FatalErrorIn
        (
            "inline const point& polyMeshGenModifier::movePoint"
            "(const label, const point&) const"
        ) << "Point " << pointI
          << " is locked and cannot be moved" << abort(FatalError);
    }

    return mesh_.points_[pointI];
}

inline void polyMeshGenModifier::backupPoints()
{
    Info << "Backing up points" << endl;
    # ifdef USE_OMP
    # pragma omp critical(backupPoints)
    # endif
    {
        mesh_.origPoints_.clear();
        forAll(mesh_.points_, pointI)
            mesh_.origPoints_[pointI] = mesh_.points_[pointI];
    }
}

template<class labelListType>
inline void polyMeshGenModifier::backupPoints(const labelListType& l)
{
    # ifdef USE_OMP
    # pragma omp critical(backupPoints)
    # endif
    {
        forAll(l, i)
        {
            const label pointI = l[i];
            mesh_.origPoints_[pointI] = mesh_.points_[pointI];
        }
    }
}

inline void polyMeshGenModifier::backupPoint(const label pointI)
{
    # ifdef USE_OMP
    # pragma omp critical(backupPoints)
    # endif
    {
        mesh_.origPoints_[pointI] = mesh_.points_[pointI];
    }
}

inline void polyMeshGenModifier::setBackupPoint
(
    const label pointI,
    const point& pOrig
)
{
    # ifdef USE_OMP
    # pragma omp critical(backupPoints)
    # endif
    {
        mesh_.origPoints_[pointI] = pOrig;
    }
}

inline void polyMeshGenModifier::revertPoints()
{
    for
    (
        std::map<label, point>::const_iterator it = mesh_.origPoints_.begin();
        it!=mesh_.origPoints_.end();
        ++it
    )
        mesh_.points_[it->first] = it->second;
}

inline void polyMeshGenModifier::revertPoint(const label pointI)
{
    std::map<label, point>::const_iterator it = mesh_.origPoints_.find(pointI);

    if( it != mesh_.origPoints_.end() )
        mesh_.points_[it->first] = it->second;
}

inline void polyMeshGenModifier::clearPointBackup(const label pointI)
{
    # ifdef USE_OMP
    # pragma omp critical(clearPointBackup)
    # endif
    {
        mesh_.origPoints_.erase(pointI);
    }
}

inline void polyMeshGenModifier::clearPointBackup()
{
    # ifdef USE_OMP
    # pragma omp critical(clearPointBackup)
    # endif
    {
        mesh_.origPoints_.clear();
    }
}

inline void polyMeshGenModifier::lockPoints()
{
    # ifdef USE_OMP
    # pragma omp critical(lockPoints)
    # endif
    {
        mesh_.lockedPoints_.clear();
        forAll(mesh_.points_, pointI)
            mesh_.lockedPoints_.insert(pointI);
    }
}

inline void polyMeshGenModifier::lockPoint(const label pointI)
{
    # ifdef USE_OMP
    # pragma omp critical(lockPoints)
    # endif
    {
        mesh_.lockedPoints_.insert(pointI);
    }
}

inline void polyMeshGenModifier::unlockPoints()
{
    # ifdef USE_OMP
    # pragma omp critical(lockPoints)
    # endif
    {
        mesh_.lockedPoints_.clear();
    }
}

inline void polyMeshGenModifier::unlockPoint(const label pointI)
{
    # ifdef USE_OMP
    # pragma omp critical(lockPoints)
    # endif
    {
        mesh_.lockedPoints_.erase(pointI);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
