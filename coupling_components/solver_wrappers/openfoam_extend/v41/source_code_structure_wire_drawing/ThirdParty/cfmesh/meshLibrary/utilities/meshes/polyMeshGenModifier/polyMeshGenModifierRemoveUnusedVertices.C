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

#include "polyMeshGenModifier.H"
#include "demandDrivenData.H"
#include "labelList.H"

# ifdef USE_OMP
#include <omp.h>
# endif

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::removeUnusedVertices()
{
    faceListPMG& faces = mesh_.faces_;
    pointFieldPMG& points = mesh_.points_;

    boolList usePoint(points.size());
    labelLongList newLabel(points.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(usePoint, pointI)
        {
            usePoint[pointI] = false;
            newLabel[pointI] = -1;
        }

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(faces, faceI)
        {
            const face& f = faces[faceI];

            forAll(f, pI)
                usePoint[f[pI]] = true;
        }

        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            label nPoints(0);
            forAll(points, pI)
                if( usePoint[pI] )
                    newLabel[pI] = nPoints++;

            //- remove unused points from the list
            forAll(newLabel, pI)
                if( (newLabel[pI] != -1) && (newLabel[pI] < pI) )
                {
                    points[newLabel[pI]] = points[pI];
                }

            points.setSize(nPoints);
        }

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(faces, faceI)
        {
            face& f = faces[faceI];

            forAll(f, pI)
                f[pI] = newLabel[f[pI]];
        }
    }

    //- update point subsets
    mesh_.updatePointSubsets(newLabel);

    //- update locked points
    std::set<label> newLockedPoints;
    forAllConstIter(std::set<label>, mesh_.lockedPoints_, it)
    {
        if( newLabel[*it] >= 0 )
            newLockedPoints.insert(newLabel[*it]);
    }
    mesh_.lockedPoints_ = newLockedPoints;
    newLockedPoints.clear();

    //- update point backup
    std::map<label, point> newOrigPoints;
    for
    (
        std::map<label, point>::const_iterator it = mesh_.origPoints_.begin();
        it!=mesh_.origPoints_.end();
        ++it
    )
    {
        const label newPointI = newLabel[it->first];

        if( newPointI >= 0 )
            newOrigPoints[newPointI] = it->second;
    }

    mesh_.origPoints_ = newOrigPoints;
    newOrigPoints.clear();

    //- delete all allocated data
    mesh_.clearOut();
    this->clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
