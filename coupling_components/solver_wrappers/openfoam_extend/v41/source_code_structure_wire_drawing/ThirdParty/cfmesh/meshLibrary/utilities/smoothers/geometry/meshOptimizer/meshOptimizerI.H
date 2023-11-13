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

#include "meshOptimizer.H"
#include "helperFunctions.H"
#include "polyMeshGenAddressing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class labelListType>
inline void meshOptimizer::lockCells(const labelListType& l)
{
    boolList lockedFace(mesh_.faces().size(), false);
    const cellListPMG& cells = mesh_.cells();
    forAll(l, lcI)
    {
        const cell& c = cells[l[lcI]];

        forAll(c, fI)
            lockedFace[c[fI]] = true;
    }

    //- Finally, mark locked points and faces
    const faceListPMG& faces = mesh_.faces();
    forAll(lockedFace, faceI)
    {
        if( lockedFace[faceI] )
        {
            lockedFaces_.append(faceI);

            const face& f = faces[faceI];

            forAll(f, pI)
                vertexLocation_[f[pI]] |= LOCKED;
        }
    }

    //- reduce information accross all processors
    lockPoints(labelLongList());

    # ifdef DEBUGSmoothing
    const label lockedFacesI = mesh_.addFaceSubset("lockedFaces");
    forAll(lockedFaces_, lfI)
        mesh_.addFaceToSubset(lockedFacesI, lockedFaces_[lfI]);

    const label lockPointsI = mesh_.addPointSubset("lockedPoints");
    forAll(vertexLocation_, pointI)
        if( vertexLocation_[pointI] & LOCKED )
            mesh_.addPointToSubset(lockPointsI, pointI);
    # endif
}

template<class labelListType>
inline void meshOptimizer::lockFaces(const labelListType& lf)
{
    boolList lockedFace(mesh_.faces().size(), false);
    forAll(lf, lfI)
    {
        lockedFace[lf[lfI]] = true;
    }

    //- Finally, mark locked points and faces
    const faceListPMG& faces = mesh_.faces();
    forAll(lockedFace, faceI)
    {
        if( lockedFace[faceI] )
        {
            lockedFaces_.append(faceI);

            const face& f = faces[faceI];

            forAll(f, pI)
                vertexLocation_[f[pI]] |= LOCKED;
        }
    }

    //- reduce information accross all processors
    lockPoints(labelLongList());

    # ifdef DEBUGSmoothing
    const label lockedFacesI = mesh_.addFaceSubset("lockedFaces");
    forAll(lockedFaces_, lfI)
        mesh_.addFaceToSubset(lockedFacesI, lockedFaces_[lfI]);

    const label lockPointsI = mesh_.addPointSubset("lockedPoints");
    forAll(vertexLocation_, pointI)
        if( vertexLocation_[pointI] & LOCKED )
            mesh_.addPointToSubset(lockPointsI, pointI);
    # endif
}

template<class labelListType>
inline void meshOptimizer::lockPoints(const labelListType& lp)
{
    forAll(lp, lpI)
        vertexLocation_[lp[lpI]] |= LOCKED;

    if( Pstream::parRun() )
    {
        //- reduce the information accross all processors
        const Map<label>& globalToLocal =
            mesh_.addressingData().globalToLocalPointAddressing();
        const VRWGraph& pointAtProcs = mesh_.addressingData().pointAtProcs();
        const DynList<label>& neiProcs = mesh_.addressingData().pointNeiProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(neiProcs, i)
            exchangeData[neiProcs[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label pointI = it();

            if( vertexLocation_[pointI] & LOCKED )
            {
                forAllRow(pointAtProcs, pointI, i)
                {
                    const label neiProc = pointAtProcs(pointI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append(it.key());
                }
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            vertexLocation_[globalToLocal[receivedData[i]]] |= LOCKED;
        }
    }

    # ifdef DEBUGSmoothing
    const label lockPointsI = mesh_.addPointSubset("lockedPoints");
    forAll(vertexLocation_, pointI)
        if( vertexLocation_[pointI] & LOCKED )
            mesh_.addPointToSubset(lockPointsI, pointI);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
