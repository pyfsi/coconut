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

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::generateCellZonesFromSubsets
(
    const DynList<label>& indices
)
{
    forAll(indices, i)
    {
        const word sName = mesh_.cellSubsetName(indices[i]);

        labelLongList elmtsInSubset;
        mesh_.cellsInSubset(indices[i], elmtsInSubset);

        const label zoneId = mesh_.addCellZone(sName);
        forAll(elmtsInSubset, elI)
            mesh_.addCellToZone(zoneId, elmtsInSubset[elI]);
    }
}

void polyMeshGenModifier::generateCellZonesFromSubsets()
{
    DynList<label> indices;
    mesh_.cellSubsetIndices(indices);

    forAll(indices, i)
    {
        const word sName = mesh_.cellSubsetName(indices[i]);

        labelLongList elmtsInSubset;
        mesh_.cellsInSubset(indices[i], elmtsInSubset);

        const label zoneId = mesh_.addCellZone(sName);
        forAll(elmtsInSubset, elI)
            mesh_.addCellToZone(zoneId, elmtsInSubset[elI]);
    }
}

void polyMeshGenModifier::generateFaceZonesFromSubsets
(
    const DynList<label>& indices
)
{
    forAll(indices, i)
    {
        const word sName = mesh_.faceSubsetName(indices[i]);

        labelLongList elmtsInSubset;
        mesh_.facesInSubset(indices[i], elmtsInSubset);

        const label zoneId = mesh_.addFaceZone(sName);
        forAll(elmtsInSubset, elI)
            mesh_.addFaceToZone(zoneId, elmtsInSubset[elI]);
    }
}

void polyMeshGenModifier::generateFaceZonesFromSubsets()
{
    DynList<label> indices;
    mesh_.faceSubsetIndices(indices);

    forAll(indices, i)
    {
        const word sName = mesh_.faceSubsetName(indices[i]);

        labelLongList elmtsInSubset;
        mesh_.facesInSubset(indices[i], elmtsInSubset);

        const label zoneId = mesh_.addFaceZone(sName);
        forAll(elmtsInSubset, elI)
            mesh_.addFaceToZone(zoneId, elmtsInSubset[elI]);
    }
}

void polyMeshGenModifier::generatePointZonesFromSubsets
(
    const DynList<label>& indices
)
{
    forAll(indices, i)
    {
        const word sName = mesh_.pointSubsetName(indices[i]);

        labelLongList elmtsInSubset;
        mesh_.pointsInSubset(indices[i], elmtsInSubset);

        const label zoneId = mesh_.addPointZone(sName);
        forAll(elmtsInSubset, elI)
            mesh_.addPointToZone(zoneId, elmtsInSubset[elI]);
    }
}

void polyMeshGenModifier::generatePointZonesFromSubsets()
{
    DynList<label> indices;
    mesh_.pointSubsetIndices(indices);

    forAll(indices, i)
    {
        const word sName = mesh_.pointSubsetName(indices[i]);

        labelLongList elmtsInSubset;
        mesh_.pointsInSubset(indices[i], elmtsInSubset);

        const label zoneId = mesh_.addPointZone(sName);
        forAll(elmtsInSubset, elI)
            mesh_.addPointToZone(zoneId, elmtsInSubset[elI]);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
