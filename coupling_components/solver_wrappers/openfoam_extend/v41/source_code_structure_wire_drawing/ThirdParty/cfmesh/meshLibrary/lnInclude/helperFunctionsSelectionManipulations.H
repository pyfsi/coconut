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

InNamespace
    Foam::help

Description
    String conversion functions

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef helperFunctionsSelectionManipulations_H
#define helperFunctionsSelectionManipulations_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "word.H"
#include "labelLongList.H"
#include "scalar.H"

namespace Foam
{

// Forward declarations
class polyMeshGen;

/*---------------------------------------------------------------------------*\
                    Namespace help functions Declaration
\*---------------------------------------------------------------------------*/

namespace help
{

    //- add additional layers for cells into a selection
    label selectAdditionalLayersOfCells
    (
        const polyMeshGen& mesh,
        const label nLayers,
        const labelLongList& activeCells,
        labelLongList& selectedCells
    );

    //- add additional layers fo cells into a subset
    label selectAdditionalLayersOfCells
    (
        polyMeshGen&,
        const label nLayers,
        const label subsetId
    );

    label selectAdditionalLayersOfCells
    (
        polyMeshGen&,
        const label nLayers,
        const word subsetName
    );

    //- create a subset of faces attached to selected points
    label selectedFacesAttachedToPoints
    (
        const polyMeshGen& mesh,
        const labelLongList& activePoints,
        labelLongList& selectedFaces
    );

    label createSubsetOfFacesAttachedToPoints
    (
        polyMeshGen& mesh,
        const labelLongList& activePoints,
        const word subsetName
    );

    //- create a subset of cells attacghed to selected points
    label selectCellsAttachedToPoints
    (
        const polyMeshGen& mesh,
        const labelLongList& activePoints,
        labelLongList& selectedCells
    );

    label createSubsetOfCellsAttachedToPoints
    (
        polyMeshGen& mesh,
        const labelLongList& activePoints,
        const word subsetName
    );

    //- create a subset of cells attached to selected faces
    label selectCellsAttachedToFaces
    (
        const polyMeshGen& mesh,
        const labelLongList& activeFaces,
        labelLongList& selectedCells
    );

    label createSubsetOfCellsAttachedToFaces
    (
        polyMeshGen&,
        const labelLongList& activeFaces,
        const word subsetName
    );

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
