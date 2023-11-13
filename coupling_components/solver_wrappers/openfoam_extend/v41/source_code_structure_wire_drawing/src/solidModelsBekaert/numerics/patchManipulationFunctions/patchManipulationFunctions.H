/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Functions to clean, re-order and manipulate patches.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef patchManipulationFunctions_H
#define patchManipulationFunctions_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace patchManipulationFunctions
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Function that given two meshes with same number of patches and the
    // patches have the same name, this utility will reorder the patches in the
    //  second mesh to be in the same order as in the first mesh.
    void reorderPatchesConsistently
    (
        const polyMesh& meshSource,
        polyMesh& meshTarget
    );

    // Remove all patches that have no faces
    void removePatchesWithNoFaces(polyMesh& mesh);

    // Remove patch, if it has no faces
    void removePatchWithNoFaces(polyMesh& mesh, const word& patchName);

    // Merge two patches
    // This moves all the faces from patch2 to patch1 and then removes the
    // patch2 patch
    void mergeTwoPatches
    (
        polyMesh& mesh,
        const word& patchName1,
        const word& patchName2
    );

    //- Add all faces in patch to a new face zone
    void addFaceZoneFromPatch
    (
        polyMesh& mesh,
        const word& patchName,
        const word& faceZoneName,
        const bool flipMap
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchManipulationFunctions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
