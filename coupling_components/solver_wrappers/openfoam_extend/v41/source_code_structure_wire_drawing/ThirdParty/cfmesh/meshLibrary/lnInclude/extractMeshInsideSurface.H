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

Class
    extractMeshInsideSurface

Description
    Extracts edges and corners to cornform the mesh to the geometry

SourceFiles
    extractMeshInsideSurface.C

\*---------------------------------------------------------------------------*/

#ifndef extractMeshInsideSurface_H
#define extractMeshInsideSurface_H

#include "polyMeshGenModifier.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class extractMeshInsideSurface Declaration
\*---------------------------------------------------------------------------*/

class extractMeshInsideSurface
{
    // Private data
        //- mesh
        polyMeshGen& mesh_;

        //- octree
        const meshOctree& octree_;

        //- assing each mesh cells to a region that can be walked without
        //- crossing the cell intersected by the surface mesh
        labelLongList cellInGroup_;

        //- number of cell regions
        label nRegions_;

        //- face candidates to become interfaces
        boolList interfaceFace_;


    // Private member functions
        //- find cells intersected by the surface mesh
        void findIntersectedCells();

        //- detect regions of cells in the mesh
        void assignCellsToRegions();

        //- remove unused cells from the mesh
        void cleanupUnusedCells();

        //- Disallow default bitwise copy construct
        extractMeshInsideSurface(const extractMeshInsideSurface&);

        //- Disallow default bitwise assignment
        void operator=(const extractMeshInsideSurface&);

public:

    // Constructors

        //- Construct from mesh and octree
        extractMeshInsideSurface(polyMeshGen& mesh, const meshOctree& octree);

    // Destructor

        ~extractMeshInsideSurface();

    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
