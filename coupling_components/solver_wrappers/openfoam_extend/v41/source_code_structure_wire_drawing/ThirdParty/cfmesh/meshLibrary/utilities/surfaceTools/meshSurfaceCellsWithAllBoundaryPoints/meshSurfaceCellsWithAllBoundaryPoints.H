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
    meshSurfaceCellsWithAllBoundaryPoints

Description
    Finds cells with all vertices at the boundary

SourceFiles
    meshSurfaceCellsWithAllBoundaryPoints.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceCellsWithAllBoundaryPoints_H
#define meshSurfaceCellsWithAllBoundaryPoints_H

#include "polyMeshGenModifier.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class meshSurfaceCellsWithAllBoundaryPoints Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceCellsWithAllBoundaryPoints
{
    // Private data
        //- const reference to the mesh
        const polyMeshGen& mesh_;

        //- cells aith all points at the boundary
        labelLongList problematicCells_;

    // Private member functions
        //- mnark cells with all points at the boundary
        void markProblematicCells();

        //- Disallow default bitwise copy construct
        meshSurfaceCellsWithAllBoundaryPoints
        (
            const meshSurfaceCellsWithAllBoundaryPoints&
        );

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceCellsWithAllBoundaryPoints&);

public:

    // Constructors

        //- Construct from mesh
        meshSurfaceCellsWithAllBoundaryPoints(const polyMeshGen& mesh);

    // Destructor

        ~meshSurfaceCellsWithAllBoundaryPoints();

    // Member Functions
        //- mark vertices in regions where all cell vertices are at the surface
        //- this indicates that the cells are too large
        void findProblematicPoints(labelHashSet& problematicPoints);

        //- create backup of problematic points
        void createBackupOfProblematicPoints();

        //- return cells with all vertices at the boundary
        const labelLongList& cellsWithAllPointsAtTheBoundary() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
