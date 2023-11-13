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
    partTetMeshSimplex

Description
    A simplex which is used for smoothing purposes

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef partTetMeshSimplex_H
#define partTetMeshSimplex_H

#include "partTetMesh.H"
#include "parPartTet.H"
#include "labelledTri.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class partTetMesh;

/*---------------------------------------------------------------------------*\
                    Class partTetMeshSimplex Declaration
\*---------------------------------------------------------------------------*/

class partTetMeshSimplex
{
    // Private data
        //- points making the simplex
        DynList<point, 128> pts_;

        //- tets making the simplex
        DynList<partTet, 256> tets_;

        //- triangles at th boundary of the mesh
        DynList<labelledTri, 32> bndTriangles_;

public:

    // Friend class
        friend class partTetMesh;

    // Constructors
        //- Construct from partTetMeshSimplex and point label
        partTetMeshSimplex(const partTetMesh& tm, const label pI);

        //- contruct from a list of points and a list of partTets
        //- ad a point label
        partTetMeshSimplex
        (
            const DynList<point, 128>& pts,
            const DynList<partTet, 256>& tets,
            const DynList<labelledTri, 32>& bndTriangles,
            const label pointI
        );

    // Destructor
        ~partTetMeshSimplex();

    // Member functions
        //- return points
        inline DynList<point, 128>& pts()
        {
            return pts_;
        }

        //- return points
        inline const DynList<point, 128>& pts() const
        {
            return pts_;
        }

        //- return tets
        inline const DynList<partTet, 256>& tets() const
        {
            return tets_;
        }

        //- return boundary triangles
        inline const DynList<labelledTri, 32>& boundaryTriangles() const
        {
            return bndTriangles_;
        }

        //- return centre point coordinates
        inline const point& centrePoint() const
        {
            return pts_[tets_[0][3]];
        }

        //- return the square of the maximum edge length
        scalar maxEdgeLengthSq() const;

        //- return the square of the minimum edge length
        scalar minEdgeLengthSq() const;

        //- write simplex to a vtk file
        void writeToVTK(const fileName& fName) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
