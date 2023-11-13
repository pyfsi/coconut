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

#include "Map.H"
#include "partTetMeshSimplex.H"
#include "OFstream.H"
#include "polyMeshGen.H"
#include "polyMeshGenAddressing.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

partTetMeshSimplex::partTetMeshSimplex
(
    const partTetMesh& tm,
    const label pointI
)
:
    pts_(),
    tets_(),
    bndTriangles_()
{
    tm.createSimplex(pointI, *this);
}

partTetMeshSimplex::partTetMeshSimplex
(
    const DynList<point, 128>& pts,
    const DynList<partTet, 256>& tets,
    const DynList<labelledTri, 32>& bndTriangles,
    const label pointI
)
:
    pts_(pts),
    tets_(tets.size()),
    bndTriangles_(bndTriangles)
{
    forAll(tets, tetI)
    {
        const partTet& tet = tets[tetI];

        const label pos = tet.whichPosition(pointI);

        switch( pos )
        {
            case 0:
            {
                tets_[tetI] =
                    partTet
                    (
                        tet.b(),
                        tet.d(),
                        tet.c(),
                        tet.a()
                    );
            } break;
            case 1:
            {
                tets_[tetI] =
                    partTet
                    (
                        tet.a(),
                        tet.c(),
                        tet.d(),
                        tet.b()
                    );
            } break;
            case 2:
            {
                tets_[tetI] =
                    partTet
                    (
                        tet.a(),
                        tet.d(),
                        tet.b(),
                        tet.c()
                    );
            } break;
            case 3:
            {
                tets_[tetI] =
                    partTet
                    (
                        tet.a(),
                        tet.b(),
                        tet.c(),
                        tet.d()
                    );
            } break;
            default:
            {
                FatalErrorIn
                (
                    "partTetMeshSimplex::partTetMeshSimplex"
                    "(const DynList<point, 128>& pts,"
                    "const DynList<partTet, 128>& tets, const label pointI)"
                ) << "Point " << pointI << " is not present in tet" << tet
                    << abort(FatalError);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

partTetMeshSimplex::~partTetMeshSimplex()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar partTetMeshSimplex::maxEdgeLengthSq() const
{
    scalar dSq(0.0);

    const point& p = centrePoint();

    forAll(tets_, tetI)
    {
        const partTet& tet = tets_[tetI];

        dSq = max(dSq, magSqr(pts_[tet.a()] - p));
        dSq = max(dSq, magSqr(pts_[tet.b()] - p));
        dSq = max(dSq, magSqr(pts_[tet.c()] - p));
    }

    return dSq;
}

scalar partTetMeshSimplex::minEdgeLengthSq() const
{
    scalar dSq(VGREAT);

    const point& p = centrePoint();

    forAll(tets_, tetI)
    {
        const partTet& tet = tets_[tetI];

        dSq = min(dSq, magSqr(pts_[tet.a()] - p));
        dSq = min(dSq, magSqr(pts_[tet.b()] - p));
        dSq = min(dSq, magSqr(pts_[tet.c()] - p));
    }

    return dSq;
}

void partTetMeshSimplex::writeToVTK(const fileName& fName) const
{
    OFstream file(fName);

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    //- write points
    file << "\nPOINTS " << pts_.size() << " float\n";
    forAll(pts_, pI)
    {
        const point& p = pts_[pI];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
    }

    //- write lines
    file << "\nCELLS " << tets_.size()
         << " " << 5*tets_.size() << nl;
    forAll(tets_, tI)
    {
        const partTet& t = tets_[tI];

        file << 4 << " " << t.a() << " " << t.b() << " "
                         << t.c() << " " << t.d() << nl;
    }

    file << "\nCELL_TYPES " << tets_.size() << nl;
    forAll(tets_, tI)
        file << 10 << nl;

    file << "\n";
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
