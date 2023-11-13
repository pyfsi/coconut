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
#include "partTriMeshSimplex.H"
#include "triSurfModifier.H"
#include "OFstream.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

partTriMeshSimplex::partTriMeshSimplex
(
    const partTriMesh& tm,
    const label pI
)
:
    pts_(),
    trias_()
{
    const pointField& points = tm.points();
    const LongList<labelledTri>& trias = tm.triangles();
    const VRWGraph& pt = tm.pointTriangles();

    trias_.setSize(pt.sizeOfRow(pI));
    label counter(0);

    Map<label> addr(2*pt.sizeOfRow(pI));
    forAllRow(pt, pI, tI)
    {
        const labelledTri& tri = trias[pt(pI, tI)];
        for(label i=0;i<3;++i)
        {
            const label tpI = tri[i];

            if( !addr.found(tpI) )
            {
                addr.insert(tpI, counter);
                pts_.append(points[tpI]);
                ++counter;
            }
        }

        # ifdef DEBUGSmooth
        Info << "Tet " << tetI << " is " << tet << endl;
        # endif

        label pos(-1);
        for(label i=0;i<3;++i)
            if( tri[i] == pI )
            {
                pos = i;
                break;
            }

        switch( pos )
        {
            case 0:
            {
                trias_[tI] = triFace(addr[tri[0]], addr[tri[1]], addr[tri[2]]);
            } break;
            case 1:
            {
                trias_[tI] = triFace(addr[tri[1]], addr[tri[2]], addr[tri[0]]);
            } break;
            case 2:
            {
                trias_[tI] = triFace(addr[tri[2]], addr[tri[0]], addr[tri[1]]);
            } break;
            default:
            {
                FatalErrorIn
                (
                    "partTriMeshSimplex::partTriMeshSimplex("
                    "(const partTriMesh& tm, const label pI)"
                ) << "Point " << pI << " is not present in triangle" << tri
                    << abort(FatalError);
            }
        }
    }

    # ifdef DEBUGSmooth
    Info << "Simplex at point " << pI << " points " << pts_ << endl;
    Info << "Simplex at point " << pI << " triangles " << trias_ << endl;
    # endif

    if( pts_.size() == 0 || trias_.size() == 0 )
        FatalError << "Simplex at point " << pI << " is not valid "
                   << pts_ << " triangles " << trias_ << abort(FatalError);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

partTriMeshSimplex::~partTriMeshSimplex()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vector partTriMeshSimplex::normal() const
{
    vector normal(vector::zero);
    scalar magN(0.0);

    forAll(trias_, tI)
    {
        const triFace& t = trias_[tI];

        vector n
        (
            0.5 * ((pts_[t[1]] - pts_[t[0]]) ^ (pts_[t[2]] - pts_[t[0]]))
        );
        const scalar magn = mag(n);

        normal += n;
        magN += magn;
    }

    return (normal / (magN + VSMALL));
}


void partTriMeshSimplex::createSurfaceMesh(triSurf& surf) const
{
    triSurfModifier sMod(surf);
    pointField& sPoints = sMod.pointsAccess();
    sPoints.setSize(pts_.size());
    forAll(sPoints, i)
        sPoints[i] = pts_[i];

    LongList<labelledTri>& sTrias = sMod.facetsAccess();
    sTrias.setSize(trias_.size());
    forAll(sTrias, i)
    {
        labelledTri& tf = sTrias[i];
        tf[0] = trias_[i][0];
        tf[1] = trias_[i][1];
        tf[2] = trias_[i][2];

        tf.region() = 0;
    }
    sMod.patchesAccess().setSize(1);
    sMod.patchesAccess()[0].name() = "bnd";
    sMod.patchesAccess()[0].geometricType() = "patch";
}

void partTriMeshSimplex::writeToVTK(const fileName& fName) const
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
    file << "\nCELLS " << trias_.size()
         << " " << 4*trias_.size() << nl;
    forAll(trias_, tI)
    {
        const triFace& t = trias_[tI];

        file << 3 << " " << t[0] << " " << t[1] << " "<< t[2] << nl;
    }

    file << "\nCELL_TYPES " << trias_.size() << nl;
    forAll(trias_, tI)
        file << 5 << nl;

    file << "\n";
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
