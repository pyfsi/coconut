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

#include "error.H"
#include "helperFunctionsGeometryQueries.H"
#include "helperFunctionsTopologyManipulation.H"
#include "edgeList.H"
#include "pointField.H"
#include "boolList.H"
#include "triSurf.H"
#include "matrix3D.H"
#include "HashSet.H"
#include "tetrahedron.H"
#include "boundBox.H"
#include "Pair.H"
#include "Map.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace help
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

template<class ListType>
inline bool isnan(const ListType& l)
{
    forAll(l, i)
        if( l[i] != l[i] )
            return true;

    return false;
}

template<class ListType>
bool isinf(const ListType& l)
{
    forAll(l, i)
        if( (l[i] < -VGREAT) || (l[i] > VGREAT) )
            return true;

    return false;
}

inline point interpolate
(
    const scalar u,
    const scalar v,
    const point& P00,
    const point& P10,
    const point& P01,
    const point& P11
)
{
    return point
    (
        (1.0 - u) * (1.0 - v) * P00 +
        u * (1.0 - v) * P10 +
        (1.0 - u) * v * P01 +
        u * v * P11
    );
}

template<class PointField, class Face>
point faceCentre(const PointField& points, const Face& f)
{
    label nPoints = f.size();

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if( nPoints == 3 )
    {
        const point& p0 = points[f[0]];
        const point& p1 = points[f[1]];
        const point& p2 = points[f[2]];

        return (p0 + p1 + p2) / 3.0;
    }

    scalar sumA = 0.0;
    vector sumAc = vector::zero;

    point fCentre = vector::zero;
    forAll(f, pI)
    {
        fCentre += points[f[pI]];
    }

    fCentre /= nPoints;

    forAll(f, pI)
    {
        const point& p = points[f[pI]];
        const point& nextPoint = points[f[(pI+1)%f.size()]];

        vector c = p + nextPoint + fCentre;
        vector n = (nextPoint - p)^(fCentre - p);
        scalar a = mag(n) + VSMALL;

        sumA += a;
        sumAc += a*c;
    }

    return (1.0/3.0) * sumAc / sumA;
}

template<class PointField, class Face>
vector faceAreaVector(const PointField& points, const Face& f)
{
    label nPoints = f.size();

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if( nPoints == 3 )
    {
        const point& p0 = points[f[0]];
        const point& p1 = points[f[1]];
        const point& p2 = points[f[2]];

        return  0.5 * ((p1 - p0) ^ (p2 - p0));
    }

    vector sumN = vector::zero;

    point fCentre = vector::zero;
    forAll(f, pI)
    {
        fCentre += points[f[pI]];
    }

    fCentre /= nPoints;

    forAll(f, pI)
    {
        const point& p = points[f[pI]];
        const point& nextPoint = points[f[(pI+1)%f.size()]];

        vector n = (nextPoint - p) ^ (fCentre - p);

        sumN += n;
    }

    return 0.5 * sumN;
}

template<class PointField, class Face>
scalar magFaceArea(const PointField& points, const Face& f)
{
    return mag(faceAreaVector(points, f));
}

template<class Face1, class Face2>
inline bool isSharedEdgeConvex
(
    const pointField& points,
    const Face1& f1,
    const Face2& f2
)
{
    DynList<label, 3> triOwn(3);
    DynList<label, 3> triNei(3);

    forAll(f1, pI)
    {
        label pos(-1);
        forAll(f2, pJ)
            if( f2[pJ] == f1[pI] )
            {
                pos = pJ;
                break;
            }

        if( pos < 0 )
            continue;

        triNei[0] = f2[pos];
        triNei[1] = f2[f2.fcIndex(pos)];
        triNei[2] = f2[f2.rcIndex(pos)];

        triOwn[0] = f1[pI];
        triOwn[1] = f1[f1.fcIndex(pI)];
        triOwn[2] = f1[f1.rcIndex(pI)];

        scalar vol(0.0);

        forAll(triOwn, pJ)
        {
            if( !triNei.contains(triOwn[pJ]) )
            {
                tetrahedron<point, point> tet
                (
                    points[triNei[0]],
                    points[triNei[1]],
                    points[triNei[2]],
                    points[triOwn[pJ]]
                );

                vol = tet.mag();
                break;
            }
        }

        if( vol > -VSMALL )
            return false;
    }

    return true;
}

template<class Face1, class Face2>
inline scalar angleBetweenFaces
(
    const pointField& points,
    const Face1& f1,
    const Face2& f2
)
{
    DynList<label, 3> triOwn(3);
    DynList<label, 3> triNei(3);

    scalar angle(0.0);
    label counter(0);

    forAll(f1, pI)
    {
        label pos(-1);
        forAll(f2, pJ)
            if( f2[pJ] == f1[pI] )
            {
                pos = pJ;
                break;
            }

        if( pos < 0 )
            continue;

        triNei[0] = f2[pos];
        triNei[1] = f2[f2.fcIndex(pos)];
        triNei[2] = f2[f2.rcIndex(pos)];

        triOwn[0] = f1[pI];
        triOwn[1] = f1[f1.fcIndex(pI)];
        triOwn[2] = f1[f1.rcIndex(pI)];

        scalar vol(0.0);

        forAll(triOwn, pJ)
        {
            if( !triNei.contains(triOwn[pJ]) )
            {
                tetrahedron<point, point> tet
                (
                    points[triNei[0]],
                    points[triNei[1]],
                    points[triNei[2]],
                    points[triOwn[pJ]]
                );

                vol = tet.mag();
                break;
            }
        }

        vector nOwn
        (
            (points[triOwn[1]] - points[triOwn[0]]) ^
            (points[triOwn[2]] - points[triOwn[0]])
        );
        nOwn /= (mag(nOwn) + VSMALL);

        vector nNei
        (
            (points[triNei[1]] - points[triNei[0]]) ^
            (points[triNei[2]] - points[triNei[0]])
        );
        nNei /= (mag(nNei) + VSMALL);

        const scalar dot = Foam::max(-1.0, Foam::min(nOwn & nNei, 1.0));

        if( vol > -VSMALL )
        {
            //- the angle is in the interval [Pi, 2Pi>
            const scalar ang = Foam::acos(dot);

            angle += ang + M_PI;
            ++counter;
        }
        else
        {
            //- the angle is in the interval [0, Pi>
            const scalar ang = Foam::acos(-dot);

            angle += ang;
            ++counter;
        }
    }

    if( counter == 0 )
    {
        FatalErrorIn
        (
            "scalar angleBetweenFaces"
            "(const pointField&, const face&, const face&)"
        ) << "Faces " << f1 << " and " << f2
          << " do no share an edge" << abort(FatalError);
    }

    return angle / counter;
}

inline faceList mergePatchFaces
(
    const List< DynList<label> >& pfcs,
    const pointField& polyPoints
)
{
    //- merge faces which share a common edge
    faceList patchFaces(pfcs.size());
    label counter(0);
    forAll(pfcs, faceI)
        if( pfcs[faceI].size() > 2 )
        {
            const DynList<label>& f = pfcs[faceI];
            face f_(f.size());
            forAll(f_, fJ)
                f_[fJ] = f[fJ];

            patchFaces[counter++] = f_;
        }

    patchFaces.setSize(counter);

    bool merged;
    do
    {
        faceList mergedFaces(patchFaces.size());
        boolList currentlyMerged(patchFaces.size(), false);

        counter = 0;
        merged = false;

        for(label nI=0;nI<(patchFaces.size()-1);nI++)
        {
            vector n0 = patchFaces[nI].normal(polyPoints);
            n0 /= mag(n0);

            for(label nJ=nI+1;nJ<patchFaces.size();nJ++)
            {
                vector n1 = patchFaces[nI].normal(polyPoints);
                n1 /= mag(n1);
                if(
                    help::shareAnEdge(patchFaces[nI], patchFaces[nJ]) &&
                    mag(n0 & n1) > 0.95
                )
                {
                    merged = true;
                    currentlyMerged[nI] = currentlyMerged[nJ] = true;
                    mergedFaces[counter++] =
                        help::mergeTwoFaces
                        (
                            patchFaces[nI],
                            patchFaces[nJ]
                        );

                    break;
                }
            }

            if( merged ) break;
        }

        forAll(patchFaces, pfI)
            if( !currentlyMerged[pfI] )
                mergedFaces.newElmt(counter++) = patchFaces[pfI];

        if( merged )
        {
            patchFaces.setSize(counter);
            for(label k=0;k<counter;k++)
                patchFaces[k] = mergedFaces[k];
        }

    } while( merged );

    return patchFaces;
}

inline bool vertexOnLine(const point& p, const edge& e, const pointField& ep)
{
    vector v = e.vec(ep);
    v /= mag(v);

    vector pv = p - ep[e.start()];
    pv /= mag(pv);

    if( mag(pv & v) > (1.0-SMALL) )
        return true;

    return false;
}

inline bool vertexInPlane(const point& p, const plane& pl)
{
    const vector& n = pl.normal();
    const point& fp = pl.refPoint();

    vector d = p - fp;
    if( mag(d) > VSMALL )
        d /= mag(d);

    if( mag(d & n) < SMALL )
        return true;

    return false;
}

inline bool planeIntersectsEdge
(
    const point& start,
    const point& end,
    const plane& pl,
    point& intersection
)
{
    const vector v = end - start;

    const vector& n = pl.normal();
    const point& fp = pl.refPoint();

    if( mag(n & (v / (mag(v) + VSMALL))) < SMALL )
        return false;

    const scalar t((n & (fp - start)) / (n & v));

    if( (t > -SMALL) && (t < (1.0+SMALL)) )
    {
        intersection = start + v * t;
        return true;
    }

    return false;
}

inline bool pointInTetrahedron
(
    const point& p,
    const tetrahedron<point, point>& tet
)
{
    const vector v0 = tet.a() - tet.d();
    const vector v1 = tet.b() - tet.d();
    const vector v2 = tet.c() - tet.d();
    const vector sp = p - tet.d();

    matrix3D mat;
    FixedList<scalar, 3> source;
    for(label i=0;i<3;++i)
    {
        mat[i][0] = v0[i];
        mat[i][1] = v1[i];
        mat[i][2] = v2[i];
        source[i] = sp[i];
    }

    //- check the determinant of the transformation
    const scalar det = mat.determinant();

    if( mag(det) < VSMALL )
        return false;

    //- get the coordinates of the point in the barycentric corrdinate system
    const scalar u0 = mat.solveFirst(source);

    if( (u0 < -SMALL) || (u0 > (1.0+SMALL)) )
        return false;

    const scalar u1 = mat.solveSecond(source);

    if( (u1 < -SMALL) || ((u0+u1) > (1.0+SMALL)) )
        return false;

    const scalar u2 = mat.solveThird(source);

    if( (u2 < -SMALL) || (u2 > (1.0+SMALL)) )
        return false;


    const scalar u3 = 1.0 - u0 - u1 - u2;

    if( (u3 < -SMALL) || (u3 > (1.0+SMALL)) )
        return false;

    return true;
}

inline bool nearestEdgePointToTheLine
(
    const point& edgePoint0,
    const point& edgePoint1,
    const point& lp0,
    const point& lp1,
    point& nearestOnEdge,
    point& nearestOnLine
)
{
    const vector v = lp1 - lp0;
    const vector d = lp0 - edgePoint0;
    const vector e = edgePoint1 - edgePoint0;

    const scalar vMag = mag(v);
    if( vMag < ROOTVSMALL )
        return false;

    const scalar eMag = mag(e);
    if( eMag < ROOTVSMALL )
    {
        nearestOnEdge = edgePoint0;
        nearestOnLine = nearestPointOnTheEdge(lp0, lp1, nearestOnEdge);
        return true;
    }

    if( mag((v/vMag) & (e/eMag)) > (1.0 - SMALL) )
        return false;

    tensor mat(tensor::zero);
    mat.xx() = (v&v);
    mat.xy() = mat.yx() = -1.0 * (v&e);
    mat.yy() = (e&e);
    mat.zz() = SMALL;

    vector source(vector::zero);
    source[0] = -1.0 * (d&v);
    source[1] = (d&e);

    const vector sol = (inv(mat) & source);

    nearestOnLine = lp0 + v * sol[0];
    if( sol[1] > 1.0 )
    {
        nearestOnEdge = edgePoint1;
    }
    else if( sol[1] < 0.0 )
    {
        nearestOnEdge = edgePoint0;
    }
    else
    {
        nearestOnEdge = edgePoint0 + e * sol[1];
    }

    return true;
}

inline point nearestPointOnTheTriangle
(
    const triangle<point, point>& tri,
    const point& p
)
{
    const vector edge0 = tri.b() - tri.a();

    if( magSqr(edge0) < VSMALL )
        return nearestPointOnTheEdge(tri.a(), tri.c(), p);

    const vector edge1 = tri.c() - tri.a();

    if( magSqr(edge1) < VSMALL )
        return nearestPointOnTheEdge(tri.a(), tri.b(), p);

    if( magSqr(tri.c() - tri.b()) < VSMALL )
        return nearestPointOnTheEdgeExact(tri.a(), tri.c(), p);

    const vector v0 = tri.a() - p;

    const scalar a = (edge0 & edge0) + VSMALL;
    const scalar b = edge0 & edge1;
    const scalar c = (edge1 & edge1) + VSMALL;
    const scalar d = edge0 & v0;
    const scalar e = edge1 & v0;

    const scalar det = mag(a * c - b * b);

    scalar u = b * e - c * d;
    scalar v = b * d - a * e;

    if( u + v <= det )
    {
        if( u < 0.0 )
        {
            if( v < 0.0 )
            {
                //- region 4
                if( d < 0.0 )
                {
                    u = Foam::min(Foam::max(-d/a, 0.0), 1.0);
                    v = 0.0;
                }
                else
                {
                    u = 0.0;
                    v = Foam::min(Foam::max(-e/c, 0.0), 1.0);
                }
            }
            else
            {
                //- region 3
                u = 0.0;
                v = Foam::min(Foam::max(-e/c, 0.0), 1.0);
            }
        }
        else if( v < 0.0 )
        {
            //- region 5
            u = Foam::min(Foam::max(-d/a, 0.0), 1.0);
            v = 0.0;
        }
        else if( mag(det) > VSMALL )
        {
            //- region 0
            const scalar invDet = 1.0 / det;
            u *= invDet;
            v *= invDet;
        }
        else
        {
            //- the area of the triangle is zero
            scalar dSqMax(0.0);
            point s, e;

            scalar dSq = Foam::magSqr(tri.b() - tri.a());
            s = tri.a();
            e = tri.b();
            dSqMax = dSq;

            dSq = Foam::magSqr(tri.c() - tri.b());
            if( dSq > dSqMax )
            {
                s = tri.b();
                e = tri.c();
                dSqMax = dSq;
            }

            dSq = Foam::magSqr(tri.a() - tri.c());
            if( dSq > dSqMax )
            {
                s = tri.c();
                e = tri.a();
                dSqMax = dSq;
            }

            return nearestPointOnTheEdgeExact(s, e, p);
        }
    }
    else
    {
        if( u < 0.0 )
        {
            //- region 2
            const scalar tmp0 = b + d;
            const scalar tmp1 = c + e;

            if( tmp1 > tmp0 )
            {
                const scalar numer = tmp1 - tmp0;
                const scalar denom = a - 2.0 * b + c + VSMALL;
                u = Foam::min(Foam::max(numer/denom, 0.0), 1.0);
                v = 1.0 - u;
            }
            else
            {
                u = 0.0;
                v = Foam::min(Foam::max(-e/c, 0.0), 1.0);
            }
        }
        else if( v < 0.0 )
        {
            //- region 6
            if( a+d > b+e )
            {
                const scalar numer = c + e - b - d;
                const scalar denom = a - 2.0 * b + c + VSMALL;

                u = Foam::min(Foam::max(numer/denom, 0.0), 1.0);
                v = 1.0 - u;
            }
            else
            {
                u = Foam::min(Foam::max(-d/a, 0.0), 1.0);
                v = 0.0;
            }
        }
        else
        {
            //- region 1
            const scalar numer = c + e - b - d;
            const scalar denom = a - 2.0 * b + c + VSMALL;

            u = Foam::min(Foam::max(numer/denom, 0.0), 1.0);
            v = 1.0 - u;
        }
    }

    return tri.a() + u * edge0 + v * edge1;
}

inline point nearestPointOnTheTriangle
(
    const label tI,
    const triSurf& surface,
    const point& p
)
{
    const triangle<point, point> tri = surfaceTriangle(surface, tI);

    return nearestPointOnTheTriangle(tri, p);
}

inline void nearestPointsOnTriangles
(
    const triangle<point, point>& tri0,
    const triangle<point, point>& tri1,
    point& nearestOnTri0,
    point& nearestOnTri1
)
{
    scalar dSqMin(VGREAT);

    //- find the distance of the points in tri1 from tri0
    point np = nearestPointOnTheTriangle(tri0, tri1.a());
    scalar dSq = magSqr(np - tri1.a());
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearestOnTri0 = np;
        nearestOnTri1 = tri1.a();
    }

    np = nearestPointOnTheTriangle(tri0, tri1.b());
    dSq = magSqr(np - tri1.b());
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearestOnTri0 = np;
        nearestOnTri1 = tri1.b();
    }

    np = nearestPointOnTheTriangle(tri0, tri1.c());
    dSq = magSqr(np - tri1.c());
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearestOnTri0 = np;
        nearestOnTri1 = tri1.c();
    }

    //- find the distance of the points in tri0 from tri1
    np = nearestPointOnTheTriangle(tri1, tri0.a());
    dSq = magSqr(np - tri0.a());
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearestOnTri0 = tri0.a();
        nearestOnTri1 = np;
    }

    np = nearestPointOnTheTriangle(tri1, tri0.b());
    dSq = magSqr(np - tri0.b());
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearestOnTri0 = tri0.b();
        nearestOnTri1 = np;
    }

    np = nearestPointOnTheTriangle(tri1, tri0.c());
    dSq = magSqr(np - tri0.c());
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearestOnTri0 = tri0.c();
        nearestOnTri1 = np;
    }

    //- find the nearest distance between the edges in tri0 from tri1
    point noe, notr;
    nearestPointsBetweenLineAndTriangle(tri1, tri0.a(), tri0.b(), notr, noe);
    dSq = magSqr(noe - notr);
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearestOnTri0 = noe;
        nearestOnTri1 = notr;
    }

    nearestPointsBetweenLineAndTriangle(tri1, tri0.b(), tri0.c(), notr, noe);
    dSq = magSqr(noe - notr);
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearestOnTri0 = noe;
        nearestOnTri1 = notr;
    }

    nearestPointsBetweenLineAndTriangle(tri1, tri0.c(), tri0.a(), notr, noe);
    dSq = magSqr(noe - notr);
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearestOnTri0 = noe;
        nearestOnTri1 = notr;
    }

    //- find the nearest distance between the edges in tri1 from tri0
    nearestPointsBetweenLineAndTriangle(tri0, tri1.a(), tri1.b(), notr, noe);
    dSq = magSqr(noe - notr);
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearestOnTri0 = notr;
        nearestOnTri1 = noe;
    }

    nearestPointsBetweenLineAndTriangle(tri0, tri1.b(), tri1.c(), notr, noe);
    dSq = magSqr(noe - notr);
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearestOnTri0 = notr;
        nearestOnTri1 = noe;
    }

    nearestPointsBetweenLineAndTriangle(tri0, tri1.c(), tri1.a(), notr, noe);
    dSq = magSqr(noe - notr);
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearestOnTri0 = notr;
        nearestOnTri1 = noe;
    }
}

inline void nearestPointsOnEdges
(
    const point& e0s,
    const point& e0e,
    const point& e1s,
    const point& e1e,
    point& nearest0,
    point& nearest1
)
{
    scalar dSqMin(VGREAT);
    point np = nearestPointOnTheEdgeExact(e1s, e1e, e0s);
    scalar dSq = magSqr(np - e0s);
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearest0 = e0s;
        nearest1 = np;
    }

    np = nearestPointOnTheEdgeExact(e1s, e1e, e0e);
    dSq = magSqr(np - e0e);
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearest0 = e0e;
        nearest1 = np;
    }

    np = nearestPointOnTheEdgeExact(e0s, e0e, e1s);
    dSq = magSqr(np - e1s);
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearest0 = np;
        nearest1 = e1s;
    }

    np = nearestPointOnTheEdgeExact(e0s, e0e, e1e);
    dSq = magSqr(np - e1e);
    if( dSq < dSqMin )
    {
        dSqMin = dSq;
        nearest0 = np;
        nearest1 = e1e;
    }

    point nl(np);
    nearestEdgePointToTheLine(e0s, e0e, e1s, e1e, np, nl);

    if( magSqr(np - nl) < dSqMin )
    {
        nearest0 = np;
        nearest1 = nl;
    }
}

inline bool findMinimizerPoint
(
    const DynList<point>& origins,
    const DynList<vector>& normals,
    point& pMin
)
{
    if( origins.size() != normals.size() )
        FatalErrorIn
        (
            "inline bool findMinimizerPoint"
            "(const DynList<point>&, const DynList<vector>&, point&)"
        ) << "Size of normals and origins do not match" << abort(FatalError);

    tensor mat(tensor::zero);
    vector source(vector::zero);

    forAll(origins, i)
    {
        //- use the normalized vector
        vector n = normals[i];
        n /= (mag(n) + VSMALL);

        const tensor t = n * n;

        mat += t;

        source += (origins[i] & n) * n;
    }

    const scalar determinant = Foam::mag(Foam::det(mat));
    if( determinant < SMALL )
        return false;

    pMin = (inv(mat, determinant) & source);

    return true;
}

inline bool triLineIntersection
(
    const triangle<point, point>& tria,
    const point& lineStart,
    const point& lineEnd,
    point& intersection
)
{
    const point& p0 = tria.a();
    const vector v(lineStart - lineEnd);
    const vector v0 = tria.b() - p0;
    const vector v1 = tria.c() - p0;
    const vector sp = lineStart - p0;

    matrix3D mat;
    FixedList<scalar, 3> source;
    for(label i=0;i<3;++i)
    {
        mat[i][0] = v0[i];
        mat[i][1] = v1[i];
        mat[i][2] = v[i];
        source[i] = sp[i];
    }

    const scalar det = mat.determinant();

    if( mag(det) < SMALL )
        return false;

    const scalar t = mat.solveThird(source);

    if( (t < -SMALL) || (t > (1.0+SMALL)) )
        return false;

    const scalar u0 = mat.solveFirst(source);

    if( u0 < -SMALL )
        return false;

    const scalar u1 = mat.solveSecond(source);

    if( (u1 < -SMALL) || ((u0+u1) > (1.0+SMALL)) )
        return false;

    intersection = lineStart - t * v;
    return true;
}

inline bool triLineIntersection
(
    const triSurf& surface,
    const label tI,
    const point& s,
    const point& e,
    point& intersection
)
{
    const triangle<point, point> tria = surfaceTriangle(surface, tI);

    return triLineIntersection(tria, s, e, intersection);
}

inline bool nearestPointsBetweenLineAndTriangle
(
    const triangle<point, point>& tria,
    const point& lineStart,
    const point& lineEnd,
    point& nearestOnTriangle,
    point& nearestOnLine
)
{
    if( magSqr(lineStart - lineEnd) < VSMALL )
    {
        //- points of the line are very close to each other
        nearestOnTriangle = help::nearestPointOnTheTriangle(tria, lineStart);
        nearestOnLine = nearestOnTriangle;
        return true;
    }
    else if( tria.mag() < ROOTVSMALL )
    {
        //- triangle has a very small area it is a line or a point
        scalar maxDist(0.0);
        point start, end;

        scalar dSq = magSqr(tria.a() - tria.b());
        start = tria.a();
        end = tria.b();
        maxDist = dSq;

        dSq = magSqr(tria.b() - tria.c());
        if( dSq > maxDist )
        {
            start = tria.b();
            end = tria.c();
            maxDist = dSq;
        }

        dSq = magSqr(tria.c() - tria.a());
        if( dSq > maxDist )
        {
            start = tria.c();
            end = tria.a();
            maxDist = dSq;
        }

        if( maxDist > VSMALL )
        {
            nearestEdgePointToTheLine
            (
                lineStart, lineEnd,
                start, end,
                nearestOnLine,
                nearestOnTriangle
            );

            return true;
        }
        else
        {
            //- triangle is a point. Find the distance from the line
            nearestOnLine =
                help::nearestPointOnTheEdgeExact(lineStart, lineEnd, start);
            nearestOnTriangle = start;

            return true;
        }
    }
    else if( triLineIntersection(tria, lineStart, lineEnd, nearestOnTriangle) )
    {
        //- the line and the triangle intersect each other
        nearestOnLine = nearestOnTriangle;
        return true;
    }
    else
    {
        point nol(lineStart), notr(tria.a());

        //- check the distance from the edge a - b
        nearestEdgePointToTheLine
        (
            lineStart, lineEnd,
            tria.a(), tria.b(),
            nol, notr
        );

        scalar dSq = magSqr(nol - notr);
        scalar dSqMin = dSq;
        nearestOnLine = nol;
        nearestOnTriangle = notr;

        //- check the distance from the edge b - c
        nearestEdgePointToTheLine
        (
            lineStart, lineEnd,
            tria.b(), tria.c(),
            nol, notr
        );

        dSq = magSqr(nol - notr);
        if( dSq < dSqMin )
        {
            dSqMin = dSq;
            nearestOnLine = nol;
            nearestOnTriangle = notr;
        }

        //- check the distance from the edge c - a
        nearestEdgePointToTheLine
        (
            lineStart, lineEnd,
            tria.c(), tria.a(),
            nol, notr
        );

        dSq = magSqr(nol - notr);
        if( dSq < dSqMin )
        {
            dSqMin = dSq;
            nearestOnLine = nol;
            nearestOnTriangle = notr;
        }

        return true;
    }

    return false;
}

inline triangle<point, point> surfaceTriangle
(
    const triSurf& surf,
    const label tI
)
{
    const pointField& pts = surf.points();
    const labelledTri& ltri = surf[tI];

    return triangle<point, point>
    (
        pts[ltri[0]],
        pts[ltri[1]],
        pts[ltri[2]]
    );
}

inline triangle<point, point> projectTriangleOntoPlane
(
    const triangle<point, point>& tri,
    const plane& pl
)
{
    const triangle<point, point> retTri
    (
        pl.nearestPoint(tri.a()),
        pl.nearestPoint(tri.b()),
        pl.nearestPoint(tri.c())
    );

    return retTri;
}

inline label pointsOnTheSameSide
(
    const triangle<point, point>& tri,
    const plane& pl
)
{
    FixedList<scalar, 3> dist;
    dist[0] = pl.distance(tri.a());
    dist[1] = pl.distance(tri.b());
    dist[2] = pl.distance(tri.c());

    //- compare tri2 to tri1
    bool hasZero(false);

    bool hasPositive(false), hasNegative(false);
    forAll(dist, pI)
    {
        if( dist[pI] > VSMALL )
        {
            hasPositive = true;
        }
        else if( dist[pI] < -VSMALL )
        {
            hasNegative = true;
        }
        else
        {
            hasZero = true;
        }
    }

    //- check if both triangles are in the same plane
    if( hasZero && !hasPositive && !hasNegative )
        return 2;

    if( hasPositive && !hasNegative )
        return 1;
    if( hasNegative && !hasPositive )
        return -1;

    return 0;
}

inline label pointsOnTheSameSide
(
    const triSurf& surf,
    const label tI,
    const plane& pl
)
{
    const triangle<point, point> tri = surfaceTriangle(surf, tI);

    return pointsOnTheSameSide(tri, pl);
}

inline label pointsOnTheSameSide
(
    const triSurf& surf,
    const label tri1,
    const label tri2
)
{
    const pointField& pts = surf.points();

    vector n1 = surf[tri1].normal(pts);
    const scalar magn1 = mag(n1);

    if( magn1 > VSMALL )
    {
        n1 /= magn1;

        plane pl(pts[surf[tri1][0]], n1);

        label side = pointsOnTheSameSide(surf, tri2, pl);

        if( side == 1 )
        {
            return 1;
        }
        else if( side == -1 )
        {
            return -1;
        }
        else if( side == 2 )
        {
            return 3;
        }
    }

    vector n2 = surf[tri2].normal(pts);
    const scalar magn2 = mag(n2);

    if( magn2 > VSMALL )
    {
        n2 /= magn2;

        plane pl(pts[surf[tri2][0]], n2);

        label side = pointsOnTheSameSide(surf, tri1, pl);

        if( side == 1 )
        {
            return 2;
        }
        else if( side == -1 )
        {
            return -2;
        }
    }

    return 0;
}

inline scalar minAngleBetweenTriangles
(
    const triangle<point, point>& tri0,
    const triangle<point, point>& tri1,
    bool& isIntersection,
    const scalar otherSignTol
)
{
    isIntersection = false;

    vector n0 = tri0.normal();
    const scalar magn0 = mag(n0) + VSMALL;
    n0 /= magn0;

    vector n1 = tri1.normal();
    const scalar magn1 = mag(n1) + VSMALL;
    n1 /= magn1;

    //- use the triangle with the greater area as a reference plane
    scalar angle = Foam::acos(max(0.0, min(mag(n0 & n1), 1.0)));

    vector n = n0 ^ n1;

    DynList<vector> normals(3);
    normals[0] = n0;
    normals[1] = n1;
    normals[2] = n;

    DynList<point> origins(3);
    origins[0] = tri0.a();
    origins[1] = tri1.a();
    origins[2] = tri0.a();

    point pOrigin;
    if( !findMinimizerPoint(origins, normals, pOrigin) )
    {
        //- there exists no intersection
        return M_PI;
    }

    //- find the max distance point in the triangle 0 from the intersection line
    point maxP0 = tri0.a();
    point maxProjP0 = nearestPointOnLine(pOrigin, n, tri0.a());
    scalar maxDistP0 = magSqr(tri0.a() - maxProjP0);

    point pProj = nearestPointOnLine(pOrigin, n, tri0.b());
    scalar d = magSqr(tri0.b() - pProj);
    if( d > maxDistP0 )
    {
        maxP0 = tri0.b();
        maxProjP0 = pProj;
        maxDistP0 = d;
    }

    pProj = nearestPointOnLine(pOrigin, n, tri0.c());
    d = magSqr(tri0.c() - pProj);
    if( d > maxDistP0 )
    {
        maxP0 = tri0.c();
        maxProjP0 = pProj;
        maxDistP0 = d;
    }

    //- find the max distance point in the triangle 1 from the intersection line
    point maxP1 = tri1.a();
    point maxProjP1 = nearestPointOnLine(pOrigin, n, tri1.a());
    scalar maxDistP1 = magSqr(tri1.a() - maxProjP1);

    pProj = nearestPointOnLine(pOrigin, n, tri1.b());
    d = magSqr(tri1.b() - pProj);
    if( d > maxDistP1 )
    {
        maxP1 = tri1.b();
        maxProjP1 = pProj;
        maxDistP1 = d;
    }

    pProj = nearestPointOnLine(pOrigin, n, tri1.c());
    d = magSqr(tri1.c() - pProj);
    if( d > maxDistP1 )
    {
        maxP1 = tri1.c();
        maxProjP1 = pProj;
        maxDistP1 = d;
    }

    //- find min and max distance of triangle 0 from the plane of triangle 1
    bool intersection01(false);
    scalar maxDist(-VGREAT), minDist(VGREAT);
    d = (tri0.a() - tri1.a()) & n1;
    maxDist = max(d, maxDist);
    minDist = min(d, minDist);

    d = (tri0.b() - tri1.a()) & n1;
    maxDist = max(d, maxDist);
    minDist = min(d, minDist);

    d = (tri0.c() - tri1.a()) & n1;
    maxDist = max(d, maxDist);
    minDist = min(d, minDist);

    if( minDist < -VSMALL && maxDist > VSMALL )
    {
        const scalar dMax = max(mag(minDist), maxDist);
        const scalar dMin = min(mag(minDist), maxDist);

        if( dMin > otherSignTol * dMax )
        {
            isIntersection = true;
            intersection01 = true;
        }
    }

    //- find min and max distance of triangle 1 from the plane of triangle 0
    bool intersection10(false);
    maxDist = -VGREAT;
    minDist = VGREAT;
    d = (tri1.a() - tri0.a()) & n0;
    maxDist = max(d, maxDist);
    minDist = min(d, minDist);

    d = (tri1.b() - tri0.a()) & n0;
    maxDist = max(d, maxDist);
    minDist = min(d, minDist);

    d = (tri1.c() - tri0.a()) & n0;
    maxDist = max(d, maxDist);
    minDist = min(d, minDist);

    if( minDist < -VSMALL && maxDist > VSMALL )
    {
        const scalar dMax = max(mag(minDist), maxDist);
        const scalar dMin = min(mag(minDist), maxDist);

        if( dMin > otherSignTol * dMax )
        {
            isIntersection = true;
            intersection10 = true;
        }
    }

    //- check if the intersection occurs
    if( intersection01 && intersection10 )
        return angle;

    //- check the dot product between the distance vectors
    if( ((maxP1 - maxProjP1) & (maxP0 - maxProjP0)) < 0.0 )
        return (M_PI - angle);

    return angle;
}

inline scalar minAngleBetweenTriangles
(
    const triSurf& surf,
    const label t0,
    const label t1,
    bool& isIntersection,
    const scalar otherSignTol
)
{
    const pointField& pts = surf.points();

    const triangle<point, point> tri0
    (
        pts[surf[t0][0]],
        pts[surf[t0][1]],
        pts[surf[t0][2]]
    );

    const triangle<point, point> tri1
    (
        pts[surf[t1][0]],
        pts[surf[t1][1]],
        pts[surf[t1][2]]
    );

    return minAngleBetweenTriangles(tri0, tri1, isIntersection, otherSignTol);
}

inline bool boundBoxLineIntersection
(
    const point& s,
    const point& e,
    const boundBox& bb
)
{
    scalar tMax(1.0+SMALL), tMin(-SMALL);

    const vector v = e - s;
    const scalar d = mag(v);

    //- check if the vector has length
    if( d < VSMALL )
    {
        if( bb.contains(s) )
            return true;

        return false;
    }

    const point& pMin = bb.min();
    const point& pMax = bb.max();

    //- check coordinates
    for(label dir=0;dir<3;++dir)
    {
        const scalar vd = v[dir];
        const scalar sd = s[dir];

        if( mag(vd) > (SMALL * d) )
        {
            if( vd >= 0.0 )
            {
                tMin = Foam::max(tMin, (pMin[dir] - sd) / vd);
                tMax = Foam::min(tMax, (pMax[dir] - sd) / vd);
            }
            else
            {
                tMin = Foam::max(tMin, (pMax[dir] - sd) / vd);
                tMax = Foam::min(tMax, (pMin[dir] - sd) / vd);
            }
        }
        else if( (sd < pMin[dir]) || (sd > pMax[dir]) )
        {
            return false;
        }
    }

    if( (tMax - tMin) > -SMALL )
        return true;

    return false;
}

inline bool lineFaceIntersection
(
    const point& sp,
    const point& ep,
    const face& f,
    const pointField& fp,
    point& intersection
)
{
    const point c = faceCentre(fp, f);

    forAll(f, pI)
    {
        const triangle<point, point> tri
        (
            fp[f[pI]],
            fp[f.nextLabel(pI)],
            c
        );

        if( triLineIntersection(tri, sp, ep, intersection) )
            return true;
    }

    return false;
}

inline bool doFaceAndTriangleIntersect
(
    const triSurf& surface,
    const label triI,
    const face& f,
    const pointField& facePoints
)
{
    const pointField& triPoints = surface.points();

    const point centre = faceCentre(facePoints, f);
    point intersection;

    //- check if any triangle edge intersects the face
    const labelledTri& tri = surface[triI];

    forAll(tri, eI)
    {
        const point& s = triPoints[tri[eI]];
        const point& e = triPoints[tri[(eI+1)%3]];

        forAll(f, pI)
        {
            const triangle<point, point> tria
            (
                facePoints[f[pI]],
                facePoints[f.nextLabel(pI)],
                centre
            );

            const bool currIntersection =
                help::triLineIntersection
                (
                    tria,
                    s,
                    e,
                    intersection
                );

            if( currIntersection )
                return true;
        }
    }

    //- check if any face edges intersect the triangle
    forAll(f, pI)
    {
        const point& s = facePoints[f[pI]];
        const point& e = facePoints[f.nextLabel(pI)];

        const bool intersected =
            help::triLineIntersection
            (
                surface,
                triI,
                s,
                e,
                intersection
            );

        if( intersected )
            return true;
    }

    return false;
}

inline bool doEdgesOverlap
(
    const point& e0p0,
    const point& e0p1,
    const point& e1p0,
    const point& e1p1,
    FixedList<point, 2>& overlappingPart,
    const scalar distTol,
    const scalar cosTol
)
{
    if( distTol < 0.0 )
    {
        WarningIn
        (
            "inline bool doEdgesOverlap(const point, const point&,"
            "const point&, const point&, const scalar, const scalar,"
            "FixedList<point, 2>& overlappingPart)"
        ) << "Distance is not specified" << endl;

        return false;
    }

    vector e0 = e0p1 - e0p0;
    const scalar de0 = mag(e0);
    e0 /= (de0 + VSMALL);

    vector e1 = e1p1 - e1p0;
    const scalar de1 = mag(e1);
    e1 /= (de1 + VSMALL);

    //- check the angle deviation between the two vectors
    if( mag(e0 & e1) < cosTol )
        return false;

    scalar t00 = (e1p0 - e0p0) & e0;
    scalar t01 = (e1p1 - e0p0) & e0;

    scalar t10 = (e0p0 - e1p0) & e1;
    scalar t11 = (e0p1 - e1p0) & e1;

    //- check if points are colinear within the tolerance
    if
    (
        (magSqr(e0p0 + t00*e0 - e1p0) <= distTol*distTol) ||
        (magSqr(e0p0 + t01*e0 - e1p1) <= distTol*distTol) ||
        (magSqr(e1p0 + t10*e1 - e0p0) <= distTol*distTol) ||
        (magSqr(e1p0 + t11*e1 - e0p1) <= distTol*distTol)
    )
    {
        vector vec = e0 + (((e0 & e1) > 0.0) ? e1 : -e1);
        vec /= (mag(vec) + VSMALL);
        const point origin = 0.25 * (e0p0 + e0p1 + e1p0 + e1p1);

        //- calculate parameters on the regression line
        t00 = (e0p0 - origin) & vec;
        t01 = (e0p1 - origin) & vec;
        t10 = (e1p0 - origin) & vec;
        t11 = (e1p1 - origin) & vec;

        //- check interval overlapping over the line
        const scalar t0Min = Foam::min(t00, t01);
        const scalar t0Max = Foam::max(t00, t01);
        const scalar t1Min = Foam::min(t10, t11);
        const scalar t1Max = Foam::max(t10, t11);

        if( t1Min < t0Max )
        {
            overlappingPart[0] = origin + t1Min * vec;
            overlappingPart[1] = origin + t0Max * vec;

            return true;
        }
        else if( t0Min < t1Max )
        {
            overlappingPart[0] = origin + t0Min * vec;
            overlappingPart[1] = origin + t1Max * vec;

            return true;
        }
    }

    return false;
}

inline bool doTrianglesOverlap
(
    const triangle<point, point>& tri0,
    const triangle<point, point>& tri1,
    DynList<point>& overlappingPolygon,
    const scalar distTol,
    const scalar cosTol
)
{
    if( distTol < 0.0 )
    {
        WarningIn
        (
            "inline bool doTrianglesOverlap(const triangle<point, point>&,"
            " const triangle<point, point>&, DynList<point>&,"
            " const scalar, const scalar)"
        ) << "Distance is not specified" << endl;

        return false;
    }

    vector n0 = tri0.normal();
    const scalar dn0 = mag(n0);
    n0 /= (dn0 + VSMALL);

    if( dn0 < VSMALL )
        return false;

    vector n1 = tri1.normal();
    const scalar dn1 = mag(n1);
    n1 /= (dn1 + VSMALL);

    if( dn1 < VSMALL )
        return false;

    //- check the angle deviation between the two vectors
    if( (mag(n0 & n1) < cosTol) && (dn0 >= VSMALL) && (dn1 >= VSMALL) )
        return false;

    //- check if the two nearest points are within tolerance
    if
    (
        (mag((tri1.a() - tri0.a()) & n0) < distTol) ||
        (mag((tri1.b() - tri0.a()) & n0) < distTol) ||
        (mag((tri1.c() - tri0.a()) & n0) < distTol) ||
        (mag((tri0.a() - tri1.a()) & n1) < distTol) ||
        (mag((tri0.b() - tri1.a()) & n1) < distTol) ||
        (mag((tri0.c() - tri1.a()) & n1) < distTol)
    )
    {
        vector vec = n0 + ((n0 & n1) >= 0.0 ? n1 : -n1);
        vec /= (mag(vec) + VSMALL);
        const point origin = 0.5 * (tri0.centre() + tri1.centre());

        overlappingPolygon.clear();

        //- calculate max distance point from the origin
        point bestPoint = tri0.a();
        scalar distSq = magSqr(tri0.a() - origin);

        scalar dSq = magSqr(tri0.b() - origin);
        if( dSq > distSq )
        {
            distSq = dSq;
            bestPoint = tri0.b();
        }

        dSq = magSqr(tri0.c() - origin);
        if( dSq > distSq )
        {
            distSq = dSq;
            bestPoint = tri0.c();
        }

        dSq = magSqr(tri1.a() - origin);
        if( dSq > distSq )
        {
            distSq = dSq;
            bestPoint = tri1.a();
        }

        dSq = magSqr(tri1.b() - origin);
        if( dSq > distSq )
        {
            distSq = dSq;
            bestPoint = tri1.b();
        }

        dSq = magSqr(tri1.c() - origin);
        if( dSq > distSq )
        {
            distSq = dSq;
            bestPoint = tri1.c();
        }

        if( distSq < VSMALL )
            return false;

        //- transform into planar coordinates
        vector x = (bestPoint - origin) - ((bestPoint - origin) & vec) * vec;
        x /= (mag(x) + VSMALL);
        vector y = vec ^ x;

        DynList<point2D, 6> poly2D(3);
        poly2D[0] = point2D((tri0.a() - origin) & x, (tri0.a() - origin) & y);
        poly2D[1] = point2D((tri0.b() - origin) & x, (tri0.b() - origin) & y);
        poly2D[2] = point2D((tri0.c() - origin) & x, (tri0.c() - origin) & y);

        FixedList<point2D, 3> t1Proj;
        t1Proj[0] = point2D((tri1.a() - origin) & x, (tri1.a() - origin) & y);
        t1Proj[1] = point2D((tri1.b() - origin) & x, (tri1.b() - origin) & y);
        t1Proj[2] = point2D((tri1.c() - origin) & x, (tri1.c() - origin) & y);

        forAll(t1Proj, eI)
        {
            const vector2D vec = t1Proj[(eI+1)%3] - t1Proj[eI];

            DynList<scalar, 6> distance(poly2D.size());
            forAll(poly2D, pI)
            {
                const vector2D pVec = poly2D[pI] - t1Proj[eI];
                distance[pI] = vec.y() * pVec.x() - vec.x() * pVec.y();
            }

            DynList<point2D, 6> newPoly2D;

            forAll(distance, pI)
            {
                if( distance[pI] >= 0.0 )
                {
                    newPoly2D.append(poly2D[pI]);

                    if( distance.fcElement(pI) < 0.0 )
                    {
                        //- this is very sensitive to floaing point tolerances
                        const point2D newP =
                            (
                                mag(distance[pI]) * poly2D.fcElement(pI) +
                                mag(distance.fcElement(pI)) * poly2D[pI]
                            ) /
                            (
                                mag(distance[pI]) +
                                mag(distance.fcElement(pI)) +
                                VSMALL
                            );

                        newPoly2D.append(newP);
                    }
                }
                else if( distance.fcElement(pI) >= 0.0 )
                {
                    //- this is very sensitive to floaing point tolerances
                    const point2D newP =
                        (
                            mag(distance[pI]) * poly2D.fcElement(pI) +
                            mag(distance.fcElement(pI)) * poly2D[pI]
                        ) /
                        (
                            mag(distance[pI]) +
                            mag(distance.fcElement(pI)) +
                            VSMALL
                        );

                    newPoly2D.append(newP);
                }
            }

            poly2D = newPoly2D;
        }

        //- check if the overlapping polygon exists
        if( poly2D.size() == 0 )
        {
            overlappingPolygon.clear();
            return false;
        }

        //- fill the overlapping polygon
        overlappingPolygon.setSize(poly2D.size());
        forAll(poly2D, pI)
        {
            const point2D& pp = poly2D[pI];
            overlappingPolygon[pI] = origin + x * pp.x() + y * pp.y();
        }

        //- cleanup overlapping polygon
        DynList<bool> removePoint(overlappingPolygon.size(), false);
        forAll(overlappingPolygon, i)
        {
            const point& p = overlappingPolygon[i];

            for(label j=i+1;j<overlappingPolygon.size();++j)
            {
                const point& op = overlappingPolygon[j];

                if( magSqr(p - op) < Foam::sqr(distTol) )
                    removePoint[j] = true;
            }
        }

        DynList<point> copy;
        forAll(overlappingPolygon, i)
        {
            if( !removePoint[i] )
                copy.append(overlappingPolygon[i]);
        }

        overlappingPolygon = copy;

        return true;
    }

    overlappingPolygon.clear();
    return false;
}

inline bool doTrianglesIntersect
(
    const triangle<point, point> &tri0,
    const triangle<point, point> &tri1,
    const scalar distTol
)
{
    if( distTol < 0.0 )
    {
        WarningIn
        (
            "inline bool doTrianglesIntersect(const triangle<point, point>&,"
            " const triangle<point, point>&,"
            " const scalar, const scalar)"
        ) << "Distance is not specified" << endl;

        return false;
    }

    //- find distances of points from the second triangle
    point np = nearestPointOnTheTriangle(tri1, tri0.a());
    if( magSqr(np - tri0.a()) < distTol*distTol )
        return true;
    np = nearestPointOnTheTriangle(tri1, tri0.b());
    if( magSqr(np - tri0.b()) < distTol*distTol )
        return true;
    np = nearestPointOnTheTriangle(tri1, tri0.c());
    if( magSqr(np - tri0.c()) < distTol*distTol )
        return true;

    //- find distances of points from the first triangle
    np = nearestPointOnTheTriangle(tri0, tri1.a());
    if( magSqr(np - tri1.a()) < distTol*distTol )
        return true;
    np = nearestPointOnTheTriangle(tri0, tri1.b());
    if( magSqr(np - tri1.b()) < distTol*distTol )
        return true;
    np = nearestPointOnTheTriangle(tri0, tri1.c());
    if( magSqr(np - tri1.c()) < distTol*distTol )
        return true;

    //- find edge intersections
    point intersection;
    if( triLineIntersection(tri0, tri1.a(), tri1.b(), intersection) )
        return true;
    if( triLineIntersection(tri0, tri1.b(), tri1.c(), intersection) )
        return true;
    if( triLineIntersection(tri0, tri1.c(), tri1.a(), intersection) )
        return true;

    if( triLineIntersection(tri1, tri0.a(), tri0.b(), intersection) )
        return true;
    if( triLineIntersection(tri1, tri0.b(), tri0.c(), intersection) )
        return true;
    if( triLineIntersection(tri1, tri0.c(), tri0.a(), intersection) )
        return true;

    return false;
}

inline bool doTrianglesIntersect
(
    const triangle<point, point>& tri0,
    const triangle<point, point>& tri1,
    DynList<point> &intersectionPoints,
    const scalar distTol
)
{
    if( distTol < 0.0 )
    {
        WarningIn
        (
            "inline bool doTrianglesIntersect(const triangle<point, point>&,"
            " const triangle<point, point>&, DynList<point>&,"
            " const scalar, const scalar)"
        ) << "Distance is not specified" << endl;

        return false;
    }

    vector n0 = tri0.normal();
    const scalar dn0 = mag(n0);
    n0 /= (dn0 + VSMALL);

    vector n1 = tri1.normal();
    const scalar dn1 = mag(n1);
    n1 /= (dn1 + VSMALL);

    //- distance of the points of the first triangle from the plane
    //- of the second triangle
    FixedList<scalar, 3> distancesTri0;
    distancesTri0[0] = (tri0.a() - tri1.a()) & n1;
    distancesTri0[1] = (tri0.b() - tri1.a()) & n1;
    distancesTri0[2] = (tri0.c() - tri1.a()) & n1;

    forAll(distancesTri0, i)
        if( mag(distancesTri0[i]) < distTol )
            distancesTri0[i] = 0.0;

    bool hasPositive(false), hasNegative(false), hasZero(false);
    DynList<point, 2> intersectionLine0;
    forAll(distancesTri0, pI)
    {
        if( distancesTri0[pI] >= distTol )
        {
            hasPositive = true;
        }
        else if( distancesTri0[pI] <= -distTol )
        {
            hasNegative = true;
        }
        else
        {
            hasZero = true;
        }
    }

    //- find points on the intersection line
    if
    (
        (hasPositive && !hasNegative && !hasZero) ||
        (hasNegative && !hasPositive && !hasZero)
    )
    {
        //- there can be no intersection
        intersectionPoints.clear();
        return false;
    }
    else if
    (
        (hasPositive && (hasNegative || hasZero)) ||
        (hasNegative && (hasPositive || hasZero))
    )
    {
        //- find points on the intersection line
        if( mag(distancesTri0[0]) < distTol)
        {
            intersectionLine0.append(tri0.a());
        }
        else if( distancesTri0[0] * distancesTri0[1] < 0.0 )
        {
            intersectionLine0.append
            (
                (tri0.a() * distancesTri0[1] - tri0.b() * distancesTri0[0]) /
                (distancesTri0[0] - distancesTri0[1])
            );
        }

        if( mag(distancesTri0[1]) < distTol )
        {
            intersectionLine0.append(tri0.b());
        }
        else if( distancesTri0[1] * distancesTri0[2] < 0.0 )
        {
            intersectionLine0.append
            (
                (tri0.b() * distancesTri0[2] - tri0.c() * distancesTri0[1]) /
                (distancesTri0[1] - distancesTri0[2])
            );
        }

        if( mag(distancesTri0[2]) < distTol )
        {
            intersectionLine0.append(tri0.c());
        }
        else if( distancesTri0[2] * distancesTri0[0] < 0.0 )
        {
            intersectionLine0.append
            (
                (tri0.c() * distancesTri0[0] - tri0.a() * distancesTri0[2]) /
                (distancesTri0[2] - distancesTri0[0])
            );
        }
    }
    else
    {
        //- these triangles are in the same plane
        //- check if they overlap
        return doTrianglesOverlap(tri0, tri1, intersectionPoints, distTol);
    }

    //- distance of the points of the second triangle from the plane
    //- of the first triangle
    FixedList<scalar, 3> distancesTri1;
    distancesTri1[0] = (tri1.a() - tri0.a()) & n0;
    distancesTri1[1] = (tri1.b() - tri0.a()) & n0;
    distancesTri1[2] = (tri1.c() - tri0.a()) & n0;

    forAll(distancesTri1, i)
        if( mag(distancesTri1[i]) < distTol )
            distancesTri1[i] = 0.0;

    hasPositive = false;
    hasNegative = false;
    hasZero = false;

    DynList<point, 2> intersectionLine1;
    forAll(distancesTri1, pI)
    {
        if( distancesTri1[pI] >= distTol )
        {
            hasPositive = true;
        }
        else if( distancesTri1[pI] <= -distTol )
        {
            hasNegative = true;
        }
        else
        {
            hasZero = true;
        }
    }

    if
    (
        (hasPositive && !hasNegative && !hasZero) ||
        (hasNegative && !hasPositive && !hasZero)
    )
    {
        //- there can be no intersection
        intersectionPoints.clear();
        return false;
    }
    else if
    (
        (hasPositive && (hasNegative || hasZero)) ||
        (hasNegative && (hasPositive || hasZero))
    )
    {
        //- find points on the intersection line
        if( mag(distancesTri1[0]) < distTol)
        {
            intersectionLine1.append(tri1.a());
        }
        else if( distancesTri1[0] * distancesTri1[1] < 0.0 )
        {
            intersectionLine1.append
            (
                (tri1.a() * distancesTri1[1] - tri1.b() * distancesTri1[0]) /
                (distancesTri1[0] - distancesTri1[1])
            );
        }

        if( mag(distancesTri1[1]) < distTol)
        {
            intersectionLine1.append(tri1.b());
        }
        else if( distancesTri1[1] * distancesTri1[2] < 0.0 )
        {
            intersectionLine1.append
            (
                (tri1.b() * distancesTri1[2] - tri1.c() * distancesTri1[1]) /
                (distancesTri1[1] - distancesTri1[2])
            );
        }

        if( mag(distancesTri1[2]) < distTol)
        {
            intersectionLine1.append(tri1.c());
        }
        else if( distancesTri1[2] * distancesTri1[0] < 0.0 )
            intersectionLine1.append
            (
                (tri1.c() * distancesTri1[0] - tri1.a() * distancesTri1[2]) /
                (distancesTri1[2] - distancesTri1[0])
            );
    }
    else
    {
        //- these triangles are in the same plane
        //- check if they overlap
        return doTrianglesOverlap(tri0, tri1, intersectionPoints, distTol);
    }

    //- calculate the direction of the intersection line
    vector vec = n0 ^ n1;

    //- n0 and n1 are nearly coplanar, try overlap intersection
    if( magSqr(vec) < SMALL )
    {
        //- these triangles are in the same plane
        return doTrianglesOverlap(tri0, tri1, intersectionPoints, distTol);
    }

    //- calculate intervals on the intersection line
    scalar t0Min(VGREAT), t0Max(-VGREAT);
    point p0Min(vector::zero), p0Max(vector::zero);
    forAll(intersectionLine0, i)
    {
        const scalar t = intersectionLine0[i] & vec;

        if( t < t0Min )
        {
            t0Min = t;
            p0Min = intersectionLine0[i];
        }
        if( t > t0Max )
        {
            t0Max = t;
            p0Max = intersectionLine0[i];
        }
    }

    scalar t1Min(VGREAT), t1Max(-VGREAT);
    point p1Min(vector::zero), p1Max(vector::zero);
    forAll(intersectionLine1, i)
    {
        const scalar t = intersectionLine1[i] & vec;

        if( t < t1Min )
        {
            t1Min = t;
            p1Min = intersectionLine1[i];
        }
        if( t > t1Max )
        {
            t1Max = t;
            p1Max = intersectionLine1[i];
        }
    }

    if( (t1Min <= t0Max) && (t1Max >= t0Min) )
    {
        intersectionPoints.setSize(2);
        intersectionPoints[0] = p1Min;
        intersectionPoints[1] = p0Max;

        return true;
    }
    else if( (t0Min <= t1Max) && (t0Max >= t1Min) )
    {
        intersectionPoints.setSize(2);
        intersectionPoints[0] = p0Min;
        intersectionPoints[1] = p1Max;

        return true;
    }

    intersectionPoints.setSize(0);
    return false;
}

inline point nearestPointOnFace
(
    const face& f,
    const pointField& points,
    const point& p
)
{
    point np;
    scalar dSq(VGREAT);

    const point c = faceCentre(points, f);
    forAll(f, eI)
    {
        const triangle<point, point> tri
        (
            points[f[eI]],
            points[f.nextLabel(eI)],
            c
        );

        const point nearest = nearestPointOnTheTriangle(tri, p);

        const scalar currDistSq = magSqr(nearest - p);

        if( currDistSq < dSq )
        {
            dSq = currDistSq;
            np = nearest;
        }
    }

    return np;
}

inline bool pointInsideFace
(
    const point& p,
    const face& f,
    const vector& n,
    const pointField& fp,
    const scalar distTol
)
{
    const edgeList fe = f.edges();
    forAll(f, pI)
    {
        if( mag(p - fp[f[pI]]) < distTol )
            return true;

        vector pv = p - fp[f[pI]];
        pv /= mag(pv);

        vector lv = n ^ fe[pI].vec(fp);
        lv /= mag(lv);

        const scalar d = pv & lv;
        if( d < -distTol )
        {
            return false;
        }
    }

    return true;
}

inline bool pointInsideFace
(
    const point& p,
    const face& f,
    const pointField& fp,
    const scalar distTol
)
{
    const point c = faceCentre(fp, f);
    const scalar tolSqr = sqr(distTol);

    forAll(f, eI)
    {
        const edge fe = f.faceEdge(eI);

        triangle<point, point> tri
        (
            fp[fe[0]],
            fp[fe[1]],
            c
        );

        const point np = nearestPointOnTheTriangle(tri, p);

        if( magSqr(np - p) <= tolSqr )
            return true;
    }

    return false;
}

inline bool isFaceConvexAndOk
(
    const face& f,
    const pointField& fp,
    DynList<bool>& OkPoints
)
{
    bool valid(true);

    vector normal = f.normal(fp);
    const scalar magN = mag(normal);

    //- face has zero area. All points are inverted
    if( magN < VSMALL )
    {
        OkPoints.setSize(f.size());
        OkPoints = false;

        return false;
    }

    normal /= (magN + VSMALL);

    //- project points into the plane formed by the face centre
    //- and the face normal
    plane pl(faceCentre(fp, f), normal);

    DynList<point> projFace(f.size());
    forAll(f, pI)
        projFace[pI] = pl.nearestPoint(fp[f[pI]]);

    //- calculate normalised edge vectors
    DynList<vector> edgeVecs(f.size());
    forAll(f, eI)
    {
        vector& v = edgeVecs[eI];

        const edge e = f.faceEdge(eI);
        v = e.vec(fp);
        v /= (mag(v) + VSMALL);
    }

    //- check if the face is convex
    OkPoints.setSize(f.size());
    forAll(f, pI)
    {
        const vector& pv = edgeVecs[f.rcIndex(pI)];
        const vector& nv = edgeVecs[pI];

        if( ((pv ^ nv) & normal) >= -0.05 )
        {
            OkPoints[pI] = true;
        }
        else
        {
            OkPoints[pI] = false;
            valid = false;
        }
    }

    return valid;
}

inline scalar tetCircumRadius(const tetrahedron<point, point>& tet)
{
    try
    {
        const vector a = tet.b() - tet.a();
        const vector b = tet.c() - tet.a();
        const vector c = tet.d() - tet.a();

        const scalar lambda = magSqr(c) - (a & c);
        const scalar mu = magSqr(b) - (a & b);

        const vector ba = b ^ a;
        const vector ca = c ^ a;

        vector num = lambda*ba - mu*ca;
        scalar denom = (c & ba);

        if (Foam::mag(denom) < ROOTVSMALL * (1.0 + mag(num)) )
        {
            // Degenerate tetrahedron, returning GREAT for circumRadius.
            return GREAT;
        }

        return Foam::mag(0.5*(a + num/denom));
    }
    catch(...)
    {
        return GREAT;
    }
}

inline scalar tetQuality(const tetrahedron<point, point>& tet)
{
    return
        tet.mag() /
        (
            (8.0/(9.0*sqrt(3.0))) * pow3(min(tetCircumRadius(tet), GREAT)) +
            ROOTVSMALL
        );
}

inline point nearestPointOnLine
(
    const point& origin,
    const vector& v,
    const point& p
)
{
    const scalar dSq = magSqr(v);
    const vector k = p - origin;

    if( dSq < VSMALL )
        return origin;

    return origin + ((v / dSq) * (v & k));
}

inline point nearestPointOnTheEdge
(
    const point& edgePoint0,
    const point& edgePoint1,
    const point& p
)
{
    const vector e = edgePoint1 - edgePoint0;
    const scalar d = mag(e);
    const vector k = p - edgePoint0;

    if( d < ROOTVSMALL )
        return edgePoint0;

    return edgePoint0 + ((e / (d*d)) * (e & k));
}

inline point nearestPointOnTheEdgeExact
(
    const point& edgePoint0,
    const point& edgePoint1,
    const point& p
)
{
    const vector e = edgePoint1 - edgePoint0;
    const scalar dSq = magSqr(e);
    const vector k = p - edgePoint0;

    if( dSq < VSMALL )
        return edgePoint0;

    const scalar t = (e & k) / (dSq + VSMALL);
    if( t > 1.0 )
    {
        return edgePoint1;
    }
    else if( t < 0.0 )
    {
        return edgePoint0;
    }

    return edgePoint0 + (e * t);
}

inline scalar distanceOfPointFromTheEdge
(
    const point& edgePoint0,
    const point& edgePoint1,
    const point& p
)
{
    return mag(nearestPointOnTheEdge(edgePoint0, edgePoint1, p) - p);
}

inline bool isGapTest
(
    const point& p0,
    const vector& n0,
    const point& p1,
    const vector& n1,
    const scalar minCosBetweenNormals,
    const scalar minCosDistToNormal
)
{
    const scalar normalCos = n0 & n1;

    if( normalCos >= minCosBetweenNormals )
    {
        vector n = n0 + n1;
        n /= (mag(n) + VSMALL);

        vector dVec = p1 - p0;
        dVec /= (mag(dVec) + VSMALL);

        if( mag(dVec & n) > minCosDistToNormal )
            return true;
    }
    else if( normalCos <= -minCosBetweenNormals )
    {
        vector n = n0 - n1;
        n /= (mag(n) + VSMALL);

        vector dVec = p1 - p0;
        dVec /= (mag(dVec) + VSMALL);

        if( mag(dVec & n) > minCosDistToNormal )
            return true;
    }

    return false;
}

inline bool doTriangleMakeAGap
(
    const triangle<point, point>& tri0,
    const triangle<point, point>& tri1,
    scalar& distSq,
    const scalar minCosBetweenNormals,
    const scalar minCosDistToNormal
)
{
    point nearest0, nearest1;

    nearestPointsOnTriangles(tri0, tri1, nearest0, nearest1);

    vector n0 = tri0.normal();
    n0 /= (mag(n0) + VSMALL);
    vector n1 = tri1.normal();
    n1 /= (mag(n1) + VSMALL);

    distSq = magSqr(nearest0 - nearest1);

    return isGapTest
    (
        nearest0,
        n0,
        nearest1,
        n1,
        minCosBetweenNormals,
        minCosDistToNormal
    );
}

inline bool doTrianglesMakeAGap
(
    const triSurf& surf,
    const label t0,
    const label t1,
    scalar& distSq,
    const scalar minCosBetweenNormals,
    const scalar minCosDistToNormal
)
{
    const triangle<point, point> tri0 = surfaceTriangle(surf, t0);
    const triangle<point, point> tri1 = surfaceTriangle(surf, t1);

    return doTriangleMakeAGap
    (
        tri0,
        tri1,
        distSq,
        minCosBetweenNormals,
        minCosDistToNormal
    );
}

inline label numberOfFaceGroups
(
    const labelHashSet& containedElements,
    const point& centre,
    const scalar range,
    const triSurf& surface,
    const List<DynList<label> >* triangleTrianglesPtr
)
{
    const pointField& points = surface.points();
    const edgeLongList& edges = surface.edges();
    const VRWGraph& faceEdges = surface.facetEdges();
    const VRWGraph& edgeFaces = surface.edgeFacets();

    labelHashSet triaInRange(containedElements.size());
    const scalar rangeSq = range * range;
    forAllConstIter(labelHashSet, containedElements, it)
    {
        const point p = nearestPointOnTheTriangle(it.key(), surface, centre);
        if( magSqr(p - centre) < rangeSq )
            triaInRange.insert(it.key());
    }

    Map<label> elGroup(triaInRange.size());
    Map<label> testEdge(triaInRange.size());

    label nGroups(0);

    DynList<label> front;

    forAllConstIter(labelHashSet, triaInRange, it)
        if( !elGroup.found(it.key()) )
        {
            front.clear();
            front.append(it.key());
            elGroup.insert(it.key(), nGroups);

            while( front.size() )
            {
                const label fLabel = front.removeLastElement();

                forAllRow(faceEdges, fLabel, feI)
                {
                    const label edgeI = faceEdges(fLabel, feI);

                    //- check if the edge intersects the bounding box
                    if( testEdge.found(edgeI) )
                    {
                        if( !testEdge[edgeI] )
                            continue;
                    }
                    else
                    {
                        const point& s = points[edges[edgeI][0]];
                        const point& e = points[edges[edgeI][1]];
                        const point np = nearestPointOnTheEdge(s, e, centre);
                        if( magSqr(np - centre) < rangeSq )
                        {
                            testEdge.insert(edgeI, 1);
                        }
                        else
                        {
                            testEdge.insert(edgeI, 0);
                            continue;
                        }
                    }

                    forAllRow(edgeFaces, edgeI, efI)
                    {
                        const label nei = edgeFaces(edgeI, efI);
                        if
                        (
                            triaInRange.found(nei) &&
                            !elGroup.found(nei)
                        )
                        {
                            elGroup.insert(nei, nGroups);
                            front.append(nei);
                        }
                    }
                }

                if( triangleTrianglesPtr )
                {
                    const DynList<label>& neiTriangles =
                        triangleTrianglesPtr->operator[](fLabel);

                    forAll(neiTriangles, i)
                    {
                        const label nei = neiTriangles[i];

                        if
                        (
                            triaInRange.found(nei) &&
                            !elGroup.found(nei)
                        )
                        {
                            elGroup.insert(nei, nGroups);
                            front.append(nei);
                        }
                    }
                }
            }

            ++nGroups;
        }

    return nGroups;
}

inline label numberOfEdgeGroups
(
    const labelHashSet& containedEdges,
    const point& centre,
    const scalar range,
    const triSurf& surface,
    const List<DynList<label> >* edgeEdgesPtr
)
{
    const pointField& points = surface.points();
    const edgeLongList& edges = surface.edges();
    const VRWGraph& pointEdges = surface.pointEdges();

    const scalar rangeSq = range * range;
    labelHashSet edgesInRange(containedEdges.size());
    forAllConstIter(labelHashSet, containedEdges, it)
    {
        const edge& e = edges[it.key()];
        const point& sp = points[e[0]];
        const point& ep = points[e[1]];

        const point p = nearestPointOnTheEdgeExact(sp, ep, centre);
        if( magSqr(p - centre) < rangeSq )
            edgesInRange.insert(it.key());
    }

    Map<label> elGroup(edgesInRange.size());
    Map<label> pointTest(edgesInRange.size());

    label nGroups(0);

    DynList<label> front;

    forAllConstIter(labelHashSet, edgesInRange, it)
        if( !elGroup.found(it.key()) )
        {
            front.clear();
            front.append(it.key());
            elGroup.insert(it.key(), nGroups);

            while( front.size() )
            {
                const label eLabel = front.removeLastElement();
                const edge& e = edges[eLabel];

                for(label i=0;i<2;++i)
                {
                    if( pointTest.found(e[i]) )
                    {
                        if( !pointTest[e[i]] )
                            continue;
                    }
                    else
                    {
                        if( magSqr(points[e[i]] - centre) < rangeSq )
                        {
                            pointTest.insert(e[i], 1);
                        }
                        else
                        {
                            pointTest.insert(e[i], 0);
                            continue;
                        }
                    }

                    forAllRow(pointEdges, e[i], peI)
                    {
                        const label nei = pointEdges(e[i], peI);

                        if( edgesInRange.found(nei) && !elGroup.found(nei) )
                        {
                            elGroup.insert(nei, nGroups);
                            front.append(nei);
                        }
                    }
                }

                if( edgeEdgesPtr )
                {
                    const DynList<label>& neiEdges =
                        edgeEdgesPtr->operator[](eLabel);

                    forAll(neiEdges, i)
                    {
                        const label nei = neiEdges[i];

                        if( edgesInRange.found(nei) && !elGroup.found(nei) )
                        {
                            elGroup.insert(nei, nGroups);
                            front.append(nei);
                        }
                    }
                }
            }

            ++nGroups;
        }

    return nGroups;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
