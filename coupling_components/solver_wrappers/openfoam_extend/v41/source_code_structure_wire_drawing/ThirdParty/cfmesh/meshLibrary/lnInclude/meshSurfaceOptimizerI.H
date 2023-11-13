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

#include "demandDrivenData.H"
#include "meshSurfaceOptimizer.H"
#include "plane.H"
#include "Map.H"
#include "surfaceOptimizer.H"
#include "surfaceOptimizerHeight.H"
#include "helperFunctions.H"
#include "partTriMeshSimplex.H"
#include "quadricFitting.H"

//#define DEBUGSmooth

# ifdef DEBUGSmooth
#include "triSurf.H"
#include "triSurfModifier.H"
# endif

#define transformToPlane

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const partTriMesh& meshSurfaceOptimizer::triMesh() const
{
    if( !triMeshPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
        {
            FatalErrorIn
            (
                "inline const partTriMesh& meshSurfaceOptimizer"
                "::triMesh() const"
            ) << "Contructing a triangle mesh inside a parallel region"
              << ". This is not thread-safe" << abort(FatalError);
        }
        # endif

        calculateTrianglesAndAddressing();
    }

    return *triMeshPtr_;
}

inline void meshSurfaceOptimizer::updateTriMesh(const labelLongList& selPoints)
{
    if( !triMeshPtr_ )
        FatalErrorIn
        (
            "inline void meshSurfaceOptimizer::updateTriMesh"
            "(const labelLongList&)"
        ) << "triMeshPtr_ is not allocated " << abort(FatalError);

    triMeshPtr_->updateVertices(selPoints);
}

inline void meshSurfaceOptimizer::updateTriMesh()
{
    if( !triMeshPtr_ )
        FatalErrorIn
        (
            "inline void meshSurfaceOptimizer::updateTriMesh()"
        ) << "triMeshPtr_ is not allocated " << abort(FatalError);

    triMeshPtr_->updateVertices();
}

inline bool meshSurfaceOptimizer::transformIntoPlane
(
    const label bpI,
    const plane& pl,
    vector& vecX,
    vector& vecY,
    scalar& minCosSq,
    DynList<point>& pts,
    DynList<triFace>& trias
) const
{
    # ifdef DEBUGSmooth
    Pout << "Transforming boundary node " << bpI << endl;
    # endif

    const partTriMesh& triMesh = this->triMesh();
    const label triPointI = triMesh.meshSurfacePointLabelInTriMesh()[bpI];
    if( triPointI < 0 )
        return false;

    partTriMeshSimplex simplex(triMesh, triPointI);

    const DynList<point, 32>& sPts = simplex.pts();

    minCosSq = 1.0;

    //- create vecX and vecY
    const point& p = simplex.centrePoint();
    pts.setSize(0);
    bool found(false);
    scalar dMax(0.0);
    forAll(sPts, pI)
    {
        const point sp = pl.nearestPoint(sPts[pI]);
        const scalar dSqOrig = magSqr(sPts[pI] - p);
        const scalar d = mag(sp - p);

        if( dSqOrig > VSMALL )
            minCosSq = min(minCosSq, sqr(d) / (dSqOrig + VSMALL));

        if( (d > VSMALL) && (d > dMax) )
        {
            dMax = d;
            vecX = sp - p;
            vecX /= d;
            vecY = pl.normal() ^ vecX;
            vecY /= mag(vecY);
            found = true;
        }
    }

    if( !found )
        return false;

    trias = simplex.triangles();

    # ifdef DEBUGSmooth
    Pout << "VecX " << vecX << endl;
    Pout << "vecY " << vecY << endl;
    Pout << "minCosSq " << minCosSq << endl;
    # endif

    //- transform the vertices
    pts.setSize(sPts.size());

    forAll(sPts, pI)
    {
        const point pt
        (
            ((sPts[pI] - p) & vecX),
            ((sPts[pI] - p) & vecY),
            0.0
        );

        pts[pI] = pt;
    }

    # ifdef DEBUGSmooth
    Pout << "Original triangles " << endl;
    forAll(simplex.triangles(), triI)
        Info << "Tri " << triI << " is " << simplex.triangles()[triI] << endl;
    Pout << "Transformed triangles are " << trias << endl;
    Pout << "Transformed vertices " << pts << endl;

    triSurf surf;
    triSurfModifier sMod(surf);
    pointField& sPoints = sMod.pointsAccess();
    sPoints.setSize(pts.size());
    forAll(sPoints, i)
        sPoints[i] = pts[i];
    LongList<labelledTri>& sTrias = sMod.facetsAccess();
    sTrias.setSize(trias.size());
    forAll(sTrias, i)
    {
        labelledTri& tf = sTrias[i];
        tf[0] = trias[i][0];
        tf[1] = trias[i][1];
        tf[2] = trias[i][2];

        tf.region() = 0;
    }
    sMod.patchesAccess().setSize(1);
    sMod.patchesAccess()[0].name() = "bnd";
    sMod.patchesAccess()[0].geometricType() = "patch";

    fileName sName("bndSimplex_");
    sName += help::scalarToText(bpI);
    sName += ".stl";
    surf.writeSurface(sName);
    # endif

    return found;
}

inline void meshSurfaceOptimizer::localDataFaceFlatness
(
    const label bpI,
    vector& disp,
    scalar& sumW,
    const bool weighted
) const
{
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();
    const pointFieldPMG& points = surfaceEngine_.points();
    const vectorLongList& faceCentres = surfaceEngine_.faceCentres();
    const vectorLongList& faceNormals = surfaceEngine_.faceNormals();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    const label pointI = bPoints[bpI];
    const point& p = points[pointI];

    disp = vector::zero;
    sumW = 0.0;
    scalar w = 1.0;

    forAllRow(pointFaces, bpI, pfI)
    {
        const label bfI = pointFaces(bpI, pfI);

        const point& c = faceCentres[bfI];
        const vector& fn = faceNormals[bfI];
        const scalar fnSq = magSqr(fn) + VSMALL;

        const vector d = ((c - p) & fn) * fn / fnSq;

        if( weighted )
            w = 1.0 / fnSq;

        disp += w * d;
        sumW += w;
    }
}

inline point meshSurfaceOptimizer::newPositionFaceFlatness
(
    const label bpI,
    const vector& disp,
    const scalar& sumW
) const
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    return points[bPoints[bpI]] + disp / sumW;
}

inline point meshSurfaceOptimizer::newPositionFaceFlatness
(
    const label bpI,
    const bool weighted
) const
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    if( vertexType_[bpI] & LOCKED )
        return points[bPoints[bpI]];

    vector disp(vector::zero);
    scalar sumW(0.0);

    localDataFaceFlatness(bpI, disp, sumW, weighted);

    return newPositionFaceFlatness(bpI, disp, sumW);
}

inline point meshSurfaceOptimizer::newPositionLaplacian
(
    const label bpI,
    const DynList<point>& neiPoints,
    const bool transformIntoPlane,
    const bool weighted
) const
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const vectorLongList& pNormals = surfaceEngine_.pointNormals();

    if( vertexType_[bpI] & LOCKED )
        return points[bPoints[bpI]];

    const point& p = points[bPoints[bpI]];

    vector newP(vector::zero);
    if( transformIntoPlane )
    {
        const vector& pNormal = pNormals[bpI];

        if( magSqr(pNormal) > VSMALL )
        {
            # ifdef transformToPlane
            plane pl(p, pNormal);

            scalar sumW(0.0);
            scalar minCosSq(1.0);
            scalar w = 1.0;

            forAll(neiPoints, i)
            {
                const point& neiP = neiPoints[i];
                const point pProj = pl.nearestPoint(neiP);

                const scalar dSq = magSqr(neiP - p);
                const scalar hSq = magSqr(pProj - neiP);

                if( weighted )
                    w = (1.0 / (sqrt(dSq) + VSMALL));

                newP += w * pProj;
                sumW += w;

                minCosSq = min(minCosSq, mag(1.0 - (hSq / (dSq + VSMALL))));
            }

            newP /= sumW;

            minCosSq = 1.0;
            newP = p + minCosSq * (newP - p);

            scalar avgEdge(0.0), maxEdge(0.0), minEdge(VGREAT);
            forAll(neiPoints, i)
            {
                const scalar d = mag(neiPoints[i] - p);
                avgEdge += d;
                minEdge = min(d, minEdge);
                maxEdge = max(d, maxEdge);
            }
            avgEdge /= neiPoints.size();

            return newP;
            # else
            //- quadric patch fitting onto the plane
            DynList<point, 32> patchPoints;
            forAll(neiPoints, i)
                patchPoints.append(neiPoints[i]);

            Info << "Patch points " << patchPoints << endl;
            const quadricFitting qFit(p, pNormal, patchPoints);

            if( qFit.isValid() )
            {
                const DynList<point, 32>& transformedPoints =
                    qFit.transformedPoints();

                scalar sumW(0.0);

                if( weighted )
                {
                    forAll(transformedPoints, pI)
                    {
                        const scalar w =
                            (1.0 / (mag(transformedPoints[pI]) + VSMALL));
                        newP += w * transformedPoints[pI];
                        sumW += w;
                    }
                }
                else
                {
                    forAll(transformedPoints, pI)
                    {
                        newP += transformedPoints[pI];
                        sumW += 1.0;
                    }
                }
                newP /= sumW;

                newP = qFit.pointCoordinates(newP.x(), newP.y());

                return newP;
            }
            # endif
        }
    }

    //- points are not transformed into a surface coordinate system
    scalar sumW(0.0);

    if( weighted )
    {
        forAll(neiPoints, i)
        {
            const point& neiP = neiPoints[i];

            const scalar w = (1.0 / (mag(neiP - p) + VSMALL));
            newP += w * neiP;
            sumW += w;
        }
    }
    else
    {
        forAll(neiPoints, i)
        {
            const point& neiP = neiPoints[i];

            newP += neiP;
            sumW += 1.0;
        }
    }

    newP /= sumW;

    return newP;
}

inline point meshSurfaceOptimizer::newPositionLaplacian
(
    const label bpI,
    const bool transformIntoPlane,
    const bool weighted
) const
{
    const VRWGraph& pPoints = surfaceEngine_.pointPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    if( vertexType_[bpI] & LOCKED )
        return points[bPoints[bpI]];

    DynList<point> neiPoints;
    forAllRow(pPoints, bpI, i)
        neiPoints.append(points[bPoints[pPoints(bpI, i)]]);

    return newPositionLaplacian(bpI, neiPoints, transformIntoPlane, weighted);
}

inline point meshSurfaceOptimizer::newPositionLaplacianFC
(
    const label bpI,
    const bool transformIntoPlane,
    const bool weighted
) const
{
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();
    const pointFieldPMG& points = surfaceEngine_.points();
    const vectorLongList& faceCentres = surfaceEngine_.faceCentres();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    if(vertexType_[bpI] & LOCKED )
        return points[bPoints[bpI]];

    DynList<point> neiPoints;
    forAllRow(pointFaces, bpI, pfI)
        neiPoints.append(faceCentres[pointFaces(bpI, pfI)]);

    return newPositionLaplacian(bpI, neiPoints, transformIntoPlane, weighted);
}

inline point meshSurfaceOptimizer::newPositionSurfaceOptimizer
(
    const label bpI,
    const scalar tol
) const
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const vectorLongList& pNormals = surfaceEngine_.pointNormals();

    if( vertexType_[bpI] & LOCKED )
        return points[bPoints[bpI]];

    const point& p = points[bPoints[bpI]];

    const vector& pNormal = pNormals[bpI];
    if( magSqr(pNormal) < VSMALL )
        return points[bPoints[bpI]];

    # ifdef DEBUGSmooth
    Pout << "Smoothing boundary node " << bpI << endl;
    Pout << "Node label in the mesh is " << bPoints[bpI] << endl;
    Pout << "Point coordinates " << p << endl;
    # endif

    # ifdef transformToPlane
    //- project vertices onto the plane
    const plane pl(p, pNormal);

    DynList<point> pts;
    DynList<triFace> trias;
    vector vecX(vector::zero), vecY(vector::zero);
    scalar minCosSq(1.0);
    bool success =
        transformIntoPlane(bpI, pl, vecX, vecY, minCosSq, pts, trias);

    if( !success )
    {
        return p;
    }

    minCosSq = max(minCosSq, 0.05);

    //surfaceOptimizer so(pts, trias);
    surfaceOptimizerHeight so(pts, trias);
    const point newPoint = so.optimizePoint(tol);

    const vector disp = vecX * newPoint.x() + vecY * newPoint.y();

    const point newP = p + minCosSq * disp;
    # else

    const partTriMesh& triMesh = this->triMesh();
    const label triPointI = triMesh.meshSurfacePointLabelInTriMesh()[bpI];
    if( triPointI < 0 )
        return p;

    partTriMeshSimplex simplex(triMesh, triPointI);

    vector n = simplex.normal();
    const scalar magN = mag(n);
    n /= (magN + VSMALL);

    const DynList<point, 32>& sPts = simplex.pts();
    const DynList<triFace, 32>& trias = simplex.triangles();

    //- quadric patch fitting onto the plane
    const quadricFitting qFit(p, n, sPts);

    if( !qFit.isValid() )
    {
        Warning << "Cannot move point " << bpI << endl;
        return p;
    }

    DynList<point> transformedPoints = qFit.transformedPoints();
    forAll(transformedPoints, i)
        transformedPoints[i].z() = 0.0;

    //surfaceOptimizer so(transformedPoints, trias);
    surfaceOptimizerHeight so(transformedPoints, trias);
    const point newPoint = so.optimizePoint(tol);

    const point newP = qFit.pointCoordinates(newPoint.x(), newPoint.y());
    # endif

    if( help::isnan(newP) || help::isinf(newP) )
    {
        WarningIn
        (
            "inline point meshSurfaceOptimizer::newPositionSurfaceOptimizer"
            "(const label, const scalar) const"
        ) << "Cannot move point " << bpI << endl;

        return points[bPoints[bpI]];
    }

    return newP;
}

inline point meshSurfaceOptimizer::newEdgePositionLaplacian
(
    const point& p0,
    const point& p,
    const point& p1
)
{
    const vector d0 = p - p0;
    const vector d1 = p1 - p;

    if( (d0 & d1) > VSMALL )
    {
        //- return the point in the middle of the cubic spline
        const vector a = 2.0 * (p0 - p1) + d0 + d1;
        const vector b = -3.0 * p0 + 3.0 * p1 - 2.0 * d0 - d1;

        return 0.125 * a + 0.25 * b + 0.5 * d0 + p0;
    }
    else
    {
        //- the mesh is tangled here. Revert to simple linear average
        //- that untangles the feature edge
        return 0.5 * (p0 + p1);
    }
}

inline point meshSurfaceOptimizer::newEdgePositionLaplacian
(
    const label bpI
) const
{
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const edgeLongList& edges = surfaceEngine_.edges();
    const VRWGraph& bpEdges = surfaceEngine_.boundaryPointEdges();
    const pointFieldPMG& points = surfaceEngine_.points();

    if( vertexType_[bpI] & LOCKED )
        return points[bPoints[bpI]];

    const labelHashSet& featureEdges = partitionerPtr_->featureEdges();

    DynList<label> edgePoints;

    forAllRow(bpEdges, bpI, i)
    {
        const label beI = bpEdges(bpI, i);

        if( featureEdges.found(beI) )
        {
            edgePoints.append(edges[beI].otherVertex(bPoints[bpI]));
        }
    }

    if( edgePoints.size() != 2 )
        return points[bPoints[bpI]];

    # ifdef DEBUGSearch
    Info << "Edge points " << edgePoints << endl;
    # endif

    return newEdgePositionLaplacian
    (
        points[edgePoints[0]],
        points[bPoints[bpI]],
        points[edgePoints[1]]
    );
}

template<class labelListType>
inline void meshSurfaceOptimizer::lockBoundaryFaces(const labelListType& l)
{
    lockedSurfaceFaces_ = l;

    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
    const labelLongList& bp = surfaceEngine_.bp();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(lockedSurfaceFaces_, lfI)
    {
        const face& bf = bFaces[lockedSurfaceFaces_[lfI]];

        forAll(bf, pI)
            vertexType_[bp[bf[pI]]] |= LOCKED;
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
        const DynList<label>& bpNeiProcs = surfaceEngine_.bpNeiProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(bpNeiProcs, i)
            exchangeData[bpNeiProcs[i]].clear();

        //- prepare data which will be sent to other processors
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( vertexType_[bpI] & LOCKED )
            {
                forAllRow(bpAtProcs, bpI, i)
                {
                    const label neiProc = bpAtProcs(bpI, i);

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
            const label bpI = globalToLocal[receivedData[i]];

            vertexType_[bpI] |= LOCKED;
        }
    }
}

template<class labelListType>
inline void meshSurfaceOptimizer::lockBoundaryPoints(const labelListType& l)
{
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(l, i)
    {
        const label bpI = l[i];

        vertexType_[bpI] |= LOCKED;
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
        const DynList<label>& bpNeiProcs = surfaceEngine_.bpNeiProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(bpNeiProcs, i)
            exchangeData[bpNeiProcs[i]].clear();

        //- prepare data which will be sent to other processors
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( vertexType_[bpI] & LOCKED )
            {
                forAllRow(bpAtProcs, bpI, i)
                {
                    const label neiProc = bpAtProcs(bpI, i);

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
            const label bpI = globalToLocal[receivedData[i]];

            vertexType_[bpI] |= LOCKED;
        }
    }
}

//- lock edge points
inline void meshSurfaceOptimizer::lockFeatureEdges()
{
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(vertexType_, bpI)
        if( vertexType_[bpI] & (EDGE | CORNER) )
            vertexType_[bpI] |= LOCKED;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
