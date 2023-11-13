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
#include "boundaryLayerOptimisation.H"
#include "meshSurfacePartitioner.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "helperFunctions.H"
#include "refLabelledPoint.H"
#include "refLabelledPointScalar.H"
#include "polyMeshGenAddressing.H"
#include "partTriMesh.H"
#include "partTetMeshSimplex.H"
#include "volumeOptimizer.H"

//#define DEBUGLayer

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayerOptimisation::calculateNormalVectors
(
    const direction eType,
    pointNormalsType& pointPatchNormal
) const
{
    const meshSurfaceEngine& mse = meshSurface();
    const labelLongList& facePatch = mse.boundaryFacePatches();
    const labelLongList& bp = mse.bp();
    const VRWGraph& pointFaces = mse.pointFaces();
    const vectorLongList& fNormals = mse.faceNormals();

    //- calculate point normals with respect to all patches at a point
    pointPatchNormal.clear();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        if( !(hairEdgeType_[hairEdgeI] & eType) )
            continue;

        const label bpI = bp[hairEdges_[hairEdgeI][0]];

        //- create an entry in a map
        patchNormalType* patchNormalPtr(NULL);
        # ifdef USE_OMP
        # pragma omp critical
            patchNormalPtr = &pointPatchNormal[bpI];
        # else
        patchNormalPtr = &pointPatchNormal[bpI];
        # endif
        patchNormalType& patchNormal = *patchNormalPtr;

        //- sum normals of faces attached to a point
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            const label patchI = facePatch[bfI];

            if( patchNormal.find(patchI) == patchNormal.end() )
            {
                patchNormal[patchI].first = fNormals[bfI];
                patchNormal[patchI].second = mag(fNormals[bfI]);
            }
            else
            {
                patchNormal[patchI].first += fNormals[bfI];
                patchNormal[patchI].second += mag(fNormals[bfI]);
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- gather information about face normals on other processors
        const Map<label>& globalToLocal =
            mse.globalToLocalBndPointAddressing();
        const DynList<label>& neiProcs = mse.bpNeiProcs();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();

        std::map<label, LongList<refLabelledPointScalar> > exchangeData;
        forAll(neiProcs, i)
            exchangeData[neiProcs[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( pointPatchNormal.find(bpI) != pointPatchNormal.end() )
            {
                const patchNormalType& patchNormal = pointPatchNormal[bpI];

                forAllRow(bpAtProcs, bpI, i)
                {
                    const label neiProc = bpAtProcs(bpI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    forAllConstIter(patchNormalType, patchNormal, pIt)
                        exchangeData[neiProc].append
                        (
                            refLabelledPointScalar
                            (
                                it.key(),
                                labelledPointScalar
                                (
                                    pIt->first,
                                    pIt->second.first,
                                    pIt->second.second
                                )
                            )
                        );
                }
            }
        }

        LongList<refLabelledPointScalar> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const refLabelledPointScalar& rlps = receivedData[i];
            const label bpI = globalToLocal[rlps.objectLabel()];

            patchNormalType& patchNormal = pointPatchNormal[bpI];

            const labelledPointScalar& lps = rlps.lps();
            patchNormal[lps.pointLabel()].first += lps.coordinates();
            patchNormal[lps.pointLabel()].second += lps.scalarValue();
        }
    }

    //- finally, calculate normal vectors
    # ifdef USE_OMP
    # pragma omp parallel
    # pragma omp single nowait
    # endif
    forAllIter(pointNormalsType, pointPatchNormal, it)
    {
        # ifdef USE_OMP
        # pragma omp task firstprivate(it)
        # endif
        {
            patchNormalType& patchNormal = it->second;

            forAllIter(patchNormalType, patchNormal, pIt)
            {
                pIt->second.first /= (pIt->second.second + VSMALL);
            }
        }
    }
}

void boundaryLayerOptimisation::calculateNormalVectorsSmother
(
    const direction eType,
    pointNormalsType& pointPatchNormal
)
{
    const meshSurfacePartitioner& mPart = surfacePartitioner();
    const meshSurfaceEngine& mse = mPart.surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const labelLongList& bp = mse.bp();

    partTriMesh triMesh(mPart);

    const pointField& triMeshPoints = triMesh.points();
    const VRWGraph& pTriangles = triMesh.pointTriangles();
    const LongList<labelledTri>& triangles = triMesh.triangles();
    const labelLongList& triPointLabel =
        triMesh.meshSurfacePointLabelInTriMesh();
    const labelLongList& surfPointLabel = triMesh.pointLabelInMeshSurface();

    Info << "Calculating normals using smoother " << endl;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdgeType_, heI)
    {
        if( !(hairEdgeType_[heI] & eType) )
            continue;

        const edge& he = hairEdges_[heI];
        const label bpI = bp[he.start()];

        const label triPointI = triPointLabel[bpI];

        //- create an entry in a map
        patchNormalType* patchNormalPtr(NULL);
        # ifdef USE_OMP
        # pragma omp critical
            patchNormalPtr = &pointPatchNormal[bpI];
        # else
        patchNormalPtr = &pointPatchNormal[bpI];
        # endif

        patchNormalType& patchNormal = *patchNormalPtr;

        //- find patches at this point
        DynList<label> patchesAtPoint;
        forAllRow(pTriangles, triPointI, ptI)
        {
            patchesAtPoint.appendIfNotIn
            (
                triangles[pTriangles(triPointI, ptI)].region()
            );
        }

        forAll(patchesAtPoint, ptchI)
        {
            const label patchI = patchesAtPoint[ptchI];

            DynList<point, 128> pts(2);
            DynList<partTet, 256> tets;

            //- create points
            pts[0] = points[he.start()];
            pts[1] = points[he.end()];

            Map<label> bpToSimplex;
            bpToSimplex.insert(bpI, 0);

            forAllRow(pTriangles, triPointI, ptI)
            {
                const labelledTri& tri = triangles[pTriangles(triPointI, ptI)];

                if( tri.region() == patchI )
                {
                    //- create points originating from triangles
                    FixedList<label, 3> triLabels;
                    forAll(tri, pI)
                    {
                        const label spLabel = tri[pI];
                        const label bpLabel = surfPointLabel[spLabel];

                        if( bpLabel < 0 )
                        {
                            triLabels[pI] = pts.size();
                            pts.append(triMeshPoints[spLabel]);
                            continue;
                        }

                        if( !bpToSimplex.found(bpLabel) )
                        {
                            bpToSimplex.insert(bpLabel, pts.size());
                            pts.append(triMeshPoints[spLabel]);
                        }

                        triLabels[pI] = bpToSimplex[bpLabel];
                    }

                    //- create a new tetrahedron
                    tets.append
                    (
                        partTet(triLabels[2], triLabels[1], triLabels[0], 1)
                    );
                }
            }

            //- create partTeMeshSimplex
            partTetMeshSimplex simplex
            (
                pts,
                tets,
                DynList<labelledTri, 32>(),
                1
            );

            //- activate volume optimizer
            volumeOptimizer vOpt(simplex);

            vOpt.optimizeNodePosition();

            const point newP = simplex.centrePoint();

            vector n = -1.0 * (newP - pts[0]);
            const scalar magN = (mag(n) + VSMALL);

            patchNormal[patchI].first = (n / magN);
            patchNormal[patchI].second = magN;
        }
    }

    Info << "Finished calculating normals using smoother " << endl;
}

void boundaryLayerOptimisation::calculateHairVectorsAtTheBoundary
(
    vectorLongList& hairVecs
)
{
    //- set the size of hairVecs
    hairVecs.setSize(hairEdges_.size());

    const pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = meshSurface();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelLongList& bp = mse.bp();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const edgeLongList& edges = mse.edges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        const direction hairType = hairEdgeType_[hairEdgeI];

        if( hairType & BOUNDARY )
        {
            const edge& he = hairEdges_[hairEdgeI];
            vector& hv = hairVecs[hairEdgeI];
            hv = vector::zero;

            if( (hairType & FEATUREEDGE) || (nEndPointsAtPoint_[he.end()] > 1) )
            {
                //- do not modify hair vectors at feature and constrained edges
                hv = he.vec(points);
            }
            else if( hairType & (ATEDGE|ATCORNER) )
            {
                //- this is a case of O-layer at a corner or feature edge

                //- find the surface edges corresponding to the hair edge
                label beI(-1);

                const label bps = bp[he.start()];
                forAllRow(bpEdges, bps, bpeI)
                {
                    const label beJ = bpEdges(bps, bpeI);

                    if( edges[beJ] == he )
                    {
                        beI = beJ;
                        continue;
                    }
                }

                if( beI < 0 )
                    FatalErrorIn
                    (
                        "boundaryLayerOptimisation::"
                        "calculateHairVectorsAtTheBoundary(vectorField&)"
                    ) << "Cannot find hair edge "
                      << hairEdgeI << abort(FatalError);

                //- find the vector at the same angle from both feature edges
                forAllRow(edgeFaces, beI, befI)
                {
                    const face& bf = bFaces[edgeFaces(beI, befI)];
                    const vector fNormal = bf.normal(points);

                    const label pos = bf.which(he.start());

                    if( pos < 0 )
                        FatalErrorIn
                        (
                            "boundaryLayerOptimisation::"
                            "calculateHairVectorsAtTheBoundary(vectorField&)"
                        ) << "Cannot find hair edge "
                          << hairEdgeI << " in face " << bf
                          << abort(FatalError);

                    if( he.end() == bf.prevLabel(pos) )
                    {
                        const edge fe = bf.faceEdge(pos);
                        const vector ev = fe.vec(points);

                        vector hev = fNormal ^ ev;
                        hev /= (mag(hev) + VSMALL);

                        hv += hev;
                    }
                    else
                    {
                        const edge fe = bf.faceEdge(bf.rcIndex(pos));
                        const vector ev = fe.vec(points);

                        vector hev = fNormal ^ ev;
                        hev /= (mag(hev) + VSMALL);

                        hv += hev;
                    }
                }
            }
            else
            {
                FatalErrorIn
                (
                    "boundaryLayerOptimisation::"
                    "calculateHairVectorsAtTheBoundary(vectorField&)"
                ) << "Invalid hair type " << label(hairType)
                  << abort(FatalError);
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- collect data at inter-processor boundaries
        const Map<label>& globalToLocalBnd =
            mse.globalToLocalBndPointAddressing();

        const polyMeshGenAddressing& addr = mesh_.addressingData();
        const labelLongList& globalPointLabel = addr.globalPointLabel();
        const Map<label>& globalToLocal = addr.globalToLocalPointAddressing();

        //- allocate data for sending
        std::map<label, LongList<refLabelledPoint> > exchangeData;
        forAll(hairEdgeNeiProcs_, i)
            exchangeData[hairEdgeNeiProcs_[i]].clear();

        //- prepare data for sending
        forAllConstIter(Map<label>, globalToLocalBnd, it)
        {
            const label bpI = it();

            forAllRow(hairEdgesAtBndPoint_, bpI, peI)
            {
                const label hairEdgeI = hairEdgesAtBndPoint_(bpI, peI);

                if( hairEdgeType_[hairEdgeI] & BOUNDARY )
                {
                    const edge& he = hairEdges_[hairEdgeI];

                    std::map<label, DynList<label, 3> >::const_iterator eIt =
                        hairEdgeAtProcs_.find(hairEdgeI);

                    //- check if the end point is
                    //- at the inter-processor boundary
                    const label geI = globalPointLabel[he.end()];
                    if( !globalToLocal.found(geI) )
                        continue;

                    if( nEndPointsAtPoint_[he.end()] == 1 )
                    {
                        const DynList<label, 3>& neiProcs = eIt->second;

                        forAll(neiProcs, i)
                        {
                            const label neiProc = neiProcs[i];

                            if( neiProc == Pstream:: myProcNo() )
                                continue;

                            exchangeData[neiProc].append
                            (
                                refLabelledPoint
                                (
                                    it.key(),
                                    labelledPoint
                                    (
                                        globalPointLabel[he.end()],
                                        hairVecs[hairEdgeI]
                                    )
                                )
                            );
                        }
                    }
                }
            }
        }

        //- exchange information between processors
        LongList<refLabelledPoint> receiveData;
        help::exchangeMap(exchangeData, receiveData);

        //- update local information
        forAll(receiveData, i)
        {
            const refLabelledPoint& rlp = receiveData[i];

            const label bpI = globalToLocalBnd[rlp.objectLabel()];

            const label ge = globalToLocal[rlp.lPoint().pointLabel()];

            forAllRow(hairEdgesAtBndPoint_, bpI, j)
            {
                const label hairEdgeJ = hairEdgesAtBndPoint_(bpI, j);

                if( hairEdges_[hairEdgeJ].end() == ge )
                {
                    hairVecs[hairEdgeJ] += rlp.lPoint().coordinates();
                }
            }
        }
    }

    //- calculate new normal vectors
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairVecs, hairEdgeI)
    {
        if( hairEdgeType_[hairEdgeI] & BOUNDARY )
            hairVecs[hairEdgeI] /= (mag(hairVecs[hairEdgeI]) + VSMALL);
    }

    # ifdef DEBUGLayer
    Info << "Saving bnd hair vectors" << endl;
    writeHairEdges("bndHairVectors.vtk", (BOUNDARY | ATEDGE), hairVecs);
    # endif
}

void boundaryLayerOptimisation::optimiseHairsNearCorners()
{
    const labelHashSet& featureEdges = surfacePartitioner().featureEdges();
    const meshSurfaceEngine& mse = meshSurface();
    const labelLongList& bp = mse.bp();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const labelLongList& bPoints = mse.boundaryPoints();
    const edgeLongList& edges = mse.edges();

    const pointFieldPMG& points = mesh_.points();

    std::map<label, DynList<label, 2> > neiPoints;
    meshSurfaceEngineModifier surfaceModifier(mse);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100)
    # endif
    forAll(bpEdges, bpI)
    {
        if( boundaryPointType_[bpI] & FLOATINGATEDGE )
        {
            DynList<label, 2> otherPoints;

            forAllRow(bpEdges, bpI, bpeI)
            {
                const label beI = bpEdges(bpI, bpeI);

                if( featureEdges.found(beI) )
                {
                    const edge& e = edges[beI];
                    otherPoints.append(bp[e.otherVertex(bPoints[bpI])]);
                }
            }

            if( otherPoints.size() == 2 )
            {
                //- both edges are local to this processor
                const point& p = points[bPoints[bpI]];
                const point& p0 = points[bPoints[otherPoints[0]]];
                const point& p1 = points[bPoints[otherPoints[1]]];

                const vector d0 = p - p0;
                const vector d1 = p1 - p;

                if( (d0 & d1) > VSMALL )
                {
                    //- return the point in the middle of the cubic spline
                    const vector a = 2.0 * (p0 - p1) + d0 + d1;
                    const vector b = -3.0 * p0 + 3.0 * p1 - 2.0 * d0 - d1;

                    const point newP = 0.125 * a + 0.25 * b + 0.5 * d0 + p0;
                    surfaceModifier.moveBoundaryVertexNoUpdate(bpI, newP);
                }
                else
                {
                    //- the mesh is tangled here. Revert to simple
                    //- linear average that untangles the feature edge
                    const point newP = 0.5 * (p0 + p1);
                    surfaceModifier.moveBoundaryVertexNoUpdate(bpI, newP);
                }
            }
            else if( otherPoints.size() == 1 )
            {
                # ifdef USE_OMP
                # pragma omp critical(collectFloatingPoints)
                # endif
                {
                    neiPoints[bpI] = otherPoints;
                }
            }
            else
            {
                FatalErrorIn
                (
                    "void boundaryLayerOptimisation::optimiseHairsNearCorners()"
                ) << "Invalid number of edges at point" << abort(FatalError);
            }
        }
    }

    if( Pstream::parRun() )
    {
        const labelLongList& globalPointLabel = mse.globalBoundaryPointLabel();
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();

        //- allocate map for exchanging of information
        std::map<label, LongList<refLabelledPoint> > exchangeData;
        forAll(mse.bpNeiProcs(), i)
            exchangeData[mse.bpNeiProcs()[i]].clear();

        //- the following map is used for storing information about
        //- neighbouring vertices of each vertex at inter-processor boundaries
        std::map<label, std::map<label, point> > neiPointCoordinates;

        //- prepare data that shall be exchanged
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            std::map<label, DynList<label, 2> >::const_iterator pIt =
                neiPoints.find(bpI);
            if( pIt == neiPoints.end() )
                continue;

            const DynList<label, 2>& np = pIt->second;

            neiPointCoordinates[bpI][globalPointLabel[np[0]]] =
                points[bPoints[np[0]]];

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                LongList<refLabelledPoint>& dts = exchangeData[neiProc];
                dts.append
                (
                    refLabelledPoint
                    (
                        it.key(),
                        labelledPoint
                        (
                            globalPointLabel[np[0]],
                            points[bPoints[np[0]]]
                        )
                    )
                );
            }
        }

        //- exchange data
        LongList<refLabelledPoint> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        //- perform reduction
        forAll(receivedData, i)
        {
            const refLabelledPoint& rlp = receivedData[i];

            std::map<label, point>& neiPts =
                neiPointCoordinates[globalToLocal[rlp.objectLabel()]];
            neiPts[rlp.lPoint().pointLabel()] =
                rlp.lPoint().coordinates();
        }

        //- calculate new coordinates
        for
        (
            std::map<label, std::map<label, point> >::const_iterator it=
                neiPointCoordinates.begin();
            it!=neiPointCoordinates.end();
            ++it
        )
        {
            const label bpI = it->first;
            const std::map<label, point>& neiPts = it->second;

            if( neiPts.size() == 2 )
            {
                const point& p = points[bPoints[bpI]];
                std::map<label, point>::const_iterator pIt = neiPts.begin();
                const point& p0 = pIt->second;
                const point& p1 = (++pIt)->second;

                const vector d0 = p - p0;
                const vector d1 = p1 - p;

                if( (d0 & d1) > VSMALL )
                {
                    //- return the point in the middle of the cubic spline
                    const vector a = 2.0 * (p0 - p1) + d0 + d1;
                    const vector b = -3.0 * p0 + 3.0 * p1 - 2.0 * d0 - d1;

                    const point newP = 0.125 * a + 0.25 * b + 0.5 * d0 + p0;
                    surfaceModifier.moveBoundaryVertexNoUpdate(bpI, newP);
                }
                else
                {
                    //- the mesh is tangled here. Revert to simple
                    //- linear average that untangles the feature edge
                    const point newP = 0.5 * (p0 + p1);
                    surfaceModifier.moveBoundaryVertexNoUpdate(bpI, newP);
                }
            }
            else
            {
                FatalErrorIn
                (
                    "void boundaryLayerOptimisation::"
                    "optimiseHairsNearCorners()"
                ) << "The number of neighbour points is not 2"
                  << abort(FatalError);
            }
        }
    }

    surfaceModifier.updateGeometry();
}

bool boundaryLayerOptimisation::optimiseConstrainedEdges()
{
    const pointFieldPMG& points = mesh_.points();
    const meshSurfaceEngine& mse = meshSurface();
    const labelLongList& bp = mse.bp();

    bool modified(false);

    polyMeshGenModifier meshModifier(mesh_);

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            for
            (
                std::map<label, label>::const_iterator it=
                    nEndPointsAtPoint_.begin();
                it!=nEndPointsAtPoint_.end();
                ++it
            )
            {
                if( (it->second == 2) && (bp[it->first] < 0) )
                {
                    # ifdef USE_OMP
                    # pragma omp task shared(meshModifier) firstprivate(it)
                    # endif
                    {
                        //- point has been extruded from a point on an edge
                        const DynList<label, 3>& hEdges =
                            hairEndPointsAtPoint_[it->first];

                        if( hEdges.size() == 2 )
                        {
                            //- find pairs of hair edges parallel to each other
                            FixedList<labelPair, 2> parallelEdges;
                            parallelEdges = labelPair(-1, -1);

                            forAll(hEdges, i)
                            {
                                const label heI = hEdges[i];

                                const edge& refHe = hairEdges_[heI];

                                parallelEdges[i].second() = heI;

                                const label heJ =
                                    hairEndPointsAtPoint_[refHe.start()][0];
                                parallelEdges[(i+1)%2].first() = heJ;
                            }

                            //- parallel edges are calculated every time
                            //- to avoid confusion at boundary edges
                            //- where two layers meet
                            if( !constraintsCalculated_ )
                            {
                                # ifdef USE_OMP
                                # pragma omp critical(constrainedHairs)
                                # endif
                                {
                                    //- point is modified at this processor
                                    hairEndPointOwnedByProc_.insert
                                    (
                                        it->first
                                    );

                                    forAll(parallelEdges, i)
                                    {
                                        const label heI =
                                            parallelEdges[i].second();
                                        const label heJ =
                                            parallelEdges[i].first();
                                        constrainedHairs_[heI].appendIfNotIn
                                        (
                                            heJ
                                        );
                                        constrainedHairs_[heJ].appendIfNotIn
                                        (
                                            heI
                                        );
                                    }
                                }
                            }

                            //- find the shortest distance over
                            //- each pair of hairs
                            FixedList<scalar, 2> distance;
                            forAll(parallelEdges, i)
                            {
                                const label pe0 = parallelEdges[i][0];
                                const label pe1 = parallelEdges[i][1];
                                const scalar d0 = hairEdges_[pe0].mag(points);
                                const scalar d1 = hairEdges_[pe1].mag(points);

                                distance[i] = min(d0, d1);
                                if( thinnedHairEdge_[pe0] )
                                {
                                    thinnedHairEdge_[pe1] = true;
                                    modified = true;
                                }
                                if( thinnedHairEdge_[pe1] )
                                {
                                    thinnedHairEdge_[pe0] = true;
                                    modified = true;
                                }
                            }

                            //- calculate directions of hair edges
                            //- at the boundary
                            const edge& he0 = hairEdges_[parallelEdges[0][0]];
                            vector n0 = he0.vec(points);
                            n0 /= (mag(n0) + VSMALL);

                            const edge& he1 = hairEdges_[parallelEdges[1][0]];
                            vector n1 = he1.vec(points);
                            n1 /=(mag(n1) + VSMALL);

                            //- move the vertices
                            const point& origin = points[he0.start()];
                            meshModifier.movePoint
                            (
                                he0.end(),
                                origin + n0 * distance[0]
                            );
                            meshModifier.movePoint
                            (
                                he1.end(),
                                origin + n1 * distance[1]
                            );
                            meshModifier.movePoint
                            (
                                it->first,
                                origin + n0 * distance[0] + n1 * distance[1]
                            );
                        }
                    }
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            for
            (
                std::map<label, label>::const_iterator it=
                    nEndPointsAtPoint_.begin();
                it!=nEndPointsAtPoint_.end();
                ++it
            )
            {
                if( it->second == 3 && (bp[it->first] < 0) )
                {
                    if( Pstream::parRun() )
                    {
                        const DynList<label, 3>& pHairs =
                            hairEndPointsAtPoint_[it->first];

                        //- check if the current processor is the owner
                        //- of the cell extruded from a corner point
                        if( pHairs.size() < 3 )
                            continue;
                    }

                    # ifdef USE_OMP
                    # pragma omp task shared(meshModifier) firstprivate(it)
                    # endif
                    {
                        //- point has been extruded from a corner
                        const DynList<label, 3>& hEdges =
                        hairEndPointsAtPoint_[it->first];

                        FixedList<label, 12> heAtPos(-1);
                        const label s0 = hairEdges_[hEdges[0]].start();
                        heAtPos[0] = hEdges[0];
                        const label s1 = hairEdges_[hEdges[1]].start();
                        heAtPos[4] = hEdges[1];
                        const label s2 = hairEdges_[hEdges[2]].start();
                        heAtPos[8] = hEdges[2];

                        label ss01(-1);
                        if
                        (
                            (
                                hairEndPointsAtPoint_.find(s0) !=
                                hairEndPointsAtPoint_.end()
                            ) &&
                            (
                                hairEndPointsAtPoint_.find(s1) !=
                                hairEndPointsAtPoint_.end()
                            )
                        )
                        {
                            forAll(hairEndPointsAtPoint_[s0], i)
                            {
                                const label heI = hairEndPointsAtPoint_[s0][i];
                                const edge& he = hairEdges_[heI];

                                forAll(hairEndPointsAtPoint_[s1], j)
                                {
                                    const label heJ =
                                        hairEndPointsAtPoint_[s1][j];
                                    const edge& nhe = hairEdges_[heJ];

                                    if( he.start() == nhe.start() )
                                    {
                                        heAtPos[1] = heJ;
                                        heAtPos[5] = heI;
                                        ss01 = he.start();
                                    }
                                }
                            }
                        }

                        label ss12(-1);
                        if
                        (
                            (
                                hairEndPointsAtPoint_.find(s1) !=
                                hairEndPointsAtPoint_.end()
                            ) &&
                            (
                                hairEndPointsAtPoint_.find(s2) !=
                                hairEndPointsAtPoint_.end()
                            )
                        )
                        {
                            forAll(hairEndPointsAtPoint_[s1], i)
                            {
                                const label heI = hairEndPointsAtPoint_[s1][i];
                                const edge& he = hairEdges_[heI];

                                forAll(hairEndPointsAtPoint_[s2], j)
                                {
                                    const label heJ =
                                        hairEndPointsAtPoint_[s2][j];
                                    const edge& nhe = hairEdges_[heJ];

                                    if( he.start() == nhe.start() )
                                    {
                                        heAtPos[11] = heI;
                                        heAtPos[7] = heJ;
                                        ss12 = he.start();
                                    }
                                }
                            }
                        }

                        label ss20(-1);
                        if
                        (
                            (
                                hairEndPointsAtPoint_.find(s0) !=
                                hairEndPointsAtPoint_.end()
                            ) &&
                            (
                                hairEndPointsAtPoint_.find(s2) !=
                                hairEndPointsAtPoint_.end()
                            )
                        )
                        {
                            forAll(hairEndPointsAtPoint_[s0], i)
                            {
                                const label heI = hairEndPointsAtPoint_[s0][i];
                                const edge& he = hairEdges_[heI];

                                forAll(hairEndPointsAtPoint_[s2], j)
                                {
                                    const label heJ =
                                        hairEndPointsAtPoint_[s2][j];
                                    const edge& nhe = hairEdges_[heJ];

                                    if( he.start() == nhe.start() )
                                    {
                                        heAtPos[9] = heI;
                                        heAtPos[3] = heJ;
                                        ss20 = he.start();
                                    }
                                }
                            }
                        }

                        label origin(-1);
                        if
                        (
                            (
                                hairEndPointsAtPoint_.find(ss20) !=
                                hairEndPointsAtPoint_.end()
                            ) &&
                            (
                                hairEndPointsAtPoint_.find(ss01) !=
                                hairEndPointsAtPoint_.end()
                            ) &&
                            (
                                hairEndPointsAtPoint_.find(ss12) !=
                                hairEndPointsAtPoint_.end()
                            )
                        )
                        {
                            forAll(hairEndPointsAtPoint_[ss20], i)
                            {
                                const label heI =
                                    hairEndPointsAtPoint_[ss20][i];
                                const edge& he = hairEdges_[heI];

                                forAll(hairEndPointsAtPoint_[ss01], j)
                                {
                                    const label heJ =
                                        hairEndPointsAtPoint_[ss01][j];
                                    const edge& nhe = hairEdges_[heJ];

                                    if( he.start() == nhe.start() )
                                    {
                                        forAll(hairEndPointsAtPoint_[ss12], k)
                                        {
                                            const label heK =
                                                hairEndPointsAtPoint_[ss12][k];
                                            const edge& nnhe = hairEdges_[heK];

                                            if( he.start() == nnhe.start() )
                                            {
                                                heAtPos[2] = heK;
                                                heAtPos[6] = heI;
                                                heAtPos[10] = heJ;
                                                origin = he.start();
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if( origin > -1 )
                        {
                            if( !constraintsCalculated_ )
                            {
                                # ifdef USE_OMP
                                # pragma omp critical(constrainedHairs)
                                # endif
                                {
                                    //- point is modified at this processor
                                    hairEndPointOwnedByProc_.insert(it->first);

                                    forAll(heAtPos, i)
                                    {
                                        DynList<label, 3>& nhe =
                                            constrainedHairs_[heAtPos[i]];
                                        const label shift = 4 * (i / 4);
                                        for(label j=1;j<4;++j)
                                        {
                                            const label pos =
                                                (i-shift+j)%4 + shift;

                                            nhe.appendIfNotIn(heAtPos[pos]);
                                        }
                                    }
                                }
                            }

                            //- calculate distances in each direction
                            FixedList<scalar, 3> distance(VGREAT);
                            forAll(heAtPos, i)
                            {
                                const label dir = i / 4;

                                if( thinnedHairEdge_[heAtPos[i]] )
                                    modified = true;

                                const scalar d =
                                    hairEdges_[heAtPos[i]].mag(points);
                                distance[dir] = min(distance[dir], d);
                            }

                            //- calculate direction vectors in each direction
                            vector n0 = hairEdges_[heAtPos[2]].vec(points);
                            n0 /= (mag(n0)+ VSMALL);
                            vector n1 = hairEdges_[heAtPos[6]].vec(points);
                            n1 /= (mag(n1) + VSMALL);
                            vector n2 = hairEdges_[heAtPos[10]].vec(points);
                            n2 /= (mag(n2) + VSMALL);

                            const point& p = points[origin];
                            meshModifier.movePoint(ss12, p + n0 * distance[0]);
                            meshModifier.movePoint(ss20, p + n1 * distance[1]);
                            meshModifier.movePoint(ss01, p + n2 * distance[2]);
                            meshModifier.movePoint
                            (
                                s2, p + n0 * distance[0] + n1 * distance[1]
                            );
                            meshModifier.movePoint
                            (
                                s1, p + n0 * distance[0] + n2 * distance[2]
                            );
                            meshModifier.movePoint
                            (
                                s0, p + n1 * distance[1] + n2 * distance[2]
                            );
                            meshModifier.movePoint
                            (
                                it->first,
                                p + n0 * distance[0] +
                                n1 * distance[1] + n2 * distance[2]
                            );
                        }
                    }
                }
            }
        }
    }

    constraintsCalculated_ = true;

    if( Pstream::parRun() )
    {
        reduce(constraintsCalculated_, maxOp<bool>());

        const polyMeshGenAddressing& addr = mesh_.addressingData();
        const VRWGraph& pAtProcs = addr.pointAtProcs();
        const Map<label>& globalToLocal = addr.globalToLocalPointAddressing();

        std::map<label, DynList<labelledPoint, 512> > exchangeData;
        forAll(addr.pointNeiProcs(), i)
            exchangeData[addr.pointNeiProcs()[i]].clear();

        //- exchange new coordinates of points with other processors
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label pointI = it();

            std::map<label, label>::const_iterator nIt =
                nEndPointsAtPoint_.find(pointI);
            if( nIt == nEndPointsAtPoint_.end() )
                continue;

            if
            (
                hairEndPointOwnedByProc_.find(pointI) !=
                hairEndPointOwnedByProc_.end()
            )
            {
                //- the current processor owns the point
                //- point has been moved
                forAllRow(pAtProcs, pointI, i)
                {
                    const label neiProc = pAtProcs(pointI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append
                    (
                        labelledPoint(it.key(), points[pointI])
                    );
                }
            }
        }

        LongList<labelledPoint> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const labelledPoint& lp = receivedData[i];

            const label pI = globalToLocal[lp.pointLabel()];

            //- update points that were not moved locally
            if
            (
                hairEndPointOwnedByProc_.find(pI) ==
                hairEndPointOwnedByProc_.end()
            )
                meshModifier.movePoint(pI, lp.coordinates());
        }

        reduce(modified, maxOp<bool>());
    }

    # ifdef DEBUGLayer
    for(label i=0;i<Pstream::nProcs();++i)
    {
        if( i == Pstream::myProcNo() )
        {
            for
            (
                std::map<label, DynList<label, 3> >::const_iterator it =
                    constrainedHairs_.begin();
                it!=constrainedHairs_.end();
                ++it
            )
            {
                const DynList<label, 3>& cHairs = it->second;

                Pout << "Hair edge " << it->first
                     << " length " << hairEdges_[it->first].mag(points)
                     << " direction " << hairEdges_[it->first].vec(points)
                     << " constrained neighbours " << cHairs
                     << endl;
                forAll(cHairs, i)
                {
                    Pout << "Constrained edge " << cHairs[i]
                         << " length " << hairEdges_[cHairs[i]].mag(points)
                         << "direction " << hairEdges_[cHairs[i]].vec(points)
                         << endl;
                }
            }
        }

        returnReduce(1, sumOp<label>());
    }
    # endif

    return modified;
}

void boundaryLayerOptimisation::optimiseHairNormalsAtTheBoundary()
{
    const pointFieldPMG& points = mesh_.points();
    polyMeshGenModifier meshModifier(mesh_);

    //- calculate direction of hair vector based on the surface normal
    const meshSurfaceEngine& mse = meshSurface();
    const labelLongList& bp = mse.bp();
    const VRWGraph& pFaces = mse.pointFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();

    //- calculate hair vectors
    //- they point in the normal direction to the surface
    vectorLongList hairVecs(hairEdges_.size());

    if( reCalculateNormals_ )
    {
        //- calulate new normal vectors
        calculateHairVectorsAtTheBoundary(hairVecs);
    }
    else
    {
        if( nSmoothNormals_ == 0 )
            return;

        //- keep existing hair vectors
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdges_, heI)
            hairVecs[heI] = hairEdges_[heI].vec(points);
    }

    Info << "Smoothing boundary hair vectors" << endl;

    //- smooth the variation of normals to reduce the twisting of faces
    label nIter(0);
    boolList modifiedEdge(hairEdges_.size());

    while( nIter++ < nSmoothNormals_ )
    {
        vectorLongList newNormals(hairVecs.size());

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdgesNearHairEdge_, hairEdgeI)
        {
            modifiedEdge[hairEdgeI] = false;

            vector& newNormal = newNormals[hairEdgeI];
            newNormal = vector::zero;

            const direction eType = hairEdgeType_[hairEdgeI];

            if( eType & BOUNDARY )
            {
                const edge& he = hairEdges_[hairEdgeI];

                if( nEndPointsAtPoint_[he.end()] > 1 )
                {
                    //- do not adjust hairs that share a common end point
                    //- with another hair edge
                    //- their hair vectors have already been adjusted
                    newNormal = hairVecs[hairEdgeI];
                    continue;
                }

                if( eType & (FEATUREEDGE | ATCORNER) )
                {
                    //- hair vectors at feature edges must not be modified
                    newNormal += hairVecs[hairEdgeI];
                }
                else if( eType & ATEDGE )
                {
                    modifiedEdge[hairEdgeI] = true;

                    DynList<label, 2> edgeFaces;
                    const label bps = bp[he.start()];
                    forAllRow(pFaces, bps, pfI)
                    {
                        const label bfI = pFaces(bps, pfI);
                        const face& bf = bFaces[bfI];

                        forAll(bf, eI)
                        {
                            if( bf.faceEdge(eI) == he )
                            {
                                edgeFaces.append(bfI);
                            }
                        }
                    }

                    //- find the best fitting vector
                    //- at the surface of the mesh
                    forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                    {
                        const label hairEdgeJ =
                            hairEdgesNearHairEdge_(hairEdgeI, nheI);

                        //- check if a neighbour hair edge shares a boundary
                        //- face with the current hair edge
                        bool useNeighbour(false);
                        const edge& nhe = hairEdges_[hairEdgeJ];
                        forAll(edgeFaces, efI)
                        {
                            const face& bf = bFaces[edgeFaces[efI]];

                            forAll(bf, eI)
                            {
                                if( bf.faceEdge(eI) == nhe )
                                {
                                    useNeighbour = true;
                                    break;
                                }
                            }
                        }

                        if( useNeighbour )
                            newNormal += hairVecs[hairEdgeJ];
                    }
                }
                else
                {
                    FatalErrorIn
                    (
                        "void boundaryLayerOptimisation::"
                        "optimiseHairNormalsAtTheBoundary(const label)"
                    ) << "Cannot smooth hair with type " << label(eType)
                      << abort(FatalError);
                }
            }
        }

        if( Pstream::parRun() )
        {
            //- collect data at inter-processor boundaries
            const Map<label>& globalToLocalBnd =
                mse.globalToLocalBndPointAddressing();

            const polyMeshGenAddressing& addr = mesh_.addressingData();
            const labelLongList& globalPointLabel = addr.globalPointLabel();
            const Map<label>& globalToLocal =
                addr.globalToLocalPointAddressing();

            //- allocate space
            std::map<label, LongList<refLabelledPoint> > exchangeData;
            forAll(hairEdgeNeiProcs_, i)
                exchangeData[hairEdgeNeiProcs_[i]].clear();

            //- prepare data for sending
            forAllConstIter(Map<label>, globalToLocalBnd, it)
            {
                const label bpI = it();

                forAllRow(hairEdgesAtBndPoint_, bpI, i)
                {
                    const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                    if( !modifiedEdge[hairEdgeI] )
                        continue;

                    std::map<label, DynList<label, 3> >::const_iterator eIt =
                        hairEdgeAtProcs_.find(hairEdgeI);

                    //- skip hair edges that are not at inter-processor boundary
                    if( eIt == hairEdgeAtProcs_.end() )
                        continue;

                    const DynList<label, 3>& eProcs = eIt->second;

                    const edge& he = hairEdges_[hairEdgeI];

                    const label ge = globalPointLabel[he.end()];

                    forAll(eProcs, i)
                    {
                        const label neiProc = eProcs[i];

                        if( neiProc == Pstream::myProcNo() )
                            continue;

                        exchangeData[neiProc].append
                        (
                            refLabelledPoint
                            (
                                it.key(),
                                labelledPoint(ge, newNormals[hairEdgeI])
                            )
                        );
                    }
                }
            }

            LongList<refLabelledPoint> receiveData;
            help::exchangeMap(exchangeData, receiveData);

            forAll(receiveData, i)
            {
                const refLabelledPoint& rlp = receiveData[i];

                const label bpI = globalToLocalBnd[rlp.objectLabel()];

                const label ge = globalToLocal[rlp.lPoint().pointLabel()];

                forAllRow(hairEdgesAtBndPoint_, bpI, j)
                {
                    const label heJ = hairEdgesAtBndPoint_(bpI, j);
                    const edge& he = hairEdges_[heJ];

                    if( he.end() == ge )
                    {
                        newNormals[heJ] += rlp.lPoint().coordinates();
                    }
                }
            }
        }

        //- calculate new normal vectors
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(newNormals, heI)
        {
            if( modifiedEdge[heI] )
            {
                newNormals[heI] /= (mag(newNormals[heI]) + VSMALL);
                newNormals[heI] = 0.5 * (newNormals[heI] + hairVecs[heI]);
                newNormals[heI] /= (mag(newNormals[heI]) + VSMALL);
            }
        }

        //- transfer new hair vectors to the hairVecs list
        hairVecs.transfer(newNormals);

        # ifdef DEBUGLayer
        if( true )
        {
            writeHairEdges
            (
                "bndHairVectors_"+help::scalarToText(nIter)+".vtk",
                (BOUNDARY | ATEDGE),
                hairVecs
            );
        }
        # endif
    }

    //- move vertices to the new locations
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        if( hairEdgeType_[hairEdgeI] & BOUNDARY )
        {
            modifiedEdge[hairEdgeI] = true;
            const edge& he = hairEdges_[hairEdgeI];

            const vector& hv = hairVecs[hairEdgeI];

            meshModifier.movePoint
            (
                he.end(),
                points[he.start()] + hv * he.mag(points)
            );
        }
        else
        {
            modifiedEdge[hairEdgeI] = false;
        }
    }

    //- ensure the same position of points at all processors
    unifyCoordinatesParallel(modifiedEdge);

    Info << "Finished smoothing boundary hair vectors" << endl;
}

void boundaryLayerOptimisation::optimiseHairNormalsInside()
{
    const pointFieldPMG& points = mesh_.points();
    polyMeshGenModifier meshModifier(mesh_);

    //- calculate direction of hair vector based on the surface normal
    const meshSurfaceEngine& mse = meshSurface();
    const labelLongList& bp = mse.bp();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const edgeLongList& edges = mse.edges();

    //- they point in the normal direction to the surface
    vectorLongList hairVecs(hairEdges_.size());
    boolList modifyNormal(hairEdges_.size());

    if( reCalculateNormals_ )
    {
        //- calculate point normals with respect to all patches at a point
        pointNormalsType pointPatchNormal;
        calculateNormalVectors(INSIDE, pointPatchNormal);
        //calculateNormalVectorsSmother(INSIDE, pointPatchNormal);

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdges_, hairEdgeI)
        {
            const direction hairType = hairEdgeType_[hairEdgeI];

            const edge& he = hairEdges_[hairEdgeI];

            if( hairType & INSIDE )
            {
                vector& hv = hairVecs[hairEdgeI];

                if( nEndPointsAtPoint_[he.end()] > 1 )
                {
                    modifyNormal[hairEdgeI] = false;
                    hv = he.vec(points);
                    hv /= (mag(hv) + VSMALL);
                    continue;
                }

                hv = vector::zero;

                const label bpI = bp[hairEdges_[hairEdgeI].start()];

                label counter(0);
                const patchNormalType& patchNormals = pointPatchNormal[bpI];
                forAllConstIter(patchNormalType, patchNormals, pIt)
                {
                    hv -= pIt->second.first;
                    ++counter;
                }

                if( counter == 0 )
                {
                    FatalErrorIn
                    (
                        "void boundaryLayerOptimisation::"
                        "optimiseHairNormalsInside()"
                    ) << "No valid patches for boundary point "
                      << bp[hairEdges_[hairEdgeI].start()] << abort(FatalError);
                }

                hv /= (mag(hv) + VSMALL);
                modifyNormal[hairEdgeI] = true;
            }
            else
            {
                //- initialise boundary hair vectors. They influence internal
                //- hairs connected to them
                vector hvec = he.vec(points);
                hvec /= (mag(hvec) + VSMALL);
                hairVecs[hairEdgeI] = hvec;

                modifyNormal[hairEdgeI] = false;
            }
        }
    }
    else
    {
        if( nSmoothNormals_ == 0 )
            return;

        Info << "Using existing hair vectors" << endl;

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdges_, hairEdgeI)
        {
            //- initialise boundary hair vectors.
            vector hvec = hairEdges_[hairEdgeI].vec(points);
            hvec /= (mag(hvec) + VSMALL);
            hairVecs[hairEdgeI] = hvec;

            modifyNormal[hairEdgeI] = false;
        }
    }

    # ifdef DEBUGLayer
    writeHairEdges("insideHairVectors.vtk", (INSIDE|BOUNDARY), hairVecs);
    # endif

    Info << "Smoothing internal hair vectors" << endl;

    //- smooth the variation of normals to reduce twisting of faces
    label nIter(0);

    while( nIter++ < nSmoothNormals_ )
    {
        vectorLongList newNormals(hairVecs.size());

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(hairEdgesNearHairEdge_, hairEdgeI)
        {
            vector& newNormal = newNormals[hairEdgeI];
            newNormal = vector::zero;

            const direction eType = hairEdgeType_[hairEdgeI];

            const edge& he = hairEdges_[hairEdgeI];
            const vector& heVec = hairVecs[hairEdgeI];

            if( eType & INSIDE )
            {
                if( (eType & ATCORNER) || (nEndPointsAtPoint_[he.end()] > 1) )
                {
                    //- hair vectors at feature edges must not be modified
                    newNormal = hairVecs[hairEdgeI];
                    modifyNormal[hairEdgeI] = false;
                }
                else if( eType & ATEDGE )
                {
                    modifyNormal[hairEdgeI] = true;

                    //- find the normal vector deviating least from
                    //- other normals at the feature edge
                    forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                    {
                        const label hairEdgeJ =
                            hairEdgesNearHairEdge_(hairEdgeI, nheI);

                        const edge& nhe = hairEdges_[hairEdgeJ];
                        const vector& nhVec = hairVecs[hairEdgeJ];

                        vector n = nhVec ^ (points[nhe[0]] - points[he[0]]);
                        n /= (mag(n) + VSMALL);

                        vector newVec = heVec - (heVec & n) * n;
                        newVec /= (mag(newVec) + VSMALL);

                        scalar weight = 1.0;

                        if( Pstream::parRun() )
                        {
                            //- edges at inter-processor boundaries contribute
                            //- at two sides are given weight 0.5
                            const edge be(he[0], nhe[0]);
                            const label bpI = bp[he[0]];

                            forAllRow(bpEdges, bpI, bpeI)
                            {
                                const edge& bndEdge = edges[bpEdges(bpI, bpeI)];

                                if( bndEdge == be )
                                {
                                    weight = 0.5;
                                    break;
                                }
                            }
                        }

                        newNormal += weight * newVec;
                    }
                }
                else
                {
                    modifyNormal[hairEdgeI] = true;

                    //- find the best fitting vector
                    //- at the surface of the mesh
                    forAllRow(hairEdgesNearHairEdge_, hairEdgeI, nheI)
                    {
                        const label hairEdgeJ =
                            hairEdgesNearHairEdge_(hairEdgeI, nheI);

                        scalar weight = 1.0;

                        //- hairs at inter-processor boundaries contribute
                        //- twice and hence are given a weight 0.5
                        if
                        (
                            hairEdgeAtProcs_.find(hairEdgeJ) !=
                            hairEdgeAtProcs_.end()
                        )
                            weight = 0.5;

                        newNormal += weight * hairVecs[hairEdgeJ];
                    }
                }
            }
            else
            {
                //- copy the existing hair vector
                modifyNormal[hairEdgeI] = false;
            }
        }

        if( Pstream::parRun() )
        {
            //- collect data at inter-processor boundaries
            const Map<label>& globalToLocalBnd =
                mse.globalToLocalBndPointAddressing();

            const polyMeshGenAddressing& addr = mesh_.addressingData();
            const labelLongList& globalPointLabel = addr.globalPointLabel();
            const Map<label>& globalToLocal =
                addr.globalToLocalPointAddressing();

            //- allocate space
            std::map<label, LongList<refLabelledPoint> > exchangeData;
            forAll(hairEdgeNeiProcs_, i)
                exchangeData[hairEdgeNeiProcs_[i]].clear();

            //- prepare data for sending
            forAllConstIter(Map<label>, globalToLocalBnd, it)
            {
                const label bpI = it();

                forAllRow(hairEdgesAtBndPoint_, bpI, i)
                {
                    const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                    //- skip edges that are not modified
                    if( !modifyNormal[hairEdgeI] )
                        continue;

                    std::map<label, DynList<label, 3> >::const_iterator hIt =
                        hairEdgeAtProcs_.find(hairEdgeI);

                    //- skip edges that are not at inter-processor boundaries
                    if( hIt == hairEdgeAtProcs_.end() )
                        continue;

                    const DynList<label, 3>& eProcs = hIt->second;

                    const edge& he = hairEdges_[hairEdgeI];

                    const label ge = globalPointLabel[he.end()];

                    forAll(eProcs, i)
                    {
                        const label neiProc = eProcs[i];

                        if( neiProc == Pstream::myProcNo() )
                            continue;

                        exchangeData[neiProc].append
                        (
                            refLabelledPoint
                            (
                                it.key(),
                                labelledPoint(ge, newNormals[hairEdgeI])
                            )
                        );
                    }
                }
            }

            LongList<refLabelledPoint> receiveData;
            help::exchangeMap(exchangeData, receiveData);

            forAll(receiveData, i)
            {
                const refLabelledPoint& rlp = receiveData[i];

                const label bpI = globalToLocalBnd[rlp.objectLabel()];

                const label ge = globalToLocal[rlp.lPoint().pointLabel()];

                forAllRow(hairEdgesAtBndPoint_, bpI, j)
                {
                    const label heJ = hairEdgesAtBndPoint_(bpI, j);
                    const edge& he = hairEdges_[heJ];

                    if( he.end() == ge )
                    {
                        modifyNormal[heJ] = true;
                        newNormals[heJ] += rlp.lPoint().coordinates();
                    }
                }
            }
        }

        //- calculate new normal vectors
        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(newNormals, hairEdgeI)
        {
            if( modifyNormal[hairEdgeI] )
                newNormals[hairEdgeI] /= (mag(newNormals[hairEdgeI]) + VSMALL);
        }

        //- transfer new hair vectors to the hairVecs list
        hairVecs.transfer(newNormals);

        # ifdef DEBUGLayer
        if( true )
        {
            writeHairEdges
            (
                "insideHairVectors_"+help::scalarToText(nIter)+".vtk",
                INSIDE,
                hairVecs
            );
        }
        # endif
    }

    //- move vertices to the new locations
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        if( modifyNormal[hairEdgeI] )
        {
            const edge& he = hairEdges_[hairEdgeI];

            const vector& hv = hairVecs[hairEdgeI];

            meshModifier.movePoint
            (
                he.end(),
                points[he.start()] + hv * he.mag(points)
            );
        }
    }

    //- ensure that points at inter-processor boundaries remain
    //- at the same positions
    unifyCoordinatesParallel(modifyNormal);

    Info << "Finished smoothing internal hair vectors" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
