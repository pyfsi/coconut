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

#include "meshOctree.H"
#include "triSurf.H"
#include "triSurfacePartitioner.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"
#include "meshSurfacePartitioner.H"
#include "labelledScalar.H"

#include "helperFunctions.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper::findMappingDistanceSquared
(
    const labelLongList& nodesToMap,
    scalarList& mappingDistanceSq
) const
{
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const VRWGraph& pointInFace = surfaceEngine_.pointInFaces();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    const polyMeshGen& mesh = surfaceEngine_.mesh();

    //- generate search distance for corner nodes
    mappingDistanceSq.setSize(nodesToMap.size());
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(nodesToMap, i)
    {
        const label bpI = nodesToMap[i];
        const label pointI = bPoints[bpI];

        mappingDistanceSq[i] = 0.0;

        //- check if the point has its backup
        point p;
        mesh.getOrigPoint(pointI, p);

        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces(bpI, pfI);
            const face& bf = bFaces[bfI];
            const label posI = pointInFace(bpI, pfI);

            for(label pI=1;pI<bf.size();++pI)
            {
                const label fPointI = bf[(posI+pI)%bf.size()];

                //- get the coordinates
                point fp;
                mesh.getOrigPoint(fPointI, fp);

                //- calculate the squared distance between the points
                const scalar dSq = magSqr(fp - p);

                //- find the maximum distance
                mappingDistanceSq[i] = Foam::max(mappingDistanceSq[i], dSq);
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- make sure that corner nodesd at parallel boundaries
        //- have the same range in which they accept the corners
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
        const labelLongList& globalPointLabel =
            surfaceEngine_.globalBoundaryPointLabel();

        //- create the map for exchanging data
        std::map<label, DynList<labelledScalar> > exchangeData;
        const DynList<label>& neiProcs = surfaceEngine_.bpNeiProcs();
        forAll(neiProcs, i)
            exchangeData.insert
            (
                std::make_pair(neiProcs[i], DynList<labelledScalar>())
            );

        Map<label> globalToLocal(nodesToMap.size());

        forAll(nodesToMap, nI)
        {
            const label bpI = nodesToMap[nI];

            if( bpAtProcs.sizeOfRow(bpI) == 0 )
                continue;

            globalToLocal.insert(globalPointLabel[bpI], nI);

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);
                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append
                (
                    labelledScalar(globalPointLabel[bpI], mappingDistanceSq[nI])
                );
            }
        }

        //- exchange data between processors
        LongList<labelledScalar> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        //- select the maximum mapping distance for processor points
        forAll(receivedData, i)
        {
            const labelledScalar& ls = receivedData[i];

            const label nI = globalToLocal[ls.scalarLabel()];

            //- choose the maximum value for the mapping distance
            mappingDistanceSq[nI] = max(mappingDistanceSq[nI], ls.value());
        }
    }
}

void meshSurfaceMapper::limitDisplacement
(
    const label pointI,
    const scalar maxDistSq,
    point& mapPoint,
    scalar& dSq
) const
{
    const point& p = surfaceEngine_.points()[pointI];

    vector disp = mapPoint - p;
    const scalar magSqrDisp = magSqr(disp);

    if( magSqrDisp > maxDistSq )
    {
        disp *= sqrt(maxDistSq / (magSqrDisp + VSMALL));

        mapPoint = p + disp;
        dSq = magSqr(mapPoint - p);
    }
}

scalar meshSurfaceMapper::faceMetricInPatch
(
    const label bfI,
    const label patchI
) const
{
    const face& bf = surfaceEngine_.boundaryFaces()[bfI];

    const pointFieldPMG& points = surfaceEngine_.points();

    const point centre = help::faceCentre(points, bf);
    const vector area = bf.normal(points);

    point projCentre;
    scalar dSq;
    label nt;

    meshOctree_.findNearestSurfacePointInRegion
    (
        projCentre,
        dSq,
        nt,
        patchI,
        centre
    );

    DynList<point> projPoints(bf.size());
    forAll(bf, pI)
    {
        meshOctree_.findNearestSurfacePointInRegion
        (
            projPoints[pI],
            dSq,
            nt,
            patchI,
            points[bf[pI]]
        );
    }

    vector projArea(vector::zero);
    forAll(bf, pI)
    {
        projArea +=
            triPointRef
            (
                projPoints[pI],
                projPoints[bf.fcIndex(pI)],
                projCentre
            ).normal();
    }

    return magSqr(centre - projCentre) + mag(mag(projArea) - mag(area));
}

void meshSurfaceMapper::mapCorners(const labelLongList& nodesToMap)
{
    const triSurfacePartitioner& sPartitioner = surfacePartitioner();
    const labelList& surfCorners = sPartitioner.corners();
    const List<DynList<label> >& cornerPatches = sPartitioner.cornerPatches();

    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& corners = mPart.corners();
    const VRWGraph& pPatches = mPart.pointPatches();

    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    const polyMeshGen& mesh = surfaceEngine_.mesh();

    const triSurf& surf = meshOctree_.surface();
    const pointField& sPoints = surf.points();

    scalarList mappingDistanceSq;
    findMappingDistanceSquared(nodesToMap, mappingDistanceSq);

    //- for every corner in the mesh surface find the nearest corner in the
    //- triSurface
    meshSurfaceEngineModifier& sMod = meshSurfaceModifier();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(nodesToMap, cornerI)
    {
        const label bpI = nodesToMap[cornerI];
        if( !corners.found(bpI) )
            FatalErrorIn
            (
                "meshSurfaceMapper::mapCorners(const labelLongList&)"
            ) << "Trying to map a point that is not a corner"
                << abort(FatalError);

        const label pointI = bPoints[bpI];

        if( mesh.isLockedPoint(pointI) )
            continue;

        const point& p = points[pointI];
        const scalar maxDistSq = mappingDistanceSq[cornerI];

        //- find the nearest position to the given point patches
        const DynList<label> patches = pPatches[bpI];

        //- find the nearest point in all required patches
        DynList<point> nearestAtPatch(patches.size());

        point mapPointApprox(p);
        scalar distSqApprox;

        label iter(0);
        while( iter++ < 20 )
        {
            point newP(vector::zero);
            forAll(patches, patchI)
            {
                point& np = nearestAtPatch[patchI];
                label nt;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    np,
                    distSqApprox,
                    nt,
                    patches[patchI],
                    mapPointApprox
                );

                newP += np;
            }

            newP /= patches.size();

            if( magSqr(newP - mapPointApprox) < 1e-8 * maxDistSq )
                break;

            mapPointApprox = newP;
        }
        distSqApprox = magSqr(mapPointApprox - p);

        //- check whether the nearest points in all required patches
        //- are close to each other
        forAll(nearestAtPatch, i)
        {
            const point& npp = nearestAtPatch[i];

            for(label j=i+1;j<nearestAtPatch.size();++j)
            {
                if( magSqr(npp - nearestAtPatch[j]) > 0.125 * maxDistSq )
                {
                    mapPointApprox = p;
                    distSqApprox = 0.0;

                    break;
                }
            }
        }

        //- find the nearest triSurface corner for the given corner
        scalar distSq(mappingDistanceSq[cornerI]);
        point mapPoint(p);
        bool foundExact(false);
        forAll(surfCorners, scI)
        {
            const label cornerID = surfCorners[scI];
            const point& sp = sPoints[cornerID];

            if( Foam::magSqr(sp - p) < distSq )
            {
                bool store(true);
                const DynList<label>& cPatches = cornerPatches[scI];

                forAll(patches, i)
                {
                    if( !cPatches.contains(patches[i]) )
                    {
                        store = false;
                        break;
                    }
                }

                if( store )
                {
                    mapPoint = sp;
                    distSq = Foam::magSqr(sp - p);
                    foundExact = true;
                }
            }
        }

        //- use the point with the smallest distance
        if( (distSq > distSqApprox) || !foundExact )
        {
            mapPoint = mapPointApprox;
            distSq = distSqApprox;
        }

        //- check if the moving distance is within the limits
        limitDisplacement(pointI, maxDistSq, mapPoint, distSq);

        //- move the point to the nearest corner
        sMod.moveBoundaryVertexNoUpdate(bpI, mapPoint);
    }

    sMod.syncVerticesAtParallelBoundaries(nodesToMap);
    sMod.updateGeometry(nodesToMap);
}

void meshSurfaceMapper::mapEdgeNodes(const labelLongList& nodesToMap)
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    const meshSurfacePartitioner& mPart = meshPartitioner();
    const VRWGraph& pPatches = mPart.pointPatches();

    const polyMeshGen& mesh = surfaceEngine_.mesh();

    //- find mapping distance for selected vertices
    scalarList mappingDistanceSq;
    findMappingDistanceSquared(nodesToMap, mappingDistanceSq);

    const VRWGraph* bpAtProcsPtr(NULL);
    if( Pstream::parRun() )
        bpAtProcsPtr = &surfaceEngine_.bpAtProcs();

    LongList<parMapperHelper> parallelBndNodes;

    meshSurfaceEngineModifier& sMod = meshSurfaceModifier();

    //- map point to the nearest vertex on the surface mesh
    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(nodesToMap, i)
    {
        const label bpI = nodesToMap[i];
        const label pointI = bPoints[bpI];
        const point& p = points[pointI];

        if( mesh.isLockedPoint(pointI) )
            continue;

        //- find patches at this edge point
        const DynList<label> patches = pPatches[bpI];

        const scalar maxDSq = mappingDistanceSq[i];

        //- find approximate position of the vertex on the edge
        DynList<point> nearestAtPatch(patches.size());

        point mapPointApprox(p);
        scalar distSqApprox;
        label iter(0);
        while( iter++ < 20 )
        {
            point newP(vector::zero);

            forAll(patches, patchI)
            {
                point& np = nearestAtPatch[patchI];
                label nt;
                meshOctree_.findNearestSurfacePointInRegion
                (
                    np,
                    distSqApprox,
                    nt,
                    patches[patchI],
                    mapPointApprox
                );

                newP += np;
            }

            newP /= patches.size();

            if( magSqr(newP - mapPointApprox) < 1e-8 * maxDSq )
                break;

            mapPointApprox = newP;
        }
        distSqApprox = magSqr(mapPointApprox - p);

        //- check whether the nearest points in all required patches
        //- are close to each other
        forAll(nearestAtPatch, i)
        {
            const point& npp = nearestAtPatch[i];

            for(label j=i+1;j<nearestAtPatch.size();++j)
            {
                if( magSqr(npp - nearestAtPatch[j]) > 0.125 * maxDSq )
                {
                    mapPointApprox = p;
                    distSqApprox = 0.0;

                    break;
                }
            }
        }

        //- find the nearest vertex on the triSurface feature edge
        point mapPoint(mapPointApprox);
        scalar distSq(distSqApprox);
        label nse;
        const bool foundExact =
            meshOctree_.findNearestEdgePoint(mapPoint, distSq, nse, p, patches);

        //- use the vertex with the smallest mapping distance
        if( (distSq > distSqApprox) || !foundExact )
        {
            mapPoint = mapPointApprox;
            distSq = distSqApprox;
        }

        //- check if the moving distance is within the limits
        limitDisplacement(pointI, maxDSq, mapPoint, distSq);

        //- move the point to the nearest edge vertex
        sMod.moveBoundaryVertexNoUpdate(bpI, mapPoint);

        if( bpAtProcsPtr && bpAtProcsPtr->sizeOfRow(bpI) )
        {
            # ifdef USE_OMP
            # pragma omp critical
            # endif
            {
                parallelBndNodes.append
                (
                    parMapperHelper
                    (
                        mapPoint,
                        distSq,
                        bpI,
                        -1
                    )
                );
            }
        }
    }

    mapToSmallestDistance(parallelBndNodes);

    sMod.updateGeometry(nodesToMap);
}

void meshSurfaceMapper::mapCornersAndEdges()
{
    const meshSurfacePartitioner& mPart = meshPartitioner();
    const labelHashSet& cornerPoints = mPart.corners();
    labelLongList selectedPoints;
    forAllConstIter(labelHashSet, cornerPoints, it)
        selectedPoints.append(it.key());

    mapCorners(selectedPoints);

    selectedPoints.clear();
    const labelHashSet& edgePoints = mPart.edgePoints();
    forAllConstIter(labelHashSet, edgePoints, it)
        selectedPoints.append(it.key());

    mapEdgeNodes(selectedPoints);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
