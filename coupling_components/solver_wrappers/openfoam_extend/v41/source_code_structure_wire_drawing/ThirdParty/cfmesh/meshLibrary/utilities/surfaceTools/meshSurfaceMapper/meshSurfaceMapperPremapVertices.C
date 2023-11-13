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
#include "meshSurfaceEngineModifier.H"
#include "meshSurfaceMapper.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "refLabelledPoint.H"
#include "refLabelledPointScalar.H"
#include "labelledScalar.H"
#include "helperFunctions.H"
#include "meshSurfaceOptimizer.H"
#include "meshSurfaceCheckInvertedVertices.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGMapping

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper::preMapVertices(const label nIterations)
{
    Info << "Smoothing mesh surface before mapping." << endl;

    const polyMeshGen& mesh = surfaceEngine_.mesh();
    const pointFieldPMG& points = mesh.points();

    const labelLongList& boundaryPoints = surfaceEngine_.boundaryPoints();
    const vectorLongList& faceCentres = surfaceEngine_.faceCentres();
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();
    const VRWGraph& pointInFace = surfaceEngine_.pointInFaces();
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

    //- information about vertices at inter-processor boundaries
    const VRWGraph* bpAtProcsPtr = NULL;
    if( Pstream::parRun() )
    {
        bpAtProcsPtr = &surfaceEngine_.bpAtProcs();
    }

    //- calculate squared distances of face centres from the vertices
    //- These distances will be used as weighting factors for the procedure
    //- aimed at smoothing the stairs in the surface
    List<DynList<scalar, 6> > faceCentreDistances(bFaces.size());
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(bFaces, bfI)
    {
        const point& c = faceCentres[bfI];
        const face& bf = bFaces[bfI];

        faceCentreDistances[bfI].setSize(bf.size());

        forAll(bf, pI)
        {
            faceCentreDistances[bfI][pI] = magSqr(points[bf[pI]] - c);
        }
    }

    //- max displacement limit for a single iteration
    //- do not allow a point to move more than half of the shortest edge length
    scalarField maxDispSq(boundaryPoints.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(boundaryPoints, bpI)
    {
        scalar dSq(0.0);

        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            const label posI = pointInFace(bpI, pfI);
            dSq = max(dSq, faceCentreDistances[bfI][posI]);
        }

        maxDispSq[bpI] = dSq;
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();

        std::map<label, LongList<labelledScalar> > exchangeData;
        forAll(surfaceEngine_.bpNeiProcs(), i)
            exchangeData[surfaceEngine_.bpNeiProcs()[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append
                (
                    labelledScalar(it.key(), maxDispSq[bpI])
                );
            }
        }

        LongList<labelledScalar> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const label bpI = globalToLocal[receivedData[i].scalarLabel()];

            maxDispSq[bpI] = max(maxDispSq[bpI], receivedData[i].value());
        }
    }

    //- allocate space for face weights
    List<DynList<scalar> > weights(boundaryPoints.size());
    forAll(weights, bpI)
    {
        weights[bpI].setSize(pointFaces.sizeOfRow(bpI));
    }

    //- use the shrinking laplace first
    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(weights, bpI)
    {
        forAllRow(pointInFace, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            const label pI = pointInFace(bpI, pfI);

            weights[bpI][pfI] =
                (1.0 / (sqrt(faceCentreDistances[bfI][pI]) + VSMALL));
        }
    }

    //- start smoothing out the stairs at the surface of the volume mesh
    //- and moving the suface of the mesh towards the input geometry
    meshSurfaceEngineModifier& surfaceModifier = meshSurfaceModifier();

    for(label iterI=0;iterI<nIterations;++iterI)
    {
        # ifdef DEBUGMapping
        surfaceEngine_.writeSurfaceVTK
        (
            "preMapInitialIter_"+help::labelToText(iterI)+".vtk"
        );
        # endif

        //- find patches in the vicinity of a boundary face
        List<DynList<label> > boundaryPointPatches(boundaryPoints.size());
        findPointPatches(boundaryPointPatches);

        //- find displacement of surface points using a shrinking laplacian
        vectorField dispLaplace;
        calculateMovementShrinkingSurfaceLaplaceFC(weights, dispLaplace);

        //- move points as requested by the smoother
        pointField origVertices(boundaryPoints.size());
        # ifdef USE_OMP
        # pragma omp parallel for schedule(static, 1)
        # endif
        forAll(boundaryPoints, bpI)
        {
            if( surfaceEngine_.mesh().isLockedPoint(boundaryPoints[bpI]) )
                continue;

            origVertices[bpI] = points[boundaryPoints[bpI]];
            const point newP =
                points[boundaryPoints[bpI]] + 0.5 * dispLaplace[bpI];

            surfaceModifier.moveBoundaryVertexNoUpdate(bpI, newP);
        }

        //- project the points to the nearest point at the geometry
        //- and enforce mesh movement in the direction towards the surface mesh
        pointField newPoints;
        boolList isMapped;
        findMappingVertices(boundaryPointPatches, newPoints, isMapped);

        LongList<parMapperHelper> parallelBndNodes;

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(boundaryPoints, bpI)
        {
            if( surfaceEngine_.mesh().isLockedPoint(boundaryPoints[bpI]) )
                continue;

            const point& p = points[boundaryPoints[bpI]];

            vector dispProjection = vector::zero;
            if( isMapped[bpI] )
            {
                //- limit the laplacian displacement to the tangential
                //- direction, only
                dispProjection = newPoints[bpI] - p;
            }

            //- calculate new position
            point newP = (p + dispProjection);
            vector disp = 0.5 * (newP - origVertices[bpI]);
            scalar dispDSq = magSqr(disp);

            //- do not allow displacement larger than maxDisp
            if( dispDSq > maxDispSq[bpI] )
            {
                disp *= sqrt(maxDispSq[bpI] / (dispDSq + VSMALL));
                dispDSq = maxDispSq[bpI];
            }

            //- move the point to the new location
            newP = origVertices[bpI] + disp;

            surfaceModifier.moveBoundaryVertexNoUpdate(bpI, newP);

            if( bpAtProcsPtr && (bpAtProcsPtr->sizeOfRow(bpI) != 0) )
            {
                # ifdef USE_OMP
                # pragma omp critical(mapParallelBndNodes)
                # endif
                {
                    parallelBndNodes.append
                    (
                        parMapperHelper
                        (
                            newP,
                            dispDSq,
                            bpI,
                            -1
                        )
                    );
                }
            }
        }

        //- make sure that the vertices at inter-processor boundaries
        //- are mapped onto the same location
        mapToSmallestDistance(parallelBndNodes);

        //- update the surface geometry
        surfaceModifier.updateGeometry();

        # ifdef DEBUGMapping
        surfaceEngine_.writeSurfaceVTK
        (
            "preMapAfterProjectionIter_"+help::labelToText(iterI)+".vtk"
        );
        # endif

        Info << "Untangling surface vertices" << endl;
        meshSurfaceOptimizer(meshPartitioner()).untangleAndRevertSurface();

        # ifdef DEBUGMapping
        surfaceEngine_.writeSurfaceVTK
        (
            "preMapAfterUntanglingIter_"+help::labelToText(iterI)+".vtk"
        );
        # endif
    }

    meshSurfaceCheckInvertedVertices checkInverted(meshPartitioner());

    if( returnReduce(checkInverted.invertedVertices().size(), sumOp<label>()) )
    {
        polyMeshGen& mesh = const_cast<polyMeshGen&>(surfaceEngine_.mesh());
        const label pId = mesh.addPointSubset("invertedAfterPresmoothing");

        const labelHashSet& inverted = checkInverted.invertedVertices();
        forAllConstIter(labelHashSet, inverted, it)
            mesh.addPointToSubset(pId, it.key());
    }

    Info << "Finished smoothing mesh surface before mapping." << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
