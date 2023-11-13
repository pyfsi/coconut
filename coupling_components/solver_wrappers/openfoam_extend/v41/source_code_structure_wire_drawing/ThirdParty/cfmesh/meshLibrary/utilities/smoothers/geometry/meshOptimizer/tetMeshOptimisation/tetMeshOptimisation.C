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
#include "tetMeshOptimisation.H"
#include "partTetMesh.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceEngineModifier.H"

#include "partTetMeshSimplex.H"
#include "meshUntangler.H"
#include "volumeOptimizer.H"
#include "volumeOptimizerHeight.H"
#include "knuppMetric.H"

#include "helperFunctions.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

// #define DEBUGOptimisation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

point tetMeshOptimisation::remapPoint
(
    const partTetMeshSimplex& simplex,
    const point& origCentre
) const
{
    //- find boundary faces of the simplex
    const DynList<point, 128>& pts = simplex.pts();
    const DynList<labelledTri, 32>& bndTriangles =
        simplex.boundaryTriangles();

    point p = simplex.centrePoint();

    label nIter(0);
    scalar tol;
    do
    {
        DynList<label> patches;
        DynList<point> nearestPointInPatch;
        DynList<scalar> nearestDistSq;

        forAll(bndTriangles, j)
        {
            const labelledTri& t = bndTriangles[j];

            const triangle<point, point> tri
            (
                origCentre,
                pts[t[1]],
                pts[t[2]]
            );

            const point pp = help::nearestPointOnTheTriangle(tri, p);

            const label pos =
                patches.containsAtPosition(t.region());

            if( pos >= 0 )
            {
                //- check if the point is the nearest one
                const scalar dSq = mag(pp - p);
                if( dSq < nearestDistSq[pos] )
                {
                    nearestDistSq[pos] = dSq;
                    nearestPointInPatch[pos] = pp;
                }
            }
            else
            {
                //- append the nearest point
                patches.append(t.region());
                nearestPointInPatch.append(pp);
                nearestDistSq.append(magSqr(pp - p));
            }
        }

        //- new point position is the average of patch positions
        p = vector::zero;
        forAll(nearestPointInPatch, j)
            p += nearestPointInPatch[j];
        p /= nearestPointInPatch.size();

        if( patches.size() == 1 )
            break;

        scalar maxDistSq(0.0);
        forAll(nearestDistSq, j)
            maxDistSq = max(maxDistSq, nearestDistSq[j]);

        tol = maxDistSq / simplex.maxEdgeLengthSq();

    } while( tol > 1e-4 && (++nIter < 10) );

    return p;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
tetMeshOptimisation::tetMeshOptimisation(partTetMesh& mesh)
:
    tetMesh_(mesh)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

tetMeshOptimisation::~tetMeshOptimisation()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetMeshOptimisation::optimiseUsingKnuppMetric(const label nIterations)
{
    const LongList<direction>& smoothVertex = tetMesh_.smoothVertex();
    const labelLongList& movablePoints = tetMesh_.movableInternalPoints();

    boolList negativeNode(smoothVertex.size());

    label nIter(0), nNegative, nNegativeBefore;

    do
    {
        //- find the number of inverted tets
        nNegative = tetMesh_.findInvertedPoints(negativeNode);

        # ifdef DEBUGOptimisation
        const scalar startKn = omp_get_wtime();
        # endif

        //- smooth the mesh
        List<LongList<labelledPoint> > newPositions;
        # ifdef USE_OMP
        # pragma omp parallel if( smoothVertex.size() > 100 )
        # endif
        {
            # ifdef USE_OMP
            # pragma omp single
            {
                newPositions.setSize(omp_get_num_threads());
            }

            LongList<labelledPoint>& np = newPositions[omp_get_thread_num()];
            # else
            newPositions.setSize(1);
            LongList<labelledPoint>& np = newPositions[0];
            # endif

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 10)
            # endif
            forAll(movablePoints, i)
            {
                const label nodeI = movablePoints[i];

                if( !negativeNode[nodeI] )
                    continue;

                partTetMeshSimplex simplex(tetMesh_, nodeI);
                knuppMetric(simplex).optimizeNodePosition();
                np.append(labelledPoint(nodeI, simplex.centrePoint()));
            }
        }

        # ifdef DEBUGOptimisation
        const scalar startUp = omp_get_wtime();
        Info << "Knupp optimisation time " << (startUp - startKn) << endl;
        # endif

        //- update mesh vertices
        tetMesh_.updateVerticesSMP(newPositions);
        newPositions.clear();

        # ifdef DEBUGOptimisation
        Info << "Updating time " << (omp_get_wtime() - startUp) << endl;
        # endif

        //- check which tets have been repaired
        nNegativeBefore = nNegative;
        nNegative = tetMesh_.findInvertedPoints(negativeNode);

        reduce(nNegative, sumOp<label>());
        if( nNegative == 0 )
            return;

    } while( (nNegative < nNegativeBefore) && (nIter++ < nIterations) );
}

void tetMeshOptimisation::optimiseUsingMeshUntangler(const label nIterations)
{
    const labelLongList& movablePoints = tetMesh_.movableInternalPoints();

    boolList negativeNode;

    label nIter(0), nNegative, nNegativeBefore;

    do
    {
        //- find the number of inverted tets
        nNegative = tetMesh_.findInvertedPoints(negativeNode);

        //- smooth the mesh
        List<LongList<labelledPoint> > newPositions;
        # ifdef USE_OMP
        # pragma omp parallel if( movablePoints.size() > 100 )
        # endif
        {
            # ifdef USE_OMP
            # pragma omp single
            {
                newPositions.setSize(omp_get_num_threads());
            }

            LongList<labelledPoint>& np = newPositions[omp_get_thread_num()];
            # else
            newPositions.setSize(1);
            LongList<labelledPoint>& np = newPositions[0];
            # endif

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 10)
            # endif
            forAll(movablePoints, i)
            {
                const label nodeI = movablePoints[i];

                if( !negativeNode[nodeI] )
                    continue;

                partTetMeshSimplex simplex(tetMesh_, nodeI);
                meshUntangler(simplex).optimizeNodePosition();
                np.append(labelledPoint(nodeI, simplex.centrePoint()));
            }
        }

        //- update mesh vertices
        tetMesh_.updateVerticesSMP(newPositions);
        newPositions.clear();

        //- check which tets have been repaired
        nNegativeBefore = nNegative;
        nNegative = tetMesh_.findInvertedPoints(negativeNode);

        reduce(nNegative, sumOp<label>());
        if( nNegative == 0 )
            return;

    } while( (nNegative < nNegativeBefore) && (++nIter < nIterations) );

}

void tetMeshOptimisation::optimiseUsingVolumeOptimizer(const label nIterations)
{
    const labelLongList& movablePoints = tetMesh_.movableInternalPoints();

    //- use mesh optimizer to improve the result
    for(label i=0;i<nIterations;++i)
    {
        List<LongList<labelledPoint> > newPositions;

        # ifdef USE_OMP
        # pragma omp parallel if( movablePoints.size() > 100 )
        # endif
        {
            # ifdef USE_OMP
            # pragma omp single
            {
                newPositions.setSize(omp_get_num_threads());
            }

            LongList<labelledPoint>& np = newPositions[omp_get_thread_num()];
            # else
            newPositions.setSize(1);
            LongList<labelledPoint>& np = newPositions[0];
            # endif

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 10)
            # endif
            forAll(movablePoints, i)
            {
                const label nodeI = movablePoints[i];

                partTetMeshSimplex simplex(tetMesh_, nodeI);

                volumeOptimizer vOpt(simplex);
                vOpt.optimizeNodePosition(1e-5);

                np.append(labelledPoint(nodeI, simplex.centrePoint()));
            }
        }

        //- update mesh vertices
        tetMesh_.updateVerticesSMP(newPositions);
        newPositions.clear();
    }
}

void tetMeshOptimisation::optimiseUsingHeightOptimizer
(
    const label nIterations,
    const scalar relativeTolerance
)
{
    const LongList<direction>& smoothVertex = tetMesh_.smoothVertex();
    const labelLongList& movablePoints = tetMesh_.movableInternalPoints();

    scalar relaxFactor(0.5);

    boolList activeNode(smoothVertex.size(), true);

    //- use mesh optimizer to improve the result
    for(label i=0;i<nIterations;++i)
    {
        if( Pstream::parRun() )
            unifyNegativePoints(activeNode);

        List<LongList<labelledPoint> > newPositions;

        scalar maxDisp(0.0), avgDisp(0.0);
        label nMoved(0);

        # ifdef USE_OMP
        # pragma omp parallel if( movablePoints.size() > 100 )
        # endif
        {
            # ifdef USE_OMP
            # pragma omp single
            {
                newPositions.setSize(omp_get_num_threads());
            }

            LongList<labelledPoint>& np = newPositions[omp_get_thread_num()];
            # else
            newPositions.setSize(1);
            LongList<labelledPoint>& np = newPositions[0];
            # endif

            scalar localMaxDisp(0.0);

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 10) reduction(+ : nMoved,avgDisp)
            # endif
            forAll(movablePoints, i)
            {
                const label nodeI = movablePoints[i];

                if( !activeNode[nodeI] )
                    continue;

                partTetMeshSimplex simplex(tetMesh_, nodeI);

                const point p = simplex.centrePoint();

                volumeOptimizerHeight vOptH(simplex);
                vOptH.optimizeNodePosition(1e-5);

                const vector disp = simplex.centrePoint() - p;

                const scalar dSq =
                    magSqr(disp) / (simplex.minEdgeLengthSq() + VSMALL);
                localMaxDisp = max(localMaxDisp, dSq);
                avgDisp += dSq;
                ++nMoved;

                if( dSq < relativeTolerance )
                    activeNode[nodeI] = false;

                np.append(labelledPoint(nodeI, p + relaxFactor * disp));
            }

            # ifdef USE_OMP
            # pragma omp critical(maxDisp)
            # endif
            maxDisp = max(localMaxDisp, maxDisp);
        }

        reduce(avgDisp, sumOp<scalar>());
        reduce(nMoved, sumOp<label>());

        # ifdef DEBUGOptimisation
        Info << "Iteration " << i << " : Max displacement "
             << returnReduce(maxDisp, maxOp<scalar>())
             << " average displacement "
             << (avgDisp/nMoved) << endl;
        # endif

        //- update mesh vertices
        tetMesh_.updateVerticesSMP(newPositions);
        newPositions.clear();

        if( returnReduce(maxDisp, maxOp<scalar>()) < relativeTolerance )
            break;
    }
}

void tetMeshOptimisation::optimiseBoundaryVolumeOptimizer
(
    const label nIterations,
    const bool nonShrinking
)
{
    const labelLongList& movableBndPoints = tetMesh_.movableBoundaryPoints();

    for(label iterI=0;iterI<nIterations;++iterI)
    {
        List<LongList<labelledPoint> > newPositions;

        # ifdef USE_OMP
        # pragma omp parallel
        # endif
        {
            # ifdef USE_OMP
            # pragma omp single
            {
                newPositions.setSize(omp_get_num_threads());
            }

            LongList<labelledPoint>& np = newPositions[omp_get_thread_num()];
            # else
            newPositions.setSize(1);
            LongList<labelledPoint>& np = newPositions[0];
            # endif

            # ifdef USE_OMP
            # pragma omp for schedule(guided, 5)
            # endif
            forAll(movableBndPoints, i)
            {
                const label nodeI = movableBndPoints[i];

                partTetMeshSimplex simplex(tetMesh_, nodeI);

                const point origCentre = simplex.centrePoint();

                volumeOptimizer vOpt(simplex);
                vOpt.optimizeNodePosition(1e-5);

                if( nonShrinking )
                {
                    //- find the nearest point at the boundary of the simplex
                    const point p = remapPoint(simplex, origCentre);

                    //- append the new position of the vertex
                    np.append(labelledPoint(nodeI, 0.5 * (p + origCentre)));
                }
                else
                {
                    //- move the vertex without constraining it
                    np.append
                    (
                        labelledPoint
                        (
                            nodeI,
                            0.5 * (origCentre + simplex.centrePoint())
                        )
                    );
                }
            }
        }

        //- update tetMesh
        tetMesh_.updateVerticesSMP(newPositions);
        newPositions.clear();
    }
}

void tetMeshOptimisation::optimiseBoundaryHeightOptimizer
(
    const label nIterations,
    const bool nonShrinking
)
{
    const labelLongList& movableBndPoints = tetMesh_.movableBoundaryPoints();
    const labelLongList& movableIntPoints = tetMesh_.movableInternalPoints();

    const scalar relaxFactor = 0.5;

    for(label iterI=0;iterI<nIterations;++iterI)
    {
        List<LongList<labelledPoint> > newPositions;

        # ifdef USE_OMP
        # pragma omp parallel
        # endif
        {
            # ifdef USE_OMP
            # pragma omp single
            {
                newPositions.setSize(omp_get_num_threads());
            }

            LongList<labelledPoint>& np = newPositions[omp_get_thread_num()];
            # else
            newPositions.setSize(1);
            LongList<labelledPoint>& np = newPositions[0];
            # endif

            # ifdef USE_OMP
            # pragma omp for schedule(guided, 5)
            # endif
            forAll(movableBndPoints, i)
            {
                const label nodeI = movableBndPoints[i];

                partTetMeshSimplex simplex(tetMesh_, nodeI);

                const point p = simplex.centrePoint();

                volumeOptimizerHeight vOpt(simplex);
                vOpt.optimizeNodePosition(1e-5);

                if( nonShrinking )
                {
                    //- find the nearest point at the boundary of the simplex
                    const point rp = remapPoint(simplex, p);

                    //- append the new position of the vertex
                    np.append(labelledPoint(nodeI, p + relaxFactor * (rp - p)));
                }
                else
                {
                    //- move the vertex without constraining it
                    np.append
                    (
                        labelledPoint
                        (
                            nodeI,
                            p + relaxFactor * (simplex.centrePoint() - p)
                        )
                    );
                }
            }

            //- move internal points, too
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 10)
            # endif
            forAll(movableIntPoints, i)
            {
                const label nodeI = movableIntPoints[i];

                partTetMeshSimplex simplex(tetMesh_, nodeI);

                const point p = simplex.centrePoint();

                volumeOptimizerHeight vOptH(simplex);
                vOptH.optimizeNodePosition(1e-5);

                const vector disp = simplex.centrePoint() - p;

                np.append(labelledPoint(nodeI, p + relaxFactor * disp));
            }
        }

        //- update tetMesh
        tetMesh_.updateVerticesSMP(newPositions);
        newPositions.clear();
    }
}

void tetMeshOptimisation::optimiseBoundarySurfaceLaplace
(
    const label nIterations
)
{
    const labelLongList& movableBndPoints = tetMesh_.movableBoundaryPoints();

    for(label iterI=0;iterI<nIterations;++iterI)
    {
        List<LongList<labelledPoint> > newPositions;

        # ifdef USE_OMP
        # pragma omp parallel
        # endif
        {
            # ifdef USE_OMP
            # pragma omp single
            {
                newPositions.setSize(omp_get_num_threads());
            }

            LongList<labelledPoint>& np =
                newPositions[omp_get_thread_num()];
            # else
            newPositions.setSize(1);
            LongList<labelledPoint>& np = newPositions[0];
            # endif

            # ifdef USE_OMP
            # pragma omp for schedule(guided, 5)
            # endif
            forAll(movableBndPoints, i)
            {
                const label nodeI = movableBndPoints[i];

                partTetMeshSimplex simplex(tetMesh_, nodeI);

                const point origCentre = simplex.centrePoint();

                //- find boundary faces of the simplex
                const DynList<point, 128>& pts = simplex.pts();
                const DynList<labelledTri, 32>& bndTriangles =
                    simplex.boundaryTriangles();

                //- calculate new coordinates as the average of triangle centres
                point newP(vector::zero);
                forAll(bndTriangles, tI)
                {
                    const labelledTri& tri = bndTriangles[tI];

                    point c(vector::zero);
                    for(label pI=0;pI<3;++pI)
                        c += pts[tri[pI]];

                    newP += c;
                }

                if( bndTriangles.size() != 0 )
                {
                    newP /= (3.0 * bndTriangles.size());
                    np.append(labelledPoint(nodeI, 0.5 * (newP + origCentre)));
                }
            }
        }

        //- update tetMesh with new vertex positions
        tetMesh_.updateVerticesSMP(newPositions);
        newPositions.clear();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
