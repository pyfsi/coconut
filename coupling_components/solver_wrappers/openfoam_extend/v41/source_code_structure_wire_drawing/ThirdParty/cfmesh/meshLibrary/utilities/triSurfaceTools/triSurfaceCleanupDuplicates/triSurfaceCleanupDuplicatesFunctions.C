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

#include "triSurfaceCleanupDuplicates.H"
#include "triSurfModifier.H"
#include "meshOctree.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceCleanupDuplicates::calculateMergingTolerance()
{
    tolerance_ = VGREAT;

    const pointField& pts = surf_.points();

    // PhilipC, 21-Feb-19: disable as gcc-5.4.0 gives the following error on
    // compilation:
    // ‘Foam::triSurfaceCleanupDuplicates::tolerance_’ is not a variable in
    // clause ‘reduction’
    // # ifdef USE_OMP
    // # pragma omp for schedule(dynamic, 50) reduction(min:tolerance_)
    // # endif
    forAll(surf_, tI)
    {
        const labelledTri& tri = surf_[tI];

        forAll(tri, i)
        {
            const scalar d = mag(pts[tri[i]] - pts[tri[(i+1)%3]]);

            tolerance_ = min(tolerance_, 0.5 * d);
        }
    }
}

bool triSurfaceCleanupDuplicates::checkDuplicateTriangles()
{
    Info << "Checking for existence of duplicate triangles" << endl;

    labelLongList newTriangleLabel(surf_.size(), -1);

    const VRWGraph& pointTriangles = surf_.pointFacets();

    //- check if there exist duplicate triangles
    label counter(0);

    forAll(surf_, triI)
    {
        if( newTriangleLabel[triI] != -1 )
            continue;

        newTriangleLabel[triI] = counter;
        ++counter;

        const labelledTri& tri = surf_[triI];

        forAll(pointTriangles[tri[0]], ptI)
        {
            const label triJ = pointTriangles(tri[0], ptI);

            if( triJ <= triI )
                continue;

            const labelledTri& otherTri = surf_[triJ];

            if( tri == otherTri )
                newTriangleLabel[triJ] = newTriangleLabel[triI];
        }
    }

    Info << "Found " << (newTriangleLabel.size()-counter)
        << " duplicate triangles" << endl;

    //- return if there exist no duplicate triangles
    if( counter == newTriangleLabel.size() )
        return false;

    Info << "Current number of triangles " << surf_.size() << endl;
    Info << "New number of triangles " << counter << endl;

    //- create new list of triangles and store it in the surface mesh
    LongList<labelledTri> newTriangles(counter);

    forAll(newTriangleLabel, triI)
    {
        newTriangles[newTriangleLabel[triI]] = surf_[triI];
    }

    updateTriangleLabels(newTriangleLabel);

    triSurfModifier(surf_).facetsAccess().transfer(newTriangles);
    surf_.updateFacetsSubsets(newTriangleLabel);

    surf_.clearAddressing();
    surf_.clearGeometry();

    return true;
}

bool triSurfaceCleanupDuplicates::mergeDuplicatePoints()
{
    pointField& pts = const_cast<pointField&>(surf_.points());
    labelLongList newPointLabel(surf_.nPoints());
    bool foundDuplicates(false);

    const scalar sqTol = sqr(tolerance_);

    DynList<label> ct, containedLeaves;

    # ifdef USE_OMP
    # pragma omp parallel private(ct,containedLeaves)
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(newPointLabel, pI)
            newPointLabel[pI] = pI;

        //- check if there exist any vertices closer
        //- than the prescribed tolerance
        std::map<label, std::map<scalar, std::set<label> > > localNearPoints;
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 50)
        # endif
        forAll(pts, pI)
        {
            containedLeaves.clear();
            octree_.findLeavesInSphere(pts[pI], tolerance_, containedLeaves);

            const point& p = pts[pI];

            //- find points with indices smaller than the current point
            //- located in the vicinity
            std::set<label> points;

            forAll(containedLeaves, i)
            {
                ct.clear();
                octree_.containedTriangles(containedLeaves[i], ct);

                forAll(ct, ctI)
                {
                    const label triI = newTriangleLabel_[ct[ctI]];

                    if( triI < 0 )
                        continue;

                    const labelledTri& tri = surf_[triI];

                    forAll(tri, i)
                    {
                        if( tri[i] < pI )
                        {
                            points.insert(tri[i]);
                        }
                    }
                }
            }

            //- check if any points in the vicinity is closer to the current
            //- point than the prescribed tolerance
            for
            (
                std::set<label>::const_iterator it=points.begin();
                it!=points.end();
                ++it
            )
            {
                const label pointI = *it;
                const point& op = pts[pointI];

                const scalar dSq = magSqr(p - op);

                if( dSq < sqTol )
                {
                    foundDuplicates = true;
                    localNearPoints[pI][dSq].insert(pointI);
                }
            }
        }

        //- update new label of points. Merge the point with the nearest
        //- point
        for( const auto it : localNearPoints )
        {
            const std::map<scalar, std::set<label> >& nearPts =
                it.second;
            const std::set<label>& pointIndices = nearPts.begin()->second;

            //- merge with the closest point, which has the smallest label
            //- out of all points at the same distance
            newPointLabel[it.first] = *pointIndices.begin();
        }
    }

    //- find if there exist no duplicate points
    if( !foundDuplicates )
        return false;

    //- remove vertices and update node labels
    label counter(0);
    forAll(pts, pI)
        if( newPointLabel[pI] == pI )
        {
            newPointLabel[pI] = counter;
            if( counter < pI )
                pts[counter] = pts[pI];
            ++counter;
        }
        else
        {
            const label origI = newPointLabel[pI];
            newPointLabel[pI] = newPointLabel[origI];
        }

    Info << "Found " << (pts.size() - counter) << " duplicate points" << endl;

    pts.setSize(counter);

    //- remove triangles containing duplicate points
    LongList<labelledTri> newTriangles(surf_.facets());
    labelLongList newTriangleLabel(surf_.size(), -1);

    counter = 0;
    forAll(surf_, triI)
    {
        const labelledTri& tri = surf_[triI];
        const labelledTri newTri
        (
            newPointLabel[tri[0]],
            newPointLabel[tri[1]],
            newPointLabel[tri[2]],
            tri.region()
        );

        bool store(true);
        for(label i=0;i<2;++i)
            for(label j=i+1;j<3;++j)
                if( newTri[i] == newTri[j] )
                {
                    store = false;
                    break;
                }

        if( store )
        {
            newTriangles[counter] = newTri;
            newTriangleLabel[triI] = counter;
            ++counter;
        }
    }

    newTriangles.setSize(counter);

    updateTriangleLabels(newTriangleLabel);

    //- update the surface
    triSurfModifier(surf_).facetsAccess().transfer(newTriangles);
    surf_.updateFacetsSubsets(newTriangleLabel);

    //- update feature edges
    edgeLongList newEdges;
    const edgeLongList& featureEdges = surf_.featureEdges();
    forAll(featureEdges, feI)
    {
        const edge& fe = featureEdges[feI];

        if( fe.start() < 0 || fe.start() >= newPointLabel.size() )
            continue;
        if( fe.end() < 0 || fe.end() >= newPointLabel.size() )
            continue;

        const edge newEdge
        (
            newPointLabel[fe.start()],
            newPointLabel[fe.end()]
        );

        if( newEdge.start() != newEdge.end() )
            newEdges.append(newEdge);
    }

    triSurfModifier(surf_).featureEdgesAccess().transfer(newEdges);

    //- clear old addressing data
    surf_.clearAddressing();
    surf_.clearGeometry();

    return true;
}

void triSurfaceCleanupDuplicates::updateTriangleLabels
(
    const labelLongList& newTriangleLabel
)
{
    //- update addressing between the original triangles and the cleaned mesh
    forAll(newTriangleLabel_, triI)
    {
        const label origI = newTriangleLabel_[triI];

        if( origI >= 0 )
            newTriangleLabel_[triI] = newTriangleLabel[origI];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
