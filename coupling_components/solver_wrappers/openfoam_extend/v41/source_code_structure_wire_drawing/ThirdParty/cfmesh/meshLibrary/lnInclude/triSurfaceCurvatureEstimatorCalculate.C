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

#include "triSurfaceCurvatureEstimator.H"
#include "matrix3D.H"
#include "quadricFitting.H"
#include "helperFunctions.H"
#include "HashSet.H"
#include "boolList.H"
#include "Map.H"
#include "DynList.H"

#include "OFstream.H"

#include <map>
#include <set>

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGCurvatureEstimator

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeSurfaceToVTK
(
    OFstream& file,
    const triSurf& surf
)
{
    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    const pointField& points = surf.points();
    file << "POINTS " << points.size() << " float\n";
    forAll(points, pI)
    {
        const point& p = points[pI];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << ' ';

        if( pI % 5 == 0 )
            file << "\n";
    }

    //- write triangles
    file << "\n";
    file << "\nPOLYGONS " << surf.size() << " " << 4*surf.size() << endl;
    forAll(surf, triI)
    {
        const labelledTri& tri = surf[triI];
        file << 3 << " " << tri[0] << " " << tri[1] << " " << tri[2] << endl;
    }
}

void writeSurfaceToVTK
(
    OFstream& file,
    const triSurf& surf,
    const DynList<label>& triangles,
    const std::set<label>& labels
)
{
    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    std::map<label, label> newPointLabel;
    forAll(triangles, tI)
    {
        const labelledTri& tri = surf[triangles[tI]];

        for(label pI=0;pI<3;++pI)
        {
            newPointLabel[tri[pI]];
        }
    }

    forAllConstIter(std::set<label>, labels, it)
        newPointLabel[*it];

    const pointField& points = surf.points();
    file << "POINTS " << label(newPointLabel.size()) << " float\n";
    label counter(0);
    for
    (
        std::map<label, label>::iterator it=newPointLabel.begin();
        it!=newPointLabel.end();
        ++it
    )
    {
        it->second = counter++;

        const point& p = points[it->first];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << ' ';

        file << nl;
    }

    //- write triangles
    file << "\n";
    file << "\nPOLYGONS " << triangles.size()
         << " " << 4*triangles.size() << endl;
    forAll(triangles, tI)
    {
        const labelledTri& tri = surf[triangles[tI]];

        file << 3
             << " " << newPointLabel[tri[0]]
             << " " << newPointLabel[tri[1]]
             << " " << newPointLabel[tri[2]] << endl;
    }
}

void writeSurfaceToVTK
(
    const word& name,
    const triSurf& surf,
    const scalarField& data
)
{
    OFstream file(name+".vtk");

    writeSurfaceToVTK(file, surf);

    //- write curvature fields
    const pointField& points = surf.points();
    file << "\n";
    file << "\nPOINT_DATA " << points.size() << "\n";

    file << "SCALARS Values double\n";
    file << "LOOKUP_TABLE default\n";
    forAll(points, pI)
    {
        file << data[pI] << " ";

        if( pI && pI % 5 == 0 )
            file << endl;
    }

    file.flush();
}

void writeSurfaceToVTK
(
    const word& name,
    const triSurf& surf,
    const List<DynList<scalar, 1> >& data
)
{
    OFstream file(name+".vtk");

    writeSurfaceToVTK(file, surf);

    //- write curvature fields
    const pointField& points = surf.points();
    file << "\n";
    file << "\nPOINT_DATA " << points.size() << "\n";

    file << "SCALARS Values double\n";
    file << "LOOKUP_TABLE default\n";
    forAll(points, pI)
    {
        scalar avgCurv(0.0);
        forAll(data[pI], i)
            avgCurv += data[pI][i];
        avgCurv /= data[pI].size();

        file << avgCurv << " ";

        if( pI && pI % 5 == 0 )
            file << endl;
    }

    file.flush();
}

void writeSurfaceToVTK
(
    const word& name,
    const triSurf& surf,
    const List<DynList<vector, 1> >& data
)
{
    OFstream file(name+".vtk");

    writeSurfaceToVTK(file, surf);

    //- write curvature fields
    const pointField& points = surf.points();
    file << "\n";
    file << "\nPOINT_DATA " << points.size() << "\n";

    file << "VECTORS Values double\n";
    file << "LOOKUP_TABLE default\n";
    forAll(points, pI)
    {
        const vector& v = data[pI][0];

        file << v[0] << " " << v[1] << " " << v[2] << " ";

        if( pI && pI % 5 == 0 )
            file << endl;
    }

    file.flush();
}

void triSurfaceCurvatureEstimator::calculateEdgeCurvature()
{
    const pointField& points = surface_.points();
    const edgeLongList& edges = surface_.edges();
    const VRWGraph& pointEdges = surface_.pointEdges();
    const VRWGraph& edgeFaces = surface_.edgeFacets();

    edgePointCurvature_.setSize(points.size());
    boolList featureEdge(edges.size());

    List<FixedList<label, 2> > neiPoints(points.size());
    scalarField smoothEdgeCurvature(points.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(edgePointCurvature_, i)
        {
            edgePointCurvature_[i] = 0.0;
            neiPoints[i] = -1;
        }

        //- mark feature edges
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(edgeFaces, eI)
        {
            if( edgeFaces.sizeOfRow(eI) != 2 )
            {
                featureEdge[eI] = true;
                continue;
            }

            if(
                 surface_[edgeFaces(eI, 0)].region() ==
                 surface_[edgeFaces(eI, 1)].region()
            )
            {
                featureEdge[eI] = false;
            }
            else
            {
                featureEdge[eI] = true;
            }
        }

        //- loop through the points and calculate the curvature for points
        //- attached to two feature edges
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 20)
        # endif
        forAll(pointEdges, pI)
        {
            DynList<label> features;
            forAllRow(pointEdges, pI, peI)
            {
                const label edgeI = pointEdges(pI, peI);
                if( featureEdge[edgeI] )
                    features.append(edgeI);
            }

            if( features.size() == 2 )
            {
                //- find neighbour points
                neiPoints[pI][0] = edges[features[0]].otherVertex(pI);
                neiPoints[pI][1] = edges[features[1]].otherVertex(pI);

                //- coordinates of vertices
                const vector p = points[pI];
                const vector p1 = points[neiPoints[pI][0]];
                const vector p2 = points[neiPoints[pI][1]];

                //- sides
                const vector a = p1 - p;
                const vector b = p2 - p;
                const vector c = p2 - p1;

                const scalar curv =
                    mag(2.0 * mag(a ^ c) /
                    ((mag(a) * mag(b) * mag(c)) + VSMALL));

                edgePointCurvature_[pI] = Foam::mag(curv);
            }
        }

        //- smooth edge curvature
        for(label iterI=0;iterI<2;++iterI)
        {
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 20)
            # endif
            forAll(neiPoints, pI)
            {
                smoothEdgeCurvature[pI] = edgePointCurvature_[pI];

                if( neiPoints[pI][0] != -1 && neiPoints[pI][1] != -1 )
                {
                    const scalar curvMeanNei =
                        (
                            edgePointCurvature_[neiPoints[pI][0]] +
                            edgePointCurvature_[neiPoints[pI][1]]
                        ) / 2.0;

                    smoothEdgeCurvature[pI] =
                        0.5 * (edgePointCurvature_[pI] + curvMeanNei);
                }
            }

            # ifdef USE_OMP
            # pragma omp for schedule(static, 1)
            # endif
            forAll(smoothEdgeCurvature, pI)
                edgePointCurvature_[pI] = smoothEdgeCurvature[pI];
        }
    }

    # ifdef DEBUGCurvatureEstimator
    Info << "Max edge curvature " << max(edgePointCurvature_) << endl;
    Info << "min edge curvature " << min(edgePointCurvature_) << endl;

    writeSurfaceToVTK("edgeCurvature.vtk", surface_, edgePointCurvature_);
    # endif
}

void triSurfaceCurvatureEstimator::calculateSurfaceCurvatures()
{
    const VRWGraph& pointTriangles = surface_.pointFacets();

    const pointField& points = surface_.points();

    patchPositions_.setSize(surface_.size());
    gaussianCurvature_.setSize(points.size());
    meanCurvature_.setSize(points.size());
    maxCurvature_.setSize(points.size());
    minCurvature_.setSize(points.size());
    maxCurvatureVector_.setSize(points.size());
    minCurvatureVector_.setSize(points.size());

    List<DynList<label, 4> > pointPatches(points.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(pointTriangles, pointI)
    {
        typedef std::map<label, std::set<label> > labelToSetMap;
        typedef std::map<label, vector> labelToVectorMap;

        //- allocate storage
        std::map<label, DynList<label> > regionTriangles;
        labelToSetMap otherNodesInRegion;
        labelToVectorMap normalsInRegion;

        //- find points near the given point, their regions, and the approximate
        //- surface normal
        forAllRow(pointTriangles, pointI, ptI)
        {
            const label triI = pointTriangles(pointI, ptI);
            const label regionI = surface_[triI].region();

            if( normalsInRegion.find(regionI) == normalsInRegion.end() )
            {
                normalsInRegion[regionI] = vector::zero;
                regionTriangles[regionI].clear();
                otherNodesInRegion[regionI].clear();
            }

            regionTriangles[regionI].append(triI);

            forAll(surface_[triI], tpI)
            {
                const label pI = surface_[triI][tpI];

                if( pointI == pI )
                    continue;

                otherNodesInRegion[regionI].insert(pI);
                normalsInRegion[regionI] += surface_[triI].normal(points);
            }
        }

        //- select and additional layer of points
        forAllConstIter(labelToSetMap, otherNodesInRegion, it)
        {
            const std::set<label>& currLabels = it->second;

            std::set<label> additionalPoints;
            forAllConstIter(std::set<label>, currLabels, lit)
            {
                const label neiPointI = *lit;
                const constRow pTriangles = pointTriangles[neiPointI];

                forAll(pTriangles, ptI)
                {
                    const labelledTri& nTri = surface_[pTriangles[ptI]];

                    if( nTri.region() != it->first )
                        continue;

                    forAll(nTri, pI)
                        additionalPoints.insert(nTri[pI]);
                }
            }

            forAllConstIter(std::set<label>, additionalPoints, aIter)
                otherNodesInRegion[it->first].insert(*aIter);
        }

        //- normalize the normal vector for each region
        forAllIter(labelToVectorMap, normalsInRegion, nit)
            nit->second /= (Foam::mag(nit->second) + VSMALL);

        forAllConstIter(labelToSetMap, otherNodesInRegion, it)
        {
            const std::set<label>& labels = it->second;

            const label regionI = it->first;

            //- store the patches
            pointPatches[pointI].append(regionI);

            //- find the position of point in the patch triangles
            const DynList<label>& rTriangles = regionTriangles[regionI];
            forAll(rTriangles, i)
            {
                const label tI = rTriangles[i];

                label pos(-1);
                forAll(surface_[tI], j)
                    if( surface_[tI][j] == pointI )
                    {
                        pos = j;
                        break;
                    }

                patchPositions_(tI, pos) = gaussianCurvature_[pointI].size();
            }

            //- store point coordinates
            DynList<point> op;

            forAllConstIter(std::set<label>, labels, lit)
                op.append(points[*lit]);

            //- fit the quadric patch to the surface
            quadricFitting qfit(points[pointI], normalsInRegion[it->first], op);

            # ifdef DEBUGCurvatureEstimator
            Info << "Point " << pointI << " in patch " << regionI
                << " normal " << normalsInRegion[it->first]
                //<< " evalution points " << labels
                << " has max curvature " << qfit.maxCurvature()
                << " and min curvature " << qfit.minCurvature() << endl;
            OFstream file
            (
                "point_" +
                help::scalarToText(pointI) +
                "_region_" +
                help::scalarToText(regionI) +
                "_triangles.vtk"
            );
            writeSurfaceToVTK(file, surface_, rTriangles, labels);
            # endif

            //- store curvatures
            if( qfit.isValid() )
            {
                gaussianCurvature_[pointI].append(qfit.gaussianCurvature());
                meanCurvature_[pointI].append(qfit.meanCurvature());
                maxCurvature_[pointI].append(qfit.maxCurvature());
                minCurvature_[pointI].append(qfit.minCurvature());
                maxCurvatureVector_[pointI].append(qfit.maxCurvatureVector());
                minCurvatureVector_[pointI].append(qfit.minCurvatureVector());
            }
            else
            {
                gaussianCurvature_[pointI].append(0.0);
                meanCurvature_[pointI].append(0.0);
                maxCurvature_[pointI].append(0.0);
                minCurvature_[pointI].append(0.0);
                maxCurvatureVector_[pointI].append(vector::zero);
                minCurvatureVector_[pointI].append(vector::zero);
            }
        }
    }

    //- smooth curvatures using weighted Laplace
    Info << "Smoothing curvature" << endl;
    List<DynList<scalar, 1> > smoothMinCurv(points.size());
    List<DynList<scalar, 1> > smoothMaxCurv(points.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        for(label iteration=0;iteration<2;++iteration)
        {
            # ifdef USE_OMP
            # pragma omp for schedule(static, 1)
            # endif
            forAll(pointTriangles, pointI)
            {
                const constRow pTriangles = pointTriangles[pointI];

                //- find neighbouring points for each patch
                std::map<label, DynList<label> > patchNeiPoints;
                forAll(pointPatches[pointI], ppI)
                    patchNeiPoints.insert
                    (
                        std::make_pair
                        (
                            pointPatches[pointI][ppI],
                            DynList<label>()
                        )
                    );

                forAll(pTriangles, ptI)
                {
                    const labelledTri& tri = surface_[pTriangles[ptI]];

                    if
                    (
                        patchNeiPoints.find(tri.region()) ==
                        patchNeiPoints.end()
                    )
                        patchNeiPoints[tri.region()] = DynList<label>();

                    forAll(tri, pI)
                    {
                        const label neiPointI = tri[pI];

                        if( neiPointI == pointI )
                            continue;

                        patchNeiPoints[tri.region()].appendIfNotIn(neiPointI);
                    }
                }

                smoothMinCurv[pointI].setSize(minCurvature_[pointI].size());
                smoothMaxCurv[pointI].setSize(maxCurvature_[pointI].size());

                //- update min curvature for all point patches
                forAll(minCurvature_[pointI], patchI)
                {
                    const label cPatch = pointPatches[pointI][patchI];

                    scalar minCurv(0.0);
                    scalar maxCurv(0.0);

                    const DynList<label>& neiPoints = patchNeiPoints[cPatch];

                    if( neiPoints.size() == 0 )
                    {
                        smoothMinCurv[pointI][patchI] =
                            minCurvature_[pointI][patchI];

                        smoothMaxCurv[pointI][patchI] =
                            maxCurvature_[pointI][patchI];
                    }

                    forAll(neiPoints, i)
                    {
                        const label neiPointI = neiPoints[i];

                        const label pos =
                            pointPatches[neiPointI].containsAtPosition(cPatch);

                        minCurv += minCurvature_[neiPointI][pos];
                        maxCurv += maxCurvature_[neiPointI][pos];
                    }

                    minCurv /= neiPoints.size();
                    maxCurv /= neiPoints.size();

                    //- store the value
                    smoothMinCurv[pointI][patchI] =
                        0.5 * (minCurv + minCurvature_[pointI][patchI]);

                    smoothMaxCurv[pointI][patchI] =
                        0.5 * (maxCurv + maxCurvature_[pointI][patchI]);
                }
            }

            # ifdef USE_OMP
            # pragma omp for schedule(static, 1)
            # endif
            forAll(minCurvature_, pointI)
            {
                forAll(minCurvature_[pointI], i)
                {
                    minCurvature_[pointI][i] = smoothMinCurv[pointI][i];
                    maxCurvature_[pointI][i] = smoothMaxCurv[pointI][i];
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp barrier
        # endif

        //- update Gaussian and mean curvatures
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(minCurvature_, pointI)
        {
            const DynList<scalar, 1>& minCurv = minCurvature_[pointI];
            const DynList<scalar, 1>& maxCurv = maxCurvature_[pointI];

            forAll(minCurv, i)
            {
                gaussianCurvature_[pointI][i] = minCurv[i] * maxCurv[i];
                meanCurvature_[pointI][i] = 0.5 * (minCurv[i] + maxCurv[i]);
            }
        }
    }

    # ifdef DEBUGCurvatureEstimator
    Info << "Finished calculating curvature" << endl;
    word name = "surfaceMeanCurv";
    writeSurfaceToVTK(name, surface_, meanCurvature_);
    writeSurfaceToVTK("surfaceGaussianCurv", surface_, gaussianCurvature_);
    writeSurfaceToVTK("surfaceMaxCurv", surface_, maxCurvature_);
    writeSurfaceToVTK("surfaceMinCurv", surface_, minCurvature_);
    writeSurfaceToVTK("surfaceMaxCurvVec", surface_, maxCurvatureVector_);
    writeSurfaceToVTK("surfaceMinCurvVec", surface_, minCurvatureVector_);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
