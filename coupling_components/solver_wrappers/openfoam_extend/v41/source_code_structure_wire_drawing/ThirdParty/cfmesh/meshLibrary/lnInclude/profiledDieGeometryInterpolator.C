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

#include "profiledDieGeometryInterpolator.H"
#include "demandDrivenData.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceEngine.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "rollingMillMesh.H"
#include "boundBox.H"
#include "helperFunctions.H"
#include "cubicBSpline.H"
#include "surfaceOptimizerHeight.H"

#include <algorithm>

#include "OFstream.H"

//#define DEBUGDie

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

scalar profiledDieGeometryInterpolator::radiusAtPosition(const scalar x) const
{
    //- return the first element if x is smaller than the x of the first element
    if( x < radiusAtPosition_.begin()->first )
        return radiusAtPosition_.begin()->second;

    for
    (
        auto rIt=radiusAtPosition_.begin();
        rIt!=radiusAtPosition_.end();
    )
    {
        const scalar xc = rIt->first;
        const scalar rc = rIt->second;

        //- get the next element
        ++rIt;

        if( rIt==radiusAtPosition_.end() )
        {
            //- return the radius of the latest point
            return rc;
        }

        if( x >= xc && x <= rIt->first )
        {
            //- x is within the given interval
            //- calculate the radius via linear interpolation
            const scalar alpha = mag(x - xc) / (mag(rIt->first - xc) + VSMALL);

            return alpha * rIt->second + (1.0 - alpha) * rc;
        }
    }

    FatalErrorIn
    (
        "scalar profiledDieGeometryInterpolator::radiusAtPosition"
        "(const scalar x) const"
    ) << "Cannot interppolate average radius at position " << x
      << abort(FatalError);

    return -1.0;
}

scalar profiledDieGeometryInterpolator::outerRadiusAtPosition
(
    const scalar x
) const
{
    const triSurf* axialProfilePtr = dieProfiles_.axialCrossSectionSurface();

    const pointField& pts = axialProfilePtr->points();

    scalar r = 0.0;

    forAll(*axialProfilePtr, tI)
    {
        const labelledTri& tri = axialProfilePtr->operator[](tI);

        forAll(tri, pI)
        {
            const edge e(tri[pI], tri[(pI+1)%3]);

            const point& ps = pts[e.start()];
            const point& pe = pts[e.end()];

            const scalar xs = ps.x();
            const scalar xe = pe.x();

            if( xs <= (x - SMALL) && xe >= (x + SMALL) )
            {
                const point px = ((x-xs) * pe + (xe-x) * ps) / (xe-xs+VSMALL);

                r = max(r, sqrt(sqr(px.y()) + sqr(px.z())));
            }
            else if( xe <= (x - SMALL) && xs >= (x + SMALL) )
            {
                const point px = ((x-xe) * ps + (xs-x) * pe) / (xs-xe+VSMALL);

                r = max(r, sqrt(sqr(px.y()) + sqr(px.z())));
            }

            if( mag(x - xs) < SMALL )
            {
                r = max(r, sqrt(sqr(ps.y()) + sqr(ps.z())));
            }

            if( mag(x - xe) < SMALL )
            {
                r = max(r, sqrt(sqr(pe.y()) + sqr(pe.z())));
            }
        }
    }

    deleteDemandDrivenData(axialProfilePtr);

    return r;
}

void writeProfileToVTK
(
    const pointField& points,
    const edgeLongList& profileEdges,
    const fileName fName
)
{
    OFstream file(fName);

    std::map<label, label> newPointLabel;
    LongList<point> newPts;
    forAll(profileEdges, eI)
    {
        const edge& e = profileEdges[eI];

        forAll(e, i)
        {
            if( newPointLabel.find(e[i]) == newPointLabel.end() )
            {
                newPointLabel[e[i]] = newPts.size();
                newPts.append(points[e[i]]);
            }
        }
    }

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "\nPOINTS " << newPts.size() << " float\n";
    forAll(newPts, pI)
    {
        const point& p = newPts[pI];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
    }

    //- write lines
    file << "\nLINES " << profileEdges.size()
         << " " << 3*profileEdges.size() << nl;
    forAll(profileEdges, eI)
    {
        const edge& e = profileEdges[eI];

        if
        (
            newPointLabel.find(e.start()) == newPointLabel.end() ||
            newPointLabel.find(e.end()) == newPointLabel.end()
        )
        {
            FatalError << "Unknown point in the map" << exit(FatalError);
        }

        file << 2 << " " << newPointLabel[e.start()]
             << " " << newPointLabel[e.end()] << nl;
    }

    file << "\n";
    file.flush();
}

void writeCrossSectionToVTK(const List<List<point> >& points, const fileName f)
{
    const label nRows = points.size();
    const label nInColumn = points[0].size();

    //- write the header
    OFstream file(f);

    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "\nPOINTS " << (nRows * nInColumn) << " float\n";
    forAll(points, rI)
    {
        forAll(points[rI], cI)
        {
            const point& p = points[rI][cI];

            file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
        }
    }

    //- write lines
    file << "\nPOLYGONS " << (nRows-1)*nInColumn
         << " " << 5 * (nRows-1) * nInColumn << nl;
    for(label i=0;i<nRows-1;++i)
    {
        for(label j=0;j<nInColumn;++j)
        {
            file << 4
                 << " " << (i * nInColumn + j)
                 << " " << ((i+1) * nInColumn + j)
                 << " " << ((i+1) * nInColumn + ((j+1)%nInColumn))
                 << " " << (i * nInColumn + ((j+1)%nInColumn))
                 << nl;
        }
    }

    file << "\n";
    file.flush();
}

std::shared_ptr<cubicBSpline>
profiledDieGeometryInterpolator::createCubicSpline(const triSurf& surf) const
{
    const scalar xAvg = boundBox(surf.points()).midpoint().x();

    //- detect all edges with x smaller than xAvg
    edgeLongList profileEdges;

    const pointField& pts = surf.points();

    forAll(surf, tI)
    {
        const labelledTri& tri = surf[tI];

        forAll(tri, pI)
        {
            const edge e(tri[pI], tri[(pI+1)%3]);

            const scalar xMax = max(pts[e.start()], pts[e.end()]).x();
            const scalar xMin = min(pts[e.start()], pts[e.end()]).x();

            if( xMin < xAvg && xMax < xAvg )
            {
                profileEdges.append(e);
            }
        }
    }

    //- create point-edges addressing
    VRWGraph pointEdges;
    pointEdges.reverseAddressing(profileEdges);

    //- search for the point with max y coordinate
    //- it will be used as a reference
    scalar maxY(-VGREAT);
    label startPointId(-1);
    forAll(pointEdges, i)
    {
        if( pointEdges.sizeOfRow(i) )
        {
            if( pointEdges.sizeOfRow(i) != 2 )
            {
                fileName fName = "profileAtPos_"+std::to_string(xAvg)+".vtk";
                writeProfileToVTK(pts, profileEdges, fName);

                pointField ePts(2*pointEdges.sizeOfRow(i));
                edgeLongList badEdges;
                forAllRow(pointEdges, i, j)
                {
                    const edge& e = profileEdges[pointEdges(i, j)];

                    ePts[2*j] = pts[e.start()];
                    ePts[2*j+1] = pts[e.end()];

                    badEdges.append(edge(2*j, 2*j+1));

                    fName = "badEdgesAtPos_"+std::to_string(xAvg)+".vtk";
                    writeProfileToVTK(ePts, badEdges, fName);
                }

                FatalErrorIn
                (
                    "std::shared_ptr<cubicBSpline>"
                    "profiledDieGeometryInterpolator::createCubicSpline"
                    "(const triSurf&) const"
                ) << "Number of profile edges at point " << i << " is not 2."
                  << " Edges at point are " << pointEdges[i]
                  << " Edges of the profile are written to " << fName
                  << exit(FatalError);
            }

            forAllRow(pointEdges, i, j)
            {
                const scalar y = pts[i].y();

                if( y > maxY )
                {
                    maxY = y;
                    startPointId = i;
                }
            }
        }
    }

    //- sort edges into a closed loop of points
    labelLongList sortedPoints;
    sortedPoints.append(startPointId);
    bool finished = false;
    do
    {
        const label lastEl = sortedPoints.size() - 1;

        const label pointI = sortedPoints[lastEl];

        forAllRow(pointEdges, pointI, i)
        {
            const edge& e = profileEdges[pointEdges(pointI, i)];

            const label otherPointI = e.otherVertex(pointI);

            if
            (
                sortedPoints.size() > 1 &&
                sortedPoints[lastEl-1] == otherPointI
            )
                continue;

            if( sortedPoints[0] == otherPointI )
            {
                finished = true;
                continue;
            }

            //- add the point
            sortedPoints.append(otherPointI);
            break;
        }
    } while( !finished );

    if( sortedPoints.size() != profileEdges.size() )
    {
        FatalErrorIn
        (
            "std::shared_ptr<cubicBSpline>"
            "profiledDieGeometryInterpolator::createCubicSpline"
            "(const triSurf&) const"
        ) << "The number of points in the polygon " << sortedPoints.size()
          << " does not match the number of edges " << profileEdges.size()
          << abort(FatalError);
    }

    //- ensure that the normal of the polygon points in the positive x direction
    vector n = vector::zero;

    point pAvg(vector::zero);
    forAll(sortedPoints, i)
        pAvg += pts[sortedPoints[i]];
    pAvg /= sortedPoints.size();

    forAll(sortedPoints, i)
    {
        const triangle<point, point> tri
        (
            pts[sortedPoints[i]],
            pts[sortedPoints[(i+1)%sortedPoints.size()]],
            pAvg
        );

        n += tri.normal();
    }

    //- create points for the spline
    LongList<point> splinePoints;

    const plane ySymm(vector::zero, vector(0.0, 1.0, 0.0));
    const plane zSymm(vector::zero, vector(0.0, 0.0, 1.0));
    point intersection;

    if( n.x() >= 0.0 )
    {
        //- create a point exactly at the symmetry plane
        if
        (
            help::planeIntersectsEdge
            (
                pts[sortedPoints[sortedPoints.size()-1]],
                pts[sortedPoints[0]],
                zSymm,
                intersection
            )
        )
        {
            splinePoints.append(intersection);
        }

        forAll(sortedPoints, i)
        {
            splinePoints.append(pts[sortedPoints[i]]);

            //- check and add vertices at y=0 and z=0 planes
            if
            (
                help::planeIntersectsEdge
                (
                    pts[sortedPoints[i]],
                    pts[sortedPoints[(i+1)%sortedPoints.size()]],
                    zSymm,
                    intersection
                ) ||
                help::planeIntersectsEdge
                (
                    pts[sortedPoints[i]],
                    pts[sortedPoints[(i+1)%sortedPoints.size()]],
                    ySymm,
                    intersection
                )
            )
            {
                splinePoints.append(intersection);
            }
        }
    }
    else
    {
        //- create a point exactly at the symmetry plane
        if
        (
            help::planeIntersectsEdge
            (
                pts[sortedPoints[sortedPoints.size()-1]],
                pts[sortedPoints[0]],
                zSymm,
                intersection
            )
        )
        {
            splinePoints.append(intersection);
        }

        forAllReverse(sortedPoints, i)
        {
            splinePoints.append(pts[sortedPoints[i]]);

            //- check and add vertices at y=0 and z=0 planes
            const label previ = (i-1+sortedPoints.size()) % sortedPoints.size();
            if
            (
                help::planeIntersectsEdge
                (
                    pts[sortedPoints[i]],
                    pts[sortedPoints[previ]],
                    zSymm,
                    intersection
                ) ||
                help::planeIntersectsEdge
                (
                    pts[sortedPoints[i]],
                    pts[sortedPoints[previ]],
                    ySymm,
                    intersection
                )
            )
            {
                splinePoints.append(intersection);
            }
        }
    }

    label maxZPoint(-1), minZPoint(-1), maxYPoint(-1), minYPoint(-1);

    while( true )
    {
        scalar maxz(-VGREAT), minz(VGREAT), maxy(-VGREAT), miny(VGREAT);

        forAll(splinePoints, i)
        {
            const point& p = splinePoints[i];

            if( mag(p.z()) < SMALL )
            {
                if( p.y() < miny )
                {
                    miny = p.y();
                    minYPoint = i;
                }

                if( p.y() > maxy )
                {
                    maxy = p.y();
                    maxYPoint = i;
                }
            }
            else if( mag(p.y()) < SMALL )
            {
                if( p.z() > maxz )
                {
                    maxz = p.z();
                    maxZPoint = i;
                }

                if( p.z() < minz )
                {
                    minz = p.z();
                    minZPoint = i;
                }
            }
        }

        if( maxYPoint == 0 )
            break;

        LongList<point> copyPoints;

        forAll(splinePoints, i)
            copyPoints.append(splinePoints[(i+maxYPoint)%splinePoints.size()]);

        splinePoints = copyPoints;
    }

    # ifdef DEBUGDie
    Info << "Orig spline points " << splinePoints << endl;
    Info << "Orig maxYPoint " << maxYPoint << endl;
    Info << "Orig maxZPoint " << maxZPoint << endl;
    Info << "Orig minYPoint " << minYPoint << endl;
    Info << "Orig minZPoint " << minZPoint << endl;
    # endif

    if( maxYPoint >= maxZPoint )
    {
        FatalErrorIn
        (
            "std::shared_ptr<cubicBSpline>"
            "profiledDieGeometryInterpolator::createCubicSpline"
            "(const triSurf&) const"
        ) << "Invalid ordering in the posYposZ quadrant"
          << abort(FatalError);
    }

    if( maxZPoint >= minYPoint )
    {
        FatalErrorIn
        (
            "std::shared_ptr<cubicBSpline>"
            "profiledDieGeometryInterpolator::createCubicSpline"
            "(const triSurf&) const"
        ) << "Invalid ordering in the negYposZ quadrant"
          << abort(FatalError);
    }

    if( minYPoint >= minZPoint )
    {
        FatalErrorIn
        (
            "std::shared_ptr<cubicBSpline>"
            "profiledDieGeometryInterpolator::createCubicSpline"
            "(const triSurf&) const"
        ) << "Invalid ordering in the negYnegZ quadrant"
          << abort(FatalError);
    }

    //- append the first point to close the spline
    splinePoints.append(splinePoints[0]);

    # ifdef DEBUGDie
    Info << "Spline points " << splinePoints << endl;
    Info << "maxYPoint " << maxYPoint << endl;
    Info << "maxZPoint " << maxZPoint << endl;
    Info << "minYPoint " << minYPoint << endl;
    Info << "minZPoint " << minZPoint << endl;
    # endif

    //- create a parameter for each point spline point
    //- special attention is paid to the points crossing the y=0 and z=0 planes
    //- because their location shall be exactly at these planes in case
    //- symmetry is invoked
    scalarLongList parameters(splinePoints.size());

    //- calculate distance from the zero position
    scalarLongList distances(splinePoints.size());
    scalar sumLength(0.0);
    forAll(parameters, i)
    {
        distances[i] = sumLength;

        sumLength +=
            mag(splinePoints[i] - splinePoints[(i+1)%splinePoints.size()]);
    }

    sumLength = distances[distances.size()-1];

    //- parameters in the range posYposZ quadrant
    for(label i=0;i<maxZPoint;++i)
    {
        parameters[i] =
            (0.25 / (distances[maxZPoint] - distances[0])) *
            (distances[i] - distances[0]);
    }

    //- parameter in the range negYposZ quadrant
    for(label i=maxZPoint;i<minYPoint;++i)
    {
        parameters[i] =
            0.25 +
            (0.25 / (distances[minYPoint] - distances[maxZPoint])) *
            (distances[i] - distances[maxZPoint]);
    }

    //- parameter in the range negYnegZ quadrant
    for(label i=minYPoint;i<minZPoint;++i)
    {
        parameters[i] =
            0.5 +
            (0.25 / (distances[minZPoint] - distances[minYPoint])) *
            (distances[i] - distances[minYPoint]);
    }

    //- parameter in the range posYnegZ quadrant
    for(label i=minZPoint;i<splinePoints.size();++i)
    {
        parameters[i] =
            0.75 +
            (0.25 / (sumLength - distances[minZPoint])) *
            (distances[i] - distances[minZPoint]);
    }

    //- construct the spline
    std::shared_ptr<cubicBSpline> splinePtr =
        std::make_shared<cubicBSpline>
        (
            splinePoints,
            parameters,
            "cubicBSpline"
        );

    return splinePtr;
}

std::shared_ptr<triSurf> profiledDieGeometryInterpolator::interpolateProfiles
(
    const scalar x,
    const triSurf& profile1Surf,
    const std::pair<scalar, std::shared_ptr<cubicBSpline> >& profile1,
    const std::pair<scalar, std::shared_ptr<cubicBSpline> >& profile2
) const
{
    if( profile1.first > profile2.first )
    {
        FatalErrorIn
        (
            "std::shared_ptr<triSurf> profiledDieGeometryInterpolator::"
            "interpolateProfiles(const scalar, "
            "const std::pair<scalar, std::shared_ptr<triSurf> >&, "
            "const std::pair<scalar, std::shared_ptr<triSurf> >&) const"
        ) << "Axial position of profile1 is greater than of profile2."
          << " Please revert them" << abort(FatalError);
    }

    const scalar xmin = profile1.first;
    const scalar xmax = profile2.first;

    const scalar dx = (xmax - xmin) / max(label(xToId_.size()), 1);

    const scalar alpha = mag(x - xmin) / (mag(xmax - xmin) + VSMALL);

    //- create cubic splines for both profiles
    const std::shared_ptr<cubicBSpline> p1Ptr = profile1.second;
    const std::shared_ptr<cubicBSpline> p2Ptr = profile2.second;

    //- get the desired number of points in the profile and
    //- exclude the last point collocated to the first one
    const label nProfilePoints = p1Ptr->numberOfControlPoints() - 1;

    pointField newPoints(2 * nProfilePoints);

    for(label i=0;i<nProfilePoints;++i)
    {
        //- calculate parameter of the spline
        const scalar param = min(max(0.0, p1Ptr->pointParam(i)), 1.0);

        //- evaluate coordinates of a new point
        point newP =
            alpha * p2Ptr->evaluate(param) +
            (1.0 - alpha) * p1Ptr->evaluate(param);

        newP.x() = x;

        //- apply the coordinate to the
        newPoints[i] = newP;
        newP.x() += dx;
        newPoints[i+nProfilePoints] = newP;
    }

    //- create new surface as a copy of the first one
    std::shared_ptr<triSurf> retPtr = std::make_shared<triSurf>();

    //- interpolate new points at a position
    triSurfModifier sMod(*retPtr);
    sMod.pointsAccess().transfer(newPoints);

    if( (newPoints.size() % 2) != 0 )
    {
        FatalErrorIn
        (
            "std::shared_ptr<triSurf> profiledDieGeometryInterpolator::"
            "interpolateProfiles(const scalar, "
            "const std::pair<scalar, std::shared_ptr<triSurf> >&, "
            "const std::pair<scalar, std::shared_ptr<triSurf> >&) const"
        ) << "The number of points in the profile is not divisible by 2"
          << " It is not a 2D profile" << abort(FatalError);
    }

    sMod.patchesAccess() = profile1Surf.patches();

    for(label i=0;i<nProfilePoints;++i)
    {
        retPtr->appendTriangle
        (
            labelledTri
            (
                i,
                (i+1)%nProfilePoints,
                (i+1)%nProfilePoints+nProfilePoints,
                0
            )
        );
        retPtr->appendTriangle
        (
            labelledTri
            (
                i,
                (i+1)%nProfilePoints+nProfilePoints,
                i+nProfilePoints,
                0
            )
        );
    }

    return retPtr;
}

namespace edgePolylines
{

class edgeConnectionsNeighbourOperator
{
    const VRWGraph& faceEdges_;
    const VRWGraph& edgeFaces_;
    const edgeLongList& edges_;
    const VRWGraph& bpEdges_;
    const labelLongList& bp_;

public:

    edgeConnectionsNeighbourOperator(const meshSurfaceEngine& mse)
    :
        faceEdges_(mse.faceEdges()),
        edgeFaces_(mse.edgeFaces()),
        edges_(mse.edges()),
        bpEdges_(mse.boundaryPointEdges()),
        bp_(mse.bp())
    {}

    label size() const
    {
        return edgeFaces_.size();
    }

    void operator()(const label beI, DynList<label>& neiBndEdges) const
    {
        neiBndEdges.clear();

        if( edgeFaces_.sizeOfRow(beI) != 2 )
        {
            FatalErrorIn
            (
                "void edgeConnectionsNeighbourOperator::"
                "operator()(const label, DynList<label>&) const"
            ) << "Boundary edge " << beI << " is not attached"
              << " to 2 adjacent faces. This indicates an error in the "
              << "mesh or that you are running using MPI"
              << abort(FatalError);
        }

        const label bf0 = edgeFaces_(beI, 0);
        const label bf1 = edgeFaces_(beI, 1);

        const edge& be = edges_[beI];

        forAll(be, pI)
        {
            const label bpI = bp_[be[pI]];

            //- the topology of the mesh requires existence of exactly
            //- 4 edges at each surface vertex
            if( bpEdges_.sizeOfRow(bpI) != 4 )
            {
                continue;
            }

            label nNeigboursAtPoint(0);

            forAllRow(bpEdges_, bpI, bpeI)
            {
                const label beJ = bpEdges_(bpI, bpeI);

                //- skip current edge
                if( beJ == beI )
                    continue;
                //- skip edges contained in attached faces
                if
                (
                    faceEdges_.contains(bf0, beJ) ||
                    faceEdges_.contains(bf1, beJ)
                )
                {
                    continue;
                }

                //- this must be the only edge that is a neighbour along the
                //- polyline
                ++nNeigboursAtPoint;
                neiBndEdges.append(beJ);
            }

            if( nNeigboursAtPoint > 1 )
            {
                FatalErrorIn
                (
                    "void edgeConnectionsNeighbourOperator::"
                    "operator()(const label, DynList<label>&) const"
                ) << "More than one neighbour edge found at boundary point "
                  << bpI << abort(FatalError);
            }
        }
    }

    template<class labelListType>
    void collectGroups
    (
        std::map<label, DynList<label> >&,
        const labelListType&,
        const DynList<label>&
    ) const
    {
        notImplemented
        (
            "template<class labelListType>"
            "void collectGroups"
            "(std::map<label, DynList<label> >&, "
            "const labelListType&, const DynList<label>&) const"
        );
    }
};

class edgeConnectionsSelectorOperator
{
    const VRWGraph& edgeFaces_;

    const labelLongList& fPatch_;

    const labelHashSet& skipPatches_;

public:

    edgeConnectionsSelectorOperator
    (
        const meshSurfaceEngine& mse,
        const labelHashSet& skipPatches
    )
    :
        edgeFaces_(mse.edgeFaces()),
        fPatch_(mse.boundaryFacePatches()),
        skipPatches_(skipPatches)
    {}

    bool operator()(const label beI) const
    {
        if( edgeFaces_.sizeOfRow(beI) != 2 )
        {
            FatalErrorIn
            (
                "bool edgeConnectionsSelectorOperator::"
                "operator()(const label) const"
            ) << "Boundary edge " << beI << " is not attached"
              << " to 2 adjacent faces. This indicates an error in the "
              << "mesh or that you are running using MPI"
              << abort(FatalError);
        }

        const label bf0 = edgeFaces_(beI, 0);
        const label bf1 = edgeFaces_(beI, 1);

        //- do not consider edges that do not belong to active patches
        if
        (
            skipPatches_.found(fPatch_[bf0]) &&
            skipPatches_.found(fPatch_[bf1])
        )
        {
            return false;
        }

        return true;
    }
};

}

void profiledDieGeometryInterpolator::detectProfilePointsInVolMesh()
{
    std::set<word> skipNames;
    skipNames.insert
    (
        patchHandler_.patchNameForDie(rollingMillPatchNamesHandler::DIEHOUSING)
    );

    const word upstreamName =
        patchHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIEUPSTREAM
        );
    skipNames.insert(upstreamName);
    const word downstreamName =
        patchHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIEDOWNSTREAM
        );
    skipNames.insert(downstreamName);

    skipNames.insert
    (
        patchHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIESYMMY
        )
    );
    skipNames.insert
    (
        patchHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIESYMMZ
        )
    );
    skipNames.insert
    (
        patchHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIEFRONT
        )
    );
    skipNames.insert
    (
        patchHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIEBACK
        )
    );

    const word symmYName =
        patchHandler_.patchNameForDie(rollingMillPatchNamesHandler::DIESYMMY);
    const word symmZName =
        patchHandler_.patchNameForDie(rollingMillPatchNamesHandler::DIESYMMZ);

    labelHashSet skipPatches, symmYPatches, symmZPatches;
    label upstreamPatchId(-1);
    label downstreamPatchId(-1);
    forAll(mesh_.boundaries(), patchI)
    {
        const boundaryPatch& ptch = mesh_.boundaries()[patchI];

        for(auto it=skipNames.begin();it!=skipNames.end();++it)
        {
            if( ptch.patchName().find(*it) != word::npos )
            {
                skipPatches.insert(patchI);
            }

            if( ptch.patchName().find(downstreamName) != word::npos )
            {
                downstreamPatchId = patchI;
            }

            if( ptch.patchName().find(upstreamName) != word::npos )
            {
                upstreamPatchId = patchI;
            }

            if( ptch.patchName().find(symmYName) != word::npos )
            {
                symmYPatches.insert(patchI);
            }

            if( ptch.patchName().find(symmZName) != word::npos )
            {
                symmZPatches.insert(patchI);
            }
        }
    }

    //- start analysing the surface of the volume mesh with a goal
    //- to detect edges belonging to each radial cross-section
    meshSurfaceEngine mse(mesh_);

    //- start by assigning each surface edge into a group representing
    //- a polyline. These polylines intersect at vertices
    labelLongList edgeInPolyLine;
    const label nPolyLines =
        help::groupMarking
        (
            edgeInPolyLine,
            edgePolylines::edgeConnectionsNeighbourOperator(mse),
            edgePolylines::edgeConnectionsSelectorOperator(mse, skipPatches)
        );

    # ifdef DEBUGDie
    {
        labelList pIds(nPolyLines);
        forAll(pIds, i)
            pIds[i] = mesh_.addPointSubset("polyLine_"+help::labelToText(i));

        const edgeLongList& edges = mse.edges();
        forAll(edgeInPolyLine, eI)
        {
            if( edgeInPolyLine[eI] < 0 )
                continue;

            const edge& be = edges[eI];

            mesh_.addPointToSubset(pIds[edgeInPolyLine[eI]], be.start());
            mesh_.addPointToSubset(pIds[edgeInPolyLine[eI]], be.end());
        }
        mesh_.write();
        Info << "Mesh surf polylines written" << endl;
    }
    # endif

    //- detect points at symmetry and outer boundaries
    const pointFieldPMG& points = mesh_.points();
    const VRWGraph& bpEdges = mse.boundaryPointEdges();
    const VRWGraph& pFaces = mse.pointFaces();
    const labelLongList& fPatch = mse.boundaryFacePatches();
    const faceList::subList& bFaces = mse.boundaryFaces();

    forAll(fPatch, bfI)
    {
        const face& bf = bFaces[bfI];

        if( symmYPatches.found(fPatch[bfI]) )
        {
            //- point is at the y=0 symmetry plane
            forAll(bf, pI)
                pointsInYSymmetry_.insert(bf[pI]);

            //- skip checking if the point is at the outer radius of a section
            continue;
        }

        if( symmZPatches.found(fPatch[bfI]) )
        {
            //- point is at z = 0 symmetry plane
            forAll(bf, pI)
                pointsInZSymmetry_.insert(bf[pI]);

            //- skip checking if the point is at the outer radius of a section
            continue;
        }

        if
        (
            skipPatches.found(fPatch[bfI]) &&
            (
                (mag(help::faceAreaVector(points, bf).y()) > SMALL) ||
                (mag(help::faceAreaVector(points, bf).z()) > SMALL)
            )
        )
        {
            //- point is at the outer radius of its radial cross-section
            forAll(bf, pI)
                pointsAtFixedRadius_.insert(bf[pI]);
        }
    }

    # ifdef DEBUGDie
    {
        const label symmYId = mesh_.addPointSubset("pointsAtYSymm");
        const label symmZId = mesh_.addPointSubset("pointsAtZSymm");
        const label outerRId = mesh_.addPointSubset("pointsAtOuterRadius");

        forAllConstIter(labelHashSet, pointsInYSymmetry_, it)
            mesh_.addPointToSubset(symmYId, it.key());
        forAllConstIter(labelHashSet, pointsInZSymmetry_, it)
            mesh_.addPointToSubset(symmZId, it.key());

        forAllConstIter(labelHashSet, pointsAtFixedRadius_, it)
            mesh_.addPointToSubset(outerRId, it.key());

        mesh_.write();
        Info << "Mesh with point subsets at symmetry and outer radius" << endl;
    }
    # endif

    //- filter out polylines that are not closed or if they end at upstream
    //- or downstream patches
    List<labelHashSet> patchesAtPolyline(nPolyLines);

    forAll(bpEdges, bpI)
    {
        //- detect patches adjacent to this point
        DynList<label> patches;
        forAllRow(pFaces, bpI, pfI)
        {
            patches.appendIfNotIn(fPatch[pFaces(bpI, pfI)]);
        }

        //- assign patches to polylines
        forAllRow(bpEdges, bpI, peI)
        {
            const label beI = bpEdges(bpI, peI);

            const label polylineI = edgeInPolyLine[beI];

            if( polylineI >= 0 )
            {
                forAll(patches, i)
                    patchesAtPolyline[polylineI].insert(patches[i]);
            }
        }
    }

    //- detect active polylines
    Map<label> activePolylines;
    forAll(patchesAtPolyline, polylineI)
    {
        if
        (
            patchesAtPolyline[polylineI].found(upstreamPatchId) &&
            patchesAtPolyline[polylineI].found(downstreamPatchId)
        )
        {
            continue;
        }

        activePolylines.insert(polylineI, activePolylines.size());
    }

    //- detect points in each polyline
    const edgeLongList& edges = mse.edges();

    std::map<label, std::set<label> > pointsInPolyline;
    forAll(edges, beI)
    {
        const label polylineI = edgeInPolyLine[beI];

        if( polylineI >= 0 && activePolylines.found(polylineI) )
        {
            const edge& be = edges[beI];

            pointsInPolyline[polylineI].insert(be.start());
            pointsInPolyline[polylineI].insert(be.end());
        }
    }

    //- sort indices of polylines depending on their x coordinate
    std::map<scalar, label> xToPolylineId;

    for(auto it=pointsInPolyline.begin();it!=pointsInPolyline.end();++it)
    {
        const label polylineI = it->first;

        scalar xMin(VGREAT), xMax(-VGREAT);

        const std::set<label>& pointsInCurrentPolyline = it->second;

        forAllConstIter(std::set<label>, pointsInCurrentPolyline, pIt)
        {
            const point& p = points[*pIt];

            xMin = min(p.x(), xMin);
            xMax = max(p.x(), xMax);
        }

        if( mag(xMax - xMin) < SMALL )
        {
            xToPolylineId[xMin] = polylineI;
        }
        else
        {
            FatalErrorIn
            (
                "void profiledDieGeometryInterpolator::"
                "detectProfilePointsInVolMesh()"
            ) << "Selected polyline does not have a constant x coordinate."
              << " The x coordinate ranges between " << xMin << " and " << xMax
              << abort(FatalError);
        }
    }

    //- copy profile points
    sortedProfilePoints_.setSize(xToPolylineId.size());
    xToId_.clear();

    label counter(0);
    for(auto it=xToPolylineId.begin();it!=xToPolylineId.end();++it)
    {
        const label polylineI = it->second;

        const std::set<label>& pointsInCurrentPolyline =
            pointsInPolyline[polylineI];

        forAllConstIter(std::set<label>, pointsInCurrentPolyline, pIt)
            sortedProfilePoints_.append(counter, *pIt);

        xToId_[it->first] = counter;

        ++counter;
    }

    # ifdef DEBUGDie
    forAll(sortedProfilePoints_, i)
    {
        const label pId =
            mesh_.addPointSubset("pointsInProfile_"+help::labelToText(i));

        forAllRow(sortedProfilePoints_, i, j)
            mesh_.addPointToSubset(pId, sortedProfilePoints_(i, j));
    }
    mesh_.write();
    Info << "Profiles created" << endl;
    # endif
}

void profiledDieGeometryInterpolator::analyseAxialCrossSection()
{
    //- get the reference to the axial profile of a die
    const triSurf* axialCrossSectionPtr =
        dieProfiles_.axialCrossSectionSurface();
    const triSurf& axialCrossSection = *axialCrossSectionPtr;

    const word housingName =
        patchHandler_.patchNameForDie(rollingMillPatchNamesHandler::DIEHOUSING);
    const word upstreamName =
        patchHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIEUPSTREAM
        );
    const word downstreamName =
        patchHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIEDOWNSTREAM
        );

    labelHashSet skipPatches;
    forAll(axialCrossSection.patches(), patchI)
    {
        const geometricSurfacePatch& ptch = axialCrossSection.patches()[patchI];

        if
        (
            ptch.name() == housingName ||
            ptch.name() == downstreamName ||
            ptch.name() == upstreamName
        )
        {
            skipPatches.insert(patchI);
        }
    }

    //- detect average radius at axial position
    radiusAtPosition_.clear();
    forAll(axialCrossSection, tI)
    {
        const labelledTri& tri = axialCrossSection[tI];

        if( skipPatches.found(tri.region()) )
            continue;

        forAll(tri, pI)
        {
            const point& p = axialCrossSection.points()[tri[pI]];

            radiusAtPosition_[p.x()] = p.y();
        }
    }

    //- cleanup cross section
    deleteDemandDrivenData(axialCrossSectionPtr);
}

void profiledDieGeometryInterpolator::calculateProfilesAtRequiredPositions()
{
    //- const reference to radial cross-sections in the geometry
    const DynList<std::pair<scalar, std::shared_ptr<triSurf> > >& sections =
        dieProfiles_.crossSections();

    if( sections.size() == 0 )
    {
        //- this is not a profiled die
        return;
    }

    //- create cubic splines for the given profiles
    DynList<std::pair<scalar, std::shared_ptr<cubicBSpline> > > profileSplines;
    forAll(sections, i)
    {
        profileSplines.append
        (
            std::make_pair
            (
                sections[i].first,
                createCubicSpline(*sections[i].second)
            )
        );
    }

    //- scale cross section via linear interpolation all mesh positions
    profiles_.setSize(xToId_.size());
    label profileI(0);

    for(auto it=xToId_.begin();it!=xToId_.end();++it)
    {
        const scalar currX = it->first;

        if( currX <= sections[0].first )
        {
            //- find the scaling factor and scale the first profile
            //- start by calculating the ratio between the average radius
            //- at a given position
            const scalar rCurr = radiusAtPosition(currX);
            const scalar rFirst = radiusAtPosition(sections[0].first);

            const triSurf& s = *sections[0].second;

            std::shared_ptr<triSurf> scaledSurfPtr =
                std::make_shared<triSurf>(s);

            triSurfModifier sMod(*scaledSurfPtr);
            pointField& pts = sMod.pointsAccess();

            //- calculate scaling factor as the ratio of radius
            const scalar scalingFactor = (rCurr / (rFirst + VSMALL));

            forAll(pts, pI)
            {
                point& p = pts[pI];

                p.x() += currX - sections.lastElement().first;
                p.y() *= scalingFactor;
                p.z() *= scalingFactor;
            }

            //- insert data into the list of profiles
            it->second = profileI;
            profiles_[profileI] = std::make_tuple(currX, false, scaledSurfPtr);
            ++profileI;
        }
        else if( currX >= sections.lastElement().first )
        {
            //- find the scaling factor and scale the last profile

            //- start by calculating the ratio between the average radius
            //- at a given position
            const scalar rCurr = radiusAtPosition(currX);
            const scalar rLast = radiusAtPosition(sections.lastElement().first);

            const triSurf& s = *sections.lastElement().second;

            std::shared_ptr<triSurf> scaledSurfPtr =
                std::make_shared<triSurf>(s);

            triSurfModifier sMod(*scaledSurfPtr);
            pointField& pts = sMod.pointsAccess();

            const scalar scalingFactor = (rCurr / (rLast + VSMALL));

            forAll(pts, pI)
            {
                point& p = pts[pI];

                p.x() += currX - sections.lastElement().first;
                p.y() *= scalingFactor;
                p.z() *= scalingFactor;
            }

            //- insert data into the list of profiles
            it->second = profileI;
            profiles_[profileI] = std::make_tuple(currX, false, scaledSurfPtr);
            ++profileI;
        }
        else
        {
            forAll(profileSplines, sectI)
            {
                if
                (
                    currX >= profileSplines[sectI].first &&
                    currX <= profileSplines.fcElement(sectI).first
                )
                {
                    //- current x coordinate is within this interval
                    //- calculate the current cross-section as an interpolate
                    //- between the given profiles
                    std::shared_ptr<triSurf> profilePtr =
                        interpolateProfiles
                        (
                            currX,
                            *sections[sectI].second,
                            profileSplines[sectI],
                            profileSplines.fcElement(sectI)
                        );

                    //- insert data into the list of profiles
                    it->second = profileI;
                    profiles_[profileI] =
                        std::make_tuple(currX, false, profilePtr);
                    ++profileI;

                    break;
                }
            }
        }
    }

    # ifdef DEBUGDie
    Info << "Stopping " << endl;
    forAll(sections, i)
    {
        Info << "Profile at position " << sections[i].first << endl;
        sections[i].second->writeSurface
        (
            "givenProfile_"+std::to_string(i)+".stl"
        );
    }
    forAll(profiles_, i)
    {
        Info << "Writting profile " << i << " at position "
             << std::get<0>(profiles_[i]) << endl;

        std::get<2>(profiles_[i])->writeSurface
        (
            "intermediateProfile_"+std::to_string(i)+".stl"
        );
    }
    Info << "Written cross sections" << endl;
    # endif
}

void profiledDieGeometryInterpolator::initializeInterpolationMatrix
(
    List<List<List<point> > >& interpolationPoints,
    std::map<scalar, label>& indexOfAxialPosition,
    List<std::pair<scalar, scalar> >& minAndMaxRadiusAtAxialPosition
) const
{
    const pointFieldPMG& points = mesh_.points();
    const edgeLongList& edges = mesh_.addressingData().edges();

    //- calculate minimum deltaR and deltaTheta
    const scalar pi2 = 2.0 * M_PI;

    scalar deltaR(VGREAT), deltaTheta(pi2);
    scalar maxR(0.0), minR(VGREAT);

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        const point& ps = points[e.start()];
        const point& pe = points[e.end()];

        const scalar rs = sqrt(sqr(ps.y()) + sqr(ps.z()));
        const scalar re = sqrt(sqr(pe.y()) + sqr(pe.z()));

        maxR = max(maxR, rs);
        maxR = max(maxR, re);
        minR = min(minR, rs);
        minR = min(minR, re);
    }

    //- calculate the number of subdivisions in the circumferential direction
    label nTheta = sortedProfilePoints_.sizeOfRow(0);
    if( pointsInYSymmetry_.size() )
        nTheta *= 2;
    if( pointsInZSymmetry_.size() )
        nTheta *= 2;

    deltaTheta = pi2 / nTheta;

    //- calculate the number of subdivions in the radial direction
    if( dieProfiles_.dieDict().found("maxCellSize") )
        dieProfiles_.dieDict().readIfPresent("maxCellSize", deltaR);
    if( dieProfiles_.dieDict().found("contactCellSize") )
        dieProfiles_.dieDict().readIfPresent("contactCellSize", deltaR);

    const label nRadius = ceil((maxR - minR) / (deltaR + VSMALL));

    //- allocate space for points and other arrays
    indexOfAxialPosition.clear();
    minAndMaxRadiusAtAxialPosition.setSize(sortedProfilePoints_.size());
    interpolationPoints.setSize(sortedProfilePoints_.size());
    forAll(interpolationPoints, i)
    {
        interpolationPoints[i].setSize(nRadius+1);
        forAll(interpolationPoints[i], j)
            interpolationPoints[i][j].setSize(nTheta);
    }

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic)
    # endif
    forAll(interpolationPoints, i)
    {
        //- initialize point positions at the current section
        List<List<point> >& pointsAtSection = interpolationPoints[i];

        const point& p = points[sortedProfilePoints_(i, 0)];
        const point origin(p.x(), 0.0, 0.0);

        # ifdef USE_OMP
        # pragma omp critical(indexOfAxialPosition)
        # endif
        {
            indexOfAxialPosition[p.x()] = i;
        }

        const scalar rinner = mag(p - origin);
        const scalar rout = outerRadiusAtPosition(p.x());

        minAndMaxRadiusAtAxialPosition[i] = std::make_pair(rinner, rout);

        const scalar dR = (rout - rinner) / nRadius;

        const triSurf& surf = *std::get<2>(profiles_[i]);
        const std::shared_ptr<cubicBSpline> splinePtr = createCubicSpline(surf);

        //- initialize point positions
        forAll(pointsAtSection[0], k)
        {
            scalar alpha = k * deltaTheta;

            const scalar factor = min(max(alpha / pi2, 0.0), 1.0);

            const point newP = splinePtr->evaluate(factor);

            for(label j=0;j<nRadius;++j)
            {
                const scalar cr = rinner + j * dR;
                pointsAtSection[j][k] =
                    point(p.x(), cr * cos(alpha), cr * sin(alpha));
            }

            pointsAtSection[0][k] = newP;
            pointsAtSection[0][k].x() = p.x();

            const scalar cr = rinner + nRadius * dR;
            pointsAtSection[nRadius][k] =
                point(p.x(), cr * cos(alpha), cr * sin(alpha));
        }

        # ifdef DEBUGDie
        writeCrossSectionToVTK
        (
            pointsAtSection,
            "beforeOpt"+std::to_string(i)+".vtk"
        );
        # endif
    }
}

void profiledDieGeometryInterpolator::calculateInterpolationMatrixElliptic
(
    List<List<List<point> > >& sectionPoints,
    const List<std::pair<scalar, scalar> >& minAndMaxRadiusAtSection
) const
{
    Info << "Meshing radial cross-sections using the elliptic mesher" << endl;

    label numIterations(-1);
    if( dieProfiles_.dieDict().found("numSmoothingIterations") )
        dieProfiles_.dieDict().readIfPresent
        (
            "numSmoothingIterations",
            numIterations
        );

    scalar weight(1.2);
    if( dieProfiles_.dieDict().found("radialGrading") )
        dieProfiles_.dieDict().readIfPresent("radialGrading", weight);

    const scalar pi2 = 2.0 * M_PI;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic)
    # endif
    forAll(sectionPoints, sectionI)
    {
        List<List<point> >& pointsAtSection = sectionPoints[sectionI];

        const label nRadius = pointsAtSection.size() - 1;
        const label nTheta = pointsAtSection[0].size();

        const label nIterations = numIterations>=0?numIterations:(20*nRadius);

        const scalar rinner = minAndMaxRadiusAtSection[sectionI].first;
        const scalar rout = minAndMaxRadiusAtSection[sectionI].second;

        const scalar dR = (rout - rinner) / nRadius;
        const scalar dR2 = 2.0 * dR;
        const scalar dRSq = sqr(dR);

        const scalar deltaTheta = pi2 / nTheta;
        const scalar dth2 = 2.0 * deltaTheta;
        const scalar dthSq = sqr(deltaTheta);

        for(label iterI=0;iterI<nIterations;++iterI)
        {
            for(label j=1;j<nRadius;++j)
            {
                const scalar r = rinner + j * dR;

                const label jn = j+1;
                const label jp = j-1;

                for(label k=0;k<nTheta;++k)
                {
                    const point& pjk = pointsAtSection[j][k];

                    const scalar dydr =
                        (
                            pointsAtSection[jn][k].y() -
                            pointsAtSection[jp][k].y()
                        ) / dR2;
                    const scalar dzdr =
                        (
                            pointsAtSection[jn][k].z() -
                            pointsAtSection[jp][k].z()
                        ) / dR2;

                    const label kn = (k+1) % nTheta;
                    const label kp = (k-1+nTheta) % nTheta;

                    const scalar dydth =
                        (
                            pointsAtSection[j][kn].y() -
                            pointsAtSection[j][kp].y()
                        ) / dth2;
                    const scalar dzdth =
                        (
                            pointsAtSection[j][kn].z() -
                            pointsAtSection[j][kp].z()
                        ) / dth2;

                    const scalar a = sqr(r) * (sqr(dydr) + sqr(dydth));
                    const scalar b = -1.0 * sqr(r) * (dydr*dzdr + dydth*dzdth);
                    const scalar c = sqr(r) * (sqr(dzdr) + sqr(dzdth));

                    point newP =
                        (
                            (a * dthSq) *
                            (
                                pointsAtSection[jn][k] +
                                pointsAtSection[jp][k]
                            ) +
                            (c * dRSq) *
                            (
                                pointsAtSection[j][kn] +
                                pointsAtSection[j][kp]
                            ) +
                            (0.5 * b * dR * deltaTheta) *
                            (
                                pointsAtSection[jn][kn] -
                                pointsAtSection[jp][kn] -
                                pointsAtSection[jn][kp] +
                                pointsAtSection[jp][kp]
                            )
                        ) / (2.0 * (a*dthSq + c*dRSq) + VSMALL);

                    newP.x() = pjk.x();

                    //- update point coordinates
                    pointsAtSection[j][k] += 0.9 * (newP - pjk);
                }
            }

            //- enforce zero gradient in circumferential direction
            for(label k=0;k<nTheta;++k)
            {
                const point& pp = pointsAtSection[nRadius-1][k];

                const scalar rp = sqrt(sqr(pp.y()) + sqr(pp.z()));

                scalar theta = acos(min(max(pp.y() / rp, -1.0), 1.0));
                if( pp.z() < 0.0 )
                    theta = pi2 - theta;

                point newP
                (
                    pp.x(),
                    rout * cos(theta),
                    rout * sin(theta)
                );

                pointsAtSection[nRadius][k] = newP;
            }
        }

        # ifdef DEBUGDie
        writeCrossSectionToVTK
        (
            pointsAtSection,
            "afterOpt"+std::to_string(sectionI)+".vtk"
        );
        # endif
    }

    Info << "Finished meshing radial cross-sections" << endl;
}

void profiledDieGeometryInterpolator::calculateInterpolationMatrixSurfOpt
(
    List<List<List<point> > >& sectionPoints
) const
{
    Info << "Meshing radial cross-sections using the optimizer" << endl;

    label numIterations(-1);
    if( dieProfiles_.dieDict().found("numSmoothingIterations") )
        dieProfiles_.dieDict().readIfPresent
        (
            "numSmoothingIterations",
            numIterations
        );

    scalar weight(1.2);
    if( dieProfiles_.dieDict().found("radialGrading") )
        dieProfiles_.dieDict().readIfPresent("radialGrading", weight);

    const scalar pi2 = 2.0 * M_PI;

    //- create a matrix of triangles and their weights
    DynList<triFace> simplexTriangles(12);
    DynList<scalar> weights(12);
    simplexTriangles[0] = triFace(4, 5, 8);
    weights[0] = 1.0;
    simplexTriangles[1] = triFace(4, 8, 7);
    weights[1] = 1.0;
    simplexTriangles[2] = triFace(4, 5, 7);
    weights[2] = 1.0;
    simplexTriangles[3] = triFace(4, 7, 6);
    weights[3] = weight;
    simplexTriangles[4] = triFace(4, 6, 3);
    weights[4] = weight;
    simplexTriangles[5] = triFace(4, 7, 3);
    weights[5] = 1.0;
    simplexTriangles[6] = triFace(4, 3, 0);
    weights[6] = weight;
    simplexTriangles[7] = triFace(4, 0, 1);
    weights[7] = weight;
    simplexTriangles[8] = triFace(4, 3, 1);
    weights[8] = 1.0;
    simplexTriangles[9] = triFace(4, 1, 2);
    weights[9] = 1.0;
    simplexTriangles[10] = triFace(4, 2, 5);
    weights[10] = 1.0;
    simplexTriangles[11] = triFace(4, 1, 5);
    weights[11] = 1.0;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic)
    # endif
    forAll(sectionPoints, sectionI)
    {
        List<List<point> >& pointsAtSection = sectionPoints[sectionI];

        const label nRadius = pointsAtSection.size() - 1;
        const label nTheta = pointsAtSection[0].size();

        const label nIterations = numIterations>=0?numIterations:(20*nRadius);

        for(label iterI=0;iterI<nIterations;++iterI)
        {
            for(label j=1;j<nRadius;++j)
            {
                const label jn = j+1;
                const label jp = j-1;

                for(label k=0;k<nTheta;++k)
                {
                    const label kn = (k+1)%nTheta;
                    const label kp = (k-1+nTheta)%nTheta;

                    DynList<point> simplexPts(9);
                    simplexPts[0] = pointsAtSection[jp][kp];
                    simplexPts[1] = pointsAtSection[j][kp];
                    simplexPts[2] = pointsAtSection[jn][kp];
                    simplexPts[3] = pointsAtSection[jp][k];
                    simplexPts[4] = pointsAtSection[j][k];
                    simplexPts[5] = pointsAtSection[jn][k];
                    simplexPts[6] = pointsAtSection[jp][kn];
                    simplexPts[7] = pointsAtSection[j][kn];
                    simplexPts[8] = pointsAtSection[jn][kn];

                    forAll(simplexPts, i)
                    {
                        simplexPts[i].x() = simplexPts[i].y();
                        simplexPts[i].y() = simplexPts[i].z();
                        simplexPts[i].z() = 0.0;
                    }

                    surfaceOptimizerHeight hOpt
                    (
                        simplexPts,
                        simplexTriangles,
                        weights
                    );

                    const point newP = hOpt.optimizePoint(0.001);

                    pointsAtSection[j][k].y() = newP.x();
                    pointsAtSection[j][k].z() = newP.y();
                }
            }

            //- enforce zero gradient in circumferential direction
            for(label k=0;k<nTheta;++k)
            {
                const point& p = pointsAtSection[nRadius][k];
                const point& pp = pointsAtSection[nRadius-1][k];

                const scalar rp = sqrt(sqr(pp.y()) + sqr(pp.z()));

                scalar theta = acos(min(max(pp.y() / rp, -1.0), 1.0));
                if( pp.z() < 0.0 )
                    theta = pi2 - theta;

                point& currp = pointsAtSection[nRadius][k];
                const scalar currRadius = sqrt(sqr(currp.y()) + sqr(currp.z()));

                point newP
                (
                    p.x(),
                    currRadius * cos(theta),
                    currRadius * sin(theta)
                );

                pointsAtSection[nRadius][k] = newP;
            }
        }

        # ifdef DEBUGDie
        writeCrossSectionToVTK
        (
            pointsAtSection,
            "afterOpt"+std::to_string(sectionI)+".vtk"
        );
        # endif
    }

    Info << "Finished meshing radial cross-sections" << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

profiledDieGeometryInterpolator::profiledDieGeometryInterpolator
(
    const dieGeometryInfo& dieGeom,
    const rollingMillPatchNamesHandler& patchHandler,
    polyMeshGen& mesh
)
:
    dieProfiles_(dieGeom),
    patchHandler_(patchHandler),
    mesh_(mesh),
    sortedProfilePoints_(),
    pointsAtFixedRadius_(),
    pointsInYSymmetry_(),
    pointsInZSymmetry_(),
    radiusAtPosition_(),
    profiles_(),
    xToId_()
{
    detectProfilePointsInVolMesh();

    analyseAxialCrossSection();

    calculateProfilesAtRequiredPositions();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

profiledDieGeometryInterpolator::~profiledDieGeometryInterpolator()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void profiledDieGeometryInterpolator::axialPositions(scalarLongList& pos) const
{
    pos.clear();

    for(auto it=profiles_.begin();it!=profiles_.end();++it)
        pos.append(std::get<0>(*it));
}

const triSurf& profiledDieGeometryInterpolator::crossSection
(
    const scalar x
) const
{
    //- get an initial guess where is the nearest element
    std::map<scalar, label>::const_iterator lowIt =
        std::lower_bound
        (
            xToId_.begin(),
            xToId_.end(),
            x,
            [](const std::pair<scalar, label>& p1, const scalar& p2)
            {
                return p1.first < p2;
            }
        );

    //- check if the previous element is a better fit
    const label profileI = lowIt->second;
    const label prevProfileI = profileI>0?profileI-1:0;

    const scalar d = mag(x - std::get<0>(profiles_[profileI]));
    const scalar dPrev = mag(x - std::get<0>(profiles_[prevProfileI]));

    if( d < dPrev )
    {
        return *std::get<2>(profiles_[profileI]);
    }
    else
    {
        return *std::get<2>(profiles_[prevProfileI]);
    }

    FatalError << "Axial profile at position " << x << " does not exist"
        << abort(FatalError);

    return *std::get<2>(profiles_[0]);
}

void profiledDieGeometryInterpolator::calculateDisplacements
(
    std::map<label, vector>& pointDisplacements
) const
{
    pointDisplacements.clear();
    const pointFieldPMG& points = mesh_.points();

    const scalar pi2 = 2.0 * M_PI;

    forAll(sortedProfilePoints_, i)
    {
        const triSurf& surf = *std::get<2>(profiles_[i]);

        const std::shared_ptr<cubicBSpline> splinePtr = createCubicSpline(surf);

        forAllRow(sortedProfilePoints_, i, j)
        {
            const label pointI = sortedProfilePoints_(i, j);

            const point& p = points[pointI];
            const point origin(p.x(), 0.0, 0.0);

            const scalar r = mag(p - origin);

            scalar alpha = acos(min(max(p.y() / r, -1.0), 1.0));

            if( p.z() < 0.0 )
                alpha = pi2 - alpha;

            const scalar factor = min(max(alpha / pi2, 0.0), 1.0);

            const point newP = splinePtr->evaluate(factor);

            pointDisplacements[pointI] = newP - p;
            pointDisplacements[pointI].x() = 0.0;
        }
    }
}

static const vector n(-1.0, 0.0, 0.0);
std::pair<vector, scalar> calculatePointDisplacement
(
    const List<List<point> >& pts,
    const label j, const label k, // indices of the point being moved
    const label jOrigin, const label kOrigin, // centre of rotation
    const label oj, const label ok //
)
{
    const point& p = pts[j][k];
    vector v = pts[oj][ok] - pts[jOrigin][kOrigin];

//    Info << "j " << j << " k " << k << endl;
//    Info << "jOrigin " << jOrigin << " kOrigin " << kOrigin << endl;
//    Info << "oj " << oj << " ok " << ok << endl;

    scalar magv = mag(v);
    scalar dv = magv;
    if( jOrigin == j )
    {
//        Info << "1. pj " << oj << " pk " << k << endl;
        const point& pn = pts[oj][k];

        dv = mag(pn - pts[oj][ok]);
    }
    else if( kOrigin == k )
    {
//        Info << "2. pj " << j << " pk " << ok << endl;
        const point& pn = pts[j][ok];

        dv = mag(pn - pts[oj][ok]);
    }
    else
    {
        FatalError << "Invalid j and k provided" << abort(FatalError);
    }

    if( dv > VSMALL && magv > VSMALL )
    {
//        Info << "dv / magv " << (dv / magv) << endl;
//        Info << "1. v " << v << endl;
        v *= min((dv / magv), 1.0);
//        Info << "2. v " << v << endl;
    }

    const label djo = oj - j;
    label dko = ok - k;
    if( dko > 1 )
        dko = -1;
    if( dko < -1 )
        dko = 1;

    const label djOrigin = jOrigin - j;
    label dkOrigin = kOrigin - k;
    if( dkOrigin > 1 )
        dkOrigin = -1;
    if( dkOrigin < -1 )
        dkOrigin = 1;

    const scalar s = -1.0 * sign(djOrigin * dko - dkOrigin * djo);

    const vector vopt = s * (n ^ v);

    scalar w = 1.0;

//    if( (jOrigin == oj) && (oj < j) )
//    {
//        w = 1.4;
//        w = max(1.0, 6.0 - 1.0 * j);
//    }
//    else if( jOrigin < j || oj < j )
//    {
//        w = 1.1;
//    }
//    else if( jOrigin != oj && ((jOrigin > j) || (oj > j)) )
//    {
//        w = 1.1;
//    }

    return std::make_pair(w * ((vopt + pts[jOrigin][kOrigin]) - p), w);
}

vector calculatePointDisplacementElement
(
    const List<List<point> >& pts,
    const label j, const label k,
    const label oj, const label ok
)
{
    const point& p = pts[j][k];

    const point& pDiag = pts[oj][ok];

    const point& p1 = pts[j][ok];
    const point& p2 = pts[oj][k];

    const vector v1 = p1 - pDiag;
    const vector v2 = p2 - pDiag;

    //- find the point orthogonal to v1 and v2
    DynList<point> origins(3);
    origins[0] = pDiag;
    origins[1] = p1;
    origins[2] = p2;

    DynList<vector> normals(3);
    normals[0] = n;
    normals[1] = v1;
    normals[2] = v2;

    point pOrtho;
    const bool orthoExist = help::findMinimizerPoint(origins, normals, pOrtho);

    //- find the point of the parallelogram
    const point pPar = pDiag + v1 + v2;

    if( orthoExist /*&& magSqr(pOrtho - pDiag) < magSqr(pPar - pDiag)*/ )
    {
        return pOrtho - p;
    }

    return (pPar - p);
}

vector calculatePointDisplacement
(
    const List<List<point> >& pts,
    const label j, const label k
)
{
    const label jp = j - 1;
    const label jn = j + 1;
    const label kn = (k + 1) % pts[j].size();
    const label kp = (k - 1 + pts[j].size()) % pts[j].size();

//    vector dispVec = calculatePointDisplacementElement(pts, j, k, jp, kp);
//    dispVec += calculatePointDisplacementElement(pts, j, k, jn, kp);
//    dispVec += calculatePointDisplacementElement(pts, j, k, jn, kn);
//    dispVec += calculatePointDisplacementElement(pts, j, k, jp, kn);

//    return 0.25 * dispVec;

    //- displacements due to bottom vectors
    std::pair<vector, scalar> disp(vector::zero, 0.0), d;
    d = calculatePointDisplacement(pts, j, k,jp, k, jp, kn);
    disp.first += d.first;
    disp.second += d.second;

    d = calculatePointDisplacement(pts, j, k, jp, k, jp, kp);
    disp.first += d.first;
    disp.second += d.second;

    //- displacements due to the right side
    d = calculatePointDisplacement(pts, j, k, j, kn, jp, kn);
    disp.first += d.first;
    disp.second += d.second;

    d = calculatePointDisplacement(pts, j, k, j, kn, jn, kn);
    disp.first += d.first;
    disp.second += d.second;

    //- displacements due to the top side
    d = calculatePointDisplacement(pts, j, k, jn, k, jn, kn);
    disp.first += d.first;
    disp.second += d.second;

    d = calculatePointDisplacement(pts, j, k, jn, k, jn, kp);
    disp.first += d.first;
    disp.second += d.second;

    //- displacements due to the left side
    d = calculatePointDisplacement(pts, j, k, j, kp, jn, kp);
    disp.first += d.first;
    disp.second += d.second;

    d = calculatePointDisplacement(pts, j, k, j, kp, jp, kp);
    disp.first += d.first;
    disp.second += d.second;

    return (disp.first / (disp.second + VSMALL));
}

void profiledDieGeometryInterpolator::calculateDisplacementsAll
(
    vectorLongList& displacements
) const
{
    const pointFieldPMG& points = mesh_.points();

    displacements.setSize(points.size());
    displacements = vector::zero;

    //- initialize data
    List<List<List<point> > > interpolationPoints;
    std::map<scalar, label> indexOfAxialPosition;
    List<std::pair<scalar, scalar> > minAndMaxRadiusAtAxialPosition;

    initializeInterpolationMatrix
    (
        interpolationPoints,
        indexOfAxialPosition,
        minAndMaxRadiusAtAxialPosition
    );

    //- construct cross-sections
    word mesherType = "elliptic";
    dieProfiles_.dieDict().readIfPresent("mesherType", mesherType);

    if( mesherType == "elliptic" )
    {
        calculateInterpolationMatrixElliptic
        (
            interpolationPoints,
            minAndMaxRadiusAtAxialPosition
        );
    }
    else if( mesherType == "optimizer" )
    {
        calculateInterpolationMatrixSurfOpt(interpolationPoints);
    }
    else
    {
        FatalError << "Unknown meshType given " << mesherType
            << ". Available mesher types are elliptic and optimizer"
            << exit(FatalError);
    }

    //- calculate and apply displacements of points at the inner profile
    std::map<label, vector> bndDisplacements;
    calculateDisplacements(bndDisplacements);

    for(auto it=bndDisplacements.begin();it!=bndDisplacements.end();++it)
    {
        point disp = it->second;

        //- ensure consistency with symmetry planes
        if( pointsInYSymmetry_.found(it->first) )
            disp.y() = 0.0;
        if( pointsInZSymmetry_.found(it->first) )
            disp.z() = 0.0;

        displacements[it->first] = disp;
    }

    //- calculate displacement of other points
    const label nRadius = interpolationPoints[0].size() - 1;

    const scalar pi2 = 2.0 * M_PI;
    const label nTheta = interpolationPoints[0][0].size();
    const scalar deltaTheta = pi2 / nTheta;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(displacements, pointI)
    {
        //- check if the displacement is already set
        if( bndDisplacements.find(pointI) != bndDisplacements.end() )
            continue;

        //- calculate position in the polar coordinate system
        const point& p = points[pointI];

        const scalar r = sqrt(sqr(p.y()) + sqr(p.z()));

        scalar theta = acos(min(max(p.y() / r, -1.0), 1.0));
        if( p.z() < 0.0 )
            theta = pi2 - theta;

        //- calculate position in the interpolation space
        std::map<scalar, label>::const_iterator lowIt =
            std::lower_bound
            (
                indexOfAxialPosition.begin(),
                indexOfAxialPosition.end(),
                p.x(),
                [](const std::pair<scalar, label>& p1, const scalar& p2)
                {
                    return p1.first < p2;
                }
            );

        if( lowIt != indexOfAxialPosition.begin() )
            --lowIt;

        const label i = max(0, lowIt->second - 1);

        const std::pair<scalar, scalar>& minAndMaxR =
            minAndMaxRadiusAtAxialPosition[i];
        const scalar dr = (minAndMaxR.second - minAndMaxR.first) / nRadius;

        const scalar rpos = (r - minAndMaxR.first) / dr;
        const label j = max(min(floor(rpos), nRadius-1), 0);
        const label k = max(min(floor(theta / deltaTheta), nTheta-1), 0);

        const scalar lowx = lowIt->first;
        ++lowIt;
        scalar upperx = lowx;
        if( lowIt != indexOfAxialPosition.end() )
            upperx = lowIt->first;

        //- calculate parametric coordinates in the cylindrical space
        const scalar u = max(0.0, min(1.0, (p.x() - lowx) / (upperx - lowx)));

        const scalar lowr = minAndMaxR.first + j * dr;
        const scalar v = max(0.0, min(1.0, (r - lowr) / dr));

        const scalar lowtheta = k * deltaTheta;
        const scalar s = max(0.0, min(1.0, (theta - lowtheta) / deltaTheta));

        //- calculate new position of the point in the plane at min x
        const List<List<point> >& minxPts = interpolationPoints[i];
        const point newPXMin =
            help::interpolate
            (
                v, s,
                minxPts[j][k], minxPts[j+1][k],
                minxPts[j][(k+1)%nTheta], minxPts[j+1][(k+1)%nTheta]
            );

        //- calculate new position of the point in the plane at max x
        const List<List<point> >& maxxPts =
            interpolationPoints[min(i+1, interpolationPoints.size()-1)];
        const point newPXMax =
            help::interpolate
            (
                v, s,
                maxxPts[j][k], maxxPts[j+1][k],
                maxxPts[j][(k+1)%nTheta], maxxPts[j+1][(k+1)%nTheta]
            );

        //- calculate new position of point at the current x coordinate
        point newP = (1.0 - u) * newPXMin + u * newPXMax;

        //- check if the point is at the outer surface and shall keep its radius
        if( pointsAtFixedRadius_.found(pointI) )
        {
            vector dv = newP - point(newP.x(), 0.0, 0.0);

            const scalar rn = sqrt(sqr(newP.y()) + sqr(newP.z()));

            newP = point(newP.x(), 0.0, 0.0) + dv * (r / (rn + VSMALL));
        }

        //- calculate point displacement
        displacements[pointI] = newP - p;
        displacements[pointI].x() = 0.0;

        if( pointsInYSymmetry_.found(pointI) )
            displacements[pointI].y() = 0.0;
        if( pointsInZSymmetry_.found(pointI) )
            displacements[pointI].z() = 0.0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
