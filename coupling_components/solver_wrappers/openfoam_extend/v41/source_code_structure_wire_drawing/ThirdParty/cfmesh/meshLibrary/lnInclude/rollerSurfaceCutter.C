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

#include "rollerSurfaceCutter.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "triSurfaceRemoveFacets.H"
#include "boundBox.H"
#include "helperFunctions.H"

//#define DEBUGRollerCutter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool rollerSurfaceCutter::classifyVertices()
{
    //- check the position of surface vertices according to the plane
    //- Find the vertices that shall be removed from the geometry
    const pointField& points = rollSurf_.points();

    boundBox bb(points);

    const scalar tol = SMALL * bb.mag();

    bool hasNegative(false), hasPositive(false);

    pointType_.setSize(points.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(pointType_, pI)
    {
        const point& p = points[pI];

        pointType_[pI] = NONE;

        const scalar d = ((p - plane_.refPoint()) & plane_.normal());

        if( d < tol )
        {
            if( d >= -tol )
            {
                hasPositive = true;

                //- point lies in the plane
                pointType_[pI] |= ATPLANE;
            }
            else
            {
                hasNegative = true;

                //- point is not used
                pointType_[pI] |= OUTSIDE;
            }
        }
        else
        {
            hasPositive = true;

            //- point is on the positive side of the plane
            pointType_[pI] |= INSIDE;
        }
    }

    if( !hasPositive )
    {
        FatalError << "The geometry does not contain any vertices on"
             << " the positive side of the plane " << plane_
             << ". All vertices will be removed" << exit(FatalError);
    }

    # ifdef DEBUGRollerCutter
    const label i0 = rollSurf_.addPointSubset("INSIDE");
    const label o0 = rollSurf_.addPointSubset("OUTSIDE");
    const label p0 = rollSurf_.addPointSubset("ATPLANE");

    forAll(pointType_, pI)
    {
        if( pointType_[pI] & OUTSIDE )
            rollSurf_.addPointToSubset(o0, pI);

        if( pointType_[pI] & INSIDE )
            rollSurf_.addPointToSubset(i0, pI);

        if( pointType_[pI] & ATPLANE )
            rollSurf_.addPointToSubset(p0, pI);
    }

    rollSurf_.writeSurface("surfaceMeshes/currRoll.fms");
    # endif

    return (hasPositive && hasNegative);
}

void rollerSurfaceCutter::newTrianglesFromCutEdges
(
    const std::map<label, label>& cutEdges
)
{
    if( cutEdges.size() == 0 || (cutEdges.size() % 2) )
        return;

    const edgeLongList& edges = rollSurf_.edges();
    const VRWGraph& edgeFacets = rollSurf_.edgeFacets();
    const VRWGraph& faceEdges = rollSurf_.facetEdges();

    for
    (
        std::map<label, label>::const_iterator it=cutEdges.begin();
        it!=cutEdges.end();
        ++it
    )
    {
        const label edgeI = it->first;

        //- find another triangle with the intersected edge
        if( edgeFacets.sizeOfRow(edgeI) != 1 )
        {
            FatalErrorIn
            (
                "void rollerSurfaceCutter::newTrianglesFromCutEdges"
                "(const std::map<label, label>&)"
            ) << "Invalid surface intersection" << abort(FatalError);
        }

        const label tI = edgeFacets(edgeI, 0);

        forAllRow(faceEdges, tI, teI)
        {
            const label eJ = faceEdges(tI, teI);

            forAllRow(edgeFacets, eJ, efJ)
            {
                const label tJ = edgeFacets(eJ, efJ);

                if( tJ == tI )
                    continue;

                forAllRow(faceEdges, tJ, teJ)
                {
                    const label oEdgeI = faceEdges(tJ, teJ);

                    //- filter out edges with smaller id to ensure
                    //- that triangles are modified only once
                    if( oEdgeI < edgeI )
                        continue;

                    if( cutEdges.find(oEdgeI) != cutEdges.end() )
                    {
                        //- create new triangles
                        const edge& e = edges[edgeI];
                        const edge& oe = edges[oEdgeI];

                        const label patchI = rollSurf_[tI].region();

                        label pe(-1);
                        if( !(pointType_[e.start()] & OUTSIDE) )
                        {
                            pe = e.start();
                        }
                        else if( !(pointType_[e.end()] & OUTSIDE) )
                        {
                            pe = e.end();
                        }
                        else
                        {
                            FatalErrorIn
                            (
                                "void rollerSurfaceCutter::"
                                "newTrianglesFromCutEdges"
                                "(const std::map<label, label>&)"
                            ) << "Cannot find surface vertex"
                              << abort(FatalError);
                        }

                        label poe(-1);
                        if( !(pointType_[oe.start()] & OUTSIDE) )
                        {
                            poe = oe.start();
                        }
                        else if( !(pointType_[oe.end()] & OUTSIDE) )
                        {
                            poe = oe.end();
                        }
                        else
                        {
                            FatalErrorIn
                            (
                                "void rollerSurfaceCutter::"
                                "newTrianglesFromCutEdges"
                                "(const std::map<label, label>&)"
                            ) << "Cannot find other vertex"
                              << abort(FatalError);
                        }

                        //- create new triangles
                        rollSurf_.appendTriangle
                        (
                            labelledTri(pe, it->second, poe, patchI)
                        );

                        std::map<label, label>::const_iterator oeIt =
                            cutEdges.find(oEdgeI);

                        if( oeIt == cutEdges.end() )
                        {
                            FatalErrorIn
                            (
                                "void rollerSurfaceCutter::"
                                "newTrianglesFromCutEdges"
                                "(const std::map<label, label>&)"
                            ) << "Other edge is not found"
                              << abort(FatalError);
                        }

                        rollSurf_.appendTriangle
                        (
                            labelledTri
                            (
                                it->second,
                                oeIt->second,
                                poe,
                                patchI
                            )
                        );
                    }
                }
            }
        }
    }

    rollSurf_.clearAddressing();
}

void rollerSurfaceCutter::sortPlanePoints(DynList<label>& planePoints) const
{
    const vector n = plane_.normal();

    const pointField& points = rollSurf_.points();

    //- calculate the centre of the patch
    vector c(vector::zero);
    forAll(planePoints, i)
        c += points[planePoints[i]];

    c /= planePoints.size();

    //- reference vector
    vector ref = points[planePoints[0]] - c;
    ref /= (mag(ref) + VSMALL);

    //- calculate angles accoding to the reference vector
    std::map<scalar, label> pointAngle;

    forAll(planePoints, i)
    {
        vector v = points[planePoints[i]] - c;
        v /= (mag(v) + VSMALL);

        scalar rad = Foam::acos(max(-1.0, min(1.0, (v & ref))));

        if( ((ref ^ v) & n) < 0.0 )
            rad = 2.0 * M_PI - rad;

        # ifdef DEBUGRollerCutter
        Info << "Plane point " << i << " with coordinates "
             << points[planePoints[i]] << " has angle " << rad << endl;
        # endif

        pointAngle[rad] = planePoints[i];
    }

    if( label(pointAngle.size()) != planePoints.size() )
    {
        # ifdef DEBUGRollerCutter
        Info << "Point angles " << label(pointAngle.size()) << endl;
        Info << "Plane points " << planePoints << endl;
        # endif

        Warning << "Collocated points found" << endl;
    }

    //- store data back in the correct order
    label counter(0);
    planePoints.setSize(pointAngle.size());
    for
    (
        std::map<scalar, label>::const_iterator it=pointAngle.begin();
        it!=pointAngle.end();
        ++it
    )
    {
        planePoints[counter++] = it->second;
    }
}

label rollerSurfaceCutter::patchForPlane() const
{
    const pointField& points = rollSurf_.points();

    //- the patch that shall become a symmetry plane must have a perfect
    //- alignment with plane normal and have the smallest distance
    //- from the plane
    const vector n = plane_.normal();

    std::set<label> filteredPatches;
    std::map<scalar, label> distanceToPatch;
    forAll(rollSurf_, triI)
    {
        vector tn = help::faceAreaVector(rollSurf_.points(), rollSurf_[triI]);
        tn /= (mag(tn) + VSMALL);

        if( mag(tn & n) > 0.999 )
        {
            const point fc = help::faceCentre(points, rollSurf_[triI]);
            const scalar d = (fc - plane_.refPoint()) & n;

            distanceToPatch[d] = rollSurf_[triI].region();
        }
        else
        {
            filteredPatches.insert(rollSurf_[triI].region());
        }
    }

    for
    (
        std::map<scalar, label>::const_iterator it=distanceToPatch.begin();
        it!=distanceToPatch.end();
        ++it
    )
    {
        if( filteredPatches.find(it->second) != filteredPatches.end() )
        {
            distanceToPatch.erase(it->first);

            it = distanceToPatch.begin();
        }
    }

    if( distanceToPatch.size() == 0 )
    {
        FatalErrorIn
        (
            "label rollerSurfaceCutter::patchForPlane"
            "(const DynList<label>&) const"
        ) << "No patch candidates for the symmetry plane" << exit(FatalError);
    }

    return distanceToPatch.begin()->second;
}

void rollerSurfaceCutter::createIntersections()
{
    DynList<label> planePoints;

    const pointField& points = rollSurf_.points();
    const edgeLongList& edges = rollSurf_.edges();

    std::map<label, label> newPointAtEdge;
    DynList<point> newPoints;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        const label es = e.start();
        const label ee = e.end();

        const point& sp = points[es];
        const point& ep = points[ee];

        //- do not cut edges spanning in the x direction
        if( mag(sp.x() - ep.x()) > SMALL )
            continue;

        if
        (
            ((pointType_[es] & OUTSIDE) && (pointType_[ee] & INSIDE)) ||
            ((pointType_[es] & INSIDE) && (pointType_[ee] & OUTSIDE))
        )
        {
            //- calculate intersection with the plane
            point intersection = vector::zero;
            const bool hasIntersection =
                help::planeIntersectsEdge(sp, ep, plane_, intersection);

            if( !hasIntersection )
            {
                Info << "es  is outside "
                     << bool(pointType_[es] & OUTSIDE) << endl;
                Info << " ee is outside "
                     << bool(pointType_[ee] & OUTSIDE) << endl;
                Info << "v " << (ep - sp) << " n " << plane_.normal() << endl;
                FatalErrorIn
                (
                    "void rollerSurfaceCutter::createIntersections()"
                ) << "Intersection must exist " << exit(FatalError);
            }

            # ifdef USE_OMP
            # pragma omp critical(appendVertex)
            # endif
            {
                const label npI = rollSurf_.nPoints()+newPoints.size();
                newPointAtEdge[edgeI] = npI;
                planePoints.appendIfNotIn(npI);

                newPoints.append(intersection);
            }
        }

        if( pointType_[es] & ATPLANE )
        {
            # ifdef USE_OMP
            # pragma omp critical(appendVertex)
            # endif
            planePoints.appendIfNotIn(es);
        }

        if( pointType_[ee] & ATPLANE )
        {
            # ifdef USE_OMP
            # pragma omp critical(appendVertex)
            # endif
            planePoints.appendIfNotIn(ee);
        }
    }

    //- create new triangles from cut edges
    newTrianglesFromCutEdges(newPointAtEdge);

    //- construct surface modifier
    triSurfModifier sMod(rollSurf_);
    pointField& pts = sMod.pointsAccess();
    geometricSurfacePatchList& patches = sMod.patchesAccess();

    //- append new vertices to the surface mesh
    const label nOrigPoints = pts.size();
    pts.setSize(nOrigPoints+newPoints.size());

    forAll(newPoints, i)
        pts[nOrigPoints+i] = newPoints[i];

    //- create triangles at intersected edges
    if( planePoints.size() >= 3 )
    {
        //- find the patch than shall become symmetry plane
        label patchI(-1);
        if( patchName_.empty() )
        {
            patchI = patchForPlane();
        }
        else
        {
            forAll(patches, i)
            {
                if( patches[i].name() == patchName_ )
                {
                    patchI = i;
                    break;
                }
            }

            if( patchI < 0 )
            {
                patchI = patches.size();
                patches.setSize(patchI+1);

                patches[patchI].name() = patchName_;
            }
        }

        patches[patchI].geometricType() = "symmetryPlane";

        //- sort vertices into cyclic order
        sortPlanePoints(planePoints);

        //- create and append new triangles
        for(label i=1;i<planePoints.size()-1;++i)
        {
            rollSurf_.appendTriangle
            (
                labelledTri
                (
                    planePoints[0], planePoints[i], planePoints[i+1], patchI
                )
            );
        }
    }
    else
    {
        FatalErrorIn
        (
            "void rollerCurfaceCutter::createIntersections()"
        ) << "Invalid number of cuts by a plane. "
          << "Found " << planePoints.size() << " and it shall be at least 4"
          << abort(FatalError);
    }

    # ifdef DEBUGRollerCutter
    rollSurf_.writeSurface("surfaceMeshes/currRollCut.fms");
    # endif
}

void rollerSurfaceCutter::removeOutsideTriangles()
{
    const VRWGraph& pTriangles = rollSurf_.pointFacets();
    const label removeId = rollSurf_.addFacetSubset("_remove_");

    forAll(pointType_, pI)
    {
        if( pointType_[pI] & OUTSIDE )
        {
            forAllRow(pTriangles, pI, ptI)
            {
                const label tI = pTriangles(pI, ptI);

                rollSurf_.addFacetToSubset(removeId, tI);
            }
        }
    }

    //- remove the triangles outside of the surface
    triSurfaceRemoveFacets remover(rollSurf_);
    remover.selectFacetsInSubset("_remove_");
    remover.removeFacets();

    rollSurf_.removeFacetSubset(removeId);

    # ifdef DEBUGRollerCutter
    rollSurf_.writeSurface("surfaceMeshes/currRollCleaned.fms");
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollerSurfaceCutter::rollerSurfaceCutter
(
    triSurf& rollSurf,
    const plane& cutPlane,
    const word patchName
)
:
    rollSurf_(rollSurf),
    plane_(cutPlane),
    patchName_(patchName),
    pointType_()
{
    //- classify vertices according to their position
    if( classifyVertices() )
    {
        //- create intersection points and triangles
        createIntersections();

        //- remove triangles and vertices that are not needed any more
        removeOutsideTriangles();
    }
}

rollerSurfaceCutter::~rollerSurfaceCutter()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
