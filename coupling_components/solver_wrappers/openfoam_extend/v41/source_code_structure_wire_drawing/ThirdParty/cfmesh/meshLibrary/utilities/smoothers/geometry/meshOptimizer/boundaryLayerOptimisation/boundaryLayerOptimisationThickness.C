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
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "refLabelledScalar.H"
#include "polyMeshGenAddressing.H"

#include "labelledPoint.H"

//#define DEBUGLayer

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar boundaryLayerOptimisation::calculateMaxThickness
(
    const label heI,
    DynList<label>& neiHairs
) const
{
    neiHairs.clear();

    const labelLongList& faceOwner = meshSurfacePtr_->faceOwners();

    DynList<label> bndBaseFaces;
    findInfluencedHairs(heI, bndBaseFaces, neiHairs);

    scalar maxThickness(hairEdges_[heI].mag(mesh_.points()));

    //- find max allowed thickness over cells
    forAll(bndBaseFaces, i)
    {
        const label bfI = bndBaseFaces[i];
        const label baseFaceI = mesh_.nInternalFaces() + bfI;
        const label cOwn = faceOwner[bfI];

        //- check if there exist any self-intersections
        const scalar cThickness =
            calculateThicknessOverCell
            (
                heI,
                cOwn,
                baseFaceI
            );


        maxThickness = min(maxThickness, cThickness);
    }

    forAll(neiHairs, i)
    {
        const scalar thickness = calculateThickness(heI, neiHairs[i]);

        maxThickness = min(maxThickness, thickness);
    }

    return maxThickness;
}

void boundaryLayerOptimisation::findInfluencedHairs
(
    const label heI,
    DynList<label>& bndBaseFaces,
    DynList<label>& neiHairs
) const
{
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();

    const meshSurfaceEngine& mse = *meshSurfacePtr_;
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelLongList& bp = mse.bp();
    const VRWGraph& pFaces = mse.pointFaces();
    const labelLongList& faceOwner = mse.faceOwners();

    const label bpI = bp[hairEdges_[heI].start()];

    neiHairs.clear();
    bndBaseFaces.clear();
    forAllRow(pFaces, bpI, pfI)
    {
        const label bfI = pFaces(bpI, pfI);

        const face& bf = bFaces[bfI];

        //- boundary face must not contain the end points
        if( bf.which(hairEdges_[heI].end()) >= 0 )
            continue;

        //- cell must contain the end point
        const label baseFaceI = mesh_.nInternalFaces() + bfI;
        bool foundEnd(false);
        const cell& c = cells[faceOwner[bfI]];
        forAll(c, fI)
        {
            if( c[fI] == baseFaceI )
                continue;

            const face& f = faces[c[fI]];

            if( f.which(hairEdges_[heI].end()) >= 0 )
            {
                foundEnd = true;
                break;
            }
        }

        if( !foundEnd )
            continue;

        bndBaseFaces.append(bfI);
        DynList<edge> hairsAtFace;
        hairEdgesAtBndFace(faceOwner[bfI], baseFaceI, hairsAtFace);

        forAll(hairsAtFace, i)
        {
            const label bps = bp[hairsAtFace[i].start()];

            forAllRow(hairEdgesAtBndPoint_, bps, pfJ)
            {
                const label heJ = hairEdgesAtBndPoint_(bps, pfJ);

                if( hairEdges_[heJ] == hairsAtFace[i] )
                    neiHairs.appendIfNotIn(heJ);
            }
        }
    }
}

void boundaryLayerOptimisation::hairEdgesAtBndFace
(
    const label cellI,
    const label baseFaceI,
    DynList<edge>& hairEdges
) const
{
    const faceListPMG& faces = mesh_.faces();

    const cell& c = mesh_.cells()[cellI];

        //- check cell topology
    DynList<edge, 48> edges;
    DynList<DynList<label, 2>, 48> edgeFaces;
    DynList<DynList<label, 10>, 24> faceEdges;
    faceEdges.setSize(c.size());
    label baseFace(-1);
    forAll(c, fI)
    {
        if( c[fI] == baseFaceI )
        {
            baseFace = fI;
        }

        const face& f = faces[c[fI]];
        faceEdges[fI].setSize(f.size());

        forAll(f, eI)
        {
            const edge e = f.faceEdge(eI);

            label pos = edges.containsAtPosition(e);

            if( pos < 0 )
            {
                pos = edges.size();
                edges.append(e);
                edgeFaces.setSize(pos+1);
            }

            edgeFaces[pos].append(fI);
            faceEdges[fI][eI] = pos;
        }
    }

    const face& bf = faces[c[baseFace]];
    hairEdges.setSize(bf.size());

    forAll(bf, pI)
    {
        const label nextEdge = faceEdges[baseFace][pI];
        const label prevEdge = faceEdges[baseFace][bf.rcIndex(pI)];

        if( edgeFaces[nextEdge].size() != 2 || edgeFaces[prevEdge].size() != 2 )
            break;

        //- find the face attached to the edge after the current point
        label otherNextFace = edgeFaces[nextEdge][0];
        if( otherNextFace == baseFace )
            otherNextFace = edgeFaces[nextEdge][1];

        //- find the face attached to the edge before the current point
        label otherPrevFace = edgeFaces[prevEdge][0];
        if( otherPrevFace == baseFace )
            otherPrevFace = edgeFaces[prevEdge][1];

        label commonEdge;
        for(commonEdge=0;commonEdge<edges.size();++commonEdge)
            if(
                edgeFaces[commonEdge].contains(otherNextFace) &&
                edgeFaces[commonEdge].contains(otherPrevFace)
            )
                break;

        if( commonEdge == edges.size() )
            break;

        //- there exists a common edge which shall be used as a hair
        if( edges[commonEdge].start() == bf[pI] )
        {
            hairEdges[pI] = edges[commonEdge];
        }
        else
        {
            hairEdges[pI] = edges[commonEdge].reverseEdge();
        }
    }
}

scalar boundaryLayerOptimisation::calculateThickness
(
    const label heI,
    const label heJ
) const
{
    const pointFieldPMG& points = mesh_.points();

    //- references to hair edges
    const edge& he = hairEdges_[heI];
    const edge& nhe = hairEdges_[heJ];

    //- distance vector between the surface points of hair edges
    const point& sp = points[he[0]];
    const point& ep = points[nhe[0]];
    const vector dv = ep - sp;
    const scalar magDv = mag(dv);

    //- calculate layer thickness
    const scalar currThickness = he.mag(points);
    scalar retThickness = currThickness;

    const scalar currNeiThickness = nhe.mag(points);
    scalar suggestedNeiThickness = currNeiThickness;

    //- calculate layer height at the current point
    const point npAlpha = help::nearestPointOnTheEdge(sp, ep, points[he[1]]);
    const scalar currHeight = mag(npAlpha - points[he[1]]);
    scalar retHeight = currHeight;
    const scalar cosAlpha = sign((npAlpha - sp) & dv) * mag(npAlpha - sp);
    const scalar alpha =
        Foam::acos
        (
            Foam::max
            (
                -1.0,
                Foam::min(1.0, cosAlpha / (currThickness + VSMALL))
            )
        );

    //- calculate the height of the layer at the neighbour
    //- point
    const point npBeta = help::nearestPointOnTheEdge(ep, sp, points[nhe[1]]);
    const scalar currNeiHeight = mag(npBeta - points[nhe[1]]);
    scalar suggestedNeiHeight = currNeiHeight;
    const scalar cosBeta = sign((npBeta - ep) & -dv) * mag(npBeta - ep);
    const scalar beta =
        Foam::acos
        (
            Foam::max
            (
                -1.0,
                Foam::min(1.0, cosBeta / (currNeiThickness + VSMALL))
            )
        );

    //- check if the current thickness is Ok for the local curvature
    if( (alpha + beta) < M_PI )
    {
        const scalar gamma = M_PI - (alpha + beta);
        const scalar sinGamma = Foam::max(SMALL, Foam::sin(gamma));
        const scalar sinAlpha = Foam::max(SMALL, Foam::sin(alpha));
        const scalar sinBeta = Foam::max(SMALL, Foam::sin(beta));

        //- max allowed thickness and layer height due to curvature
        retThickness =
            Foam::min
            (
                retThickness,
                featureSizeFactor_ * magDv * sinBeta / sinGamma
            );
        retHeight *= (retThickness / (currThickness + VSMALL));

        //- max allowed neighbour hair thickness
        //- and layer height due to curvature
        suggestedNeiThickness =
            Foam::min
            (
                suggestedNeiThickness,
                featureSizeFactor_ * magDv * sinAlpha / sinGamma
            );
        suggestedNeiHeight *=
            (suggestedNeiThickness / (currNeiThickness + VSMALL));
    }

    //- check the height variation
    const scalar tanVal = (retHeight - suggestedNeiHeight) / (magDv + VSMALL);

    if( tanVal > relThicknessTol_ )
    {
        retHeight = suggestedNeiHeight + relThicknessTol_ * magDv;

        retThickness = (retHeight / currHeight) * currThickness;
    }

    return retThickness;
}

scalar boundaryLayerOptimisation::calculateThicknessOverCell
(
    const label heI,
    const label cellI,
    const label baseFaceI
) const
{
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    const cell& c = mesh_.cells()[cellI];

    const face& bf = faces[baseFaceI];

    const edge& he = hairEdges_[heI];

    const point& sp = points[he[0]];
    const point& ep = points[he[1]];

    scalar maxThickness = he.mag(points);

    //- the base face must not contain the hair edge
    //- this is the case at exitting layers
    forAll(bf, eI)
        if( bf.faceEdge(eI) == he )
            return maxThickness;

    forAll(c, fI)
    {
        if( c[fI] == baseFaceI )
            continue;

        const face& f = faces[c[fI]];

        if( help::shareAnEdge(bf, f) && (f.which(he.start()) == -1) )
        {
            point intersection;

            if( !help::lineFaceIntersection(sp, ep, f, points, intersection) )
                continue;

            const scalar maxDist = featureSizeFactor_ * mag(intersection - sp);

            maxThickness =
                Foam::min(maxThickness, maxDist);
        }
    }

    return maxThickness;
}

bool boundaryLayerOptimisation::propagateThicknessOfConstrainedEdges
(
    scalarLongList& hairLength,
    boolList& modifiedEdge,
    const direction edgeType
)
{
    bool changed(false);

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp single nowait
        # endif
        {
            for
            (
                std::map<label, DynList<label, 3> >::const_iterator it =
                constrainedHairs_.begin();
                it!=constrainedHairs_.end();
                ++it
            )
            {
                # ifdef USE_OMP
                # pragma omp task shared(modifiedEdge) firstprivate(it)
                # endif
                {
                    const DynList<label, 3>& nHairs = it->second;

                    //- check if any of the edges constrained with the current
                    //- one is modified
                    label minEdge = it->first;
                    bool anyModified = modifiedEdge[it->first];
                    direction gType = hairEdgeType_[it->first];
                    forAll(nHairs, i)
                    {
                        if( modifiedEdge[nHairs[i]] )
                        {
                            anyModified = true;
                            minEdge = min(minEdge, nHairs[i]);
                            gType |= hairEdgeType_[nHairs[i]];
                        }
                    }

                    //- do not process edges that have not been modified
                    if
                    (
                        anyModified && (minEdge == it->first) &&
                        (gType & edgeType)
                    )
                    {


                        scalar minDist = hairLength[it->first];
                        forAll(nHairs, i)
                            minDist = min(minDist, hairLength[nHairs[i]]);

                        if( minDist < hairLength[it->first] )
                        {
                            hairLength[it->first] = minDist;
                            thinnedHairEdge_[it->first] = true;
                            changed = true;
                        }

                        forAll(nHairs, i)
                        {
                            if( minDist < hairLength[nHairs[i]] )
                            {
                                changed = true;
                                thinnedHairEdge_[nHairs[i]] = true;
                                modifiedEdge[nHairs[i]] = true;
                            }

                            hairLength[nHairs[i]] = minDist;
                        }

                        //- mark all edges as the modified ones
                        const label epI = hairEdges_[it->first].end();
                        const DynList<label, 3>& hairsAtEnd =
                            hairEndPointsAtPoint_[epI];
                        forAll(hairsAtEnd, j)
                        {
                            const label hairEdgeJ = hairsAtEnd[j];

                            modifiedEdge[hairEdgeJ] = true;

                            forAll(constrainedHairs_[hairEdgeJ], k)
                            {
                                const label hairEdgeK =
                                    constrainedHairs_[hairEdgeJ][k];

                                modifiedEdge[hairEdgeK] = true;
                            }
                        }
                    }
                }
            }
        }
    }

    return changed;
}

bool boundaryLayerOptimisation::unifyEdgeLengthParallel
(
    scalarLongList& hairLength,
    boolList& modifiedEdge
)
{
    if( !Pstream::parRun() )
        return false;

    bool changed(false);

    //- collect data at inter-processor boundaries
    const meshSurfaceEngine& mse = meshSurface();
    const Map<label>& globalToLocalBnd =
        mse.globalToLocalBndPointAddressing();

    const polyMeshGenAddressing& addr = mesh_.addressingData();
    const labelLongList& globalPointLabel = addr.globalPointLabel();

    //- allocate space
    std::map<label, LongList<refLabelledScalar> > exchangeData;
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

            //- skip edges that are not at inter-processor boundaries
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
                    refLabelledScalar
                    (
                        it.key(),
                        labelledScalar(ge, hairLength[hairEdgeI])
                    )
                );
            }
        }
    }

    LongList<refLabelledScalar> receiveData;
    help::exchangeMap(exchangeData, receiveData);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100) \
    if( receiveData.size() > 1000 )
    # endif
    forAll(receiveData, i)
    {
        const refLabelledScalar& rls = receiveData[i];

        const label bpI = globalToLocalBnd[rls.objectLabel()];

        forAllRow(hairEdgesAtBndPoint_, bpI, j)
        {
            const label heJ = hairEdgesAtBndPoint_(bpI, j);
            const edge& he = hairEdges_[heJ];

            if( globalPointLabel[he.end()] == rls.lScalar().scalarLabel() )
            {
                modifiedEdge[heJ] = true;
                thinnedHairEdge_[heJ] = true;

                if( hairLength[heJ] > rls.lScalar().value() * (1.0+SMALL) )
                {
                    hairLength[heJ] = rls.lScalar().value();
                    changed = true;
                }
            }
        }
    }

    return changed;
}

void boundaryLayerOptimisation::optimiseThicknessVariation
(
    const direction edgeType
)
{
    const pointFieldPMG& points = mesh_.points();
    polyMeshGenModifier meshModifier(mesh_);

    meshSurface().pointFaces();

    vectorLongList hairDirections(hairEdges_.size());
    scalarLongList hairLength(hairEdges_.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges_, hairEdgeI)
    {
        vector n = hairEdges_[hairEdgeI].vec(points);

        hairLength[hairEdgeI] = (Foam::mag(n) + VSMALL);
        hairDirections[hairEdgeI] = n / hairLength[hairEdgeI];
    }

    //- reduce thickness of the layer
    //- such that the variation of layer thickness
    //- It is an iterative process where the layer is thinned in the regions
    //- where the tangent is greater than the tolerance value or the curvature
    //- permits thicker boundary layers.
    boolList activeHairEdge(hairEdges_.size(), true);
    boolList modifiedEdge(hairEdges_.size());

    bool changed;
    label nIter(0);
    do
    {
        boolList influencedEdges(hairEdges_.size());
        changed = false;

        //- check if the hair edge intersects some other face in the cells
        //- attached to the hair edge
        # ifdef USE_OMP
        # pragma omp parallel
        # endif
        {
            # ifdef USE_OMP
            # pragma omp for schedule(static, 1)
            # endif
            forAll(influencedEdges, hairEdgeI)
            {
                influencedEdges[hairEdgeI] = false;
                modifiedEdge[hairEdgeI] = false;
            }

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(hairEdges_, hairEdgeI)
            {
                if
                (
                    (hairEdgeType_[hairEdgeI] & edgeType) &&
                    activeHairEdge[hairEdgeI]
                )
                {
                    DynList<label> influencers;
                    const scalar maxThickness =
                        calculateMaxThickness(hairEdgeI, influencers);

                    if( hairLength[hairEdgeI] > (maxThickness * (1.0+SMALL)) )
                    {
                        //- make the hair edge shorter
                        hairLength[hairEdgeI] = maxThickness;
                        modifiedEdge[hairEdgeI] = true;
                        influencedEdges[hairEdgeI] = true;
                        changed = true;

                        thinnedHairEdge_[hairEdgeI] = true;

                        forAll(influencers, i)
                            influencedEdges[influencers[i]] = true;
                    }
                }
            }
        }

        if( Pstream::parRun() )
        {
            changed |= unifyEdgeLengthParallel(hairLength, modifiedEdge);
        }

        //- modify thickness of edges constrained with the modified ones
        changed |=
            propagateThicknessOfConstrainedEdges
            (
                hairLength,
                modifiedEdge,
                edgeType
            );

        if( Pstream::parRun() )
        {
            changed |= unifyEdgeLengthParallel(hairLength, modifiedEdge);
        }

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 100)
        # endif
        forAll(hairEdges_, hairEdgeI)
        {
            if( modifiedEdge[hairEdgeI] )
            {
                DynList<label> bndBaseFaces, influencers;
                findInfluencedHairs(hairEdgeI, bndBaseFaces, influencers);

                forAll(influencers, i)
                    influencedEdges[influencers[i]] = true;
            }
        }

        //- reduce the information over all processors
        reduce(changed, maxOp<bool>());

        if( !changed )
            break;

        //- move vertices to the new positions
        # ifdef USE_OMP
        # pragma omp parallel
        # endif
        {
            //- adjust hair connected to corner points
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(hairEdgesAtBndPoint_, bpI)
            {
                if( hairEdgesAtBndPoint_.sizeOfRow(bpI) == 3 )
                {
                    forAllRow(hairEdgesAtBndPoint_, bpI, i)
                    {
                        const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                        if( !modifiedEdge[hairEdgeI] )
                            continue;

                        const edge& he = hairEdges_[hairEdgeI];
                        const point& p = points[he.start()];
                        const vector& hv = hairDirections[hairEdgeI];
                        const scalar length = hairLength[hairEdgeI];

                        meshModifier.movePoint(he.end(), p + hv * length);
                    }
                }
            }

            //- adjust hairs connected to edge points
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(hairEdgesAtBndPoint_, bpI)
            {
                if( hairEdgesAtBndPoint_.sizeOfRow(bpI) == 2 )
                {
                    forAllRow(hairEdgesAtBndPoint_, bpI, i)
                    {
                        const label hairEdgeI = hairEdgesAtBndPoint_(bpI, i);

                        if( !modifiedEdge[hairEdgeI] )
                            continue;

                        const edge& he = hairEdges_[hairEdgeI];
                        const point& p = points[he.start()];
                        const vector& hv = hairDirections[hairEdgeI];
                        const scalar length = hairLength[hairEdgeI];

                        meshModifier.movePoint(he.end(), p + hv * length);
                    }
                }
            }

            //- adjust remaining hairs
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(hairEdgesAtBndPoint_, bpI)
            {
                if( hairEdgesAtBndPoint_.sizeOfRow(bpI) == 1 )
                {
                    const label hairEdgeI = hairEdgesAtBndPoint_(bpI, 0);

                    if( !modifiedEdge[hairEdgeI] )
                        continue;

                    const edge& he = hairEdges_[hairEdgeI];
                    const point& p = points[he.start()];
                    const vector& hv = hairDirections[hairEdgeI];
                    const scalar length = hairLength[hairEdgeI];

                    meshModifier.movePoint(he.end(), p + hv * length);
                }
            }
        }

        //- ensure the same position of points at all processors
        unifyCoordinatesParallel(modifiedEdge);

        //- mark edges which may be changed
        activeHairEdge.transfer(influencedEdges);
    } while( changed && (++nIter < 10000) );

    # ifdef DEBUGLayer
    if( Pstream::parRun() )
    {
        const meshSurfaceEngine& mse = meshSurface();
        const labelLongList& bp = mse.bp();
        const pointFieldPMG& points = mesh_.points();

        const labelLongList& globalPointLabel =
            mesh_.addressingData().globalPointLabel();
        const Map<label>& globalToLocal =
            mesh_.addressingData().globalToLocalPointAddressing();
        const VRWGraph& pAtProcs = mesh_.addressingData().pointAtProcs();
        const DynList<label>& neiProcs =
            mesh_.addressingData().pointNeiProcs();

        std::map<label, LongList<labelledPoint> > exchangeData;
        forAll(neiProcs, i)
            exchangeData[neiProcs[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label pointI = it();

            forAllRow(pAtProcs, pointI, i)
            {
                const label neiProc = pAtProcs(pointI, i);

                exchangeData[neiProc].append
                (
                    labelledPoint(it.key(), points[pointI])
                );
            }
        }

        LongList<labelledPoint> receiveData;
        help::exchangeMap(exchangeData, receiveData);

        forAll(receiveData, i)
        {
            const labelledPoint& lp = receiveData[i];

            const label pointI = globalToLocal[lp.pointLabel()];

            if( mag(lp.coordinates() - points[pointI]) > SMALL )
            {
                Pout << "Problematic point " << lp.pointLabel() << " at procs "
                     << pAtProcs[pointI]
                     << " coordinates " << points[pointI]
                     << " at other proc " << lp.coordinates()
                     << " boundary index " << bp[pointI] << endl;

                std::map<label, DynList<label, 3> >::const_iterator it =
                    hairEndPointsAtPoint_.find(pointI);

                if( it != hairEndPointsAtPoint_.end() )
                {
                    const DynList<label, 3>& eHairs = it->second;
                    Pout << "Problematic point " << lp.pointLabel()
                         << " edges at end point " << eHairs << endl;
                    forAll(eHairs, j)
                    {
                        const edge& he = hairEdges_[eHairs[j]];

                        Pout << "Edge at problematic point " << lp.pointLabel()
                             << " curr length " << he.mag(points)
                             << " wanted length " << hairLength[eHairs[j]]
                             << " start " << globalPointLabel[he.start()]
                             << " start pos " << points[he.start()]
                             << " end " << globalPointLabel[he.end()]
                             << " end pos " << points[he.end()] << endl;
                    }
                }
            }
        }
    }
    # endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
