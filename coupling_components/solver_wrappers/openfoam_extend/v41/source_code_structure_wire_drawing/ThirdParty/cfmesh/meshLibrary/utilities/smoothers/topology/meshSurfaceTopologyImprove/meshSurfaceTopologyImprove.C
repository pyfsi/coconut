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

#include "meshSurfaceTopologyImprove.H"
#include "demandDrivenData.H"
#include "helperFunctions.H"
#include "meshSurfaceEngine.H"
#include "meshOptimizer.H"
#include "meshOctree.H"
#include "meshSurfacePartitioner.H"
#include "meshSurfaceMapper.H"
#include "meshSurfaceOptimizer.H"

# ifdef USE_OMP
#include <omp.h>
# endif

#define DEBUGTopoImprove

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

void meshSurfaceTopologyImprove::createMeshSurface() const
{
    meshSurfacePtr_ = new meshSurfaceEngine(mesh_);
}

void meshSurfaceTopologyImprove::clearMeshSurface()
{
    deleteDemandDrivenData(meshSurfacePtr_);
}

void meshSurfaceTopologyImprove::optimizeCells()
{
    //- untangle the mesh first, do not allow movement of surface points
    meshOptimizer mOpt(mesh_);
    mOpt.optimizeMeshNearBoundaries(2, 2);
    mOpt.untangleMeshFV(2, 50, 0);
}

label meshSurfaceTopologyImprove::detectProblematicFaces()
{
    const pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = meshSurface();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelLongList& fOwner = mse.faceOwners();
    const labelLongList& facePatches = mse.boundaryFacePatches();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    const vectorLongList& cellCentres = mesh_.addressingData().cellCentres();

    const Map<label>* otherFacePatchPtr = NULL;
    if( Pstream::parRun() )
        otherFacePatchPtr = &mse.otherEdgeFacePatch();

    problemAtFace_.setSize(bFaces.size());

    label nProblematic(0);

    # ifdef USE_OMP
    # pragma omp parallel reduction(+:nProblematic)
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(problemAtFace_, bfI)
            problemAtFace_[bfI] = FACEOK;

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 50)
        # endif
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            DynList<label> neiPatches(bf.size());
            DynList<label> otherPatches;

            forAllRow(faceEdges, bfI, feI)
            {
                const label beI = faceEdges(bfI, feI);

                if( edgeFaces.sizeOfRow(beI) == 2 )
                {
                    label otherFace = edgeFaces(beI, 0);
                    if( otherFace == bfI )
                        otherFace = edgeFaces(beI, 1);

                    neiPatches[feI] = facePatches[otherFace];
                    otherPatches.appendIfNotIn(facePatches[otherFace]);
                }
                else
                {
                    neiPatches[feI] =
                        otherFacePatchPtr->operator[](beI);
                    otherPatches.append(neiPatches[feI]);
                }
            }

            if
            (
                otherPatches.size() == 1 &&
                otherPatches[0] == facePatches[bfI]
            )
                continue;

            //- check the face twist
            const point fCentre = help::faceCentre(points, bf);

            forAll(neiPatches, eI)
            {
                if( neiPatches[eI] != facePatches[bfI] )
                {
                    const edge e = bf.faceEdge(eI);

                    vector ev = e.vec(points);
                    ev /= (mag(ev) + VSMALL);

                    const point np =
                        help::nearestPointOnTheEdgeExact
                        (
                            points[e.start()],
                            points[e.end()],
                            fCentre
                        );

                    vector dv = np - fCentre;
                    dv /= (mag(dv) + VSMALL);

                    if( mag(dv & ev) > VSMALL )
                        problemAtFace_[bfI] |= TWISTEDATEDGE;
                }
            }

            //- check if the face is visible from the cell centre
            const point& c = cellCentres[fOwner[bfI]];
            const point np = help::nearestPointOnFace(bf, points, c);

            vector dv = (np - c);
            dv /= (mag(dv) + VSMALL);

            vector n = help::faceAreaVector(points, bf);
            n /= (mag(n) + VSMALL);

            if( (dv & n) < VSMALL )
            {
                //- face is not visible from the cell centre
                //- this indicates that the cell is very twisted
                problemAtFace_[bfI] = BADVISIBILITY;
                ++nProblematic;
            }

            //- check if there exist multiple faces belonging to the same cell
            //- that are in different patches
            forAllRow(faceEdges, bfI, feI)
            {
                const label beI = faceEdges(bfI, feI);

                forAllRow(edgeFaces, beI, efI)
                {
                    const label bfJ = edgeFaces(beI, efI);

                    if( bfJ == bfI )
                        continue;

                    if( fOwner[bfI] == fOwner[bfJ] )
                    {
                        //- faces are in the same cell
                        if( facePatches[bfI] != facePatches[bfJ] )
                        {
                            const face& bfn = bFaces[bfJ];

                            if( !help::isSharedEdgeConvex(points, bf, bfn) )
                            {
                                problemAtFace_[bfI] |= CONCAVEEDGE;
                                problemAtFace_[bfJ] |= CONCAVEEDGE;
                            }
                        }
                    }
                }
            }
        }
    }

    # ifdef DEBUGTopoImprove
    const label concaveCellId = mesh_.addCellSubset("concaveCells");
    const label concaveId = mesh_.addFaceSubset("concaveEdgeFaces");
    const label badVisibilityId = mesh_.addFaceSubset("badVisibilityFaces");
    const label badVisibilityCellsId =
        mesh_.addCellSubset("badVisibilityCells");
    const label twistedEdgeId = mesh_.addFaceSubset("twistedEdgeFaces");
    const label twistedEdgeCellsId = mesh_.addCellSubset("twistedEdgeCells");

    forAll(problemAtFace_, bfI)
    {
        if( problemAtFace_[bfI] & BADVISIBILITY )
        {
            mesh_.addFaceToSubset(badVisibilityId, bfI+mesh_.nInternalFaces());
            mesh_.addCellToSubset(badVisibilityCellsId, fOwner[bfI]);
        }
        if( problemAtFace_[bfI] & TWISTEDATEDGE )
        {
            mesh_.addFaceToSubset(twistedEdgeId, bfI+mesh_.nInternalFaces());
            mesh_.addCellToSubset(twistedEdgeCellsId, fOwner[bfI]);
        }
        if( problemAtFace_[bfI] & CONCAVEEDGE )
        {
            mesh_.addFaceToSubset(concaveId, bfI+mesh_.nInternalFaces());
            mesh_.addCellToSubset(concaveCellId, fOwner[bfI]);
        }
    }
    # endif

    reduce(nProblematic, sumOp<label>());

    return nProblematic;
}

bool meshSurfaceTopologyImprove::modifyDeformedFaces()
{
    bool changed(false);

    const meshSurfaceEngine& mse = meshSurface();

    const pointFieldPMG& points = mse.points();

    const labelLongList& facePatches = mse.boundaryFacePatches();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    const Map<label>* otherFacePatchPtr = NULL;
    if( Pstream::parRun() )
        otherFacePatchPtr = &mse.otherEdgeFacePatch();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        //- find patches over edges
        DynList<label> otherPatch(bf.size());
        DynList<scalar> origAngles(bf.size()), newAngles(bf.size());

        forAllRow(faceEdges, bfI, feI)
        {
            //- calculate angle
            vector vn = points[bf.nextLabel(feI)] - points[bf[feI]];
            vn /= (mag(vn) + VSMALL);
            vector vp = points[bf.prevLabel(feI)] - points[bf[feI]];
            vp /= (mag(vp) + VSMALL);

            newAngles[feI] = (vn & vp);

            if
            (
                mesh_.hasBackup(bf[feI]) &&
                mesh_.hasBackup(bf.prevLabel(feI)) &&
                mesh_.hasBackup(bf.nextLabel(feI))
            )
            {
                point p, pp, np;
                mesh_.getOrigPoint(bf[feI], p);
                mesh_.getOrigPoint(bf.prevLabel(feI), pp);
                mesh_.getOrigPoint(bf.nextLabel(feI), np);

                vn = np - p;
                vn /= (mag(vn) + VSMALL);
                vp = pp - p;
                vp /= (mag(vp) + VSMALL);

                origAngles[feI] = (vn & vp);
            }

            const label beI = faceEdges(bfI, feI);

            if( edgeFaces.sizeOfRow(beI) == 2 )
            {
                label otherFace = edgeFaces(beI, 0);
                if( otherFace == bfI )
                    otherFace = edgeFaces(beI, 1);

                otherPatch[feI] = facePatches[otherFace];
            }
            else
            {
                otherPatch[feI] = otherFacePatchPtr->operator[](beI);
            }
        }


    }

    reduce(changed, maxOp<bool>());

    return changed;
}

bool meshSurfaceTopologyImprove::updateFacePatches()
{
    bool modified;

    const pointFieldPMG& points = mesh_.points();
    const cellListPMG& cells = mesh_.cells();

    modified = false;

    const meshSurfaceEngine& mse = meshSurface();
    const labelLongList& facePatches = mse.boundaryFacePatches();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();
    const labelLongList& fOwner = mse.faceOwners();
    const labelLongList& bPoints = mse.boundaryPoints();

    pointField bndPointsBefore(bPoints.size());
    forAll(bPoints, bpI)
        bndPointsBefore[bpI] = points[bPoints[bpI]];

    const vectorLongList& cellCentres = mesh_.addressingData().cellCentres();

    const Map<label>* otherFacePatchPtr = NULL;
    if( Pstream::parRun() )
        otherFacePatchPtr = &mse.otherEdgeFacePatch();

    //- try combinations of face patches at problematic regions
    newFacePatch_.setSize(bFaces.size());
    forAll(facePatches, bfI)
        newFacePatch_[bfI] = facePatches[bfI];

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(bFaces, bfI)
    {
        if( problemAtFace_[bfI] == FACEOK )
            continue;

        const face& bf = bFaces[bfI];

        if( problemAtFace_[bfI] & CONCAVEEDGE )
        {
            //- find boundary faces in the cell
            const cell& c = cells[fOwner[bfI]];

            const point& cent = cellCentres[fOwner[bfI]];

            DynList<label> bndFaces;

            forAll(c, fI)
            {
                const label bfJ = c[fI] - mesh_.nInternalFaces();

                if( bfJ < 0 || bfJ >= bFaces.size() )
                    continue;

                bndFaces.append(bfJ);
            }

            //- find a face with a concave angle that shares an edge
            //- with the current bnd face and is not in the same patch
            const point np = help::nearestPointOnFace(bf, points, cent);
            vector dv = np - cent;
            dv /= (mag(dv) + VSMALL);

            vector fn = help::faceAreaVector(points, bf);
            fn /= (mag(fn) + VSMALL);

            const scalar fAlign = (fn & dv);

            forAll(bndFaces, j)
            {
                const label bfJ = bndFaces[j];

                if( facePatches[bfJ] == facePatches[bfI] )
                    continue;

                const face& obf = bFaces[bfJ];

                if
                (
                    help::shareAnEdge(bf, obf) &&
                    !help::isSharedEdgeConvex(points, bf, obf)
                )
                {

                    const point onp =
                        help::nearestPointOnFace(obf, points, cent);

                    vector odv = onp - cent;
                    odv /= (mag(odv) + VSMALL);

                    vector ofn = help::faceAreaVector(points, obf);
                    ofn /= (mag(ofn) + VSMALL);

                    const scalar ofAlign = (odv & ofn);

                    if( ofAlign > (fAlign+SMALL) )
                    {
                        //- change the patch
                        newFacePatch_[bfI] = facePatches[bfJ];

                        if( newFacePatch_[bfI] != facePatches[bfI] )
                            modified = true;
                    }
                }
            }
        }
        else if( problemAtFace_[bfI] & BADVISIBILITY )
        {
            const point& cent = cellCentres[fOwner[bfI]];

            //- find the nearest edge
            DynList<label> neiPatch(bf.size());
            label nearestEdge(-1);
            scalar distSq(VGREAT);

            forAllRow(faceEdges, bfI, feI)
            {
                const edge e = bf.faceEdge(feI);
                const label beI = faceEdges(bfI, feI);

                if( edgeFaces.sizeOfRow(beI) == 2 )
                {
                    if( edgeFaces(beI, 0) == bfI )
                    {
                        neiPatch[feI] = facePatches[edgeFaces(beI, 1)];
                    }
                    else
                    {
                        neiPatch[feI] = facePatches[edgeFaces(beI, 0)];
                    }
                }
                else if( edgeFaces.sizeOfRow(beI) == 1 )
                {
                    neiPatch[feI] = otherFacePatchPtr->operator[](beI);
                }

                const point np =
                    help::nearestPointOnTheEdgeExact
                    (
                        points[e.start()],
                        points[e.end()],
                        cent
                    );
                const scalar dSq = magSqr(np - cent);

                if( dSq < distSq )
                {
                    nearestEdge = feI;
                    distSq = dSq;
                }
            }

            //- set the patch of the
            newFacePatch_[bfI] = neiPatch[nearestEdge];
            Info << "Orig patch " << facePatches[bfI] << " new patch "
                 << newFacePatch_[bfI] << endl;
            if( newFacePatch_[bfI] != facePatches[bfI] )
                modified = true;
        }
    }

    reduce(modified, maxOp<label>());

    return modified;
}

bool meshSurfaceTopologyImprove::analyseUpdatedPatches()
{
    const meshSurfaceEngine& mse = meshSurface();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& faceEdges = mse.faceEdges();
    const VRWGraph& edgeFaces = mse.edgeFaces();

    const pointFieldPMG& points = mesh_.points();

    Info << "Analysing updated patches" << endl;

    Map<label> newOtherPatch;
    if( Pstream::parRun() )
    {
        //- find newly-assigned patch across the inter-processor boundaries
        const Map<label>& otherFaceProc = mse.otherEdgeFaceAtProc();
        const Map<label>& globalToLocal = mse.globalToLocalBndEdgeAddressing();

        //- create a map containing messages
        std::map<label, labelLongList> exchangeData;
        forAll(mse.beNeiProcs(), i)
            exchangeData[mse.beNeiProcs()[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label beI = it();

            if( edgeFaces.sizeOfRow(beI) == 1 )
            {
                const label fPatch = newFacePatch_[edgeFaces(beI, 0)];

                //- send:
                //- 1. global edge label
                //- 2. face patch
                labelLongList& dts = exchangeData[otherFaceProc[beI]];
                dts.append(it.key());
                dts.append(fPatch);
            }
        }

        //- exchange data among processors
        labelLongList receiveData;
        help::exchangeMap(exchangeData, receiveData);

        for(label i=0;i<receiveData.size();)
        {
            const label beI = globalToLocal[receiveData[i++]];
            const label fPatch = receiveData[i++];

            newOtherPatch.insert(beI, fPatch);
        }
    }

    //- analyse changes in the assigned patches over the local neighbourhood
    //- pay attention to the requested movement of points towards the
    //- feature edges
    bool changed(false);

    //- correct island faces
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(newFacePatch_, bfI)
    {
        const face& bf = bFaces[bfI];

        //- check if some faces need to be changed to improve
        DynList<label> neiPatches(bf.size());
        DynList<label> newNeiPatches, nNeiInPatch;

        forAllRow(faceEdges, bfI, feI)
        {
            const label beI = faceEdges(bfI, feI);

            if( edgeFaces.sizeOfRow(beI) == 2 )
            {
                label neiFace = edgeFaces(beI, 0);
                if( neiFace == bfI )
                    neiFace = edgeFaces(beI, 1);

                neiPatches[feI] = newFacePatch_[neiFace];
                label pos = newNeiPatches.contains(newFacePatch_[neiFace]);
                if( pos < 0 )
                {
                    newNeiPatches.append(newFacePatch_[neiFace]);
                    nNeiInPatch.append(1);
                }
                else
                {
                    ++nNeiInPatch[pos];
                }
            }
            else if( edgeFaces.sizeOfRow(beI) == 1 )
            {
                neiPatches[feI] = newOtherPatch[beI];
                label pos = newNeiPatches.contains(newOtherPatch[beI]);

                if( pos < 0 )
                {
                    newNeiPatches.append(newOtherPatch[beI]);
                    nNeiInPatch.append(1);
                }
                else
                {
                    ++nNeiInPatch[pos];
                }
            }
        }

        if
        (
            (nNeiInPatch.size() == 1) &&
            (nNeiInPatch[0] == bf.size()) &&
            (newNeiPatches[0] != newFacePatch_[bfI])
        )
        {
            newFacePatch_[bfI] = newNeiPatches[0];
            changed = true;

            continue;
        }

        //- find new point coordinates with respect
        //- to the newly selected patches
        DynList<point> newFacePoints(bf.size());
        DynList<scalar> distSq(bf.size());

        forAll(bf, pI)
        {
            DynList<label> patches;
            patches.append(newFacePatch_[bfI]);
            patches.appendIfNotIn(neiPatches[pI]);
            patches.appendIfNotIn(neiPatches[neiPatches.rcIndex(pI)]);

            octree_.findNearestPointToPatches
            (
                newFacePoints[pI],
                distSq[pI],
                points[bf[pI]],
                patches
            );
        }
    }

    reduce(changed, maxOp<bool>());

    if( changed )
    {
        labelLongList facePatchCopy(newFacePatch_.size());
        forAll(newFacePatch_, bfI)
            facePatchCopy[bfI] = newFacePatch_[bfI];

        Info << "facePatchCopy " << facePatchCopy << endl;

        meshSurfacePartitioner mPart(mse, facePatchCopy);
        meshSurfaceMapper(mPart, octree_).mapVerticesOntoSurfacePatches();
        meshSurfaceOptimizer(mPart, octree_).optimizeSurface();

        wordList patchNames(mesh_.boundaries().size());
        forAll(patchNames, patchI)
            patchNames[patchI] = mesh_.boundaries()[patchI].patchName();

        //- replace boundary
        const labelLongList& fOwner = mse.faceOwners();

        VRWGraph newBndFaces;
        labelLongList newBndOwner;
        forAll(bFaces, bfI)
        {
            newBndFaces.appendList(bFaces[bfI]);
            newBndOwner.append(fOwner[bfI]);
        }

        Info << "replacing boundary" << endl;
        polyMeshGenModifier(mesh_).replaceBoundary
        (
            patchNames,
            newBndFaces,
            newBndOwner,
            newFacePatch_
        );

        Info << "Finished replacing boundary" << endl;

        clearMeshSurface();

        meshOptimizer(mesh_).untangleMeshFV();
        Info << "Writing mesh and stopping" << endl;
        mesh_.write();
        ::exit(0);
    }

    return changed;
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //
// Constructors

meshSurfaceTopologyImprove::meshSurfaceTopologyImprove
(
    polyMeshGen& mesh,
    const meshOctree& octree
)
:
    mesh_(mesh),
    octree_(octree),
    meshSurfacePtr_(NULL),
    problemAtFace_(),
    newFacePatch_()
{}

// * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * * * * //

meshSurfaceTopologyImprove::~meshSurfaceTopologyImprove()
{
    clearMeshSurface();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool meshSurfaceTopologyImprove::improveTopology()
{
    bool changed(false), modified;

    label nIter(0);

    do
    {
        modified = false;

        optimizeCells();

        if( detectProblematicFaces() )
        {
            if( updateFacePatches() )
            {
                modified = true;
                changed = true;
            }

            if( analyseUpdatedPatches() )
            {
                modified = true;
                changed = true;
            }

//        if( modifyDeformedFaces() )
//        {
//            modified = true;
//            changed = true;
//        }
        }
    } while( modified && (++nIter < 3) );

    mesh_.write();
    ::exit(0);

    return changed;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
