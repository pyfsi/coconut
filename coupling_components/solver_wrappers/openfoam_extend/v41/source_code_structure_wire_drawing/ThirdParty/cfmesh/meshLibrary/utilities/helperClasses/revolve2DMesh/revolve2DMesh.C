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

#include "revolve2DMesh.H"
#include "scalarLongList.H"
#include "helperFunctions.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "detectBoundaryLayers.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"

#include <map>

//#define DEBUGRevolver

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool revolve2DMesh::checkSetup() const
{
    //- check if the rotation axis is defined correctly
    if( mag(rotationAxis_) < VSMALL )
    {
        Warning << "Rotation axis has zero length" << endl;

        return false;
    }

    //- check if there exist any cells in the mesh
    if( mesh_.cells().size() == 0 )
    {
        Warning << "The mesh has not cells. Please extrude the patch first"
                << endl;

        return false;
    }

    //- check if the input mesh is a 1D or 2D mesh
    //- there must exist patches of type empty and the sum of their normals
    //- shall be a zero vector
    //- In addition, check the dot product between the rotation axis and
    //- the normals of empty patches. The dot product must not be 1.
    vector avgPatchNormal(vector::zero);
    scalar magPatchNormal(0.0);
    bool impossibleRotation(false);

    const PtrList<boundaryPatch>& patches = mesh_.boundaries();
    const faceListPMG& faces = mesh_.faces();
    const pointFieldPMG& points = mesh_.points();

    forAll(patches, patchI)
    {
        if( patches[patchI].patchName() == revolvingPatch_ )
        {
            //- check the dot product between the face vector and the
            //- rotation axis. Report error if the dot product is almost one
            //- because it generates a poor quality mesh
            const label start = patches[patchI].patchStart();
            const label end = start + patches[patchI].patchSize();

            for(label faceI=start;faceI<end;++faceI)
            {
                const face& f = faces[faceI];

                const vector fn = help::faceAreaVector(points, f);

                //- update the sum of face normals
                avgPatchNormal += fn;
                magPatchNormal += mag(fn);

                //- check the projection between the rotation axis and the
                //- face normal
                if( mag((fn/(mag(fn)+VSMALL)) & rotationAxis_) > 0.99 )
                    impossibleRotation = true;
            }
        }
    }

    //- check the ratio between the sum of face normals
    //- and the sum of face areas
    reduce(avgPatchNormal, sumOp<vector>());
    reduce(magPatchNormal, sumOp<scalar>());

    if( (mag(avgPatchNormal) / magPatchNormal) < 0.95 )
    {
        Warning << "The selected patch " << revolvingPatch_
                << " is not flat" << endl;

        return false;
    }

    avgPatchNormal /= (mag(avgPatchNormal) + VSMALL);

    reduce(impossibleRotation, maxOp<bool>());

    if( impossibleRotation )
    {
        Warning << "It is not possible to rotate the mesh"
          << " around the given rotation axis" << endl;

        return false;
    }

    //- check whether the rotation axis is below or above the mesh
    //- it must not intersect the mesh
    scalar maxDist(-VGREAT), minDist(VGREAT);

    forAll(patches, patchI)
    {
        if( patches[patchI].patchName() == revolvingPatch_ )
        {
            //- check the dot product between the face vector and the
            //- rotation axis. Report error if the dot product is almost one
            //- because it generates a poor quality mesh
            const label start = patches[patchI].patchStart();
            const label end = start + patches[patchI].patchSize();

            for(label faceI=start;faceI<end;++faceI)
            {
                const face& f = faces[faceI];

                forAll(f, pI)
                {
                    const point& p = points[f[pI]];

                    //- calculate the signed distance from the rotation axis
                    const point np =
                        help::nearestPointOnLine(origin_, rotationAxis_, p);

                    const scalar dst =
                        (((p - np) ^ avgPatchNormal) & rotationAxis_);

                    minDist = min(minDist, dst);
                    maxDist = max(maxDist, dst);
                }
            }
        }
    }

    //- check if both min and max distances are of the same sign
    const scalar tol = SMALL * max(mag(minDist), mag(maxDist));

    if( mag(minDist) < tol )
        minDist = 0.0;
    if( mag(maxDist) < tol )
        maxDist = 0.0;

    if( (mag(minDist) < tol) && (mag(maxDist) < tol) )
    {
        Warning << "The patch area is zero in the extrusion direction" << endl;

        return false;
    }

    if( mag(minDist + maxDist) < (mag(minDist) + mag(maxDist)) )
    {
        Warning << "Revolving of patch " << revolvingPatch_
                << " around the axis " << rotationAxis_
                << " will result in a tangled mesh!" << endl;

        return false;
    }

    return true;
}

void revolve2DMesh::analyseRevolvingAngles()
{
    if( angleIntervals_.size() == 0 )
    {
        FatalErrorIn
        (
            "void revolve2DMesh::analyseRevolvingAngles()"
        ) << "It is no possible to revolve the mesh because the resolution"
          << " in the circumferential direction is not specified."
          << " Please use specifyIntervalResoluction and setCircumResolution"
          << " functions." << exit(FatalError);
    }

    const scalar pi_2 = 2.0 * M_PI;

    scalar minAngle(VGREAT), maxAngle(VGREAT);

    for
    (
        std::map<std::pair<scalar, scalar>, label>::iterator it=
            angleIntervals_.begin();
        it!=angleIntervals_.end();
        ++it
    )
    {
        const scalar sa = it->first.first;
        const scalar se = it->first.second;

        # ifdef DEBUGRevolver
        Info << "Interval start " << sa << " end " << se << endl;
        # endif

        if( mag(mag(se-sa) - pi_2) < SMALL )
            completeRevolution_ = true;

        minAngle = min(minAngle, sa);
        maxAngle = max(maxAngle, se);
    }

    if( completeRevolution_ )
    {
        //- modify interval angle to fit within the [0, 2*M_PI]
        DynList<std::pair<scalar, scalar> > eraseIntervals;
        for
        (
            std::map<std::pair<scalar, scalar>, label>::iterator it=
                angleIntervals_.begin();
            it!=angleIntervals_.end();
            ++it
        )
        {
            const scalar sa = it->first.first;
            const scalar se = it->first.second;

            if( sa < -SMALL && se < pi_2 )
            {
                //- split into two new intervals
                angleIntervals_.insert
                (
                    std::make_pair(std::make_pair(sa+pi_2, pi_2), it->second)
                );

                angleIntervals_.insert
                (
                    std::make_pair(std::make_pair(0., se), it->second)
                );

                eraseIntervals.append(it->first);
            }
            else if( sa >= 0.0 && se > (pi_2+SMALL) )
            {
                //- split into two new intervals
                angleIntervals_.insert
                (
                    std::make_pair(std::make_pair(sa, pi_2), it->second)
                );

                angleIntervals_.insert
                (
                    std::make_pair(std::make_pair(0., se-pi_2), it->second)
                );

                eraseIntervals.append(it->first);
            }
            else if( sa < -SMALL && se > (pi_2+SMALL) )
            {
                //- create a single interval
                angleIntervals_.insert
                (
                    std::make_pair(std::make_pair(0., se-pi_2), it->second)
                );

                eraseIntervals.append(it->first);
            }
        }

        forAll(eraseIntervals, i)
            angleIntervals_.erase(eraseIntervals[i]);
    }
    else
    {
        if( angleIntervals_.size() > 1 )
        {
            FatalErrorIn
            (
                "void revolve2DMesh::analyseRevolvingAngles()"
            ) << "It is not possible to specify more than one"
              << " resolution interval when the revolution is not spciefied"
              << " for all 360 degrees" << exit(FatalError);
        }

        std::map<std::pair<scalar, scalar>, label>::iterator it =
            angleIntervals_.begin();

        startAngle_ = it->first.first;
        const scalar se = it->first.second;
        const label resolution = it->second;

        angleIntervals_.erase(it);

        angleIntervals_.insert
        (
            std::make_pair
            (
                std::make_pair(0, min(se-startAngle_, pi_2)),
                resolution
            )
        );

        //- check if the setup create a wedge
        if( isWedge_ )
        {
            startAngle_ = pi_2 - 0.5 * (se - startAngle_);
        }
    }

    # ifdef DEBUGRevolver
    Info << "Number of intervals " << label(angleIntervals_.size()) << endl;
    std::map<std::pair<scalar, scalar>, label>::iterator it =
        angleIntervals_.begin();
    for(it=angleIntervals_.begin();it!=angleIntervals_.end();++it)
        Info << "Interval " << it->first.first << " " << it->first.second
             << " discretisation " << it->second << endl;
    # endif
}

void revolve2DMesh::revolveMesh()
{
    //- ensure that all intervals are within the [0, 2*M_PI]
    analyseRevolvingAngles();

    //- calculate the number of subdivisions of the revolving patch
    bool changed;
    do
    {
        changed = false;

        for
        (
            std::map<std::pair<scalar, scalar>, label>::iterator it=
                angleIntervals_.begin();
            it!=angleIntervals_.end();
            ++it
        )
        {
            const scalar sa = it->first.first;
            const scalar se = it->first.second;

            const scalar angle = (se - sa) / it->second;

            # ifdef DEBUGRevolver
            Info << "Validating interval (" << sa << " " << se << ")" << endl;
            Info << "Angle " << angle << endl;
            # endif

            std::map<std::pair<scalar, scalar>, label>::iterator oIt = it;
            ++oIt;
            for(;oIt!=angleIntervals_.end();++oIt)
            {
                const scalar osa = oIt->first.first;
                const scalar ose = oIt->first.second;

                const scalar oAngle = (ose - osa) / oIt->second;

                # ifdef DEBUGRevolver
                Info << "Other interval (" << osa << " " << ose << ")" << endl;
                Info << "Other angle " << oAngle << endl;
                # endif

                if( se > (osa + SMALL) )
                {
                    //- these intervals overlap
                    if( ose < se )
                    {
                        //- the other interval is fully
                        //- embedded in the first one
                        if( oAngle > angle )
                        {
                            //- the other interval has coarser resolution
                            //- hence it shall be deleted
                            angleIntervals_.erase(oIt);

                            changed = true;
                            break;
                        }
                        else
                        {
                            //- split the current interval into 2 new intervals
                            //- other interval is kept intact
                            std::pair<scalar, scalar> fInt(sa, osa);
                            const label fDiv = ceil((osa - sa) / angle);
                            angleIntervals_.insert(std::make_pair(fInt, fDiv));

                            std::pair<scalar, scalar> sInt(ose, se);
                            const label sDiv = ceil((se - ose) / angle);
                            angleIntervals_.insert(std::make_pair(sInt, sDiv));

                            //- delete the current interval and leave the other
                            angleIntervals_.erase(it);

                            changed = true;
                            break;
                        }
                    }
                    else
                    {
                        //- delete both intervals
                        angleIntervals_.erase(it);
                        angleIntervals_.erase(oIt);

                        //- split into three new intervals
                        std::pair<scalar, scalar> fInt(sa, osa);
                        const label fDiv = max(ceil((osa - sa) / angle), 1);
                        if( mag(osa - sa) > SMALL )
                            angleIntervals_.insert(std::make_pair(fInt, fDiv));

                        std::pair<scalar, scalar> sInt(osa, se);
                        const label sDiv =
                            max(ceil((se - osa) / min(angle, oAngle)), 1);
                        if( mag(se - osa) > SMALL )
                            angleIntervals_.insert(std::make_pair(sInt, sDiv));

                        std::pair<scalar, scalar> tInt(se, ose);
                        const label tDiv = max(ceil((ose - se) / oAngle), 1);
                        if( mag(ose - se) > SMALL )
                            angleIntervals_.insert(std::make_pair(tInt, tDiv));

                        changed = true;
                        break;
                    }
                }
            }

            //- stop the loop if the intervals have changed
            //- the iterators are not in the consistent state any more
            if( changed )
                break;
        }
    } while( changed );

    scalarLongList intervalAngles;

    scalar currAngle(0.0);
    intervalAngles.append(currAngle);

    for
    (
        std::map<std::pair<scalar, scalar>, label>::const_iterator it=
            angleIntervals_.begin();
        it!=angleIntervals_.end();
        ++it
    )
    {
        # ifdef DEBUGRevolver
        Info << "Angle interval start " << it->first.first << " second "
             << it->first.second << endl;
        # endif

        const scalar step = (it->first.second - it->first.first) / it->second;

        for(label i=0;i<it->second;++i)
        {
            currAngle += step;
            intervalAngles.append(currAngle);
        }
    }

    if( intervalAngles.size() == 1 )
    {
        FatalError << "Number of division is zero" << exit(FatalError);
    }

    //- smooth the interval angles until the difference is smaller than the
    //- prescribed tolerance
    label iter(0);
    do
    {
        changed = false;

        //- enforce grading in the positive direction
        for(label i=0;i<intervalAngles.size()-1;++i)
        {
            scalar prev, next;

            if( i == 0 )
            {
                if( completeRevolution_ )
                {
                    const label n = intervalAngles.size() - 1;

                    prev = intervalAngles[n] - intervalAngles[n-1];
                    next = intervalAngles[1] - intervalAngles[0];
                }
                else
                {
                    continue;
                }
            }
            else
            {
                prev = intervalAngles[i] - intervalAngles[i-1];
                next = intervalAngles[i+1] - intervalAngles[i];
            }

            if( next > maxExpansionRatio_ * prev )
            {
                if( next > maxExpansionRatio_ * (1.+maxExpansionRatio_) * prev )
                {
                    //- insert an additional point
                    const label size = intervalAngles.size();
                    intervalAngles.setSize(size+1);

                    for(label j=size-1;j>=i+1;--j)
                        intervalAngles[j+1] = intervalAngles[j];
                }

                if( (i+2) < intervalAngles.size() )
                {
                    intervalAngles[i+1] =
                        intervalAngles[i] + maxExpansionRatio_ * prev;

                    changed = true;
                }
            }
        }

        //- enforce grading in the negative direction
        for(label i=intervalAngles.size()-2;i>0;)
        {
            //- calculate the lengths of the intervals
            scalar curr, prev;

            if( i < (intervalAngles.size()-2) )
            {
                curr = intervalAngles[i+1] - intervalAngles[i];
                prev = intervalAngles[i+2] - intervalAngles[i+1];
            }
            else if( completeRevolution_ )
            {
                curr = intervalAngles[i+1] - intervalAngles[i];
                prev = intervalAngles[1] - intervalAngles[0];
            }
            else
            {
                //- the last element. Decrement the index and continue
                --i;

                continue;
            }

            //- check the size ratio and correct it if needed
            if( curr > maxExpansionRatio_ * prev )
            {
                if( curr > maxExpansionRatio_ * (1.+maxExpansionRatio_) * prev )
                {
                    //- insert an additional point
                    const label size = intervalAngles.size();
                    intervalAngles.setSize(size+1);

                    for(label j=size-1;j>=i+1;--j)
                        intervalAngles[j+1] = intervalAngles[j];

                    ++i;
                }

                intervalAngles[i] =
                    intervalAngles[i+1] - maxExpansionRatio_ * prev;

                changed = true;
            }

            //- decrement the index
            --i;
        }

    } while( changed && (++iter < 3) );

    //- insert the required number of sheets into the mesh
    //- generate a layer over this patch
    boundaryLayers(mesh_).addLayerForPatch(revolvingPatch_);

    const VRWGraph* newVerticesAtEdgePtr = NULL;
    const VRWGraph* hairEdgesAtPointsPtr = NULL;
    const edgeLongList* hairEdgesPtr = NULL;

    if( intervalAngles.size() > 2 )
    {
        //- refine the layer into the required number of slices
        refineBoundaryLayers refLayer(mesh_);

        Info << "Refining in the circumferential direction into "
             << (intervalAngles.size()-1) << " slices." << endl;

        refLayer.setNumberOfLayersForPatch
        (
            revolvingPatch_,
            intervalAngles.size()-1
        );

        refLayer.refineLayers();

        newVerticesAtEdgePtr = new VRWGraph(refLayer.newPointsAtHairEdges());
        hairEdgesAtPointsPtr = new VRWGraph(refLayer.hairEdgesAtPoint());
        hairEdgesPtr = new edgeLongList(refLayer.hairEdges());
    }
    else
    {
        meshSurfaceEngine mse(mesh_);
        meshSurfacePartitioner mPart(mse);

        detectBoundaryLayers dbl(mPart);

        const edgeLongList& hairEdges = dbl.hairEdges();

        //- find points in the active patch
        const faceListPMG& faces = mesh_.faces();
        boolList pointInExtrusionPatch(mesh_.points().size(), false);

        forAll(mesh_.boundaries(), patchI)
        {
            const boundaryPatch& patch = mesh_.boundaries()[patchI];

            if( patch.patchName() == revolvingPatch_ )
            {
                const label start = patch.patchStart();
                const label end = start + patch.patchSize();

                for(label faceI=start;faceI<end;++faceI)
                {
                    const face& f = faces[faceI];

                    forAll(f, pI)
                        pointInExtrusionPatch[f[pI]] = true;
                }
            }
        }

        VRWGraph* localNewVerticesAtEdgePtr = new VRWGraph(hairEdges.size());
        VRWGraph* localHairEdgesAtPointsPtr =
            new VRWGraph(mesh_.points().size());

        forAll(hairEdges, heI)
        {
            const edge& he = hairEdges[heI];

            localHairEdgesAtPointsPtr->append(he.start(), heI);

            if
            (
                pointInExtrusionPatch[he.start()] &&
                !pointInExtrusionPatch[he.end()]
            )
            {
                //- this edge is generate as a consequence of the boundary
                //- layer extrusion
                //- add the second point twice because the filter below
                //- skips edges with two or lest points!!
                localNewVerticesAtEdgePtr->append(heI, he.start());
                localNewVerticesAtEdgePtr->append(heI, he.end());
                localNewVerticesAtEdgePtr->append(heI, he.end());
            }
            else
            {
                localNewVerticesAtEdgePtr->append(heI, he.start());
                localNewVerticesAtEdgePtr->append(heI, he.end());
            }
        }

        //- create data
        hairEdgesPtr = new edgeLongList(hairEdges);
        newVerticesAtEdgePtr = localNewVerticesAtEdgePtr;
        hairEdgesAtPointsPtr = localHairEdgesAtPointsPtr;
    }

    //- move the vertices to the correct location
    const VRWGraph& newVerticesAtEdge = *newVerticesAtEdgePtr;
    const VRWGraph& hairEdgesAtPoints = *hairEdgesAtPointsPtr;
    const edgeLongList& hairEdges = *hairEdgesPtr;

    //- set the number of points in each plane
    VRWGraph pointsInPlane(intervalAngles.size());
    forAll(pointsInPlane, i)
        pointsInPlane.setRowSize(i, hairEdges.size());

    polyMeshGenModifier meshModifier(mesh_);
    pointFieldPMG& points = meshModifier.pointsAccess();

    std::set<label> collapsedEdges;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges, hairEdgeI)
    {
        if( newVerticesAtEdge.sizeOfRow(hairEdgeI) <= 2 )
            continue;

        const edge& he = hairEdges[hairEdgeI];
        const point& s = points[he.start()];

        const point nps = help::nearestPointOnLine(origin_, rotationAxis_, s);

        const scalar r = mag(s - nps);

        vector xAxis = s - nps;
        xAxis /= (mag(r) + VSMALL);

        if( r < VSMALL )
        {
            # ifdef USE_OMP
            # pragma omp critical(collapsedPoints)
            # endif
            collapsedEdges.insert(hairEdgeI);
        }
        else
        {
            //- calculate the y axis of the local coordinate system
            const vector yAxis = rotationAxis_ ^ xAxis;

            //- move the points to their new locations
            forAll(intervalAngles, i)
            {
                const label pointI = newVerticesAtEdge(hairEdgeI, i);

                //- set the points into the graph
                if( hairEdgeI < pointsInPlane.sizeOfRow(i) )
                    pointsInPlane(i, hairEdgeI) = pointI;

                const scalar phi = intervalAngles[i] + startAngle_;

                const point newP =
                    nps + cos(phi) * xAxis * r + sin(phi) * yAxis * r;

                meshModifier.movePoint(pointI, newP);
            }
        }
    }

    if( generatePointSubsets_ )
    {
        forAll(pointsInPlane, i)
        {
            const label sId =
                mesh_.addPointSubset("pointsInPlane_"+help::labelToText(i));

            forAllRow(pointsInPlane, i, j)
                mesh_.addPointToSubset(sId, pointsInPlane(i, j));
        }
    }

    //- remove remaining cells
    List<direction> revolvedPoint(mesh_.points().size(), direction(0));
    forAll(newVerticesAtEdge, heI)
    {
        const label rowSize = newVerticesAtEdge.sizeOfRow(heI);

        if( rowSize <= 2 )
            continue;

        forAllRow(newVerticesAtEdge, heI, i)
            revolvedPoint[newVerticesAtEdge(heI, i)] |= 1;

        revolvedPoint[newVerticesAtEdge(heI, rowSize-1)] |= 2;
    }

    //- find the cells that shall be removed from the mesh
    const cellListPMG& cells = mesh_.cells();
    const faceListPMG& faces = mesh_.faces();

    //- find and remove cells from the mesh
    boolList removeCell(cells.size(), false);

    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];

        forAll(c, fI)
        {
            const face& f = faces[c[fI]];

            forAll(f, pI)
            {
                if( !revolvedPoint[f[pI]] )
                {
                    removeCell[cellI] = true;

                    break;
                }
            }

            if( removeCell[cellI] )
                break;
        }
    }

    if( completeRevolution_ )
    {
                //- update boundary patches for faces
        if( rplPatches_.size() > 0 )
        {
            meshSurfaceEngine mse(mesh_);
            const VRWGraph& pointFaces = mse.pointFaces();
            const labelLongList& bp = mse.bp();
            const faceList::subList& bFaces = mse.boundaryFaces();

            labelLongList facePatch = mse.boundaryFacePatches();

            //- set the interval for each boundary point
            labelLongList pointInInterval(pointFaces.size(), -1);

            forAll(pointsInPlane, i)
            {
                forAllRow(pointsInPlane, i, j)
                {
                    const label bpI = bp[pointsInPlane(i, j)];

                    if( bpI < 0 )
                        continue;

                    pointInInterval[bpI] = i;
                }
            }

            //- patch index for each boundary face
            forAll(bFaces, bfI)
            {
                const face& bf = bFaces[bfI];
                scalar minAngle = VGREAT;
                scalar maxAngle = -VGREAT;

                bool skipPatchChecking(false);

                forAll(bf, pI)
                {
                    const label intervalI = pointInInterval[bp[bf[pI]]];

                    if( intervalI < 0 )
                    {
                        //- this vertex shall be merged with another one
                        skipPatchChecking = true;
                        continue;
                    }

                    const scalar angle = intervalAngles[intervalI];

                    minAngle = min(minAngle, angle);
                    maxAngle = max(maxAngle, angle);
                }

                if( skipPatchChecking )
                    continue;

                forAllConstIter(patchesMap, rplPatches_, it)
                {
                    if
                    (
                        minAngle > it->first.first &&
                        maxAngle < it->first.second
                    )
                    {
                        llMap::const_iterator pIt =
                            it->second.find(facePatch[bfI]);

                        if( pIt != it->second.end() )
                        {
                            facePatch[bfI] = pIt->second;
                        }
                    }
                }
            }

            VRWGraph newFaces;
            forAll(bFaces, bfI)
                newFaces.appendList(bFaces[bfI]);
            labelLongList faceOwner = mse.faceOwners();

            wordList patchNames(mesh_.boundaries().size());
            wordList patchTypes(mesh_.boundaries().size());

            forAll(mesh_.boundaries(), patchI)
            {
               patchNames[patchI] = mesh_.boundaries()[patchI].patchName();
               patchTypes[patchI] = mesh_.boundaries()[patchI].patchType();
            }

            meshModifier.replaceBoundary
            (
                patchNames,
                newFaces,
                faceOwner,
                facePatch
            );

            forAll(meshModifier.boundariesAccess(), patchI)
            {
                meshModifier.boundariesAccess()[patchI].patchType() =
                    patchTypes[patchI];
            }
        }

        //- find pairs of boundary faces that shall become internal faces
        LongList<std::pair<label, label> > facePairs;

        VRWGraph pFaces;
        pFaces.reverseAddressing(faces);

        forAll(mesh_.boundaries(), patchI)
        {
            const boundaryPatch& patch = mesh_.boundaries()[patchI];

            if( patch.patchName() != revolvingPatch_ )
                continue;

            const label start = patch.patchStart();
            const label end = start + patch.patchSize();

            for(label faceI=start;faceI<end;++faceI)
            {
                const face& bf = faces[faceI];

                DynList<label> revolvedFace;

                forAll(bf, pI)
                {
                    const label pointI = bf[pI];

                    forAllRow(hairEdgesAtPoints, pointI, i)
                    {
                        const label heI = hairEdgesAtPoints(pointI, i);
                        const edge& he = hairEdges[heI];

                        if( revolvedPoint[he.end()] & 2 )
                        {
                            revolvedFace.append(he.end());
                        }
                    }
                }

                if( revolvedFace.size() == bf.size() )
                {
                    //- find the label of the revolved face
                    forAllRow(pFaces, revolvedFace[0], pfI)
                    {
                        const label faceJ = pFaces(revolvedFace[0], pfI);
                        const face& of = faces[faceJ];

                        if( help::areFacesEqual(of, revolvedFace) )
                        {
                            facePairs.append
                            (
                                std::make_pair(faceI, faceJ)
                            );
                        }
                    }
                }
            }
        }

        //- update cells
        const labelLongList& owner = mesh_.owner();
        const labelLongList& neighbour = mesh_.neighbour();
        cellListPMG& cells = meshModifier.cellsAccess();

        forAll(facePairs, fpI)
        {
            const label fI = facePairs[fpI].first;
            const label ofI = facePairs[fpI].second;

            const label nei = neighbour[ofI];
            const label own = owner[ofI];

            //- swap the face in the cell that shall be removed
            label newNeighbour(-1);
            bool isNewNeiOldOwner(false);

            if( removeCell[nei] )
            {
                cell& cNei = cells[nei];

                forAll(cNei, i)
                    if( cNei[i] == ofI )
                        cNei[i] = fI;

                isNewNeiOldOwner = false;
                newNeighbour = own;
            }
            else
            {
                cell& cOwn = cells[own];

                forAll(cOwn, i)
                    if( cOwn[i] == ofI )
                        cOwn[i] = fI;

                isNewNeiOldOwner = true;
                newNeighbour = nei;
            }

            //- update the cell containing the face fI
            cell& c = cells[owner[fI]];
            forAll(c, i)
                if( c[i] == fI )
                    c[i] = ofI;

            //- reverse the faces
            if
            (
                ((owner[fI] < newNeighbour) && isNewNeiOldOwner) ||
                ((owner[fI] > newNeighbour) && !isNewNeiOldOwner)
            )
                meshModifier.facesAccess()[ofI] = faces[ofI].reverseFace();
        }

        //- reconnect faces
        faceListPMG& faces = meshModifier.facesAccess();

        forAll(revolvedPoint, pointI)
        {
            if( revolvedPoint[pointI] & 2 )
            {
                //- find the new point label that shall be replaced
                label pLabel(-1);

                forAllRow(hairEdgesAtPoints, pointI, i)
                {
                    const edge& he = hairEdges[hairEdgesAtPoints(pointI, i)];

                    if( he.end() == pointI )
                        pLabel = he.start();
                }

                if( pLabel < 0 )
                    FatalError << "Cannot find requested point"
                         << abort(FatalError);

                //- update vertices
                forAllRow(pFaces, pointI, pfI)
                {
                    face& f = faces[pFaces(pointI, pfI)];

                    const label pos = f.which(pointI);

                    f[pos] = pLabel;
                }
            }
        }
    }
    else if( collapsedEdges.size() )
    {
        //- some edges are at the rotation axis
        //- remove edges and faces with duplicates points
        faceListPMG& faces = meshModifier.facesAccess();

        VRWGraph pointFaces;
        pointFaces.reverseAddressing(points.size(), faces);

        forAllConstIter(std::set<label>, collapsedEdges, it)
        {
            const label heI = *it;

            const label s = hairEdges[heI].start();

            for(label i=1;i<newVerticesAtEdge.sizeOfRow(heI);++i)
            {
                const label npI = newVerticesAtEdge(heI, i);

                forAllRow(pointFaces, npI, pfI)
                {
                    face& f = faces[pointFaces(npI, pfI)];

                    const label pos = f.which(npI);

                    f[pos] = s;
                }
            }
        }

        //- remove faces collapsed to edges
        boolList removeFace(faces.size(), false);

        forAll(faces, faceI)
        {
            faces[faceI].collapse();

            if( faces[faceI].size() < 3 )
                removeFace[faceI] = true;
        }

        meshModifier.removeFaces(removeFace);
    }

    //- delete allocated data
    deleteDemandDrivenData(newVerticesAtEdgePtr);
    deleteDemandDrivenData(hairEdgesPtr);
    deleteDemandDrivenData(hairEdgesAtPointsPtr);

    meshModifier.removeCells(removeCell);

    //- remove vetices that are not used any more
    meshModifier.removeUnusedVertices();

    //- cleanup empty patches
    LongList<boundaryPatch*> newPatches;
    forAll(mesh_.boundaries(), patchI)
    {
        if( mesh_.boundaries()[patchI].patchSize() )
        {
            const boundaryPatch& patch = mesh_.boundaries()[patchI];

            newPatches.append
            (
                new boundaryPatch
                (
                    patch.patchName(),
                    patch.patchType()=="empty"?"patch":patch.patchType(),
                    patch.patchSize(),
                    patch.patchStart()
                )
            );
        }
    }

    //- update new patches
    meshModifier.boundariesAccess().setSize(newPatches.size());
    forAll(newPatches, patchI)
        meshModifier.boundariesAccess().set(patchI, newPatches[patchI]);

    //- is patch a wedge
    if( isWedge_ )
    {
        //- create a wedge patch
        VRWGraph newBndFaces;
        labelLongList newBndOwner;
        labelLongList newFacePatch;

        //- calculate the total number of boundary faces
        std::map<word, word> patchToType;
        std::map<label, label> patchToNewPatchId;
        label defaultFacesId(-1);
        label revolvePatchId(-1);
        label nBoundaryFaces(0);
        label newPatchId(1);

        //- patch names after this modification
        wordList newPatchNames;

        forAll(mesh_.boundaries(), patchI)
        {
            const boundaryPatch& patch = mesh_.boundaries()[patchI];

            patchToType[patch.patchName()] = patch.patchType();

            nBoundaryFaces += patch.patchSize();

            if( patch.patchName() == revolvingPatch_ )
            {
                revolvePatchId = patchI;
                patchToNewPatchId[patchI] = 0;
            }
            else if( patch.patchName() == "defaultFaces" )
            {
                defaultFacesId = patchI;
                patchToNewPatchId[patchI] = 0;
            }
            else
            {
                patchToNewPatchId[patchI] = newPatchId++;
            }
        }

        if( defaultFacesId != -1 )
        {
            newPatchNames.setSize(newPatchId);
            newPatchNames[0] = "wedge";
        }
        else
        {
            //- something went wrong
            FatalErrorIn
            (
                "void revolve2DMesh::revolveMesh()"
            ) << "Could not find the defaultFaces patch" << exit(FatalError);
        }

        //- calculate the number of faces in new patches
        labelList nFacesInPatch(newPatchNames.size(), 0);
        for
        (
            std::map<label, label>::const_iterator it=patchToNewPatchId.begin();
            it!=patchToNewPatchId.end();
            ++it
        )
        {
            if( it->second != 0 )
            {
                newPatchNames[it->second] =
                    mesh_.boundaries()[it->first].patchName();
            }

            nFacesInPatch[it->second] +=
                mesh_.boundaries()[it->first].patchSize();
        }

        //- calcute starting faces of patches
        labelList startFaceInPatch(newPatchNames.size());
        startFaceInPatch[0] = 0;
        for(label patchI=1;patchI<nFacesInPatch.size();++patchI)
        {
            startFaceInPatch[patchI] =
                startFaceInPatch[patchI-1] + nFacesInPatch[patchI-1];
        }

        //- fill the faces
        const labelLongList& owner = mesh_.owner();
        const cellListPMG& cells = mesh_.cells();

        newBndFaces.setSize(nBoundaryFaces);
        newBndOwner.setSize(nBoundaryFaces);
        newFacePatch.setSize(nBoundaryFaces);

        forAll(mesh_.boundaries(), patchI)
        {
            //- skip the defaultFaces, it is hadled with the revolving patch
            if( patchI == defaultFacesId )
                continue;

            const label newPatchId = patchToNewPatchId[patchI];
            const label newStart = startFaceInPatch[newPatchId];

            const boundaryPatch& patch = mesh_.boundaries()[patchI];
            const label start = patch.patchStart();
            const label size =  patch.patchSize();

            for(label fI=0;fI<size;++fI)
            {
                const face& f = faces[start+fI];

                newBndFaces[newStart+fI] = f;
                newBndOwner[newStart+fI] = owner[start+fI];
                newFacePatch[newStart+fI] = newPatchId;

                if( patchI == revolvePatchId )
                {
                    //- find the opposite face of the cyclic boundary
                    const cell& c = cells[owner[start+fI]];

                    label oppositeFace(-1);
                    forAll(c, cfI)
                    {
                        //- skip the current face
                        if( c[cfI] == (start+fI) )
                            continue;

                        //- the face shall not share an edge
                        //- with the current face
                        if( !help::shareAnEdge(f, faces[c[cfI]]) )
                        {
                            oppositeFace = c[cfI];
                            break;
                        }
                    }

                    newBndFaces[newStart+size+fI] = faces[oppositeFace];
                    newBndOwner[newStart+size+fI] = owner[oppositeFace];
                    newFacePatch[newStart+size+fI] = newPatchId;
                }
            }
        }

        //- replace the boundary
        meshModifier.replaceBoundary
        (
            newPatchNames,
            newBndFaces,
            newBndOwner,
            newFacePatch
        );

        //- set the patch types correctly
        meshModifier.boundariesAccess()[0].patchType() = "wedge";

        for(label patchI=1;patchI<newPatchNames.size();++patchI)
            meshModifier.boundariesAccess()[patchI].patchType() =
                patchToType[newPatchNames[patchI]];
    }

    //- clear all obsolete data
    meshModifier.clearAll();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

revolve2DMesh::revolve2DMesh(polyMeshGen& mesh)
:
    mesh_(mesh),
    revolvingPatch_(),
    rotationAxis_(vector::zero),
    origin_(vector::zero),
    maxExpansionRatio_(1.5),
    angleIntervals_(),
    rplPatches_(),
    startAngle_(0.0),
    completeRevolution_(false),
    generatePointSubsets_(false),
    isWedge_(false)
{}

revolve2DMesh::~revolve2DMesh()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

void revolve2DMesh::setRotationAxis(const vector& axis)
{
    rotationAxis_ = (axis / (mag(axis)+VSMALL));
}

void revolve2DMesh::setRevolvingPatch(const word& pName)
{
    revolvingPatch_ = pName;
}

void revolve2DMesh::setOrigin(const point& origin)
{
    origin_ = origin;
}

void revolve2DMesh::setMaxRatio(const scalar maxExpansionRatio)
{
    maxExpansionRatio_ = maxExpansionRatio;
}

void revolve2DMesh::clearAngleIntervals()
{
    angleIntervals_.clear();
}

void revolve2DMesh::setCircumResolution(const label nDivisions)
{
    if( nDivisions < 2 )
    {
        Warning << "Number of subdivisions in the circumferential direction is "
                << nDivisions << ". Please consider a larger number" << endl;

        return;
    }

    angleIntervals_.insert
    (
        std::make_pair(std::make_pair(0, 2.0*M_PI), nDivisions)
    );
}

void revolve2DMesh::setCircumResolution(const scalar resolutionAngle)
{
    const label nDivisions = ceil(2.0 * M_PI / resolutionAngle);

    angleIntervals_.insert
    (
        std::make_pair(std::make_pair(0, 2.0*M_PI), nDivisions)
    );
}

void revolve2DMesh::setIntervalResolution
(
    const scalar startAngle,
    const scalar endAngle,
    const label nDivisions
)
{
    if( nDivisions < 1 )
    {
        Warning << "Number of subdivisions in the circumferential direction is "
                << nDivisions << ". Please consider a larger number" << endl;

        return;
    }

    if( endAngle < startAngle )
    {
        const std::pair<scalar, scalar> pp
        (
            startAngle,
            endAngle + 2.0 * M_PI
        );

        std::map<std::pair<scalar, scalar>, label>::iterator it =
            angleIntervals_.find(pp);

        if( it == angleIntervals_.end() )
        {
            angleIntervals_.insert(std::make_pair(pp, nDivisions));
        }
        else
        {
            it->second = max(it->second, nDivisions);
        }
    }
    else
    {
        const std::pair<scalar, scalar> pp(startAngle, endAngle);

        std::map<std::pair<scalar, scalar>, label>::iterator it =
            angleIntervals_.find(pp);

        if( it == angleIntervals_.end() )
        {
            angleIntervals_.insert(std::make_pair(pp, nDivisions));
        }
        else
        {
            it->second = max(it->second, nDivisions);
        }
    }
}

void revolve2DMesh::setIntervalResolution
(
    const scalar startAngle,
    const scalar endAngle,
    const scalar resolutionAngle
)
{
    if( endAngle < startAngle )
    {
        const label nDivisions =
            ceil(((endAngle+2.0*M_PI) - startAngle) / resolutionAngle);

        const std::pair<scalar, scalar> pp
        (
            startAngle,
            endAngle + 2.0 * M_PI
        );

        std::map<std::pair<scalar, scalar>, label>::iterator it =
            angleIntervals_.find(pp);

        if( it == angleIntervals_.end() )
        {
            angleIntervals_.insert(std::make_pair(pp, nDivisions));
        }
        else
        {
            it->second = max(it->second, nDivisions);
        }
    }
    else
    {
        const label nDivisions =
            ceil((endAngle - startAngle) / resolutionAngle);

        const std::pair<scalar, scalar> pp(startAngle, endAngle);

        std::map<std::pair<scalar, scalar>, label>::iterator it =
            angleIntervals_.find(pp);

        if( it == angleIntervals_.end() )
        {
            angleIntervals_.insert(std::make_pair(pp, nDivisions));
        }
        else
        {
            it->second = max(it->second, nDivisions);
        }
    }
}

void revolve2DMesh::setIntervalPatches
(
    const scalar startAngle,
    const scalar endAngle,
    const word origPatchName,
    const word newPatchName,
    const word newPatchType
)
{
    //- get patch indices first
    label origPatchId(-1), newPatchId(-1);

    forAll(mesh_.boundaries(), patchI)
    {
        const word& pName = mesh_.boundaries()[patchI].patchName();

        if( pName == origPatchName )
            origPatchId = patchI;

        if( pName == newPatchName )
            newPatchId = patchI;
    }

    if( origPatchId < 0 )
    {
        Warning << "Cannot rename non-existing patch " << origPatchName << endl;
        return;
    }

    if( newPatchId < 0 )
    {
        newPatchId = mesh_.boundaries().size();

        polyMeshGenModifier mModifier(mesh_);
        PtrList<boundaryPatch>& patches = mModifier.boundariesAccess();

        label startFace = mesh_.faces().size();
        if( Pstream::parRun() )
            startFace = mesh_.procBoundaries()[0].patchStart();

        patches.setSize(newPatchId+1);

        patches.set
        (
            newPatchId,
            new boundaryPatch(newPatchName, newPatchType, 0, startFace)
        );
    }

    //- set the patch interval
    rplPatches_[std::make_pair(startAngle, endAngle)][origPatchId] = newPatchId;
}

void revolve2DMesh::setWedgeBoundaryConditions()
{
    isWedge_ = true;
}

void revolve2DMesh::createPointSubsets()
{
    generatePointSubsets_ = true;
}

void revolve2DMesh::generateRevolvedMesh()
{
    //- check if the problem is well posed
    if( !checkSetup() )
    {
        FatalError << "Cannot revolve the mesh. Exitting.." << endl;

        return;
    }

    //- generate a 3D mesh from the 2D mesh
    revolveMesh();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
