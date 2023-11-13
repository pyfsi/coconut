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

#include "extrude2DMesh.H"
#include "vectorLongList.H"
#include "helperFunctions.H"
#include "boundaryLayers.H"
#include "extrudeLayer.H"
#include "refineBoundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "edgeMesh.H"

#include <map>

#include "ops.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool extrude2DMesh::checkSetup() const
{
    if( extrudeInPatchNormalDirection_ )
        return true;

    //- check if the extrusion axis is defined correctly
    if( mag(extrusionVector_) < VSMALL )
    {
        Warning << "Extrusion axis has zero length" << endl;

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
    bool impossibleExtrusion(false);

    const PtrList<boundaryPatch>& patches = mesh_.boundaries();
    const faceListPMG& faces = mesh_.faces();
    const pointFieldPMG& points = mesh_.points();

    forAll(patches, patchI)
    {
        if( patches[patchI].patchName() == extrusionPatch_ )
        {
            //- check the dot product between the face vector and the
            //- extrusion axis. Report error if the dot product is almost zero
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
                if( mag((fn/(mag(fn)+VSMALL)) & extrusionVector_) < 0.01 )
                    impossibleExtrusion = true;
            }
        }
    }

    //- check the ratio between the sum of face normals
    //- and the sum of face areas
    reduce(avgPatchNormal, sumOp<vector>());
    reduce(magPatchNormal, sumOp<scalar>());

    if( (mag(avgPatchNormal) / magPatchNormal) < 0.95 )
    {
        Warning << "The selected patch " << extrusionPatch_
                << " is not flat" << endl;

        return false;
    }

    avgPatchNormal /= (mag(avgPatchNormal) + VSMALL);

    reduce(impossibleExtrusion, maxOp<bool>());

    if( impossibleExtrusion )
    {
        Warning << "It is not possible to extrude the mesh at patch "
            << extrusionPatch_ << " in the direction "
            << extrusionVector_ << endl;

        return false;
    }

    return true;
}

void extrude2DMesh::extrudeMesh()
{
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();

    const label origNumCells = cells.size();

    //- insert the required number of sheets into the mesh
    //- generate a layer over this patch
    boundaryLayers(mesh_).addLayerForPatch(extrusionPatch_);

    if( !cellSubsetName_.empty() )
    {
        //- add new cells to a zone
        const label zId = mesh_.addCellSubset(cellSubsetName_);

        for(label cellI=origNumCells;cellI<cells.size();++cellI)
            mesh_.addCellToSubset(zId, cellI);
    }

    //- refine the layer into the required number of slices
    refineBoundaryLayers refLayer(mesh_);
    refLayer.setNumberOfLayersForPatch(extrusionPatch_, nSubdivisions_);
    refLayer.refineLayers();

    //- move the vertices to their correct locations
    VRWGraph newVerticesAtEdge = refLayer.newPointsAtHairEdges();
    edgeLongList hairEdges = refLayer.hairEdges();

    if( hairEdges.size() == 0 )
    {
        //- generate hair edges
        const labelList matchedIDs = mesh_.findPatches(extrusionPatch_);
        labelHashSet patchIds;
        forAll(matchedIDs, i)
        {
            patchIds.insert(matchedIDs[i]);
        }

        meshSurfaceEngine mse(mesh_);
        const labelLongList& facePatch = mse.boundaryFacePatches();
        const labelLongList& fOwner = mse.faceOwners();
        const labelLongList& bPoints = mse.boundaryPoints();
        const VRWGraph& pFaces = mse.pointFaces();

        forAll(pFaces, bpI)
        {
            bool hasExtrudedPatch(false);
            forAllRow(pFaces, bpI, pfI)
            {
                if( patchIds.found(facePatch[pFaces(bpI, pfI)]) )
                    hasExtrudedPatch = true;
            }

            if( !hasExtrudedPatch )
                continue;

            //- find edge candidates
            const label pointI = bPoints[bpI];

            DynList<edge> edgeCandidates;
            forAllRow(pFaces, bpI, pfI)
            {
                const cell& c = cells[fOwner[pFaces(bpI, pfI)]];

                forAll(c, fI)
                {
                    const face& f = faces[c[fI]];

                    const label pos = f.which(pointI);

                    if( pos < 0 )
                        continue;

                    edgeCandidates.appendIfNotIn(f.faceEdge(pos));
                    edgeCandidates.appendIfNotIn
                    (
                        f.faceEdge(f.rcIndex(pos)).reverseEdge()
                    );
                }
            }

            //- find edges in base faces
            DynList<edge> baseEdges;
            forAllRow(pFaces, bpI, pfI)
            {
                const label bfI = pFaces(bpI, pfI);
                if( patchIds.found(facePatch[bfI]) )
                {
                    const face& bf = faces[mesh_.nInternalFaces()+bfI];

                    const label pos = bf.which(pointI);

                    baseEdges.appendIfNotIn(bf.faceEdge(pos));
                    baseEdges.appendIfNotIn(bf.faceEdge(bf.rcIndex(pos)));
                }
            }

            //- eliminate base edges from edge candidates
            forAllReverse(edgeCandidates, i)
            {
                if( baseEdges.contains(edgeCandidates[i]) )
                    edgeCandidates.removeElement(i);
            }

            if( edgeCandidates.size() == 1 )
            {
                hairEdges.append(edgeCandidates[0]);
                newVerticesAtEdge.appendList(edgeCandidates[0]);
            }
        }
    }

    //- Contruct mesh modifier
    polyMeshGenModifier meshModifier(mesh_);

    //- calculate translation vector for each edge
    vectorLongList translationVector(hairEdges.size(), vector::zero);
    boolList validTranslationVector(hairEdges.size(), false);

    if( extrudeInPatchNormalDirection_ )
    {
        typedef std::map<label, vector> pointNormalsMap;
        pointNormalsMap pointNormal;
        std::map<label, pointNormalsMap> symmNormal;

        forAll(mesh_.boundaries(), patchI)
        {
            const boundaryPatch& patch = mesh_.boundaries()[patchI];

            if
            (
                (patch.patchType().find("symmetry") != word::npos) ||
                (patch.patchType().find("wedge") != word::npos)
            )
            {
                const label start = patch.patchStart();
                const label size = patch.patchSize();

                //- calculate face normal and update point normals
                for(label fI=0;fI<size;++fI)
                {
                    const face& f = faces[start+fI];

                    const vector n = help::faceAreaVector(points, f);

                    forAll(f, pI)
                    {
                        auto pIt = symmNormal.find(f[pI]);

                        if( pIt != symmNormal.end() )
                        {
                            std::map<label, vector>& pMap = pIt->second;
                            auto it = pMap.find(patchI);

                            if( it != pMap.end() )
                            {
                                it->second += n;
                            }
                            else
                            {
                                pMap[patchI] = n;
                            }
                        }
                        else
                        {
                            symmNormal[f[pI]][patchI] = n;
                        }
                    }
                }
            }

            if( patch.patchName() == extrusionPatch_ )
            {
                const label start = patch.patchStart();
                const label size = patch.patchSize();

                //- calculate face normal and update point normals
                for(label fI=0;fI<size;++fI)
                {
                    const face& f = faces[start+fI];

                    const vector n = help::faceAreaVector(points, f);

                    forAll(f, pI)
                    {
                        auto pIt = pointNormal.find(f[pI]);

                        if( pIt != pointNormal.end() )
                        {
                            pIt->second += n;
                        }
                        else
                        {
                            pointNormal[f[pI]] = n;
                        }
                    }
                }
            }
        }

        //- assign the translation vector to each hair edge
        forAll(hairEdges, heI)
        {
            const label s = hairEdges[heI].start();
            const label e = hairEdges[heI].end();

            pointNormalsMap::const_iterator it = pointNormal.find(s);

            //- skip starting points that cannot be detect
            //- this usually happens when the layer is not refined
            if( it == pointNormal.end() )
            {
                continue;
            }

            //- end point must not be in the extruded patch
            if( pointNormal.find(e) != pointNormal.end() )
                continue;

            //- make translation vector lie in the symmetry plane
            vector n = it->second;

            const auto symmIt = symmNormal.find(s);
            if( symmIt != symmNormal.end() )
            {
                const std::map<label, vector>& pMap = symmIt->second;

                if( pMap.size() == 1 )
                {
                    for(auto sIt=pMap.begin();sIt!=pMap.end();++sIt)
                    {
                        vector symmNormal = sIt->second;
                        symmNormal /= (mag(symmNormal) + VSMALL);

                        n -= (n & symmNormal) * symmNormal;
                    }
                }
                else if( pMap.size() == 2 )
                {
                    DynList<vector> nVecs;

                    for(auto sIt=pMap.begin();sIt!=pMap.end();++sIt)
                    {
                        nVecs.append(sIt->second);
                    }

                    vector symmNormal = nVecs[0] ^ nVecs[1];

                    if( nVecs.size() != 2 )
                    {
                        FatalError << "Invalid input map" << abort(FatalError);
                    }

                    n -= (n & symmNormal) * symmNormal;
                }
                else
                {
                    FatalError << "Cannot extrude mesh at corners"
                        << " with more than two symmetry planes"
                        << exit(FatalError);
                }
            }

            n /= (mag(n) + VSMALL);

            translationVector[heI] = n;
            validTranslationVector[heI] = true;
        }

        //- move start and end points of hair edges to their correct location
        forAll(hairEdges, heI)
        {
            if( !validTranslationVector[heI] )
                continue;

            const edge& he = hairEdges[heI];

            meshModifier.movePoint(he.end(), points[he.start()]);

            meshModifier.movePoint
            (
                he.start(),
                points[he.start()] + translationVector[heI] * extrusionDistance_
            );

            translationVector[heI] *= -1.0;
        }

        //- invert grading factor
        gradingFactor_ = 1.0 / (gradingFactor_ + VSMALL);
    }
    else
    {
        translationVector = extrusionVector_;
        validTranslationVector = true;
    }

    //- move the extruded points
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(hairEdges, hairEdgeI)
    {
        if( !validTranslationVector[hairEdgeI] )
            continue;
        if( newVerticesAtEdge.sizeOfRow(hairEdgeI) <= 2 )
        {
            continue;
        }

        const edge& he = hairEdges[hairEdgeI];
        const point& s = points[he.start()];

        const vector ev = translationVector[hairEdgeI] * extrusionDistance_;

        vector dv = ev / nSubdivisions_;
        if( mag(gradingFactor_ - 1.0) > SMALL )
        {
            dv =
                ev /
                (
                    (1.0 - pow(gradingFactor_, nSubdivisions_)) /
                    (1.0 - gradingFactor_)
                );
        }

        //- move the points to their new locations
        point currPoint = s;
        vector currDisp = vector::zero;
        forAllRow(newVerticesAtEdge, hairEdgeI, i)
        {
            const label pointI = newVerticesAtEdge(hairEdgeI, i);

            currPoint += currDisp;

            currDisp = dv;
            dv *= gradingFactor_;

            meshModifier.movePoint(pointI, currPoint);
        }
    }

    if( removeOriginalCells_ )
    {
        //- find and remove cells from the mesh
        boolList removeCell(cells.size(), false);

        for(label cellI=0;cellI<origNumCells;++cellI)
            removeCell[cellI] = true;

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
                        patch.patchName()=="defaultFaces"?
                        word(extrusionPatch_+"Copy"):patch.patchName(),
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
    }

    //- clear all obsolete data
    meshModifier.clearAll();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

extrude2DMesh::extrude2DMesh(polyMeshGen& mesh)
:
    mesh_(mesh),
    extrusionPatch_(),
    cellSubsetName_(),
    extrusionVector_(vector::zero),
    extrusionDistance_(0.0),
    gradingFactor_(1.0),
    nSubdivisions_(0),
    extrudeInPatchNormalDirection_(false),
    removeOriginalCells_(true)
{}

extrude2DMesh::~extrude2DMesh()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

void extrude2DMesh::setExtrusionVector(const vector v)
{
    extrusionDistance_ = mag(v);
    extrusionVector_ = v / (extrusionDistance_ + VSMALL);
}

void extrude2DMesh::extrusionLength(const scalar distance)
{
    extrusionDistance_ = distance;
}

void extrude2DMesh::setExtrusionPatch(const word& extrusionPatch)
{
    extrusionPatch_ = extrusionPatch;
}

void extrude2DMesh::setNumberOfSubdivisions(const label nSubdivisions)
{
    nSubdivisions_ = nSubdivisions;
}

void extrude2DMesh::setGradingFactor(const scalar gradingFactor)
{
    gradingFactor_ = gradingFactor;
}

void extrude2DMesh::setCellSubsetName(const word& cellSubsetName)
{
    cellSubsetName_ = cellSubsetName;
}

void extrude2DMesh::extrudeInPatchNormalDirection
(
    const bool extrudeInPatchNormalDirection
)
{
    extrudeInPatchNormalDirection_ = extrudeInPatchNormalDirection;
}

void extrude2DMesh::removeOriginalCells(const bool removeOrigCells)
{
    removeOriginalCells_ = removeOrigCells;
}

void extrude2DMesh::generateExtrudedMesh()
{
    //- check if the problem is well posed
    if( !checkSetup() )
    {
        FatalError << "Cannot extrude the mesh. Exitting.." << endl;

        return;
    }

    //- generate a 3D mesh from the 2D mesh
    extrudeMesh();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// ************************************************************************* //
