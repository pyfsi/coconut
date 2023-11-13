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

#include "triSurfaceExtrude2DEdges.H"
#include "triSurfModifier.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceExtrude2DEdges::ensureCorrectNormalOrientation() const
{
    point pAvg(vector::zero);

    const pointField& pts = surf_.points();

    if( pts.size() == 0 )
        FatalErrorIn
        (
            "void triSurfaceExtrude2DEdges::"
            "ensureCorrectNormalOrientation() const"
        ) << "Surface does not contains a single point" << exit(FatalError);

    //- calculate the bounding box and check its size in the z direction
    boundBox bb(pts);
    if( (bb.max().z() - bb.min().z()) < VSMALL )
    {
        Info << "Edges are in the x-y plane" << endl;
        return;
    }

    const scalar dMax = bb.mag();

    //- calculate the average of all points
    forAll(pts, pI)
        pAvg += pts[pI];
    pAvg /= pts.size();

    //- calculate the covariance matrix
    tensor covarianceMatrix(tensor::zero);

    forAll(pts, pI)
    {
        const vector d = pts[pI] - pAvg / dMax;
        covarianceMatrix += d * d;
    }

    covarianceMatrix.xx() += SMALL;
    covarianceMatrix.yy() += SMALL;
    covarianceMatrix.zz() += SMALL;

    const vector eVal = eigenValues(symm(covarianceMatrix));
    vector x = eigenVector(covarianceMatrix, eVal[2]);
    vector y = eigenVector(covarianceMatrix, eVal[1]);

    vector n = eigenVector(covarianceMatrix, eVal[0]);
    n /= (mag(n) + VSMALL);

    if( magSqr(n - vector(0, 0, 1)) > SMALL )
    {
        Info << "Surface is not in the x-y plane. Transforming.." << endl;

        triSurfModifier sMod(const_cast<triSurf&>(surf_));

        pointField& newPts = sMod.pointsAccess();

        forAll(newPts, pI)
        {
            const point origP = newPts[pI];
            newPts[pI] = point((origP - pAvg) & x, (origP - pAvg) & y, 0.0);
        }

        Info << "Finished transforming into x-y plane" << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceExtrude2DEdges::triSurfaceExtrude2DEdges(const triSurf& surface)
:
    surf_(surface)
{
    ensureCorrectNormalOrientation();
}

triSurfaceExtrude2DEdges::~triSurfaceExtrude2DEdges()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceExtrude2DEdges::extrudeSurface(triSurf& newSurf) const
{
    triSurfModifier sMod(newSurf);

    //- subsets of feature edges shall become patches in the extruded surface
    DynList<label> subsetIds;
    surf_.edgeSubsetIndices(subsetIds);

    geometricSurfacePatchList& patches = sMod.patchesAccess();
    patches.setSize(subsetIds.size()+1);

    std::map<label, label> subsetToPatch;
    forAll(subsetIds, i)
    {
        const word sName = surf_.edgeSubsetName(subsetIds[i]);

        subsetToPatch[subsetIds[i]] = i;

        patches[i].name() = sName;
        patches[i].geometricType() = "patch";
    }

    //- set last patch name
    const label defaultId = patches.size() - 1;
    patches[defaultId].name() = "defaultPatch";
    patches[defaultId].geometricType() = "patch";


    //- check if the edges are in the x-y plane
    const pointField& sPoints = surf_.points();
    const boundBox bb(sPoints);

    if( Foam::mag(bb.max().z() - bb.min().z()) > SMALL )
        FatalErrorIn
        (
            "void triSurfaceExtrude2DEdges::extrudeSurface(triSurf&) const"
        ) << "Cannot extrude edges which are not in the x-y plane!"
          << exit(FatalError);

    //- copy points
    pointField& pts = sMod.pointsAccess();
    pts.setSize(2 * sPoints.size());

    const label nOffset = sPoints.size();
    const scalar zOffset = 0.01 * bb.mag();

    forAll(sPoints, pI)
    {
        pts[pI] = pts[pI+nOffset] = sPoints[pI];
        pts[pI+sPoints.size()].z() += zOffset;
    }

    //- create triangles from feature edges
    LongList<labelledTri>& triangles = sMod.facetsAccess();
    const edgeLongList& edges = surf_.featureEdges();

    triangles.setSize(2 * edges.size());
    forAll(edges, eI)
    {
        const edge& e = edges[eI];
        const label tI = 2 * eI;
        triangles[tI] = labelledTri(e[0], e[1], e[1]+nOffset, defaultId);
        triangles[tI + 1] =
            labelledTri(e[0], e[1]+nOffset, e[0]+nOffset, defaultId);
    }

    //- assign patches
    forAll(subsetIds, i)
    {
        labelLongList edgesInSubset;
        surf_.edgesInSubset(subsetIds[i], edgesInSubset);

        const label patchI = subsetToPatch[subsetIds[i]];

        forAll(edgesInSubset, eI)
        {
            const label tI = 2 * edgesInSubset[eI];

            triangles[tI].region() = patchI;
            triangles[tI+1].region() = patchI;
        }
    }

    //- check if there exists any faces in the default patch
    bool hasDefault(false);
    forAll(triangles, triI)
    {
        if( triangles[triI].region() == defaultId )
        {
            hasDefault = true;
            break;
        }
    }

    if( !hasDefault )
        patches.setSize(subsetIds.size());

    //- check for existence of _corners_ point subset
    const label psId = surf_.pointSubsetIndex("_corners_");
    if( psId >= 0 )
    {
        labelLongList pointsInSubset;
        surf_.pointsInSubset(psId, pointsInSubset);

        forAll(pointsInSubset, i)
        {
            newSurf.appendFeatureEdge
            (
                edge
                (
                    pointsInSubset[i],
                    pointsInSubset[i] + sPoints.size()
                )
            );
        }
    }
}

const triSurf* triSurfaceExtrude2DEdges::extrudeSurface() const
{
    triSurf* sPtr = new triSurf();

    extrudeSurface(*sPtr);

    return sPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
