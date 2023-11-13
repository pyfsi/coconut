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

#include "triSurfaceImportSurfaceAsSubset.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceImportSurfaceAsSubset::triSurfaceImportSurfaceAsSubset(triSurf& surface)
:
    surf_(surface)
{}

triSurfaceImportSurfaceAsSubset::~triSurfaceImportSurfaceAsSubset()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceImportSurfaceAsSubset::addSurfaceAsSubset
(
    const triSurf& importSurf,
    const word& subsetName,
    const scalar angleTol,
    const scalar distanceTol
)
{
    const scalar normalTolerance = Foam::cos(angleTol*M_PI/180.0);

    const pointField& points = surf_.points();
    const vectorField& fNormals = surf_.facetNormals();
    const vectorField& fCentres = surf_.facetCentres();

    const vectorField& importFaceNormals = importSurf.facetNormals();

    boolList inSubset(surf_.size(), false);

    const scalar surfDiagonalSq =
         min(sqr(distanceTol), 0.0001 * sqr(boundBox(surf_.points()).mag()));

    meshOctree otherSurfOctree(importSurf);
    meshOctreeCreator(otherSurfOctree).createOctreeWithRefinedBoundary(20, 15);

    //- search for nearest facets in the import surface
    DynList<label> containedTriangles;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40) private(containedTriangles)
    # endif
    forAll(surf_, triI)
    {
        const labelledTri& tri = surf_[triI];

        vector tNormal = fNormals[triI];
        const scalar magTNormal = mag(tNormal);

        //- ignore sliver triangles
        if( magTNormal < VSMALL )
            continue;

        tNormal /= (magTNormal + VSMALL);

        //- find the search range
        scalar rSq(0.0);
        const point& fc = fCentres[triI];

        forAll(tri, pI)
        {
            const point& p = points[tri[pI]];

            rSq = max(rSq, magSqr(p - fc));
        }

        rSq = 0.01 * Foam::min(rSq, surfDiagonalSq);
        const scalar r = sqrt(rSq);

        const boundBox bb(fc-point(r, r, r), fc+point(r, r, r));

        //- find the nearest triangle in the surface which shall be imported
        containedTriangles.clear();
        otherSurfOctree.findTrianglesInBox(bb, containedTriangles);

        label nt(-1);
        scalar dSq(VGREAT);

        //- find the nearest triangle to the current one
        forAll(containedTriangles, ctI)
        {
            const label oTriI = containedTriangles[ctI];

            vector ctNormal = importFaceNormals[oTriI];
            const scalar magCtNormal = mag(ctNormal);

            //- ignore sliver triangles
            if( magCtNormal < VSMALL )
                continue;

            ctNormal /= (magCtNormal + VSMALL);

            //- ignore triangles with normal deviation above the max threshold
            if( mag(tNormal & ctNormal) < normalTolerance )
                continue;

            const point p =
                help::nearestPointOnTheTriangle(oTriI, importSurf, fc);

            const scalar distSq  = magSqr(p - fc);

            if( distSq < dSq )
            {
                nt = oTriI;
                dSq = distSq;
            }
        }

        //- check if the triangle has been found
        if( (nt < 0) || (dSq > rSq) )
            continue;

        //- select the triangle
        inSubset[triI] = true;
    }

    //- create a facet subset in the surface mesh and add the facets into it
    const label subsetId = surf_.addFacetSubset(subsetName);

    forAll(inSubset, triI)
    {
        if( !inSubset[triI] )
            continue;

        surf_.addFacetToSubset(subsetId, triI);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
