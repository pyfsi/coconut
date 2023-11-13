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
#include "meshSurfacePartitioner.H"
#include "partTriMesh.H"
#include "meshSurfaceEngine.H"
#include "triSurfModifier.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void partTriMesh::createPointsAndTrias
(
    const List<direction>& useFace
)
{
    const labelLongList& facePatch = mPart_.boundaryFacePatches();
    const meshSurfaceEngine& meshSurface = mPart_.surfaceEngine();
    const pointFieldPMG& points = meshSurface.points();
    const VRWGraph& pointFaces = meshSurface.pointFaces();
    const vectorLongList& faceCentres = meshSurface.faceCentres();
    const labelLongList& bPoints = meshSurface.boundaryPoints();
    const labelLongList& bp = meshSurface.bp();
    const faceList::subList& bFaces = meshSurface.boundaryFaces();

    meshSurfacePointLabelInTriMesh_.setSize(bPoints.size());
    meshSurfacePointLabelInTriMesh_ = -1;
    labelLongList nodeLabelForFace(bFaces.size(), -1);

    label nTriPoints(0);
    forAll(bFaces, bfI)
    {
        if( useFace[bfI] )
        {
            const face& bf = bFaces[bfI];

            //- create a point in the face centre
            if( bf.size() > 3 )
                nodeLabelForFace[bfI] = nTriPoints++;

            //- create points at face points
            forAll(bf, pI)
            {
                const label bpI = bp[bf[pI]];

                if( meshSurfacePointLabelInTriMesh_[bpI] == -1 )
                    meshSurfacePointLabelInTriMesh_[bpI] = nTriPoints++;
            }

            //- create triangles
            if( bf.size() > 3 )
            {
                forAll(bf, eI)
                {
                    //- add a triangle connected to face centre
                    labelledTri tri
                    (
                        meshSurfacePointLabelInTriMesh_[bp[bf[eI]]],
                        meshSurfacePointLabelInTriMesh_[bp[bf.nextLabel(eI)]],
                        nodeLabelForFace[bfI],
                        facePatch[bfI]
                    );

                    surf_.appendTriangle(tri);

                    //- add a triangle for shape
                    labelledTri secondTri
                    (
                        meshSurfacePointLabelInTriMesh_[bp[bf[eI]]],
                        meshSurfacePointLabelInTriMesh_[bp[bf.nextLabel(eI)]],
                        meshSurfacePointLabelInTriMesh_[bp[bf.prevLabel(eI)]],
                        facePatch[bfI]
                    );

                    surf_.appendTriangle(secondTri);

                    //- add remaining triangles to become more sensitive
                    //- to rotation
                    const label nTri = bf.size() - 2;
                    for(label i=0;i<nTri;++i)
                    {
                        labelledTri shapeTri
                        (
                            meshSurfacePointLabelInTriMesh_[bp[bf[eI]]],
                            meshSurfacePointLabelInTriMesh_[bp[bf[(eI+i+1)%bf.size()]]],
                            meshSurfacePointLabelInTriMesh_[bp[bf[(eI+i+2)%bf.size()]]],
                            facePatch[bfI]
                        );

                        surf_.appendTriangle(shapeTri);
                    }
                }
            }
            else
            {
                //- face is a triangle
                labelledTri tri
                (
                    meshSurfacePointLabelInTriMesh_[bp[bf[0]]],
                    meshSurfacePointLabelInTriMesh_[bp[bf[1]]],
                    meshSurfacePointLabelInTriMesh_[bp[bf[2]]],
                    facePatch[bfI]
                );

                surf_.appendTriangle(tri);

                //- add a triangle for shape
                forAll(bf, eI)
                {
                    labelledTri secondTri
                    (
                        meshSurfacePointLabelInTriMesh_[bp[bf[eI]]],
                        meshSurfacePointLabelInTriMesh_[bp[bf.nextLabel(eI)]],
                        meshSurfacePointLabelInTriMesh_[bp[bf.prevLabel(eI)]],
                        facePatch[bfI]
                    );

                    surf_.appendTriangle(secondTri);
                }
            }
        }
    }

    //- add points
    triSurfModifier sMod(surf_);
    pointField& pts = sMod.pointsAccess();
    pts.setSize(nTriPoints);

    pointType_.setSize(nTriPoints);
    pointLabelInMeshSurface_.setSize(pts.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1) nowait
        # endif
        forAll(pointType_, pI)
            pointType_[pI] = NONE;

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(pointLabelInMeshSurface_, bpI)
            pointLabelInMeshSurface_[bpI] = -1;

        # ifdef USE_OMP
        # pragma omp for schedule(guided, 100)
        # endif
        forAll(meshSurfacePointLabelInTriMesh_, bpI)
        {
            const label npI = meshSurfacePointLabelInTriMesh_[bpI];

            if( npI != -1 )
            {
                pointLabelInMeshSurface_[npI] = bpI;
                pts[npI] = points[bPoints[bpI]];
                pointType_[npI] |= SMOOTH;

                forAllRow(pointFaces, bpI, pfI)
                    if( !useFace[pointFaces(bpI, pfI)] )
                    {
                        pointType_[npI] = BOUNDARY;
                        break;
                    }
            }
        }

        # ifdef USE_OMP
        # pragma omp for schedule(guided, 100) nowait
        # endif
        forAll(nodeLabelForFace, bfI)
        {
            const label npI = nodeLabelForFace[bfI];

            if( npI != -1 )
            {
                pts[npI] = faceCentres[bfI];
                pointType_[npI] |= FACECENTRE;
            }
        }

        # ifdef USE_OMP
        # pragma omp single nowait
        # endif
        {
            //- set CORNER and FEATUREEDGE flags to surface points
            forAllConstIter(labelHashSet, mPart_.corners(), it)
            {
                # ifdef USE_OMP
                # pragma omp task firstprivate(it)
                # endif
                {
                    const label pI = meshSurfacePointLabelInTriMesh_[it.key()];
                    if( pI != -1 )
                        pointType_[pI] |= CORNER;
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp single nowait
        # endif
        {
            forAllConstIter(labelHashSet, mPart_.edgePoints(), it)
            {
                # ifdef USE_OMP
                # pragma omp task firstprivate(it)
                # endif
                {
                    const label pI = meshSurfacePointLabelInTriMesh_[it.key()];
                    if( pI != -1 )
                        pointType_[pI] |= FEATUREEDGE;
                }
            }
        }
    }

    //- create addressing for parallel runs
    if( Pstream::parRun() )
    {
        createParallelAddressing();

        createBufferLayers();
    }

    //- calculate point facets addressing
    surf_.pointFacets();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
