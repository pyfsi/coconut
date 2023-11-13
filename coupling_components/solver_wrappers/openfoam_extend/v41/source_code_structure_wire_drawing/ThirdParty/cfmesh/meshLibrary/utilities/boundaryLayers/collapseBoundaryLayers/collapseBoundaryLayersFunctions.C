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

#include "collapseBoundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "meshSurfacePartitioner.H"
#include "detectBoundaryLayers.H"
#include "helperFunctions.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGLayer

# ifdef DEBUGLayer
#include "OFstream.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void collapseBoundaryLayers::detectLowQualityHairsAndCollapseFaces()
{
    Info << "Detecting low quality faces" << endl;
    const meshSurfaceEngine& mse = surfaceEngine();
    const VRWGraph& pFaces = mse.pointFaces();
    const VRWGraph& eFaces = mse.edgeFaces();
    const labelLongList& bp = mse.bp();
    const labelLongList& facePatches = mse.boundaryFacePatches();
    const faceList::subList& bFaces = mse.boundaryFaces();

    meshSurfacePartitioner mPart(mse);

    detectBoundaryLayers dbl(mPart);
    const edgeLongList& hairEdges = dbl.hairEdges();
    const VRWGraph& hairEdgesAtBndPoint = dbl.hairEdgesAtBndPoint();

    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();
    pointFieldPMG& points = meshModifier.pointsAccess();

    boolList collapseEdge(pFaces.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(collapseEdge, bpI)
            collapseEdge[bpI] = false;

        //- collapse hairs at features that are at high angle
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 50)
        # endif
        forAll(eFaces, beI)
        {
            if( eFaces.sizeOfRow(beI) == 2 )
            {
                //- check if this is a feature edge
                if( facePatches[eFaces(beI, 0)] != facePatches[eFaces(beI, 1)] )
                {
                    //- find the shared vertices
                    const edge e =
                        help::sharedEdge
                        (
                            bFaces[eFaces(beI, 0)],
                            bFaces[eFaces(beI, 1)]
                        );

                    if( hairEdgesAtBndPoint.sizeOfRow(bp[e.start()]) == 1 )
                    {
                        const label heI = hairEdgesAtBndPoint(bp[e.start()], 0);
                        collapseEdge[heI] = true;
                    }
                    if( hairEdgesAtBndPoint.sizeOfRow((bp[e.end()])) == 1 )
                    {
                        const label heI = hairEdgesAtBndPoint(bp[e.end()], 0);
                        collapseEdge[heI] = true;
                    }
                }
            }
            else if( eFaces.sizeOfRow(beI) == 1 )
            {
                //- inter-processor boundary
            }
        }

        //- modify faces with collapsed layers
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 50)
        # endif
        forAll(faces, faceI)
        {
            face& f = faces[faceI];

            forAll(f, eI)
            {
                const edge fe = f.faceEdge(eI);

                const label bps = bp[fe.start()];
                const label bpe = bp[fe.end()];

                if( (bps >= 0) && (hairEdgesAtBndPoint.sizeOfRow(bps) == 1) )
                {
                    forAllRow(hairEdgesAtBndPoint, bps, pheI)
                    {
                        const label heI = hairEdgesAtBndPoint(bps, pheI);

                        if( !collapseEdge[heI] )
                            continue;

                        const edge& he = hairEdges[heI];

                        if( he == fe )
                        {
                            //- collapse the edge
                            f[eI] = he.end();
                            points[he.end()] = points[he.start()];

                            if( mesh_.hasBackup(he.start()) )
                            {
                                point pOrig;
                                mesh_.getOrigPoint(he.start(), pOrig);
                                meshModifier.setBackupPoint(he.end(), pOrig);
                            }
                        }
                    }
                }
                else if
                (
                    (bpe >= 0) && (hairEdgesAtBndPoint.sizeOfRow(bpe) == 1)
                )
                {
                    forAllRow(hairEdgesAtBndPoint, bpe, pheI)
                    {
                        const label heI = hairEdgesAtBndPoint(bpe, pheI);

                        if( !collapseEdge[heI] )
                            continue;

                        const edge& he = hairEdges[heI];

                        if( he == fe )
                        {
                            //- collapse the edge
                            f[f.fcIndex(eI)] = he.end();
                            points[he.end()] = points[he.start()];

                            if( mesh_.hasBackup(he.start()) )
                            {
                                point pOrig;
                                mesh_.getOrigPoint(he.start(), pOrig);
                                meshModifier.setBackupPoint(he.end(), pOrig);
                            }
                        }
                    }
                }
            }
        }

        //- modify boundary faces attached to collapsed layers
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 50)
        # endif
        forAll(hairEdges, heI)
        {
            if( !collapseEdge[heI] )
                continue;

            const edge& he = hairEdges[heI];

            const label pointI = he.start();
            const label bps = bp[pointI];

            if( hairEdgesAtBndPoint.sizeOfRow(bps) != 1 )
                continue;

            forAllRow(pFaces, bps, pfI)
            {
                face& f = faces[mesh_.nInternalFaces()+pFaces(bps, pfI)];

                const label pos = f.which(pointI);

                if( pos < 0 )
                {
                    continue;
                }

                f[pos] = he.end();
            }
        }
    }
}

void collapseBoundaryLayers::modifyLayerCells()
{
    Info << "Cleaning collapsed faces from the mesh" << endl;

    //- remove faces that are collapsed into a line
    //- and purge duplicate vertices
    polyMeshGenModifier meshModifier(mesh_);
    faceListPMG& faces = meshModifier.facesAccess();
    cellListPMG& cells = meshModifier.cellsAccess();

    const labelLongList& owner = mesh_.owner();
    const labelLongList& neighbour = mesh_.neighbour();

    boolList removeFace(faces.size());
    boolList removeCell(cells.size());
    bool invalidCells(false);

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        //- mark faces for removal
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 50)
        # endif
        forAll(faces, faceI)
        {
            face& f = faces[faceI];

            DynList<label> copyFace;
            forAll(f, pI)
                copyFace.appendIfNotIn(f[pI]);

            if( copyFace.size() < f.size() )
            {
                //- face has been modified
                if( copyFace.size() > 2 )
                {
                    f.setSize(copyFace.size());
                    forAll(f, pI)
                        f[pI] = copyFace[pI];

                    removeFace[faceI] = false;
                }
                else
                {
                    //- delete the face
                    removeFace[faceI] = true;
                }
            }
            else
            {
                //- face is not modified
                removeFace[faceI] = false;
            }
        }

        //- cleanup cells
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 50)
        # endif
        forAll(cells, cellI)
        {
            cell& c = cells[cellI];

            DynList<label> cCopy;
            forAll(c, fI)
            {
                if( removeFace[c[fI]] )
                    continue;

                cCopy.append(c[fI]);
            }

            if( cCopy.size() < c.size() )
            {
                c.setSize(cCopy.size());
                forAll(c, cfI)
                    c[cfI] = cCopy[cfI];

                if( c.size() > 2 )
                {
                    removeCell[cellI] = false;
                }
                else
                {
                    //- cell is collapsed into a face
                    invalidCells = true;
                    removeCell[cellI] = true;

                    label internalFace(-1), bndFace(-1);
                    forAll(c, fI)
                    {
                        if( c[fI] < mesh_.nInternalFaces() )
                        {
                            if( internalFace >= 0 )
                                FatalError << "Internal face already found"
                                    << abort(FatalError);
                            internalFace = c[fI];
                        }
                        else
                        {
                            if( bndFace >= 0 )
                                FatalError << "Boundary face already found"
                                    << abort(FatalError);

                            bndFace = c[fI];
                        }
                    }

                    //- find the neighbour cell over the internal face
                    //- and use the boundary face in that cell
                    label neiCell = owner[internalFace];
                    if( neiCell == cellI )
                        neiCell = neighbour[internalFace];

                    cell& nc = cells[neiCell];

                    # ifdef USE_OMP
                    # pragma omp critical(modifyInternalCell)
                    # endif
                    {
                        forAll(nc, fI)
                            if( nc[fI] == internalFace )
                                nc[fI] = bndFace;
                    }

                    c.setSize(0);
                }
            }
            else
            {
                removeCell[cellI] = false;
            }
        }
    }

    reduce(invalidCells, maxOp<bool>());

    //- remove collapsed faces from the mesh
    meshModifier.removeFaces(removeFace);

    if( invalidCells )
        meshModifier.removeCells(removeCell);

    meshModifier.removeUnusedVertices();

    # ifdef DEBUGLayer
    forAll(cells, cellI)
    {
        DynList<edge> edges;
        DynList<label> nFacesAtEdge;

        const cell& c = cells[cellI];

        forAll(c, fI)
        {
            const face& f = faces[c[fI]];

            forAll(f, eI)
            {
                const edge e = f.faceEdge(eI);

                const label pos = edges.containsAtPosition(e);

                if( pos < 0 )
                {
                    edges.append(e);
                    nFacesAtEdge.append(1);
                }
                else
                {
                    ++nFacesAtEdge[pos];
                }
            }
        }

        bool closed = true;
        forAll(nFacesAtEdge, eI)
            if( nFacesAtEdge[eI] != 2 )
            {
                closed = false;
            }

        if( !closed )
        {
            Info << "Cell " << cellI << " is not closed" << endl;
            Info << "Cell faces " << c << endl;
            Info << "Cell edges " << edges << endl;
            Info << "Faces at edge " << nFacesAtEdge << endl;

            forAll(c, fI)
                Info << "Face " << c[fI] << " is " << faces[c[fI]] << endl;
        }
    }
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
