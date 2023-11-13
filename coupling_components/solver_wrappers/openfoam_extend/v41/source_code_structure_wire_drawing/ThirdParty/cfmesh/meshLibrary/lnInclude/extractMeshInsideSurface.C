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

#include "extractMeshInsideSurface.H"
#include "polyMeshGenAddressing.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "findCellsIntersectingSurface.H"
#include "helperFunctions.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// #define DEBUGInside

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void extractMeshInsideSurface::findIntersectedCells()
{
    const triSurf& surf = octree_.surface();

    const cellListPMG& cells = mesh_.cells();
    const faceListPMG& faces = mesh_.faces();
    const pointFieldPMG& points = mesh_.points();

    const vectorLongList& cellCentres = mesh_.addressingData().cellCentres();
    const vectorLongList& faceCentres = mesh_.addressingData().faceCentres();

    findCellsIntersectingSurface cis(mesh_, octree_);

    const boolList& intersectedCells = cis.intersectedCells();
    const VRWGraph& facetsInCell = cis.facetsIntersectingCells();

    # ifdef DEBUGInside
    const label csId = mesh_.addCellSubset("intersectedCells");
    forAll(intersectedCells, cellI)
    {
        if( intersectedCells[cellI] )
            mesh_.addCellToSubset(csId, cellI);
    }
    # endif

    # ifdef USE_OMP
    # pragma omp parallel if( facetsInCell.size() > 1000 )
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(interfaceFace_, faceI)
            interfaceFace_[faceI] = false;

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(intersectedCells, cI)
        {
            if( intersectedCells[cI] )
            {
                //- find face candidates
                const cell& c = cells[cI];

                //- find the length of the shortest edge
                scalar shortestEdgeSq(VGREAT);
                forAll(c, fI)
                {
                    const face& f = faces[c[fI]];

                    forAll(f, eI)
                    {
                        const edge e = f.faceEdge(eI);

                        const scalar dSq = magSqr(points[e[1]] - points[e[0]]);

                        shortestEdgeSq = min(dSq, shortestEdgeSq);
                    }
                }

                const point& cCentre = cellCentres[cI];

                //- find the centre gradient vector
                point nearestToCentre;
                scalar minDistSq(VGREAT);
                forAllRow(facetsInCell, cI, cfI)
                {
                    const label triI = facetsInCell(cI, cfI);

                    const point np =
                        help::nearestPointOnTheTriangle(triI, surf, cCentre);

                    const scalar dSq = magSqr(np - cCentre);

                    if( dSq < minDistSq )
                    {
                        nearestToCentre = np;
                        minDistSq = dSq;
                    }
                }

                bool foundIntersection(false);

                forAll(c, fI)
                {
                    const point& fc = faceCentres[c[fI]];

                    forAllRow(facetsInCell, cI, cfI)
                    {
                        const label triI = facetsInCell(cI, cfI);

                        point intersection;
                        const bool hasIntersection =
                            help::triLineIntersection
                            (
                                surf,
                                triI,
                                cCentre,
                                fc,
                                intersection
                            );

                        if( hasIntersection )
                        {
                            interfaceFace_[c[fI]] = true;
                            foundIntersection = true;
                            break;
                        }
                    }
                }

                if( !foundIntersection || (minDistSq / shortestEdgeSq < 1e-6) )
                {
                    //- mark all faces as interfaces
                    forAll(c, fI)
                        interfaceFace_[c[fI]] = true;
                }
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- unify status across the inter-processor boundaries
        const PtrList<processorBoundaryPatch>& procBnd = mesh_.procBoundaries();

        forAll(procBnd, patchI)
        {
            //- detect interface faces at this patch
            labelLongList toSend;

            const label start = procBnd[patchI].patchStart();
            const label end = start + procBnd[patchI].patchSize();

            for(label faceI=start;faceI<end;++faceI)
            {
                if( interfaceFace_[faceI] )
                    toSend.append(faceI-start);
            }

            //- send data to the other patch
            OPstream toOtherProc
            (
                Pstream::blocking,
                procBnd[patchI].neiProcNo(),
                toSend.byteSize()
            );

            toOtherProc << toSend;
        }

        forAll(procBnd, patchI)
        {
            const label start = procBnd[patchI].patchStart();

            //- receive data from other processor
            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBnd[patchI].neiProcNo()
            );

            labelList receivedData;

            fromOtherProc >> receivedData;

            forAll(receivedData, i)
                interfaceFace_[start+receivedData[i]] = true;
        }
    }

    # ifdef DEBUGInside
    const label sId = mesh_.addFaceSubset("interfaces");
    forAll(interfaceFace_, faceI)
    {
        if( interfaceFace_[faceI] )
            mesh_.addFaceToSubset(sId, faceI);
    }

    Info << "Save-am" << endl;
    mesh_.write();
    ::exit(0);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace extractMeshInsideHelper
{

class neighbourOp
{
    const polyMeshGen& mesh_;

    const boolList& atInterface_;

public:

    inline neighbourOp
    (
        const polyMeshGen& mesh,
        const boolList& atInterface
    )
    :
        mesh_(mesh),
        atInterface_(atInterface)
    {}

    inline label size() const
    {
        return mesh_.cells().size();
    }

    inline void operator()
    (
        const label cellI,
        DynList<label>& neighbourCells
    ) const
    {
        //- find neighbour cells of the current cell
        neighbourCells.clear();

        const labelLongList& owner = mesh_.owner();
        const labelLongList& neighbour = mesh_.neighbour();

        const cell& c = mesh_.cells()[cellI];

        forAll(c, fI)
        {
            //- do not cross selected interfaces
            if( atInterface_[c[fI]] )
                continue;

            label nei = owner[c[fI]];

            if( nei == cellI )
                nei = neighbour[c[fI]];

            if( nei >= 0 )
                neighbourCells.append(nei);
        }
    }

    template<class labelListType>
    void collectGroups
    (
        std::map<label, DynList<label> >& neiGroups,
        const labelListType& elementInGroup,
        const DynList<label>& localGroupLabel
    ) const
    {
        const PtrList<processorBoundaryPatch>& procBoundaries =
            mesh_.procBoundaries();
        const labelLongList& owner = mesh_.owner();

        //- send the data to other processors
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();
            const label size = procBoundaries[patchI].patchSize();

            labelList groupOwner(procBoundaries[patchI].patchSize());
            for(label faceI=0;faceI<size;++faceI)
            {
                const label groupI = elementInGroup[owner[start+faceI]];

                if( groupI < 0 )
                {
                    groupOwner[faceI] = -1;
                    continue;
                }

                groupOwner[faceI] = localGroupLabel[groupI];
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo(),
                groupOwner.byteSize()
            );

            toOtherProc << groupOwner;
        }

        //- receive data from other processors
        forAll(procBoundaries, patchI)
        {
            const label start = procBoundaries[patchI].patchStart();

            labelList receivedData;

            IPstream fromOtherProc
            (
                Pstream::blocking,
                procBoundaries[patchI].neiProcNo()
            );

            fromOtherProc >> receivedData;

            forAll(receivedData, faceI)
            {
                if( receivedData[faceI] < 0 )
                    continue;

                const label groupI = elementInGroup[owner[start+faceI]];

                if( groupI < 0 )
                    continue;

                DynList<label>& ng = neiGroups[localGroupLabel[groupI]];

                //- store the connection over the inter-processor boundary
                ng.appendIfNotIn(receivedData[faceI]);
            }
        }
    }
};

class selectorOp
{

public:

    selectorOp()
    {}

    bool operator()(const label) const
    {
        return true;
    }
};

class binaryOrEqOp
{
public:

    inline binaryOrEqOp()
    {}

    inline labelList operator()(const labelList& l1, const labelList& l2) const
    {
        if( l1.size() != l2.size() )
        {
            FatalErrorIn
            (
                "inline const labelList& operator()"
                "(labelList&, const labelList&)"
            ) << "Lists are not of the same size" << abort(FatalError);
        }

        labelList ret(l1);

        forAll(l1, i)
            ret[i] |= l2[i];

        return ret;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extractMeshInsideHelper

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void extractMeshInsideSurface::assignCellsToRegions()
{
    nRegions_ =
        help::groupMarking
        (
            cellInGroup_,
            extractMeshInsideHelper::neighbourOp(mesh_, interfaceFace_),
            extractMeshInsideHelper::selectorOp()
        );
}

void extractMeshInsideSurface::cleanupUnusedCells()
{
    const vectorLongList& cCentres = mesh_.addressingData().cellCentres();

    labelList regionType(nRegions_, 0);

    boolList removeCell(cCentres.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(cellInGroup_, cI)
        {
            //- initialize
            removeCell[cI] = false;

            const point& c = cCentres[cI];

            const label cLabel = octree_.findLeafContainingVertex(c);

            if( cLabel < 0 )
            {
                //- cell is outside of the octree
                regionType[cellInGroup_[cI]] |= meshOctreeCube::OUTSIDE;
            }
            else
            {
                const direction type = octree_.returnLeaf(cLabel).cubeType();

                if( type & meshOctreeCube::INSIDE )
                {
                    regionType[cellInGroup_[cI]] |= meshOctreeCube::INSIDE;
                }
                else if( type & meshOctreeCube::OUTSIDE )
                {
                    regionType[cellInGroup_[cI]] |= meshOctreeCube::OUTSIDE;
                }
                else if( type & meshOctreeCube::UNKNOWN )
                {
                    regionType[cellInGroup_[cI]] |= meshOctreeCube::UNKNOWN;
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            returnReduce(regionType, extractMeshInsideHelper::binaryOrEqOp());
        }

        # ifdef USE
        # pragma omp for schedule(static, 1)
        # endif
        forAll(cellInGroup_, cI)
        {
            if
            (
                regionType[cellInGroup_[cI]] &
                (meshOctreeCube::OUTSIDE|meshOctreeCube::UNKNOWN)
            )
            {
                //- select cell for removal
                removeCell[cI] = true;
            }
        }
    }

    # ifdef DEBUGInside
    const label cId = mesh_.addCellSubset("cellsToRemove");
    forAll(removeCell, cI)
    {
        if( removeCell[cI] )
            mesh_.addCellToSubset(cId, cI);
    }
    mesh_.write();
    Info << "Stopping" << endl;
    ::exit(0);
    # endif

    //- remove unnecessary cells from the mesh
    polyMeshGenModifier(mesh_).removeCells(removeCell);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

extractMeshInsideSurface::extractMeshInsideSurface
(
    polyMeshGen& mesh,
    const meshOctree& octree
)
:
    mesh_(mesh),
    octree_(octree),
    cellInGroup_(mesh.cells().size()),
    nRegions_(),
    interfaceFace_(mesh.faces().size())
{
    //- detect cells intersected by the surface mesh
    findIntersectedCells();

    //- assign cells to region that may be walked without crossing
    //- over intersected cells
    assignCellsToRegions();

    //- cleanup cells that are not inside the domain
    cleanupUnusedCells();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

extractMeshInsideSurface::~extractMeshInsideSurface()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
