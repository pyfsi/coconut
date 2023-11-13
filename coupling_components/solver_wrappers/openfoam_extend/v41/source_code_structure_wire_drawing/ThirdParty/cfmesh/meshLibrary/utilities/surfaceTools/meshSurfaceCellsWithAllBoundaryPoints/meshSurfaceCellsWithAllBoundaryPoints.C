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

#include "meshSurfaceCellsWithAllBoundaryPoints.H"
#include "polyMeshGenAddressing.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceCellsWithAllBoundaryPoints::markProblematicCells()
{
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();

    boolList surfacePoint(mesh_.points().size());

    problematicCells_.clear();

    const label start = mesh_.nInternalFaces();
    label end = mesh_.faces().size();
    if( Pstream::parRun() )
        end = mesh_.procBoundaries()[0].patchStart();

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(surfacePoint, pI)
            surfacePoint[pI] = false;

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 50)
        # endif
        for(label faceI=start;faceI<end;++faceI)
        {
            const face& f = faces[faceI];

            forAll(f, pI)
                surfacePoint[f[pI]] = true;
        }


        if( Pstream::parRun() )
        {
            # ifdef USE_OMP
            # pragma omp single
            # endif
            {
                const Map<label>& globalToLocal =
                    mesh_.addressingData().globalToLocalPointAddressing();
                const VRWGraph& pointAtProcs =
                    mesh_.addressingData().pointAtProcs();
                const DynList<label>& neiProcs =
                    mesh_.addressingData().pointNeiProcs();

                std::map<label, labelLongList> exchangeData;
                forAll(neiProcs, i)
                    exchangeData[neiProcs[i]].clear();

                forAllConstIter(Map<label>, globalToLocal, it)
                {
                    const label pointI = it();

                    if( surfacePoint[pointI] )
                    {
                        forAllRow(pointAtProcs, pointI, i)
                        {
                            const label neiProc = pointAtProcs(pointI, i);

                            if( neiProc == Pstream::myProcNo() )
                                continue;

                            exchangeData[neiProc].append(it.key());
                        }
                    }
                }

                labelLongList receivedata;
                help::exchangeMap(exchangeData, receivedata);

                forAll(receivedata, i)
                    surfacePoint[globalToLocal[receivedata[i]]] = true;
            }
        }

        //- mark problematic cells
        std::set<label> localProblematicCells;

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 50)
        # endif
        forAll(cells, cellI)
        {
            bool allAtBnd(true);

            const cell& c = cells[cellI];

            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                forAll(f, pI)
                {
                    if( !surfacePoint[f[pI]] )
                    {
                        allAtBnd = false;
                        break;
                    }
                }

                if( !allAtBnd )
                    break;
            }

            if( allAtBnd )
                localProblematicCells.insert(cellI);
        }

        # ifdef USE_OMP
        # pragma omp critical(joinProblematicCells)
        # endif
        {
            forAllConstIter(std::set<label>, localProblematicCells, it)
                problematicCells_.append(*it);
        }
    }

    Info << "Found " << returnReduce(problematicCells_.size(), sumOp<label>())
         << " coarse cells!" << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
meshSurfaceCellsWithAllBoundaryPoints::meshSurfaceCellsWithAllBoundaryPoints
(
    const polyMeshGen& mesh
)
:
    mesh_(mesh),
    problematicCells_()
{
    markProblematicCells();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceCellsWithAllBoundaryPoints::~meshSurfaceCellsWithAllBoundaryPoints()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceCellsWithAllBoundaryPoints::findProblematicPoints
(
    labelHashSet& problematicPoints
)
{
    problematicPoints.clear();

    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();

    forAll(problematicCells_, cellI)
    {
        const cell& c = cells[problematicCells_[cellI]];

        forAll(c, fI)
        {
            const face& f = faces[c[fI]];

            forAll(f, pI)
                problematicPoints.insert(f[pI]);
        }
    }
}

void meshSurfaceCellsWithAllBoundaryPoints::createBackupOfProblematicPoints()
{
    polyMeshGen& mesh = const_cast<polyMeshGen&>(mesh_);

    polyMeshGenModifier meshModifier(mesh);

    labelHashSet badPoints(1000);

    findProblematicPoints(badPoints);

    forAllConstIter(labelHashSet, badPoints, it)
        meshModifier.backupPoint(it.key());
}

const labelLongList&
meshSurfaceCellsWithAllBoundaryPoints::cellsWithAllPointsAtTheBoundary() const
{
    return problematicCells_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
