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
#include "meshOptimizer.H"
#include "polyMeshGenAddressing.H"
#include "meshSurfaceEngine.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOptimizer::laplaceSmoother::laplacian
(
    const labelLongList& smoothPoints,
    const label nIterations
)
{
    const VRWGraph& pPoints = mesh_.addressingData().pointPoints();
    const pointFieldPMG& points = mesh_.points();
    polyMeshGenModifier meshModifier(mesh_);

    for(label iterationI=0;iterationI<nIterations;++iterationI)
    {
        labelLongList procPoints;

        forAll(smoothPoints, i)
        {
            const label pointI = smoothPoints[i];

            if( vertexLocation_[pointI] & LOCKED )
                continue;

            if( vertexLocation_[pointI] & PARALLELBOUNDARY )
            {
                procPoints.append(pointI);

                continue;
            }

            vector newP(vector::zero);

            const label nPointPoints = pPoints.sizeOfRow(pointI);

            if( nPointPoints == 0 )
                return;

            for(label pI=0;pI<nPointPoints;++pI)
                newP += points[pPoints(pointI, pI)];

            newP /= pPoints.sizeOfRow(pointI);
            meshModifier.movePoint(pointI, newP);
        }

        laplacianParallel(procPoints, false);
    }

    updateMeshGeometry(smoothPoints);
}

void meshOptimizer::laplaceSmoother::laplacianSurface
(
    const labelLongList& smoothPoints,
    const label nIterations
)
{
    const VRWGraph& pPoints = mesh_.addressingData().pointPoints();
    const pointFieldPMG& points = mesh_.points();
    polyMeshGenModifier meshModifier(mesh_);

    for(label iterationI=0;iterationI<nIterations;++iterationI)
    {
        labelLongList procPoints;

        forAll(smoothPoints, i)
        {
            const label pointI = smoothPoints[i];

            if( vertexLocation_[pointI] & LOCKED )
                continue;

            if( vertexLocation_[pointI] & PARALLELBOUNDARY )
            {
                procPoints.append(pointI);

                continue;
            }

            vector newP(vector::zero);

            label counter(0);
            forAllRow(pPoints, pointI, pI)
            {
                const label pLabel = pPoints(pointI, pI);
                if( vertexLocation_[pLabel] & INSIDE )
                    continue;

                newP += points[pLabel];
                ++counter;
            }

            if( counter != 0 )
            {
                newP /= counter;
                meshModifier.movePoint(pointI, newP);
            }
        }

        laplacianParallel(smoothPoints, true);
    }

    updateMeshGeometry(smoothPoints);
}

void meshOptimizer::laplaceSmoother::laplacianPC
(
    const labelLongList& smoothPoints,
    const label nIterations
)
{
    const VRWGraph& pointCells = mesh_.addressingData().pointCells();
    const vectorLongList& centres = mesh_.addressingData().cellCentres();
    polyMeshGenModifier meshModifier(mesh_);

    for(label iterationI=0;iterationI<nIterations;++iterationI)
    {
        labelLongList procPoints;

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 20)
        # endif
        forAll(smoothPoints, i)
        {
            const label pointI = smoothPoints[i];

            if( vertexLocation_[pointI] & LOCKED )
                continue;

            if( pointCells.sizeOfRow(pointI) == 0 )
            {
                FatalErrorIn
                (
                    "void meshOptimizer::laplaceSmoother::laplacianPC"
                    "(const labelLongList&, const label)"
                ) << "Point " << pointI << " has not cells attached to it!!"
                  << abort(FatalError);
            }

            if( vertexLocation_[pointI] & PARALLELBOUNDARY )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                procPoints.append(pointI);

                continue;
            }

            point newP(vector::zero);
            forAllRow(pointCells, pointI, pcI)
                newP += centres[pointCells(pointI, pcI)];

            newP /= pointCells.sizeOfRow(pointI);

            meshModifier.movePoint(pointI, newP);
        }

        laplacianPCParallel(procPoints);

        updateMeshGeometry(smoothPoints);
    }
}

void meshOptimizer::laplaceSmoother::laplacianWPC
(
    const labelLongList& smoothPoints,
    const label nIterations
)
{
    const VRWGraph& pointCells = mesh_.addressingData().pointCells();
    const vectorLongList& centres = mesh_.addressingData().cellCentres();
    const scalarLongList& volumes = mesh_.addressingData().cellVolumes();

    polyMeshGenModifier meshModifier(mesh_);

    for(label iterationI=0;iterationI<nIterations;++iterationI)
    {
        labelLongList procPoints;

        # ifdef USE_OMP
        # pragma omp parallel for schedule(guided, 20)
        # endif
        forAll(smoothPoints, i)
        {
            const label pointI = smoothPoints[i];

            if( vertexLocation_[pointI] & LOCKED )
                continue;

            if( pointCells.sizeOfRow(pointI) == 0 )
            {
                FatalErrorIn
                (
                    "void meshOptimizer::laplaceSmoother::laplacianWPC"
                    "(const labelLongList&, const label)"
                ) << "Point " << pointI << " has not cells attached to it!!"
                  << abort(FatalError);
            }

            if( vertexLocation_[pointI] & PARALLELBOUNDARY )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                procPoints.append(pointI);

                continue;
            }

            point newP(vector::zero);
            scalar sumWeights(0.0);
            forAllRow(pointCells, pointI, pcI)
            {
                const label cellI = pointCells(pointI, pcI);
                const scalar w = Foam::mag(volumes[cellI]) + VSMALL;
                newP += w * centres[cellI];
                sumWeights += w;
            }

            newP /= sumWeights;
            meshModifier.movePoint(pointI, newP);
        }

        laplacianWPCParallel(procPoints);

        updateMeshGeometry(smoothPoints);
    }
}

void meshOptimizer::laplaceSmoother::laplacianWHPC
(
    const labelLongList& smoothPoints,
    const label nIterations
)
{
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();
    const labelLongList& owner = mesh_.owner();
    const labelLongList& neighbour = mesh_.neighbour();
    const VRWGraph& pointCells = mesh_.addressingData().pointCells();
    const vectorLongList& centres = mesh_.addressingData().cellCentres();
    const scalarLongList& volumes = mesh_.addressingData().cellVolumes();

    polyMeshGenModifier meshModifier(mesh_);
    scalarField maxVol(volumes.size()), minVol(volumes.size());

    boolList lockedCells(cells.size());
    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 100)
    # endif
    forAll(cells, cellI)
    {
        const cell& c = cells[cellI];

        bool lockedCell(true);

        forAll(c, fI)
        {
            const face& f = faces[c[fI]];

            if( !lockedCell )
                break;

            forAll(f, pI)
            {
                if( !(vertexLocation_[f[pI]] & LOCKED) )
                    lockedCell = false;
            }
        }

        lockedCells[cellI] = lockedCell;
    }

    for(label iterationI=0;iterationI<nIterations;++iterationI)
    {
        labelLongList procPoints;

        # ifdef USE_OMP
        # pragma omp parallel
        # endif
        {
            # ifdef USE_OMP
            # pragma omp for schedule(guided, 20)
            # endif
            forAll(maxVol, cellI)
            {
                maxVol[cellI] = -mag(volumes[cellI]);
                minVol[cellI] = mag(volumes[cellI]);

                const cell& c = cells[cellI];

                scalar weight = 1.0;

                forAll(c, fI)
                {
                    const label faceI = c[fI];

                    if( neighbour[faceI] >= 0 )
                    {
                        label nei = neighbour[faceI];
                        if( nei == cellI )
                            nei = owner[faceI];

                        //- additional weighting when the locked cells
                        //- are smaller than the current cell
                        const scalar neiVol = volumes[nei];
                        if
                        (
                            lockedCells[nei] &&
                            (mag(neiVol) < mag(volumes[cellI]))
                        )
                        {
                            weight =
                                max
                                (
                                    weight,
                                    3.0 * mag(volumes[cellI] / (neiVol+VSMALL))
                                );
                        }

                        maxVol[cellI] = max(maxVol[cellI], neiVol);
                        minVol[cellI] = min(minVol[cellI], neiVol);
                    }
                }

                maxVol[cellI] *= weight;
            }

            # ifdef USE_OMP
            # pragma omp for schedule(guided, 20)
            # endif
            forAll(smoothPoints, i)
            {
                const label pointI = smoothPoints[i];

                if( vertexLocation_[pointI] & LOCKED )
                    continue;

                if( pointCells.sizeOfRow(pointI) == 0 )
                    continue;

                if( vertexLocation_[pointI] & PARALLELBOUNDARY )
                {
                    # ifdef USE_OMP
                    # pragma omp critical
                    # endif
                    procPoints.append(pointI);

                    continue;
                }

                point newP(vector::zero);
                scalar sumWeights(0.0);
                forAllRow(pointCells, pointI, pcI)
                {
                    const label cellI = pointCells(pointI, pcI);
                    scalar w = mag(maxVol[cellI] / (minVol[cellI] + VSMALL));
                    if( volumes[cellI] < 0.0 )
                        w = VSMALL;

                    newP += w * centres[cellI];
                    sumWeights += w;
                }

                newP /= sumWeights;

                meshModifier.movePoint(pointI, newP);
            }
        }

        laplacianWHPCParallel(procPoints);

        updateMeshGeometry(smoothPoints);
    }
}

void meshOptimizer::laplaceSmoother::updateMeshGeometry
(
    const labelLongList& smoothPoints
)
{
    const cellListPMG& cells = mesh_.cells();
    const VRWGraph& pointCells = mesh_.addressingData().pointCells();

    boolList chF(mesh_.faces().size());

    # ifdef USE_OMP
    # pragma omp parallel if( smoothPoints.size() > 100 )
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(chF, faceI)
            chF[faceI] = false;

        # ifdef USE_OMP
        # pragma omp for schedule(guided, 20)
        # endif
        forAll(smoothPoints, i)
        {
            const label pointI = smoothPoints[i];

            if( vertexLocation_[pointI] & LOCKED )
                continue;

            forAllRow(pointCells, pointI, pcI)
            {
                const cell& c = cells[pointCells(pointI, pcI)];

                forAll(c, fI)
                    chF[c[fI]] = true;
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- make sure that neighbouring processors get the same information
        const PtrList<processorBoundaryPatch>& pBnd = mesh_.procBoundaries();
        forAll(pBnd, patchI)
        {
            const label start = pBnd[patchI].patchStart();
            const label size = pBnd[patchI].patchSize();

            labelLongList sendData;
            for(label faceI=0;faceI<size;++faceI)
            {
                if( chF[start+faceI] )
                    sendData.append(faceI);
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                pBnd[patchI].neiProcNo(),
                sendData.byteSize()
            );

            toOtherProc << sendData;
        }

        forAll(pBnd, patchI)
        {
            labelList receivedData;

            IPstream fromOtherProc
            (
                Pstream::blocking,
                pBnd[patchI].neiProcNo()
            );

            fromOtherProc >> receivedData;

            const label start = pBnd[patchI].patchStart();
            forAll(receivedData, i)
                chF[start+receivedData[i]] = true;
        }
    }

    //- update geometry information
    const_cast<polyMeshGenAddressing&>
    (
        mesh_.addressingData()
    ).updateGeometry(chF);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Constructor of laplaceSmoother

meshOptimizer::laplaceSmoother::laplaceSmoother
(
    polyMeshGen& mesh,
    const List<direction>& vertexLocation
)
:
    mesh_(mesh),
    vertexLocation_(vertexLocation)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

meshOptimizer::laplaceSmoother::~laplaceSmoother()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Member Functions

void meshOptimizer::laplaceSmoother::optimizeLaplacian(const label nIterations)
{
    labelLongList smoothPoints;

    forAll(vertexLocation_, pointI)
    {
        if( vertexLocation_[pointI] & INSIDE )
            smoothPoints.append(pointI);
    }

    laplacian(smoothPoints, nIterations);
}

void meshOptimizer::laplaceSmoother::optimizeLaplacian
(
    const labelHashSet& badFaces,
    const label nIterations
)
{
    const faceListPMG& faces = mesh_.faces();

    std::set<label> activePoints;
    forAllConstIter(labelHashSet, badFaces, it)
    {
        const face& f = faces[it.key()];

        forAll(f, pI)
            activePoints.insert(f[pI]);
    }

    labelLongList smoothPoints;
    forAllConstIter(std::set<label>, activePoints, it)
        smoothPoints.append(*it);

    laplacian(smoothPoints, nIterations);
}

void meshOptimizer::laplaceSmoother::optimizeSurfaceLaplacian
(
    const labelHashSet& /*badFaces*/,
    const label /*nIterations*/
)
{
    FatalError << "Not implemented " << exit(FatalError);
}

void meshOptimizer::laplaceSmoother::optimizeLaplacianPC
(
    const label nIterations
)
{
    labelLongList smoothPoints;

    forAll(vertexLocation_, pointI)
    {
        if( vertexLocation_[pointI] & INSIDE )
            smoothPoints.append(pointI);
    }

    laplacianPC(smoothPoints, nIterations);
}

void meshOptimizer::laplaceSmoother::optimizeLaplacianPC
(
    const labelHashSet& badFaces,
    const label nIterations
)
{
    const faceListPMG& faces = mesh_.faces();

    std::set<label> activePoints;
    forAllConstIter(labelHashSet, badFaces, it)
    {
        const face& f = faces[it.key()];

        forAll(f, pI)
            activePoints.insert(f[pI]);
    }

    labelLongList smoothPoints;
    forAllConstIter(std::set<label>, activePoints, it)
        smoothPoints.append(*it);

    laplacianPC(smoothPoints, nIterations);
}

void meshOptimizer::laplaceSmoother::optimizeLaplacianWPC
(
    const label nIterations
)
{
    labelLongList smoothPoints;

    forAll(vertexLocation_, pointI)
    {
        if( vertexLocation_[pointI] & INSIDE )
            smoothPoints.append(pointI);
    }

    laplacianWPC(smoothPoints, nIterations);
}

void meshOptimizer::laplaceSmoother::optimizeLaplacianWPC
(
    const labelHashSet& badFaces,
    const label nIterations
)
{
    const faceListPMG& faces = mesh_.faces();

    std::set<label> activePoints;
    forAllConstIter(labelHashSet, badFaces, it)
    {
        const face& f = faces[it.key()];

        forAll(f, pI)
            activePoints.insert(f[pI]);
    }

    labelLongList smoothPoints;
    forAllConstIter(std::set<label>, activePoints, it)
        smoothPoints.append(*it);

    laplacianWPC(smoothPoints, nIterations);
}

void meshOptimizer::laplaceSmoother::optimizeLaplacianWHPC
(
    const label nIterations
)
{
    labelLongList smoothPoints;

    forAll(vertexLocation_, pointI)
    {
        if( vertexLocation_[pointI] & INSIDE )
            smoothPoints.append(pointI);
    }

    laplacianWHPC(smoothPoints, nIterations);
}

void meshOptimizer::laplaceSmoother::optimizeLaplacianWHPC
(
    const labelHashSet& badFaces,
    const label nIterations
)
{
    const faceListPMG& faces = mesh_.faces();

    std::set<label> activePoints;
    forAllConstIter(labelHashSet, badFaces, it)
    {
        const face& f = faces[it.key()];

        forAll(f, pI)
            activePoints.insert(f[pI]);
    }

    labelLongList smoothPoints;
    forAllConstIter(std::set<label>, activePoints, it)
        smoothPoints.append(*it);

    laplacianWHPC(smoothPoints, nIterations);
}

void meshOptimizer::faceFlatnessSmoother::improveFlatness
(
    const labelLongList& smoothPoints,
    const label nIterations
)
{
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();
    const VRWGraph& pointFaces = mesh_.addressingData().pointFaces();

    polyMeshGenModifier meshModifier(mesh_);

    for(label iterationI=0;iterationI<nIterations;++iterationI)
    {
        labelLongList procPoints;

        # ifdef USE_OMP
        # pragma omp parallel for schedule(guided, 20)
        # endif
        forAll(smoothPoints, i)
        {
            const label pointI = smoothPoints[i];

            if( vertexLocation_[pointI] & LOCKED )
                continue;

            if( pointFaces.sizeOfRow(pointI) == 0 )
                continue;

            if( vertexLocation_[pointI] & PARALLELBOUNDARY )
            {
                # ifdef USE_OMP
                # pragma omp critical
                # endif
                procPoints.append(pointI);

                continue;
            }

            const point& p = points[pointI];

            vector disp(vector::zero);
            scalar sumWeights(0.0);
            forAllRow(pointFaces, pointI, pfI)
            {
                const label faceI = pointFaces(pointI, pfI);

                const face& f = faces[faceI];
                const point c = help::faceCentre(points, f);
                const vector fn = help::faceAreaVector(points, f);
                const scalar fnSq = magSqr(fn) + VSMALL;

                vector d = ((c - p) & fn) * fn / fnSq;

                const scalar w = 1.0 / fnSq;

                disp += w * d;
                sumWeights += w;
            }

            disp /= sumWeights;
            meshModifier.movePoint(pointI, p + disp);
        }

        improveFlatnessParallel(procPoints);

        updateMeshGeometry(smoothPoints);
    }
}

void meshOptimizer::faceFlatnessSmoother::updateMeshGeometry
(
    const labelLongList& smoothPoints
)
{
    const cellListPMG& cells = mesh_.cells();
    const VRWGraph& pointCells = mesh_.addressingData().pointCells();

    boolList chF(mesh_.faces().size(), false);

    # ifdef USE_OMP
    # pragma omp parallel for if( smoothPoints.size() > 100 ) \
    schedule(guided, 20)
    # endif
    forAll(smoothPoints, i)
    {
        const label pointI = smoothPoints[i];

        if( vertexLocation_[pointI] & LOCKED )
            continue;

        forAllRow(pointCells, pointI, pcI)
        {
            const cell& c = cells[pointCells(pointI, pcI)];

            forAll(c, fI)
                chF[c[fI]] = true;
        }
    }

    if( Pstream::parRun() )
    {
        //- make sure that neighbouring processors get the same information
        const PtrList<processorBoundaryPatch>& pBnd = mesh_.procBoundaries();
        forAll(pBnd, patchI)
        {
            const label start = pBnd[patchI].patchStart();
            const label size = pBnd[patchI].patchSize();

            labelLongList sendData;
            for(label faceI=0;faceI<size;++faceI)
            {
                if( chF[start+faceI] )
                    sendData.append(faceI);
            }

            OPstream toOtherProc
            (
                Pstream::blocking,
                pBnd[patchI].neiProcNo(),
                sendData.byteSize()
            );

            toOtherProc << sendData;
        }

        forAll(pBnd, patchI)
        {
            labelList receivedData;

            IPstream fromOtherProc
            (
                Pstream::blocking,
                pBnd[patchI].neiProcNo()
            );

            fromOtherProc >> receivedData;

            const label start = pBnd[patchI].patchStart();
            forAll(receivedData, i)
                chF[start+receivedData[i]] = true;
        }
    }

    //- update geometry information
    const_cast<polyMeshGenAddressing&>
    (
        mesh_.addressingData()
    ).updateGeometry(chF);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Constructor of faceFlatnessSmoother

meshOptimizer::faceFlatnessSmoother::faceFlatnessSmoother
(
    polyMeshGen& mesh,
    const List<direction>& vertexLocation
)
:
    mesh_(mesh),
    vertexLocation_(vertexLocation)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

meshOptimizer::faceFlatnessSmoother::~faceFlatnessSmoother()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Member Functions

void meshOptimizer::faceFlatnessSmoother::optimizeFaceFlatness
(
    const label nIterations
)
{
    labelLongList smoothPoints;

    forAll(vertexLocation_, pointI)
    {
        if( vertexLocation_[pointI] & INSIDE )
            smoothPoints.append(pointI);
    }

    improveFlatness(smoothPoints, nIterations);
}

void meshOptimizer::faceFlatnessSmoother::optimizeFaceFlatness
(
    const labelHashSet& badFaces,
    const label nIterations
)
{
    const faceListPMG& faces = mesh_.faces();

    std::set<label> activePoints;
    forAllConstIter(labelHashSet, badFaces, it)
    {
        const face& f = faces[it.key()];

        forAll(f, pI)
        {
            if( vertexLocation_[f[pI]] & INSIDE )
                activePoints.insert(f[pI]);
        }
    }

    labelLongList smoothPoints;
    forAllConstIter(std::set<label>, activePoints, it)
        smoothPoints.append(*it);

    improveFlatness(smoothPoints, nIterations);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
