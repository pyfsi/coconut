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

#include "meshSurfaceEngineModifier.H"
#include "polyMeshGenModifier.H"
#include "demandDrivenData.H"

#include "labelledPoint.H"
#include "refLabelledPoint.H"
#include "helperFunctions.H"

#include "meshOctree.H"
#include "triSurf.H"

//#define DEBUGEngineModifier

# ifdef DEBUGEngineModifier
#include "OFstream.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceEngineModifier::updatePointNormals
(
    const boolList& updateBndPoint
)
{
    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const vectorLongList& faceNormals = surfaceEngine_.faceNormals();
    const labelLongList& facePatches = surfaceEngine_.boundaryFacePatches();

    vectorLongList& pn = *surfaceEngine_.pointNormalsPtr_;

    const VRWGraph* bpAtProcsPtr_ = NULL;

    if( Pstream::parRun() )
        bpAtProcsPtr_ = &surfaceEngine_.bpAtProcs();

    //- calculate the sum of local face normals
    std::map<label, std::map<label, vector> > pointPatchNormals;

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(updateBndPoint, bpI)
        {
            if( !updateBndPoint[bpI] )
                continue;

            //- initialise a map of normals for each patch
            std::map<label, vector> patchVector;

            forAllRow(pFaces, bpI, pfI)
            {
                const label bfI = pFaces(bpI, pfI);

                patchVector[facePatches[bfI]] = vector::zero;
            }

            //- sum up normals with respect to each patch
            forAllRow(pFaces, bpI, pfI)
            {
                const label bfI = pFaces(bpI, pfI);

                patchVector[facePatches[bfI]] += faceNormals[bfI];
            }

            if( bpAtProcsPtr_ && bpAtProcsPtr_->sizeOfRow(bpI) )
            {
                //- store data needed to exchange information among processors
                # ifdef USE_OMP
                # pragma omp critical(insertPatchVectors)
                # endif
                {
                    pointPatchNormals[bpI] = patchVector;
                }
            }
            else
            {
                //- point is available at this proc only
                vector normal(vector::zero);
                for
                (
                    std::map<label, vector>::iterator it=patchVector.begin();
                    it!=patchVector.end();
                    ++it
                )
                {
                    it->second /= (mag(it->second) + VSMALL);

                    normal += it->second;
                }

                pn[bpI] = normal;
            }
        }

        if( Pstream::parRun() )
        {
            # ifdef USE_OMP
            # pragma omp single
            # endif
            {
                //- update point normals at inter-processor boundaries
                const Map<label>& globalToLocal =
                    surfaceEngine_.globalToLocalBndPointAddressing();
                const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();

                //- start updating point normals
                std::map<label, DynList<refLabelledPoint, 64> > exchangeData;
                forAll(surfaceEngine_.bpNeiProcs(), i)
                    exchangeData[surfaceEngine_.bpNeiProcs()[i]].clear();

                //- prepare data for sending
                forAllConstIter(Map<label>, globalToLocal, it)
                {
                    const label bpI = it();

                    if( pointPatchNormals.find(bpI) != pointPatchNormals.end() )
                    {
                        const std::map<label, vector>& patchVector =
                            pointPatchNormals[bpI];

                        for
                        (
                            std::map<label, vector>::const_iterator pIt=
                                patchVector.begin();
                            pIt!=patchVector.end();
                            ++pIt
                        )
                        {
                            forAllRow(bpAtProcs, bpI, procI)
                            {
                                const label neiProc = bpAtProcs(bpI, procI);
                                if( neiProc == Pstream::myProcNo() )
                                    continue;

                                exchangeData[neiProc].append
                                (
                                    refLabelledPoint
                                    (
                                        it.key(),
                                        labelledPoint(pIt->first, pIt->second)
                                    )
                                );
                            }
                        }
                    }
                }

                //- exchange data with other procs
                LongList<refLabelledPoint> receivedData;
                help::exchangeMap(exchangeData, receivedData);

                forAll(receivedData, i)
                {
                    const refLabelledPoint& rlp = receivedData[i];
                    const label bpI = globalToLocal[rlp.objectLabel()];

                    std::map<label, vector>& patchVector =
                        pointPatchNormals[bpI];

                    //- update patch vectors
                    std::map<label, vector>::iterator pIt =
                        patchVector.find(rlp.lPoint().pointLabel());
                    if( pIt != patchVector.end() )
                    {
                        pIt->second += rlp.lPoint().coordinates();
                    }
                    else
                    {
                        patchVector[rlp.lPoint().pointLabel()] =
                            rlp.lPoint().coordinates();
                    }
                }

                //- update point normals at inter-processor boundaries
                std::map<label, std::map<label, point> >::iterator mIt;
                for
                (
                    mIt=pointPatchNormals.begin();
                    mIt!=pointPatchNormals.end();
                    ++mIt
                )
                {
                    # ifdef USE_OMP
                    # pragma omp task default(shared) firstprivate(mIt)
                    # endif
                    {
                        vector n = vector::zero;

                        std::map<label, vector>& patchVector =
                            mIt->second;

                        for
                        (
                            std::map<label, vector>::iterator pIt=
                            patchVector.begin();
                            pIt!=patchVector.end();
                            ++pIt
                        )
                        {
                            pIt->second /= (mag(pIt->second) + VSMALL);
                            n += pIt->second;
                        }

                        pn[mIt->first] = n;
                    }
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(updateBndPoint, bpI)
        {
            if( !updateBndPoint[bpI] )
                continue;

            pn[bpI] /= (mag(pn[bpI]) + VSMALL);
        }
    }

    # ifdef DEBUGEngineModifier
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

    OFstream file
    (
        "meshSurfaceNormalVectors" +
        help::labelToText(Pstream::myProcNo()) +
        ".vtk"
    );

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "POINTS " << 2*pn.size() << " float\n";
    forAll(pn, bpI)
    {
        const point& p = points[bPoints[bpI]];

        scalar maxDist(0.0);
        forAllRow(pFaces, bpI, pfI)
        {
            const point fc =
                help::faceCentre(points, bFaces[pFaces(bpI, pfI)]);
            maxDist = max(maxDist, mag(p - fc));
        }

        file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;

        const point op = p + maxDist * pn[bpI];

        file << op.x() << ' ' << op.y() << ' ' << op.z() << nl;
    }

    //- write lines
    file << "\nLINES " << pn.size()
         << " " << 3*pn.size() << nl;
    forAll(pn, eI)
    {
        file << 2 << " " << 2*eI << " " << (2*eI+1) << nl;
    }

    file << "\n";
    file.flush();

    if( Pstream::parRun() )
    {
        //- check point normals at inter-processor boundaries
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();

        std::map<label, LongList<labelledPoint> > exchangeData;
        forAll(surfaceEngine_.bpNeiProcs(), i)
            exchangeData[surfaceEngine_.bpNeiProcs()[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            forAllRow(bpAtProcs, bpI, i)
            {
                const label neiProc = bpAtProcs(bpI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append(labelledPoint(it.key(), pn[bpI]));
            }
        }

        LongList<labelledPoint> receiveData;
        help::exchangeMap(exchangeData, receiveData);

        Pout << "Received size " << receiveData.size() << endl;

        forAll(receiveData, i)
        {
            const label bpI = globalToLocal[receiveData[i].pointLabel()];

            if( mag(pn[bpI] - receiveData[i].coordinates()) > SMALL )
            {
                FatalError << "Normal at point " << bpI << " local "
                    << pn[bpI] << " other " << receiveData[i].coordinates()
                    << abort(FatalError);
            }
        }
    }
    # endif
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceEngineModifier::meshSurfaceEngineModifier
(
    meshSurfaceEngine& surfaceEngine
)
:
    surfaceEngine_(surfaceEngine),
    exchangeLabelsPtr_(NULL),
    exchangeLabelledPointsPtr_(NULL)
{
    if( Pstream::parRun() )
    {
        labelMap();
        labelledPointMap();
    }
}

meshSurfaceEngineModifier::meshSurfaceEngineModifier
(
    const meshSurfaceEngine& surfaceEngine
)
:
    surfaceEngine_(const_cast<meshSurfaceEngine&>(surfaceEngine)),
    exchangeLabelsPtr_(NULL),
    exchangeLabelledPointsPtr_(NULL)
{
    if( Pstream::parRun() )
    {
        labelMap();
        labelledPointMap();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceEngineModifier::~meshSurfaceEngineModifier()
{
    deleteDemandDrivenData(exchangeLabelsPtr_);
    deleteDemandDrivenData(exchangeLabelledPointsPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceEngineModifier::moveBoundaryVertexNoUpdate
(
    const label bpI,
    const point& newP
)
{
    const_cast<point&>
    (
        surfaceEngine_.mesh_.points()[surfaceEngine_.boundaryPoints()[bpI]]
    ) = newP;
}

void meshSurfaceEngineModifier::moveBoundaryVertex
(
    const label bpI,
    const point& newP
)
{
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.mesh_.points();
    const_cast<point&>(points[bPoints[bpI]]) = newP;

    if( surfaceEngine_.faceCentresPtr_ )
    {
        vectorLongList& faceCentres = *surfaceEngine_.faceCentresPtr_;
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces(bpI, pfI);

            faceCentres[bfI] = help::faceCentre(points, bFaces[bfI]);
        }
    }

    if( surfaceEngine_.faceNormalsPtr_ )
    {
        vectorLongList& faceNormals = *surfaceEngine_.faceNormalsPtr_;
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

        forAllRow(pFaces, bpI, pfI)
        {
            const label bfI = pFaces(bpI, pfI);

            faceNormals[bfI] = bFaces[bfI].normal(points);
        }
    }

    if( surfaceEngine_.pointNormalsPtr_ )
    {
        const vectorLongList& faceNormals = *surfaceEngine_.faceNormalsPtr_;
        const VRWGraph& pFaces = surfaceEngine_.pointFaces();
        const VRWGraph& pPoints = surfaceEngine_.pointPoints();

        vectorLongList& pn = *surfaceEngine_.pointNormalsPtr_;
        vector n(vector::zero);
        forAllRow(pFaces, bpI, pfI)
            n += faceNormals[pFaces(bpI, pfI)];

        const scalar l = mag(n);
        if( l > VSMALL )
        {
            n /= l;
        }
        else
        {
            n = vector::zero;
        }

        pn[bpI] = n;

        //- change normal of vertices connected to bpI
        forAllRow(pPoints, bpI, ppI)
        {
            const label bpJ = pPoints(bpI, ppI);
            n = vector::zero;
            forAllRow(pFaces, bpJ, pfI)
                n += faceNormals[pFaces(bpJ, pfI)];

            const scalar d = mag(n);
            if( d > VSMALL )
            {
                n /= d;
            }
            else
            {
                n = vector::zero;
            }

            pn[bpJ] = n;
        }
    }
}

void meshSurfaceEngineModifier::syncVerticesAtParallelBoundaries()
{
    if( !Pstream::parRun() )
        return;

    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();
    labelLongList syncNodes;
    forAllConstIter(Map<label>, globalToLocal, it)
        syncNodes.append(it());

    syncVerticesAtParallelBoundaries(syncNodes);
}

void meshSurfaceEngineModifier::syncVerticesAtParallelBoundaries
(
    const labelLongList& syncNodes
)
{
    if( !Pstream::parRun() )
        return;

    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
    const labelLongList& globalLabel =
        surfaceEngine_.globalBoundaryPointLabel();
    const Map<label>& globalToLocal =
        surfaceEngine_.globalToLocalBndPointAddressing();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const pointFieldPMG& points = surfaceEngine_.mesh().points();

    labelledPointMapType& exchangeData = labelledPointMap();

    //- construct the map
    forAll(syncNodes, snI)
    {
        const label bpI = syncNodes[snI];

        if( bpAtProcs.sizeOfRow(bpI) == 0 )
            continue;

        point p = points[bPoints[bpI]] / bpAtProcs.sizeOfRow(bpI);
        moveBoundaryVertexNoUpdate(bpI, p);

        forAllRow(bpAtProcs, bpI, i)
        {
            const label neiProc = bpAtProcs(bpI, i);
            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeData[neiProc].append(labelledPoint(globalLabel[bpI], p));
        }
    }

    //- exchange the data with other processors
    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    //- adjust the coordinates
    forAll(receivedData, i)
    {
        const labelledPoint& lp = receivedData[i];
        const label bpI = globalToLocal[lp.pointLabel()];
        const point newP = points[bPoints[bpI]] + lp.coordinates();
        moveBoundaryVertexNoUpdate(bpI, newP);
    }
}

void meshSurfaceEngineModifier::updateGeometry
(
    const labelLongList& updateBndNodes
)
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
    const VRWGraph& pFaces = surfaceEngine_.pointFaces();
    const labelLongList& bp = surfaceEngine_.bp();

    boolList updateFaces(bFaces.size());
    boolList updateBndPoint(pFaces.size());
    bool updateBndPointNormals(false);

    # ifdef USE_OMP
    # pragma omp parallel if( updateBndNodes.size() > 1000 )
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(updateFaces, bfI)
            updateFaces[bfI] = false;

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(updateBndNodes, i)
        {
            const label bpI = updateBndNodes[i];

            forAllRow(pFaces, bpI, j)
            {
                const label bfI = pFaces(bpI, j);

                updateFaces[bfI] = true;
            }
        }

        if( surfaceEngine_.faceCentresPtr_ )
        {
            vectorLongList& faceCentres = *surfaceEngine_.faceCentresPtr_;

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100) nowait
            # endif
            forAll(updateFaces, bfI)
            {
                if( updateFaces[bfI] )
                    faceCentres[bfI] = help::faceCentre(points, bFaces[bfI]);
            }
        }

        if( surfaceEngine_.faceNormalsPtr_ )
        {
            vectorLongList& faceNormals = *surfaceEngine_.faceNormalsPtr_;

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(updateFaces, bfI)
            {
                if( updateFaces[bfI] )
                    faceNormals[bfI] = bFaces[bfI].normal(points);
            }
        }

        if( surfaceEngine_.pointNormalsPtr_ )
        {
            updateBndPointNormals = true;

            # ifdef USE_OMP
            # pragma omp for schedule(static, 1)
            # endif
            forAll(updateBndPoint, bpI)
                updateBndPoint[bpI] = false;

            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(updateBndNodes, i)
            {
                const label bpI = updateBndNodes[i];

                forAllRow(pFaces, bpI, pfI)
                {
                    const face& bf = bFaces[pFaces(bpI, pfI)];

                    forAll(bf, pI)
                    {
                        const label bpJ = bp[bf[pI]];

                        updateBndPoint[bpJ] = true;
                    }
                }
            }

            if( Pstream::parRun() )
            {
                //- ensure that the same points are selected over all processors
                # ifdef USE_OMP
                # pragma omp single
                # endif
                {
                    //- update point normals at inter-processor boundaries
                    const Map<label>& globalToLocal =
                        surfaceEngine_.globalToLocalBndPointAddressing();
                    const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();

                    //- make sure that the points ar updated on all processors
                    labelMapType& exchangeNodeLabels = labelMap();

                    forAllConstIter(Map<label>, globalToLocal, it)
                    {
                        const label bpI = it();

                        if( updateBndPoint[bpI] )
                        {
                            forAllRow(bpAtProcs, bpI, i)
                            {
                                const label neiProc = bpAtProcs(bpI, i);

                                if( neiProc == Pstream::myProcNo() )
                                    continue;

                                exchangeNodeLabels[neiProc].append(it.key());
                            }
                        }
                    }

                    labelLongList receivedNodes;
                    help::exchangeMap(exchangeNodeLabels, receivedNodes);

                    forAll(receivedNodes, i)
                    {
                        const label bpI = globalToLocal[receivedNodes[i]];
                        updateBndPoint[bpI] = true;
                    }
                }
            }
        }
    }

    if( updateBndPointNormals )
        updatePointNormals(updateBndPoint);
}

void meshSurfaceEngineModifier::updateGeometry()
{
    labelLongList updateBndNodes(surfaceEngine_.boundaryPoints().size());

    # ifdef USE_OMP
    # pragma omp parallel for if( updateBndNodes.size() > 10000 )
    # endif
    forAll(updateBndNodes, bpI)
        updateBndNodes[bpI] = bpI;

    updateGeometry(updateBndNodes);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
