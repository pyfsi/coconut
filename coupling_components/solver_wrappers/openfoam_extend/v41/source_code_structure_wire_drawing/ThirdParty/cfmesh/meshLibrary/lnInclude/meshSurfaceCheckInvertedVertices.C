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

#include "meshSurfaceCheckInvertedVertices.H"
#include "meshSurfacePartitioner.H"
#include "boolList.H"
#include "demandDrivenData.H"
#include "refLabelledPoint.H"
#include "helperFunctions.H"

#include <map>

# ifdef USE_OMP
#include <omp.h>
# endif

//# define DEBUGInverted

# ifdef DEBUGInverted
# include "OFstream.H"
# endif

//#define timingIntverted

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# ifdef DEBUGInverted
void writeFaceToVTK
(
    const fileName& fName,
    const meshSurfaceEngine& mse,
    const label bfI
)
{
    const pointFieldPMG& points = mse.points();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const face& bf = bFaces[bfI];

    OFstream file(fName);

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "POINTS " << bf.size() << " float\n";
    forAll(bf, pI)
    {
        const point& p = points[bf[pI]];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
    }

    //- write faces
    file << "\nPOLYGONS " << 1
         << " " << (bf.size()+1) << nl;
    file << bf.size();
    forAll(bf, pI)
    {
        file << " " << pI;
    }

    file << nl;
    file << flush;
}
# endif

void meshSurfaceCheckInvertedVertices::checkVertices()
{
    # ifdef timingIntverted
    const scalar startCheckingVertices = omp_get_wtime();
    # endif

    const labelLongList& facePatch = surfacePartitioner_.boundaryFacePatches();
    const meshSurfaceEngine& mse = surfacePartitioner_.surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const labelLongList& bp = mse.bp();
    const VRWGraph& pointFaces = mse.pointFaces();
    const VRWGraph& pointInFaces = mse.pointInFaces();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const vectorLongList& pNormals = mse.pointNormals();
    const vectorLongList& fCentres = mse.faceCentres();
    const vectorLongList& fNormals = mse.faceNormals();

    const labelHashSet& corners = surfacePartitioner_.corners();
    const labelHashSet& edgePoints = surfacePartitioner_.edgePoints();

    # ifdef timingIntverted
    const scalar startNormalCalculation = omp_get_wtime();
    Info << "Starting calculation of normals "
         << (startNormalCalculation-startCheckingVertices) << endl;
    # endif

    typedef std::map<label, vector> ltvMap;
    typedef std::map<label, ltvMap> lltvMap;
    lltvMap pointPatchNormal;

    forAllConstIter(labelHashSet, corners, it)
    {
        const label bpI = it.key();

        if( activePointsPtr_ && !activePointsPtr_->operator[](bpI))
            continue;

        ltvMap& patchNormal = pointPatchNormal[bpI];
        patchNormal.clear();

        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            const label patchI = facePatch[bfI];

            if( patchNormal.find(patchI) == patchNormal.end() )
            {
                patchNormal[patchI] = fNormals[bfI];
            }
            else
            {
                patchNormal[patchI] += fNormals[bfI];
            }
        }
    }

    forAllConstIter(labelHashSet, edgePoints, it)
    {
        const label bpI = it.key();

        if( activePointsPtr_ && !activePointsPtr_->operator[](bpI))
            continue;

        ltvMap& patchNormal = pointPatchNormal[bpI];

        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            const label patchI = facePatch[bfI];

            if( patchNormal.find(patchI) == patchNormal.end() )
            {
                patchNormal[patchI] = fNormals[bfI];
            }
            else
            {
                patchNormal[patchI] += fNormals[bfI];
            }
        }
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal = mse.globalToLocalBndPointAddressing();
        const DynList<label>& neiProcs = mse.bpNeiProcs();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();

        std::map<label, LongList<refLabelledPoint> > exchangeData;
        forAll(neiProcs, i)
            exchangeData[neiProcs[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( pointPatchNormal.find(bpI) != pointPatchNormal.end() )
            {
                const ltvMap& patchNormal = pointPatchNormal[bpI];

                forAllRow(bpAtProcs, bpI, i)
                {
                    const label neiProc = bpAtProcs(bpI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    forAllConstIter(ltvMap, patchNormal, pIt)
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

        LongList<refLabelledPoint> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const refLabelledPoint& rlp = receivedData[i];
            const label bpI = globalToLocal[rlp.objectLabel()];

            lltvMap::iterator pIt = pointPatchNormal.find(bpI);
            if( pIt == pointPatchNormal.end() )
                FatalError << "point is not active at this processor"
                    << abort(FatalError);

            ltvMap& patchNormal = pIt->second;

            const labelledPoint& lp = rlp.lPoint();
            if( patchNormal.find(lp.pointLabel()) != patchNormal.end() )
            {
                patchNormal[lp.pointLabel()] += lp.coordinates();
            }
            else
            {
                patchNormal[lp.pointLabel()] = lp.coordinates();
            }
        }
    }

    forAllIter(lltvMap, pointPatchNormal, it)
    {
        ltvMap& patchNormal = it->second;

        forAllIter(ltvMap, patchNormal, pIt)
        {
            const scalar magv = mag(pIt->second) + VSMALL;

            pIt->second /= magv;
        }
    }

    # ifdef timingIntverted
    const scalar startPointChecks = omp_get_wtime();
    Info << "Time for calculation of normal vectors "
         << (startPointChecks-startNormalCalculation) << endl;
    # endif

    //- find out the problematic points at the surface
    invertedVertices_.clear();
    boolList activeFace(bFaces.size(), false);

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(activeFace, bfI)
            activeFace[bfI] = false;

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 20)
        # endif
        forAll(pointFaces, bpI)
        {
            if( activePointsPtr_ && !activePointsPtr_->operator[](bpI) )
                continue;

            forAllRow(pointFaces, bpI, pfI)
            {
                const label pI = pointInFaces(bpI, pfI);
                const label bfI = pointFaces(bpI, pfI);

                activeFace[bfI] = true;

                vector pNormal = pNormals[bpI];

                if( pointPatchNormal.find(bpI) != pointPatchNormal.end() )
                    pNormal = pointPatchNormal[bpI][facePatch[bfI]];

                const face& bf = bFaces[bfI];

                //- chech the first triangle (with the next node)
                const triangle<point, point> triNext
                (
                    points[bf[pI]],
                    points[bf.nextLabel(pI)],
                    fCentres[bfI]
                );

                vector nNext = triNext.normal();
                scalar mNext = mag(nNext);

                //- face has zero area
                if( mNext < VSMALL )
                {
                    # ifdef DEBUGInverted
                    Info << "1. inverted " << bfI << endl;
                    writeFaceToVTK
                    (
                        "1.inverted_"+help::labelToText(bfI)+".vtk",
                        mse,
                        bfI
                    );
                    # endif

                    # ifdef USE_OMP
                    # pragma omp critical(invertedVertices)
                    # endif
                    invertedVertices_.insert(bf[pI]);

                    break;
                }
                else
                {
                    nNext /= mNext;
                }

                //- collocated points
                if( magSqr(triNext.a() - triNext.b()) < VSMALL )
                {
                    # ifdef DEBUGInverted
                    Info << "2. inverted " << bfI << endl;
                    writeFaceToVTK
                    (
                        "2.inverted_"+help::labelToText(bfI)+".vtk",
                        mse,
                        bfI
                    );
                    # endif

                    # ifdef USE_OMP
                    # pragma omp critical(invertedVertices)
                    # endif
                    invertedVertices_.insert(bf[pI]);

                    break;
                }
                if( magSqr(triNext.c() - triNext.a()) < VSMALL )
                {
                    # ifdef DEBUGInverted
                    Info << "3. inverted " << bfI << endl;
                    writeFaceToVTK
                    (
                        "3.inverted_"+help::labelToText(bfI)+".vtk",
                        mse,
                        bfI
                    );
                    # endif

                    # ifdef USE_OMP
                    # pragma omp critical(invertedVertices)
                    # endif
                    invertedVertices_.insert(bf[pI]);

                    break;
                }

                //- normal vector is not visible
                if( (nNext & pNormal) < 0.0 )
                {
                    # ifdef DEBUGInverted
                    Info << "4. inverted " << bfI << endl;
                    writeFaceToVTK
                    (
                        "4.inverted_"+help::labelToText(bfI)+".vtk",
                        mse,
                        bfI
                    );
                    # endif

                    # ifdef USE_OMP
                    # pragma omp critical(invertedVertices)
                    # endif
                    invertedVertices_.insert(bf[pI]);

                    break;
                }

                //- check the second triangle (with previous node)
                const triangle<point, point> triPrev
                (
                    points[bf[pI]],
                    fCentres[bfI],
                    points[bf.prevLabel(pI)]
                );

                vector nPrev = triPrev.normal();
                scalar mPrev = mag(nPrev);

                //- face has zero area
                if( mPrev < VSMALL )
                {
                    # ifdef DEBUGInverted
                    Info << "5. inverted " << bfI << endl;
                    writeFaceToVTK
                    (
                        "5.inverted_"+help::labelToText(bfI)+".vtk",
                        mse,
                        bfI
                    );
                    # endif

                    # ifdef USE_OMP
                    # pragma omp critical(invertedVertices)
                    # endif
                    invertedVertices_.insert(bf[pI]);

                    break;
                }
                else
                {
                    nPrev /= mPrev;
                }

                //- collocated points
                if( magSqr(triPrev.a() - triPrev.b()) < VSMALL )
                {
                    # ifdef DEBUGInverted
                    Info << "6. inverted " << bfI << endl;
                    writeFaceToVTK
                    (
                        "6.inverted_"+help::labelToText(bfI)+".vtk",
                        mse,
                        bfI
                    );
                    # endif

                    # ifdef USE_OMP
                    # pragma omp critical(invertedVertices)
                    # endif
                    invertedVertices_.insert(bf[pI]);

                    break;
                }
                if( magSqr(triPrev.c() - triPrev.a()) < VSMALL )
                {
                    # ifdef DEBUGInverted
                    Info << "7. inverted " << bfI << endl;
                    writeFaceToVTK
                    (
                        "7.inverted_"+help::labelToText(bfI)+".vtk",
                        mse,
                        bfI
                    );
                    # endif

                    # ifdef USE_OMP
                    # pragma omp critical(invertedVertices)
                    # endif
                    invertedVertices_.insert(bf[pI]);

                    break;
                }

                //- normal vector is not visible
                if( (nPrev & pNormal) < 0.0 )
                {
                    # ifdef DEBUGInverted
                    Info << "8. inverted " << bfI << endl;
                    writeFaceToVTK
                    (
                        "8.inverted_"+help::labelToText(bfI)+".vtk",
                        mse,
                        bfI
                    );
                    # endif

                    # ifdef USE_OMP
                    # pragma omp critical(invertedVertices)
                    # endif
                    invertedVertices_.insert(bf[pI]);

                    break;
                }

                //- check whether the normals of both triangles
                //- point in the same direction
                if( (nNext & nPrev) < 0.0 )
                {
                    # ifdef DEBUGInverted
                    Info << "9. inverted "
                         << " nNext " << nNext
                         << " nPrev " << nPrev << (nNext & nPrev) << endl;
                    writeFaceToVTK
                    (
                        "9.inverted_"+help::labelToText(bfI)+".vtk",
                        mse,
                        bfI
                    );
                    # endif

                    # ifdef USE_OMP
                    # pragma omp critical(invertedVertices)
                    # endif
                    invertedVertices_.insert(bf[pI]);

                    break;
                }
            }
        }

        //- check if there exist concave faces
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 50)
        # endif
        forAll(bFaces, bfI)
        {
            if( !activeFace[bfI] )
                continue;

            const face& bf = bFaces[bfI];

            //- check face flatness
            const scalar magN = mag(fNormals[bfI]);
            scalar sumArea(0.0);
            forAll(bf, pI)
            {
                const triangle<point, point> tri
                (
                    points[bf[pI]],
                    points[bf.nextLabel(pI)],
                    fCentres[bfI]
                );

                sumArea += tri.mag();
            }

            const scalar flatness = magN / (sumArea+VSMALL);

            if( flatness < 0.8 )
            {
                //- set all vertices as non-flat
                forAll(bf, pI)
                {
                    if( activePointsPtr_ && !(*activePointsPtr_)[bp[bf[pI]]] )
                        continue;

                    # ifdef USE_OMP
                    # pragma omp critical(invertedVertices)
                    # endif
                    invertedVertices_.insert(bf[pI]);
                }
            }

            //- check if the face is convex
            DynList<bool> OkPoints;
            if( !help::isFaceConvexAndOk(bf, points, OkPoints) )
            {
                forAll(OkPoints, pI)
                {
                    if( activePointsPtr_ && !(*activePointsPtr_)[bp[bf[pI]]] )
                        continue;

                    if( !OkPoints[pI] )
                    {
                        # ifdef DEBUGInverted
                        Info << "Concave face " << bfI << endl;
                        writeFaceToVTK
                        (
                            "concaveFace_"+help::labelToText(bfI)+".vtk",
                            mse,
                            bfI
                        );
                        # endif

                        # ifdef USE_OMP
                        # pragma omp critical(invertedVertices)
                        # endif
                        {
                            invertedVertices_.insert(bf[pI]);
                        }
                    }
                }
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- exchange global labels of inverted points
        const labelLongList& bPoints = mse.boundaryPoints();
        const Map<label>& globalToLocal =
            mse.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = mse.bpAtProcs();
        const DynList<label>& neiProcs = mse.bpNeiProcs();

        std::map<label, labelLongList> shareData;
        forAll(neiProcs, i)
            shareData.insert(std::make_pair(neiProcs[i], labelLongList()));

        forAllConstIter(Map<label>, globalToLocal, iter)
        {
            const label bpI = iter();

            if( !invertedVertices_.found(bPoints[bpI]) )
                continue;

            forAllRow(bpAtProcs, bpI, procI)
            {
                const label neiProc = bpAtProcs(bpI, procI);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                shareData[neiProc].append(iter.key());
            }
        }

        //- exchange data with other processors
        labelLongList receivedData;
        help::exchangeMap(shareData, receivedData);

        forAll(receivedData, i)
        {
            const label bpI = globalToLocal[receivedData[i]];
            invertedVertices_.insert(bPoints[bpI]);
        }
    }

    # ifdef timingIntverted
    const scalar finishCheckingInverted = omp_get_wtime();
    Info << "Time for checking inverted points "
         << (finishCheckingInverted-startPointChecks) << endl;
    # endif
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceCheckInvertedVertices::meshSurfaceCheckInvertedVertices
(
    const meshSurfacePartitioner& mpart
)
:
    surfacePartitioner_(mpart),
    activePointsPtr_(NULL),
    invertedVertices_(mpart.surfaceEngine().boundaryPoints().size())
{
    checkVertices();
}

meshSurfaceCheckInvertedVertices::meshSurfaceCheckInvertedVertices
(
    const meshSurfacePartitioner& mpart,
    const boolList& activePoints
)
:
    surfacePartitioner_(mpart),
    activePointsPtr_(&activePoints),
    invertedVertices_(mpart.surfaceEngine().boundaryPoints().size())
{
    checkVertices();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceCheckInvertedVertices::~meshSurfaceCheckInvertedVertices()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
