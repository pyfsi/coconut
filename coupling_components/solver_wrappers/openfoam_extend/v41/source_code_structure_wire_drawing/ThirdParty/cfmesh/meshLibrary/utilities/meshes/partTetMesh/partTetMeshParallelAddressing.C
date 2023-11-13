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
#include "polyMeshGenModifier.H"
#include "partTetMesh.H"
#include "polyMeshGenAddressing.H"
#include "helperFunctions.H"
#include "parPartTet.H"

#include <map>

//#define DEBUGSmooth

# ifdef DEBUGSmooth
#include "partTetMeshSimplex.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void partTetMesh::createParallelAddressing()
{
    const pointFieldPMG& points = origMesh_.points();
    const PtrList<processorBoundaryPatch>& procBoundaries =
        origMesh_.procBoundaries();

    const polyMeshGenAddressing& addr = origMesh_.addressingData();
    const Map<label>& globalToLocal = addr.globalToLocalPointAddressing();
    const VRWGraph& pointAtProcs = addr.pointAtProcs();

    //- ensure smoothVertex_ has the sme values across all processors
    std::map<label, labelLongList> exchangeData;
    forAll(addr.pointNeiProcs(), i)
        exchangeData[addr.pointNeiProcs()[i]].clear();

    forAllConstIter(Map<label>, globalToLocal, it)
    {
        const label pointI = it();

        forAllRow(pointAtProcs, pointI, i)
        {
            const label neiProc = pointAtProcs(pointI, i);

            if( neiProc == Pstream::myProcNo() )
                continue;

            labelLongList& dts = exchangeData[neiProc];

            //- send:
            //- 1. global point label
            //- 2. type flags
            dts.append(it.key());
            dts.append(smoothVertex_[pointI]);
        }
    }

    labelLongList receiveData;
    help::exchangeMap(exchangeData, receiveData);

    for(label i=0;i<receiveData.size();)
    {
        const label pointI= globalToLocal[receiveData[i++]];
        const direction type = receiveData[i++];

        smoothVertex_[pointI] |= type;
    }

    //- create tets local to this processor
    forAllConstIter(Map<label>, globalToLocal, it)
    {
        const label pointI = it();

        //- exclude points at internal boundaries
        if( smoothVertex_[pointI] & (LOCKED|INTERNALBOUNDARY|PARALLELBOUNDARY) )
            continue;

        //- points must be smooth
        if( smoothVertex_[pointI] & (SMOOTH|BOUNDARY) )
        {
            std::map<label, label> pToPts, fcToPts, ccToPts, ptsToParPoints;

            DynList<point, 128> pts;
            DynList<partTet, 256>& tets = tetsAtPoint_[pointI];
            DynList<labelledTri, 32>& bndTriangles =
                bndTrianglesAtPoint_[pointI];
            tets.clear();
            bndTriangles.clear();

            //- create tets located at this processor
            forAllRow(pointFaces_, pointI, pfI)
            {
                createTetsAtFace
                (
                    pointFaces_(pointI, pfI),
                    pointI,
                    pToPts,
                    fcToPts,
                    ccToPts,
                    pts,
                    tets,
                    bndTriangles
                );
            }

            std::map<label, label>::const_iterator mIt;

            //- copy points
            for(mIt=pToPts.begin();mIt!=pToPts.end();++mIt)
            {
                const label pointI = mIt->first;

                if
                (
                    pointToParPoint_.find(pointI) ==
                    pointToParPoint_.end()
                )
                {
                    pointToParPoint_[pointI] = parPoints_.size();

                    parPoints_.append(pts[mIt->second]);
                }

                ptsToParPoints[mIt->second] = pointToParPoint_[pointI];
            }

            //- copy face centres
            for(mIt=fcToPts.begin();mIt!=fcToPts.end();++mIt)
            {
                const label faceI = mIt->first;

                if
                (
                    faceCentreToParPoint_.find(faceI) ==
                    faceCentreToParPoint_.end()
                )
                {
                    faceCentreToParPoint_[faceI] = parPoints_.size();

                    parPoints_.append(pts[mIt->second]);
                }

                ptsToParPoints[mIt->second] = faceCentreToParPoint_[faceI];
            }

            //- copy cell centres
            for(mIt=ccToPts.begin();mIt!=ccToPts.end();++mIt)
            {
                const label cellI = mIt->first;

                if
                (
                    cellCentreToParPoint_.find(cellI) ==
                    cellCentreToParPoint_.end()
                )
                {
                    cellCentreToParPoint_[cellI] = parPoints_.size();

                    parPoints_.append(pts[mIt->second]);
                }

                ptsToParPoints[mIt->second] = cellCentreToParPoint_[cellI];
            }

            //- update tets
            forAll(tets, tetI)
            {
                partTet& tet = tets[tetI];

                partTet tCopy
                (
                    ptsToParPoints[tet.a()],
                    ptsToParPoints[tet.b()],
                    ptsToParPoints[tet.c()],
                    ptsToParPoints[tet.d()]
                );

                tet = tCopy;
            }

            //- update boundary triangles
            forAll(bndTriangles, triI)
            {
                labelledTri& tri = bndTriangles[triI];

                labelledTri triCopy
                (
                    ptsToParPoints[tri[0]],
                    ptsToParPoints[tri[1]],
                    ptsToParPoints[tri[2]],
                    tri.region()
                );

                tri = triCopy;
            }
        }
    }

    //- select points at all processors
    exchangeData.clear();
    forAll(addr.pointNeiProcs(), i)
        exchangeData[addr.pointNeiProcs()[i]].clear();

    forAllConstIter(Map<label>, globalToLocal, it)
    {
        if( pointToParPoint_.find(it()) != pointToParPoint_.end() )
        {
            const label pointI = it();

            forAllRow(pointAtProcs, pointI, i)
            {
                const label neiProc = pointAtProcs(pointI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                exchangeData[neiProc].append(it.key());
            }
        }
    }

    receiveData.clear();
    help::exchangeMap(exchangeData, receiveData);

    forAll(receiveData, i)
    {
        const label pointI = globalToLocal[receiveData[i]];

        if( pointToParPoint_.find(pointI) == pointToParPoint_.end() )
        {
            pointToParPoint_[pointI] = parPoints_.size();
            parPoints_.append(points[pointI]);
        }
    }

    //- select face centres on all processors
    forAll(procBoundaries, patchI)
    {
        labelLongList activeFaces;

        const label start = procBoundaries[patchI].patchStart();
        const label end = start + procBoundaries[patchI].patchSize();

        //- find active faces
        for(label faceI=start;faceI<end;++faceI)
        {
            if
            (
                faceCentreToParPoint_.find(faceI) !=
                faceCentreToParPoint_.end()
            )
            {
                activeFaces.append(faceI-start);
            }
        }

        //- send the message to the other side of the patch
        OPstream toOtherProc
        (
            Pstream::blocking,
            procBoundaries[patchI].neiProcNo(),
            activeFaces.byteSize()
        );

        toOtherProc << activeFaces;
    }

    forAll(procBoundaries, patchI)
    {
        //- read the message from the other processor
        labelList activeFaces;

        IPstream fromOtherProc
        (
            Pstream::blocking,
            procBoundaries[patchI].neiProcNo()
        );

        fromOtherProc >> activeFaces;

        const label start = procBoundaries[patchI].patchStart();

        //- check active faces
        forAll(activeFaces, i)
        {
            const label faceI = start + activeFaces[i];
            if
            (
                faceCentreToParPoint_.find(faceI) ==
                faceCentreToParPoint_.end()
            )
            {
                //- activate this face
                faceCentreToParPoint_[faceI] = parPoints_.size();
                parPoints_.append(faceCentres_[faceI]);
            }
        }
    }

    //- calculate global labels
    //- global labels of vertices originating from existing points
    labelList nObjectsAtProc(Pstream::nProcs());
    label& counter = nObjectsAtProc[Pstream::myProcNo()];
    counter = 0;
    for
    (
        std::map<label, label>::const_iterator it = pointToParPoint_.begin();
        it!=pointToParPoint_.end();
        ++it
    )
    {
        //- count the points local to this proc
        //- or at te processor with the smallest label
        const label pointI = it->first;

        label minProc(Pstream::myProcNo());
        forAllRow(pointAtProcs, pointI, i)
            minProc = min(minProc, pointAtProcs(pointI, i));

        if( minProc == Pstream::myProcNo() )
            ++counter;
    }

    //- vertices originating proc face centres
    for
    (
        std::map<label, label>::const_iterator it =
            faceCentreToParPoint_.begin();
        it!=faceCentreToParPoint_.end();
        ++it
    )
    {
        const label procPatchI = origMesh_.faceIsInProcPatch(it->first);

        if( procPatchI < 0 || procBoundaries[procPatchI].owner() )
        {
            ++counter;
        }
    }

    //- vertices originating from cell centres
    counter += cellCentreToParPoint_.size();

    //- exchange data across the team of processors
    Pstream::gatherList(nObjectsAtProc);
    Pstream::scatterList(nObjectsAtProc);

    //- calculate starting label at each proc
    label startLabel(0);
    for(label i=0;i<Pstream::myProcNo();++i)
        startLabel += nObjectsAtProc[i];

    //- assign global labels to points
    globalParPointLabel_.setSize(parPoints_.size());
    globalParPointLabel_ = -1;
    globalToLocalParPointLabel_.clear();

    for
    (
        std::map<label, label>::const_iterator it = pointToParPoint_.begin();
        it!=pointToParPoint_.end();
        ++it
    )
    {
        //- count the points local to this proc
        //- or at te processor with the smallest label
        const label pointI = it->first;

        label minProc(Pstream::myProcNo());
        forAllRow(pointAtProcs, pointI, i)
            minProc = min(minProc, pointAtProcs(pointI, i));

        if( minProc == Pstream::myProcNo() )
        {
            globalToLocalParPointLabel_[startLabel] = it->second;
            globalParPointLabel_[it->second] = startLabel++;
        }
    }

    //- vertices originating from proc face centres
    for
    (
        std::map<label, label>::const_iterator it =
            faceCentreToParPoint_.begin();
        it!=faceCentreToParPoint_.end();
        ++it
    )
    {
        const label procPatchI = origMesh_.faceIsInProcPatch(it->first);

        //- set the label for all local faces (procPatchI < 0)
        //- faces at inter-processor boundaries - set labels at owner processor
        if( (procPatchI < 0) || origMesh_.procBoundaries()[procPatchI].owner() )
        {
            globalToLocalParPointLabel_[startLabel] = it->second;
            globalParPointLabel_[it->second] = startLabel++;
        }
    }

    //- vertices originating from cell centres
    for
    (
        std::map<label, label>::const_iterator it =
            cellCentreToParPoint_.begin();
        it!=cellCentreToParPoint_.end();
        ++it
    )
    {
        globalToLocalParPointLabel_[startLabel] = it->second;
        globalParPointLabel_[it->second] = startLabel++;
    }

    //- assign global label to vertices at slave processors
    exchangeData.clear();
    forAll(addr.pointNeiProcs(), i)
        exchangeData[addr.pointNeiProcs()[i]].clear();

    forAllConstIter(Map<label>, globalToLocal, it)
    {
        const label pointI = it();

        if( pointToParPoint_.find(pointI) != pointToParPoint_.end() )
        {
            const label parPointI = pointToParPoint_[pointI];

            if( globalParPointLabel_[parPointI] < 0 )
                continue;

            forAllRow(pointAtProcs, pointI, i)
            {
                const label neiProc = pointAtProcs(pointI, i);

                if( neiProc == Pstream::myProcNo() )
                    continue;

                labelLongList& dts = exchangeData[neiProc];
                //- send:
                //- 1. global point label
                //- 2. global par point label
                dts.append(it.key());
                dts.append(globalParPointLabel_[parPointI]);
            }
        }
    }

    receiveData.clear();
    help::exchangeMap(exchangeData, receiveData);

    for(label i=0;i<receiveData.size();)
    {
        const label pointI = globalToLocal[receiveData[i++]];
        const label globalParPointI = receiveData[i++];

        if( pointToParPoint_.find(pointI) == pointToParPoint_.end() )
        {
            pointToParPoint_[pointI] = parPoints_.size();
            parPoints_.append(points[pointI]);
        }

        const label parPointI = pointToParPoint_[pointI];

        globalParPointLabel_[parPointI] = globalParPointI;
        globalToLocalParPointLabel_[globalParPointI] = parPointI;
    }

    exchangeData.clear();

    //- assign global labels of face centres
    forAll(procBoundaries, patchI)
    {
        if( !procBoundaries[patchI].owner() )
            continue;

        labelLongList dataToSend;

        const label start = procBoundaries[patchI].patchStart();
        const label size = procBoundaries[patchI].patchSize();

        for(label fI=0;fI<size;++fI)
        {
            if
            (
                faceCentreToParPoint_.find(start+fI) !=
                faceCentreToParPoint_.end()
            )
            {
                const label parPointI = faceCentreToParPoint_[start+fI];

                dataToSend.append(fI);
                dataToSend.append(globalParPointLabel_[parPointI]);
            }
        }

        OPstream toOtherProc
        (
            Pstream::blocking,
            procBoundaries[patchI].neiProcNo(),
            dataToSend.byteSize()
        );

        toOtherProc << dataToSend;
    }

    forAll(procBoundaries, patchI)
    {
        if( procBoundaries[patchI].owner() )
            continue;

        //- get the data sent from the neighbour
        labelList receiveData;
        IPstream fromOtherProc
        (
            Pstream::blocking,
            procBoundaries[patchI].neiProcNo()
        );

        fromOtherProc >> receiveData;

        const label start = procBoundaries[patchI].patchStart();

        for(label i=0;i<receiveData.size();)
        {
            const label faceI = start + receiveData[i++];
            const label globalParPointI = receiveData[i++];

            if
            (
                faceCentreToParPoint_.find(faceI) ==
                faceCentreToParPoint_.end()
            )
            {
                faceCentreToParPoint_[faceI] = parPoints_.size();
                parPoints_.append(faceCentres_[faceI]);
            }

            //- set the global label for a point at face centre
            const label parPointI = faceCentreToParPoint_[faceI];
            globalParPointLabel_[parPointI] = globalParPointI;
            globalToLocalParPointLabel_[globalParPointI] = parPointI;
        }
    }

    # ifdef DEBUGSmooth
    Pout << "Finished tets at inter-processor points" << endl;
    returnReduce(1, sumOp<label>());

    const labelLongList& globalPointLabel = addr.globalPointLabel();

    for(label i=0;i<Pstream::nProcs();++i)
    {
        if( i == Pstream::myProcNo() )
        {
            Pout << "1. Points       "
                 << label(pointToParPoint_.size()) << endl;
            Pout << "1. Face centres "
                 << label(faceCentreToParPoint_.size()) << endl;
            Pout << "1. Cell centres "
                 << label(cellCentreToParPoint_.size()) << endl;

            //Pout << "Global labels " << globalParPointLabel_ << endl;

            label nInvalid(0);
            for
            (
                std::map<label, label>::const_iterator it =
                    pointToParPoint_.begin();
                it!=pointToParPoint_.end();
                ++it
            )
                if( globalParPointLabel_[it->second] < 0 )
                {
                    nInvalid++;
                    Pout << "Par point " << it->second
                         << " from point " << it->first
                         << " with global label " << globalPointLabel[it->first]
                         << " has no assigned global label. Point at procs "
                         << pointAtProcs[it->first] << endl;
                }

            Pout << "Invalid points " << nInvalid << endl;

            forAll(globalParPointLabel_, i)
                if( globalParPointLabel_[i] < 0 )
                    FatalError << "Par point " << i << " has negative global"
                         << " label" << abort(FatalError);
        }

        returnReduce(1, sumOp<label>());
    }

    for
    (
        std::map<label, DynList<partTet, 256> >::const_iterator tIt =
            tetsAtPoint_.begin();
        tIt!=tetsAtPoint_.end();
        ++tIt
    )
    {
        partTetMeshSimplex simplex(*this, tIt->first);

        simplex.writeToVTK
        (
            "simplexLocal_"+help::labelToText(Pstream::myProcNo())+"_"+
            help::labelToText(globalPointLabel[tIt->first])+".vtk"
        );
    }

    Pout << "Finished simplices after tets at inter-processor points" << endl;
    returnReduce(1, sumOp<label>());
    # endif
}

void partTetMesh::createBufferLayers()
{
    //- the goal of this function is to reconstruct the simplices
    //- all points shared by more than one processor
    const labelLongList& globalPointLabel =
        origMesh_.addressingData().globalPointLabel();
    const Map<label>& globalToLocal =
        origMesh_.addressingData().globalToLocalPointAddressing();
    const VRWGraph& pointAtProcs = origMesh_.addressingData().pointAtProcs();
    const DynList<label>& pNeiProcs =
        origMesh_.addressingData().pointNeiProcs();

    //- find the points that shall be transferred to the neighbouring processor
    std::map<label, std::set<label> > pointsToExchange;
    std::map<label, LongList<labelledPoint> > exchangePoints;
    std::map<label, labelLongList> exchangeTets;
    std::map<label, labelLongList> exchangeBndTriangles;

    forAll(pNeiProcs, i)
    {
        pointsToExchange[pNeiProcs[i]].clear();
        exchangePoints[pNeiProcs[i]].clear();
        exchangeTets[pNeiProcs[i]].clear();
        exchangeBndTriangles[pNeiProcs[i]].clear();
    }

    //- prepare data that shall be sent to other processors
    for
    (
        std::map<label, DynList<partTet, 256> >::const_iterator it =
            tetsAtPoint_.begin();
        it!=tetsAtPoint_.end();
        ++it
    )
    {
        const label pointI = it->first;
        const DynList<partTet, 256>& tets = it->second;

        forAllRow(pointAtProcs, pointI, i)
        {
            const label neiProc = pointAtProcs(pointI, i);

            if( neiProc == Pstream::myProcNo() )
                continue;

            std::set<label>& ptsToExchange = pointsToExchange[neiProc];
            LongList<labelledPoint>& pts = exchangePoints[neiProc];
            labelLongList& sendTets = exchangeTets[neiProc];
            labelLongList& sendTrias = exchangeBndTriangles[neiProc];

            //- tets are sent without the last point as it is always the same
            //- the last point in a tet and the number of tets are added first
            sendTets.append(globalPointLabel[pointI]);
            sendTets.append(globalParPointLabel_[pointToParPoint_[pointI]]);
            sendTets.append(tets.size());

            forAll(tets, tetI)
            {
                const partTet& tet = tets[tetI];

                //- send the first point
                if( ptsToExchange.find(tet.a()) == ptsToExchange.end() )
                {
                    ptsToExchange.insert(tet.a());
                    pts.append
                    (
                        labelledPoint
                        (
                            globalParPointLabel_[tet.a()],
                            parPoints_[tet.a()]
                        )
                    );

                    //- other procesors for this point
                    parPointAtOtherProcs_[tet.a()].insert(neiProc);
                }

                //- send the second point
                if( ptsToExchange.find(tet.b()) == ptsToExchange.end() )
                {
                    ptsToExchange.insert(tet.b());
                    pts.append
                    (
                        labelledPoint
                        (
                            globalParPointLabel_[tet.b()],
                            parPoints_[tet.b()]
                        )
                    );

                    //- other procesors for this point
                    parPointAtOtherProcs_[tet.b()].insert(neiProc);
                }

                //- send the third point
                if( ptsToExchange.find(tet.c()) == ptsToExchange.end() )
                {
                    ptsToExchange.insert(tet.c());
                    pts.append
                    (
                        labelledPoint
                        (
                            globalParPointLabel_[tet.c()],
                            parPoints_[tet.c()]
                        )
                    );

                    //- other procesors for this point
                    parPointAtOtherProcs_[tet.c()].insert(neiProc);
                }

                //- send the first three points of the tet
                sendTets.append(globalParPointLabel_[tet.a()]);
                sendTets.append(globalParPointLabel_[tet.b()]);
                sendTets.append(globalParPointLabel_[tet.c()]);
            }

            std::map<label, DynList<labelledTri, 32> >::const_iterator bIt =
                bndTrianglesAtPoint_.find(pointI);
            if( bIt != bndTrianglesAtPoint_.end() )
            {
                const DynList<labelledTri, 32>& bndTriangles = bIt->second;

                sendTrias.append(globalPointLabel[pointI]);
                sendTrias.append(bndTriangles.size());

                forAll(bndTriangles, triI)
                {
                    const labelledTri& tri = bndTriangles[triI];

                    //- send the first point
                    if( ptsToExchange.find(tri[0]) == ptsToExchange.end() )
                    {
                        ptsToExchange.insert(tri[0]);
                        pts.append
                        (
                            labelledPoint
                            (
                                globalParPointLabel_[tri[0]],
                                parPoints_[tri[0]]
                            )
                        );

                        //- other procesors for this point
                        parPointAtOtherProcs_[tri[0]].insert(neiProc);
                    }

                    //- send the second point
                    if( ptsToExchange.find(tri[1]) == ptsToExchange.end() )
                    {
                        ptsToExchange.insert(tri[1]);
                        pts.append
                        (
                            labelledPoint
                            (
                                globalParPointLabel_[tri[1]],
                                parPoints_[tri[1]]
                            )
                        );

                        //- other procesors for this point
                        parPointAtOtherProcs_[tri[1]].insert(neiProc);
                    }

                    //- send the third point
                    if( ptsToExchange.find(tri[2]) == ptsToExchange.end() )
                    {
                        ptsToExchange.insert(tri[2]);
                        pts.append
                        (
                            labelledPoint
                            (
                                globalParPointLabel_[tri[2]],
                                parPoints_[tri[2]]
                            )
                        );

                        //- other procesors for this point
                        parPointAtOtherProcs_[tri[2]].insert(neiProc);
                    }

                    //- send the first three points of the tet
                    sendTrias.append(globalParPointLabel_[tri[0]]);
                    sendTrias.append(globalParPointLabel_[tri[1]]);
                    sendTrias.append(globalParPointLabel_[tri[2]]);
                    sendTrias.append(tri.region());
                }
            }
        }
    }

    # ifdef DEBUGSmooth
    for(label i=0;i<Pstream::nProcs();++i)
    {
        if( i == Pstream::myProcNo() )
        {
            Pout << "2. Points       "
                 << label(pointToParPoint_.size()) << endl;
            Pout << "2. Face centres "
                 << label(faceCentreToParPoint_.size()) << endl;
            Pout << "2. Cell centres "
                 << label(cellCentreToParPoint_.size()) << endl;

            for
            (
                std::map<label, LongList<labelledPoint> >::const_iterator it =
                    exchangePoints.begin();
                it!=exchangePoints.end();
                ++it
            )
                Pout << "To proc " << it->first
                     << " points " << it->second << endl;
        }

        returnReduce(1, sumOp<label>());
    }
    # endif

    //- exchange points among processors
    LongList<labelledPoint> receivePoints;
    help::exchangeMap(exchangePoints, receivePoints);
    exchangePoints.clear();
    pointsToExchange.clear();

    forAll(receivePoints, i)
    {
        const labelledPoint& lp = receivePoints[i];

        if
        (
            globalToLocalParPointLabel_.find(lp.pointLabel()) ==
            globalToLocalParPointLabel_.end()
        )
        {
            globalParPointLabel_.append(lp.pointLabel());
            globalToLocalParPointLabel_[lp.pointLabel()] = parPoints_.size();
            parPoints_.append(lp.coordinates());
        }
    }

    //- avoid updating points that are available locally
    forAllConstIter(Map<label>, globalToLocal, it)
    {
        const label pointI = it();

        if( pointToParPoint_.find(pointI) != pointToParPoint_.end() )
        {
            const label parPointI = pointToParPoint_[pointI];
            std::set<label>& neiProcs = parPointAtOtherProcs_[parPointI];

            std::set<label> removeElement;
            forAllIter(std::set<label>, neiProcs, npIt)
                if( pointAtProcs.contains(pointI, *npIt) )
                    removeElement.insert(*npIt);

            forAllConstIter(std::set<label>, removeElement, rIt)
                neiProcs.erase(*rIt);
        }
    }

    //- avoid updating face centres available locally
    forAll(origMesh_.procBoundaries(), patchI)
    {
        const processorBoundaryPatch& patch =
            origMesh_.procBoundaries()[patchI];
        const label start = patch.patchStart();
        const label end = start + patch.patchSize();

        for(label faceI=start;faceI<end;++faceI)
        {
            if
            (
                faceCentreToParPoint_.find(faceI) !=
                faceCentreToParPoint_.end()
            )
            {
                const label parPointI = faceCentreToParPoint_[faceI];

                std::set<label>& otherProcs = parPointAtOtherProcs_[parPointI];

                if( otherProcs.find(patch.neiProcNo()) != otherProcs.end() )
                {
                    otherProcs.erase(patch.neiProcNo());
                }
            }
        }
    }

    //- exchange tets among processors
    labelLongList receiveTets;
    help::exchangeMap(exchangeTets, receiveTets);
    exchangeTets.clear();

    for(label i=0;i<receiveTets.size();)
    {
        const label pointI = globalToLocal[receiveTets[i++]];

        if( tetsAtPoint_.find(pointI) == tetsAtPoint_.end() )
            FatalErrorIn
            (
                "void partTetMesh::createBufferLayers()"
            ) << "Point " << pointI << " is not at the interface"
              << abort(FatalError);

        DynList<partTet, 256>& tets = tetsAtPoint_[pointI];

        const label d = globalToLocalParPointLabel_[receiveTets[i++]];
        const label nTets = receiveTets[i++];

        //- read tets and append them to the list
        for(label tetI=0;tetI<nTets;++tetI)
        {
            const label a = globalToLocalParPointLabel_[receiveTets[i++]];
            const label b = globalToLocalParPointLabel_[receiveTets[i++]];
            const label c = globalToLocalParPointLabel_[receiveTets[i++]];

            tets.append(partTet(a, b, c, d));
        }
    }

    //- exchange boundary triangles
    receiveTets.clear();
    help::exchangeMap(exchangeBndTriangles, receiveTets);
    exchangeBndTriangles.clear();

    for(label i=0;i<receiveTets.size();)
    {
        const label pointI = globalToLocal[receiveTets[i++]];

        if( bndTrianglesAtPoint_.find(pointI) == bndTrianglesAtPoint_.end() )
            FatalErrorIn
            (
                "void partTetMesh::createBufferLayers()"
            ) << "Point " << pointI << " is not at the interface"
              << abort(FatalError);

        DynList<labelledTri, 32>& trias = bndTrianglesAtPoint_[pointI];

        const label nTets = receiveTets[i++];

        //- read tets and append them to the list
        for(label tetI=0;tetI<nTets;++tetI)
        {
            const label a = globalToLocalParPointLabel_[receiveTets[i++]];
            const label b = globalToLocalParPointLabel_[receiveTets[i++]];
            const label c = globalToLocalParPointLabel_[receiveTets[i++]];

            const label patchI = receiveTets[i++];

            trias.append(labelledTri(a, b, c, patchI));
        }
    }

    # ifdef DEBUGSmooth
    Pout << "Finished generating buffer layers" << endl;
    returnReduce(1, sumOp<label>());

    Pout << "3. Points       " << label(pointToParPoint_.size()) << endl;
    Pout << "3. Face centres " << label(faceCentreToParPoint_.size()) << endl;
    Pout << "3. Cell centres " << label(cellCentreToParPoint_.size()) << endl;

    for
    (
        std::map<label, DynList<partTet, 256> >::const_iterator tIt =
            tetsAtPoint_.begin();
        tIt!=tetsAtPoint_.end();
        ++tIt
    )
    {
        partTetMeshSimplex simplex(*this, tIt->first);

        simplex.writeToVTK
        (
            "simplexAll_"+help::labelToText(Pstream::myProcNo())+"_"+
            help::labelToText(globalPointLabel[tIt->first])+".vtk"
        );
    }

    Pout << "Finished simplices after buffer layers" << endl;
    returnReduce(1, sumOp<label>());
    # endif
}

void partTetMesh::unifyCoordinatesAtInterProcessorBoundaries()
{
    if( !Pstream::parRun() )
        return;

    const polyMeshGenAddressing& addr = origMesh_.addressingData();
    const Map<label>& globalToLocal = addr.globalToLocalPointAddressing();
    const VRWGraph& pointAtProcs = addr.pointAtProcs();
    const DynList<label>& neiProcs = addr.pointNeiProcs();

    //- create the map
    std::map<label, LongList<labelledPoint> > exchangeData;
    forAll(neiProcs, procI)
        exchangeData[neiProcs[procI]].clear();

    //- fill in the data
    std::map<label, labelledPoint> parallelBndPoints;

    forAllConstIter(Map<label>, globalToLocal, it)
    {
        const label pointI = it();

        if( smoothVertex_[pointI] & (SMOOTH|BOUNDARY) )
            continue;

        parallelBndPoints[pointI] = labelledPoint(1, points_[pointI]);

        forAllRow(pointAtProcs, pointI, procI)
        {
            const label neiProc = pointAtProcs(pointI, procI);

            if( neiProc == Pstream::myProcNo() )
                continue;

            exchangeData[neiProc].append
            (
                labelledPoint(it.key(), points_[pointI])
            );
        }
    }

    //- send points to other processors
    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    //- gather the data
    forAll(receivedData, i)
    {
        const labelledPoint& lp = receivedData[i];

        std::map<label, labelledPoint>::iterator iter =
            parallelBndPoints.find(globalToLocal[lp.pointLabel()]);

        if( iter == parallelBndPoints.end() )
            continue;

        ++iter->second.pointLabel();
        iter->second.coordinates() += lp.coordinates();
    }

    //- move the point to the averaged position
    polyMeshGenModifier meshModifier(origMesh_);

    boolList updateFaces(faces_.size(), false);

    for
    (
        std::map<label, labelledPoint>::iterator it=parallelBndPoints.begin();
        it!=parallelBndPoints.end();
        ++it
    )
    {
        const label pI = it->first;

        forAllRow(pointFaces_, pI, i)
            updateFaces[pointFaces_(pI, i)] = true;

        const point newP = it->second.coordinates() / it->second.pointLabel();
        meshModifier.movePoint(pI, newP);
    }

    //- update geometry
    origMesh_.addressingData().updateGeometry(updateFaces);
}

void partTetMesh::updateBufferLayerPoints()
{
    # ifdef DEBUGSmooth
    Pout << "Updating buffer layers " << parPoints_.size() << endl;
    returnReduce(1, sumOp<label>());
    # endif

    //- update local points in the buffer layer
    for
    (
        std::map<label, label>::const_iterator it=pointToParPoint_.begin();
        it!=pointToParPoint_.end();
        ++it
    )
    {
        parPoints_[it->second] = points_[it->first];
    }

    //- update local face centres
    for
    (
        std::map<label, label>::const_iterator it=faceCentreToParPoint_.begin();
        it!=faceCentreToParPoint_.end();
        ++it
    )
    {
        parPoints_[it->second] = faceCentres_[it->first];
    }

    //- update local cell centres
    for
    (
        std::map<label, label>::const_iterator it=cellCentreToParPoint_.begin();
        it!=cellCentreToParPoint_.end();
        ++it
    )
    {
        parPoints_[it->second] = cellCentres_[it->first];
    }

    //- exchange the data among processors
    const DynList<label>& neiProcs = origMesh_.addressingData().pointNeiProcs();

    //- create the map
    std::map<label, LongList<labelledPoint> > exchangeData;
    forAll(neiProcs, i)
        exchangeData[neiProcs[i]].clear();

    //- add points into the map
    for
    (
        std::map<label, std::set<label> >::const_iterator it =
            parPointAtOtherProcs_.begin();
        it!=parPointAtOtherProcs_.end();
        ++it
    )
    {
        const label parPointI = it->first;
        const std::set<label>& neiProcs = it->second;

        forAllConstIter(std::set<label>, neiProcs, it)
        {
            const label neiProc = *it;

            exchangeData[neiProc].append
            (
                labelledPoint
                (
                    globalParPointLabel_[parPointI],
                    parPoints_[parPointI]
                )
            );
        }
    }

    LongList<labelledPoint> receivedData;
    help::exchangeMap(exchangeData, receivedData);

    //- update the coordinates
    forAll(receivedData, i)
    {
        const labelledPoint& lp = receivedData[i];

        parPoints_[globalToLocalParPointLabel_[lp.pointLabel()]] =
            lp.coordinates();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
