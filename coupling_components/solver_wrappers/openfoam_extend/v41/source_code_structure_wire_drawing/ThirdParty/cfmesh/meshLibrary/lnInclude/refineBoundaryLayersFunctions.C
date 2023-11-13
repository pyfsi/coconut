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

#include "refineBoundaryLayers.H"
#include "meshSurfaceEngine.H"
#include "helperFunctions.H"
#include "polyMeshGenAddressing.H"
#include "polyMeshGen2DEngine.H"
#include "VRWGraphList.H"
#include "meshSurfacePartitioner.H"
#include "detectBoundaryLayers.H"
#include "splitEdgeAttributes.H"

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



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool refineBoundaryLayers::analyseLayers()
{
    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelLongList& facePatch = mse.boundaryFacePatches();

    meshSurfacePartitioner mPart(mse);
    detectBoundaryLayers dbl(mPart, is2DMesh_);

    const label nGroups = dbl.nDistinctLayers();
    const labelLongList& faceInLayer = dbl.faceInLayer();
    faceInLayer_ = faceInLayer;

    //- get the hair edges
    splitEdges_ = dbl.hairEdges();

    # ifdef DEBUGLayer
    OFstream file("hairEdges.vtk");

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    //- write points
    file << "POINTS " << 2*splitEdges_.size() << " float\n";
    forAll(splitEdges_, seI)
    {
        const point& p = mse.mesh().points()[splitEdges_[seI].start()];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;

        const point op = mse.mesh().points()[splitEdges_[seI].end()];

        file << op.x() << ' ' << op.y() << ' ' << op.z() << nl;
    }

    //- write lines
    file << "\nLINES " << splitEdges_.size()
         << " " << 3*splitEdges_.size() << nl;
    forAll(splitEdges_, eI)
    {
        file << 2 << " " << 2*eI << " " << (2*eI+1) << nl;
    }

    file << "\n";
    # endif

    //- create point to split edges addressing
    splitEdgesAtPoint_.reverseAddressing(splitEdges_);

    //- check if the layer is valid
    bool validLayer(true);
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(faceInLayer, bfI)
    {
        if( faceInLayer[bfI] < 0 )
            continue;

        const face& bf = bFaces[bfI];

        forAll(bf, pI)
            if( splitEdgesAtPoint_.sizeOfRow(bf[pI]) == 0 )
                validLayer = false;
    }

    # ifdef DEBUGLayer
    Info << "Number of independent layers in the mesh is " << nGroups << endl;
    Info << "Is valid layer " << validLayer << endl;
    # endif

    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();

    //- create patch name to index addressing
    std::map<word, label> patchNameToIndex;
    forAll(boundaries, patchI)
        patchNameToIndex[boundaries[patchI].patchName()] = patchI;

    //- check layer labels over a patch
    layerAtPatch_.setSize(boundaries.size());
    forAll(layerAtPatch_, i)
        layerAtPatch_[i].clear();
    List<DynList<label> > groupsAtPatch(boundaries.size());
    forAll(faceInLayer, bfI)
        groupsAtPatch[facePatch[bfI]].appendIfNotIn(faceInLayer[bfI]);

    //- set the information which patches have an extruded layer
    forAll(groupsAtPatch, patchI)
    {
        const DynList<label>& layers = groupsAtPatch[patchI];

        forAll(layers, i)
        {
            if( layers[i] < 0 )
            {
                layerAtPatch_[patchI].clear();
                break;
            }
            else
            {
                layerAtPatch_[patchI].append(layers[i]);
            }
        }
    }

    # ifdef DEBUGLayer
    Info << "Layer at patch " << layerAtPatch_ << endl;
    # endif

    //- set the information which patches are a single boundary layer face
    patchesInLayer_.setSize(nGroups);
    forAll(layerAtPatch_, patchI)
    {
        const DynList<label>& layers = layerAtPatch_[patchI];

        forAll(layers, i)
            patchesInLayer_[layers[i]].append
            (
                boundaries[patchI].patchName()
            );
    }

    # ifdef DEBUGLayer
    Info << "Patches in layer " << patchesInLayer_ << endl;
    # endif

    //- set the number of boundary layers for each patch
    labelList nLayersAtPatch(layerAtPatch_.size(), -1);
    boolList protectedValue(layerAtPatch_.size(), false);

    forAll(patchesInLayer_, layerI)
    {
        const DynList<word>& layerPatches = patchesInLayer_[layerI];

        label maxNumLayers(1);
        bool hasLocalValue(false);

        //- find the maximum requested number of layers over the layer
        forAll(layerPatches, lpI)
        {
            const word pName = layerPatches[lpI];

            std::map<word, label>::const_iterator it =
                numLayersForPatch_.find(pName);

            if( it != numLayersForPatch_.end() )
            {
                //- check if the layer is interrupted at this patch
                if(
                    discontinuousLayersForPatch_.find(pName) !=
                    discontinuousLayersForPatch_.end()
                )
                {
                    //- set the number of layers and lock this location
                    nLayersAtPatch[patchNameToIndex[pName]] = it->second;
                    protectedValue[patchNameToIndex[pName]] = true;
                    hasLocalValue = true;
                }
                else
                {
                    //- take the maximum number of layers
                    maxNumLayers = Foam::max(maxNumLayers, it->second);
                    hasLocalValue = true;
                }
            }
        }

        //- apply the global value if no local values exist
        if( !hasLocalValue )
            maxNumLayers = globalNumLayers_;

        //- apply the maximum number of ayer of all unprotected patches
        forAll(layerPatches, lpI)
        {
            const label ptchI = patchNameToIndex[layerPatches[lpI]];

            if( !protectedValue[ptchI] )
                nLayersAtPatch[ptchI] = maxNumLayers;
        }
    }

    if( is2DMesh_ )
    {
        polyMeshGen2DEngine mesh2DEngine(mesh_);
        const boolList& zMinPoint = mesh2DEngine.zMinPoints();
        const boolList& zMaxPoint = mesh2DEngine.zMaxPoints();

        const faceList::subList& bFaces = mse.boundaryFaces();

        boolList allZMax(mesh_.boundaries().size(), true);
        boolList allZMin(mesh_.boundaries().size(), true);

        # ifdef USE_OMP
        # pragma omp parallel for schedule(dynamic, 50)
        # endif
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            forAll(bf, pI)
            {
                if( !zMinPoint[bf[pI]] )
                    allZMin[facePatch[bfI]] = false;
                if( !zMaxPoint[bf[pI]] )
                    allZMax[facePatch[bfI]] = false;
            }
        }

        //- mark empty patches as already used
        forAll(allZMin, patchI)
        {
            if( allZMin[patchI] ^ allZMax[patchI] )
            {
                nLayersAtPatch[patchI] = -1;
                layerAtPatch_[patchI].clear();
            }
        }
    }

    //- perform reduction over all processors
    reduce(nLayersAtPatch, maxOp<labelList>());

    # ifdef DEBUGLayer
    Pout << "nLayersAtPatch " << nLayersAtPatch << endl;
    # endif

    //- set the number of boundary layers which shall be generated above
    //- each boundary face
    nLayersAtBndFace_.setSize(facePatch.size());
    nLayersAtBndFace_ = globalNumLayers_;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(nLayersAtBndFace_, bfI)
    {
        const label patchI = facePatch[bfI];

        if( nLayersAtPatch[patchI] < 0 )
        {
            nLayersAtBndFace_[bfI] = 1;
        }
        else
        {
            nLayersAtBndFace_[bfI] = nLayersAtPatch[patchI];

            if( specialMode_ )
            {
                ++nLayersAtBndFace_[bfI];
            }
        }
    }

    reduce(validLayer, minOp<bool>());

    # ifdef DEBUGLayer
    forAll(nLayersAtBndFace_, bfI)
    Pout << "Boundary face " << bfI << " in patch "
        << facePatch[bfI] << " num layers " << nLayersAtBndFace_[bfI] << endl;
    //::exit(1);
    # endif

    return validLayer;
}

void refineBoundaryLayers::skipRefinementAtBadQualityCells()
{
    if( specialMode_ || !skipRefinementOfPoorQualityCells_ )
        return;


    const meshSurfaceEngine& mse = surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelLongList& bp = mse.bp();
    const labelLongList& facePatches = mse.boundaryFacePatches();
    const VRWGraph& pFaces = mse.pointFaces();

    const PtrList<boundaryPatch>& patches = mesh_.boundaries();

    const scalar angleCosTol = cos(maxOrthogonalityDeviationAngle_);

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(bFaces, bfI)
    {
        const label nLayers = nLayersAtBndFace_[bfI];
        if( nLayers <= 1 )
            continue;

        //- check the cell quality and disallow layer refinement
        //- at faces with low quality cells
        //- get the patch name of the face first
        const word pName = patches[facePatches[bfI]].patchName();

        //- max face deviation from flatness
        const face& bf = bFaces[bfI];

        //- calculate face centre
        const point fc = help::faceCentre(points, bf);

        //- calculate face normal
        vector n = help::faceAreaVector(points, bf);
        n /= (mag(n) + VSMALL);

        //- check face flatness
        scalar maxHeight(0.0);
        forAll(bf, pI)
        {
            const point& p = points[bf[pI]];
            const triangle<point, point> tri
            (
                p,
                points[bf.nextLabel(pI)],
                fc
            );

            vector n = tri.normal();
            n /= (mag(n) + VSMALL);

            forAll(bf, pJ)
            {
                const scalar height = mag((points[bf[pJ]] - p) & n);

                maxHeight = max(maxHeight, height);
            }
        }

        //- find hair edges at each face vertex
        forAll(bf, pI)
        {
            const label pointI = bf[pI];

            forAllRow(splitEdgesAtPoint_, pointI, i)
            {
                const edge& se = splitEdges_[splitEdgesAtPoint_(pointI, i)];

                //- calculate edge length
                vector ev = se.vec(points);
                const scalar magv = mag(ev);
                ev /= (magv + VSMALL);

                //- check if the face contains this split edge
                if( bf.which(se.end()) >= 0 )
                    continue;

                //- check the thickness ratio of the first layer
                scalar thicknessRatio(globalThicknessRatio_);
                if
                (
                    thicknessRatioForPatch_.find(pName) !=
                    thicknessRatioForPatch_.end()
                )
                    thicknessRatio = thicknessRatioForPatch_[pName];

                //- get the desired thickness of the first layer
                scalar firstThickness =
                    magv /
                    (
                        (1 - Foam::pow(thicknessRatio, nLayers))
                        / (1.0 - thicknessRatio)
                    );

                if
                (
                    maxThicknessForPatch_.find(pName) !=
                    maxThicknessForPatch_.end()
                )
                {
                    firstThickness =
                        min(firstThickness, maxThicknessForPatch_[pName]);
                }
                else
                {
                    firstThickness =
                        min(firstThickness, globalMaxThicknessFirstLayer_);
                }

                //- check if the first layer thickness is greater than the
                //- deviation of the face from flatness
                if
                (
                    (firstThickness / maxHeight) <
                    minThicknessToFlatnessDeviationHeightRatio_
                )
                {
                    //- thickness of the first layer is smaller than the
                    //- deviation from flatness. Cannot refine the layer
                    nLayersAtBndFace_[bfI] = 1;
                }

                //- check the alignment of the hair edge vector
                //- with the face normal
                if( (ev & n) > -angleCosTol )
                {
                    nLayersAtBndFace_[bfI] = 1;
                }

                if( nLayersAtBndFace_[bfI] == 1 )
                    break;
            }

            if( nLayersAtBndFace_[bfI] == 1 )
                break;
        }
    }

    //- propagate changes in the neighbourhood to avoid many layers to connect
    //- to a single one
    label layerI(1);
    bool changed;

    labelLongList nLayersAtBndFaceCopy(bFaces.size());

    if( Pstream::parRun() )
    {
        //- calculate addressing outside of a parallel region
        mse.globalToLocalBndPointAddressing();
        mse.bpAtProcs();
    }

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        forAll(nLayersAtBndFaceCopy, bfI)
            nLayersAtBndFaceCopy[bfI] = nLayersAtBndFace_[bfI];

        const label nTasks = 5 * omp_get_num_threads();
        # else
        forAll(nLayersAtBndFaceCopy, bfI)
            nLayersAtBndFaceCopy[bfI] = nLayersAtBndFace_[bfI];

        const label nTasks = 1;
        # endif
        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            do
            {
                changed = false;

                //- modify number of layers at neighbouring faces
                for(label taskI=0;taskI<nTasks;++taskI)
                {
                    # ifdef USE_OMP
                    # pragma omp task default(shared) firstprivate(taskI)
                    # endif
                    {
                        const label maxL = layerI + 1;
                        for(label bfI=taskI;bfI<bFaces.size();bfI+=nTasks)
                        {
                            if( nLayersAtBndFace_[bfI] == layerI )
                            {
                                const label layerIndex = faceInLayer_[bfI];

                                if( layerIndex < 0 )
                                    continue;

                                const face& bf = bFaces[bfI];

                                forAll(bf, pI)
                                {
                                    const label bpI = bp[bf[pI]];

                                    forAllRow(pFaces, bpI, pfI)
                                    {
                                        const label bfJ = pFaces(bpI, pfI);

                                        //- faces must be in the same layer
                                        if( faceInLayer_[bfJ] != layerIndex )
                                            continue;

                                        //- limit the number of layers
                                        if( maxL < nLayersAtBndFace_[bfJ] )
                                        {
                                            nLayersAtBndFaceCopy[bfJ] = maxL;
                                            changed = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                # ifdef USE_OMP
                # pragma omp taskwait
                # pragma omp flush(changed)
                # endif

                //- exchange data across processors
                if( Pstream::parRun() )
                {
                    const Map<label>& globalToLocal =
                        mse.globalToLocalBndPointAddressing();
                    const VRWGraph& bpAtProcs = mse.bpAtProcs();

                    std::map<label, labelLongList> exchangeData;
                    forAll(mse.bpNeiProcs(), i)
                        exchangeData[mse.bpNeiProcs()[i]].clear();

                    //- create the data that hall be sent
                    forAllConstIter(Map<label>, globalToLocal, it)
                    {
                        const label bpI = it();

                        forAllRow(pFaces, bpI, pfI)
                        {
                            const label bfI = pFaces(bpI, pfI);

                            if(  nLayersAtBndFace_[bfI] == layerI )
                            {
                                const label layerIndex = faceInLayer_[bfI];

                                if( layerIndex < 0 )
                                    continue;

                                forAllRow(bpAtProcs, bpI, i)
                                {
                                    const label neiProc = bpAtProcs(bpI, i);

                                    if( neiProc == Pstream::myProcNo() )
                                        continue;

                                    //- send:
                                    //- 1. Global point label
                                    //- 2. layer index
                                    labelLongList& dataToSend =
                                        exchangeData[neiProc];
                                    dataToSend.append(it.key());
                                    dataToSend.append(layerIndex);
                                }
                            }
                        }
                    }

                    //- exchange data among processors
                    labelLongList receiveData;
                    help::exchangeMap(exchangeData, receiveData);

                    //- apply data locally
                    for(label i=0;i<receiveData.size();)
                    {
                        const label bpI = globalToLocal[receiveData[i++]];
                        const label layerIndex = receiveData[i++];

                        forAllRow(pFaces, bpI, pfI)
                        {
                            const label bfI = pFaces(bpI, pfI);

                            if( faceInLayer_[bfI] != layerIndex )
                                continue;

                            if( (layerI + 1) < nLayersAtBndFace_[bfI] )
                            {
                                nLayersAtBndFaceCopy[bfI] = layerI + 1;
                                changed = true;
                            }
                        }
                    }

                    reduce(changed, maxOp<bool>());
                }

                if( changed )
                {
                    //- copy the content of the copy list back to the original
                    for(label taskI=0;taskI<nTasks;++taskI)
                    {
                        # ifdef USE_OMP
                        # pragma omp task default(shared) firstprivate(taskI)
                        # endif
                        {
                            for(label bfI=taskI;bfI<bFaces.size();bfI+=nTasks)
                                nLayersAtBndFace_[bfI] =
                                    nLayersAtBndFaceCopy[bfI];
                        }
                    }

                    # ifdef USE_OMP
                    # pragma omp taskwait
                    # endif
                }

                //- increment the layer index
                ++layerI;

            } while( changed );
        }
    }

    # ifdef DEBUGLayer
    std::map<label, label> nLayersToId;
    forAll(nLayersAtBndFace_, bfI)
    {
        if( nLayersToId.find(nLayersAtBndFace_[bfI]) == nLayersToId.end() )
            nLayersToId[nLayersAtBndFace_[bfI]] =
                mesh_.addFaceSubset
                (
                    "numLayers_"+help::labelToText(nLayersAtBndFace_[bfI])
                );

        mesh_.addFaceToSubset
        (
            nLayersToId[nLayersAtBndFace_[bfI]],
            mesh_.nInternalFaces()+bfI
        );
    }

    mesh_.write();
    # endif
}

void refineBoundaryLayers::generateNewVertices()
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    const pointFieldPMG& points = mesh_.points();

    const meshSurfaceEngine& mse = surfaceEngine();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const VRWGraph& pointFaces = mse.pointFaces();
    const labelLongList& facePatch = mse.boundaryFacePatches();
    const labelLongList& bp = mse.bp();

    //- allocate the data from storing parameters applying to a split edge
    LongList<scalar> firstLayerThickness(splitEdges_.size());
    LongList<scalar> thicknessRatio(splitEdges_.size());
    labelLongList nNodesAtEdge(splitEdges_.size());
    labelLongList nLayersAtEdge(splitEdges_.size());

    //- count the number of vertices for each split edge
    label nThreads(1);
    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        //- start counting vertices at each thread
        # ifdef USE_OMP
        # pragma omp single
        {
            nThreads = omp_get_num_threads();
        }

        # pragma omp for schedule(static, 1)
        # endif
        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            //- get the requested number of boundary layers
            label nLayers(1);
            scalar ratio(globalThicknessRatio_);
            scalar thickness(globalMaxThicknessFirstLayer_);
            bool overridenThickness(false);

            const label bpI = bp[e.start()];

            forAllRow(pointFaces, bpI, pfI)
            {
                const label bfI = pointFaces(bpI, pfI);
                const label pos = help::positionOfEdgeInFace(e, bFaces[bfI]);
                if( pos >= 0 )
                    continue;

                const word& patchName =
                    boundaries[facePatch[bfI]].patchName();

                //- overrride the global value with the maximum number of layers
                //- at this edge
                nLayers = Foam::max(nLayers, nLayersAtBndFace_[bfI]);

                //- override with the maximum ratio
                const std::map<word, scalar>::const_iterator rIt =
                    thicknessRatioForPatch_.find(patchName);
                if( rIt != thicknessRatioForPatch_.end() )
                {
                    ratio = rIt->second;
                }

                //- override with the minimum thickness set for this edge
                const std::map<word, scalar>::const_iterator tIt =
                    maxThicknessForPatch_.find(patchName);
                if( tIt != maxThicknessForPatch_.end() )
                {
                    if( overridenThickness )
                    {
                        thickness = Foam::min(thickness, tIt->second);
                    }
                    else
                    {
                        thickness = tIt->second;
                        overridenThickness = true;
                    }
                }
            }

            //- store the information
            firstLayerThickness[seI] = thickness;
            thicknessRatio[seI] = ratio;
            nLayersAtEdge[seI] = nLayers;

            if( !specialMode_ )
            {
                nNodesAtEdge[seI] = nLayers + 1;
            }
            else
            {
                nNodesAtEdge[seI] = 3;
            }
        }
    }

    if( Pstream::parRun() )
    {
        //- transfer the information over all processor for edges
        //- at inter-processor boundaries
        const polyMeshGenAddressing& addr = mesh_.addressingData();
        const labelLongList& globalPointLabel = addr.globalPointLabel();
        const Map<label>& globalToLocal = addr.globalToLocalPointAddressing();
        const VRWGraph& pAtProcs = addr.pointAtProcs();
        const DynList<label>& neiProcs = addr.pointNeiProcs();

        std::map<label, LongList<splitEdgeAttributes> > exchangeData;
        forAll(neiProcs, i)
            exchangeData[neiProcs[i]].clear();

        forAll(splitEdges_, seI)
        {
            const edge& se = splitEdges_[seI];

            //- find split edges located at inter-processor boundaries
            if
            (
                globalToLocal.found(globalPointLabel[se.start()]) &&
                globalToLocal.found(globalPointLabel[se.end()])
            )
            {
                const label sp = se.start();

                const splitEdgeAttributes attr
                (
                    globalPointLabel[sp],
                    globalPointLabel[se.end()],
                    nLayersAtEdge[seI],
                    nNodesAtEdge[seI],
                    firstLayerThickness[seI],
                    thicknessRatio[seI]
                );

                forAllRow(pAtProcs, sp, i)
                {
                    const label neiProc = pAtProcs(sp, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append(attr);
                }
            }
        }

        //- exchange data between processors
        LongList<splitEdgeAttributes> receivedData;
        help::exchangeMap(exchangeData, receivedData);

        # ifdef USE_OMP
        # pragma omp parallel for if( receivedData.size() > 1000 ) \
        schedule(static, 1)
        # endif
        forAll(receivedData, i)
        {
            const splitEdgeAttributes& attr = receivedData[i];

            const label sp = globalToLocal[attr.start()];
            const label ep = globalToLocal[attr.end()];

            //- find the local edge index
            label seI(-1);
            forAllRow(splitEdgesAtPoint_, sp, j)
            {
                const label seJ = splitEdgesAtPoint_(sp, j);

                if( splitEdges_[seJ].end() == ep )
                {
                    seI = seJ;
                    break;
                }
            }

            if( seI < 0 )
            {
                FatalErrorIn
                (
                    "void refineBoundaryLayers::generateNewVertices()"
                ) << "Could not recognise split edge "
                  << attr << abort(FatalError);
            }

            //- reduction of parameters
            nNodesAtEdge[seI] = max(nNodesAtEdge[seI], attr.nNodesAtEdge());
            nLayersAtEdge[seI] = max(nLayersAtEdge[seI], attr.nSplits());
            thicknessRatio[seI] =
                max(thicknessRatio[seI], attr.thicknessRatio());
            firstLayerThickness[seI] =
                min(firstLayerThickness[seI], attr.firstLayerThickness());
        }
    }

    //- calculate the number of additional vertices which will be generated
    //- on edges of the mesh
    DynList<label> numPointsAtThread;
    numPointsAtThread.setSize(nThreads);
    numPointsAtThread = 0;

    # ifdef USE_OMP
    # pragma omp parallel for schedule(static, 1)
    # endif
    forAll(nNodesAtEdge, seI)
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif

        numPointsAtThread[threadI] += nNodesAtEdge[seI] - 2;
    }

    //- allocate the space in a graph storing ids of points on a split edge
    newVerticesForSplitEdge_.setSizeAndRowSize(nNodesAtEdge);

    //- calculate the number of points which will be generated
    //- on split edges
    label numPoints = points.size();
    forAll(numPointsAtThread, threadI)
    {
        const label nPts = numPointsAtThread[threadI];
        numPointsAtThread[threadI] = numPoints;
        numPoints += nPts;
    }

    polyMeshGenModifier meshModifier(mesh_);
    meshModifier.pointsAccess().setSize(numPoints);

    # ifdef DEBUGLayer
    Info << "Generating split vertices" << endl;
    # endif

    //- generate vertices on split edges
    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        const label threadI = omp_get_thread_num();
        # else
        const label threadI(0);
        # endif

        label& nPoints = numPointsAtThread[threadI];

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(splitEdges_, seI)
        {
            const edge& e = splitEdges_[seI];

            const vector v = e.vec(points);
            const scalar magv = mag(v);

            const label nLayers = newVerticesForSplitEdge_.sizeOfRow(seI) - 1;

            scalar firstThickness = magv / nLayersAtEdge[seI];
            if( thicknessRatio[seI] > (1. + SMALL) )
            {
                firstThickness =
                    magv /
                    (
                        (1 - Foam::pow(thicknessRatio[seI], nLayersAtEdge[seI]))
                        / (1.0 - thicknessRatio[seI])
                    );

                # ifdef DEBUGLayer
                Pout << "Thread " << threadI << endl;
                Pout << "Generating vertices at split edge "
                     << " start point " << points[e.start()]
                     << " end point " << points[e.end()] << endl;
                Pout << "Edge length " << magv << endl;
                Pout << "Thickness of the first layer "
                     << firstThickness << endl;
                # endif
            }

            firstThickness =
                Foam::min
                (
                    Foam::max(firstLayerThickness[seI], SMALL),
                    firstThickness
                );

            if( specialMode_ )
            {
                scalar t = firstThickness;

                for(label i=1;i<nLayersAtEdge[seI]-1;++i)
                    t += firstThickness * Foam::pow(thicknessRatio[seI], i);

                firstThickness = t;
            }

            //- generate vertices for this edge
            newVerticesForSplitEdge_(seI, 0) = e.start();

            scalar param = firstThickness;
            const vector vec = v / (magv + VSMALL);

            for(label pI=1;pI<nLayers;++pI)
            {
                //- generate the new vertex
                const point newP = points[e.start()] + param * vec;

                # ifdef DEBUGLayer
                Pout << "Split edge " << seI << " edge points " << e
                    << " start point " << points[e.start()]
                    << " end point " << points[e.end()]
                    << " param " << param
                    << " new point " << nPoints
                    << " has coordinates " << newP << endl;
                # endif

                param += firstThickness * Foam::pow(thicknessRatio[seI], pI);

                newVerticesForSplitEdge_(seI, pI) = nPoints;
                meshModifier.movePoint(nPoints++, newP);
            }

            newVerticesForSplitEdge_(seI, nLayers) = e.end();
        }
    }

    if( specialMode_ )
    {
        //- set the number of layers to 2
        forAll(nLayersAtBndFace_, bfI)
            if( nLayersAtBndFace_[bfI] > 1 )
                nLayersAtBndFace_[bfI] = 2;
    }

    # ifdef DEBUGLayer
    for(label procI=0;procI<Pstream::nProcs();++procI)
    {
        if( procI == Pstream::myProcNo() )
        {
            forAll(splitEdges_, seI)
            {
                Pout << "\nSplit edge " << seI << " nodes " << splitEdges_[seI]
                    << " coordinates " << points[splitEdges_[seI][0]]
                    << " " << points[splitEdges_[seI][1]]
                    << " has new points "
                    << newVerticesForSplitEdge_[seI] << endl;

                forAllRow(newVerticesForSplitEdge_, seI, i)
                    Pout << "Point " << i << " on edge ha coordinates "
                         << points[newVerticesForSplitEdge_(seI, i)] << endl;
            }
        }

        returnReduce(1, sumOp<label>());
    }

    Info << "Finished generating vertices at split edges" << endl;
    //::exit(1);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
