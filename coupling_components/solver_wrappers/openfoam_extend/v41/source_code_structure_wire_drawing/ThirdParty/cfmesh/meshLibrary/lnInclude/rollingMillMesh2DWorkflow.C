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

#include "rollingMillMesh.H"
#include "triSurface2DCheck.H"
#include "polyMeshGen2DEngine.H"
#include "triSurf.H"
#include "triSurfacePatchManipulator.H"
#include "triSurfaceCleanupDuplicates.H"
#include "demandDrivenData.H"
#include "meshOctreeCreator.H"
#include "cartesianMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper2D.H"
#include "meshSurfaceEdgeExtractor2D.H"
#include "meshSurfaceOptimizer.H"
#include "topologicalCleaner.H"
#include "boundaryLayers.H"
#include "refineBoundaryLayers.H"
#include "renameBoundaryPatches.H"
#include "checkMeshDict.H"
#include "checkCellConnectionsOverFaces.H"
#include "checkIrregularSurfaceConnections.H"
#include "checkNonMappableCellConnections.H"
#include "checkBoundaryFacesSharingTwoEdges.H"
#include "triSurfaceMetaData.H"
#include "polyMeshGenGeometryModification.H"
#include "surfaceMeshGeometryModification.H"

//- include geometry creator developed for Bekaert
#include "rollerSurfaceCreator.H"
#include "meshSurfaceDistanceFromGeometry.H"
#include "meshOctreeModifier.H"
#include "meshOctreeInsideOutside.H"

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void rollingMillMesh::crossSectionMeshGenerator::cleanupSmallPatches()
{
    meshOctree octree(*surfacePtr_, true);
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(15, 30);


    triSurfaceCleanupDuplicates cleaner(octree);
    cleaner.mergeIdentities();
}

void rollingMillMesh::crossSectionMeshGenerator::createOctree()
{
    if( true )
    {
        triSurface2DCheck surfCheck(*surfacePtr_);
        if( !surfCheck.is2DSurface() )
        {
            surfCheck.createSubsets();

            Info << "Writting surface with subsets to file "
                 << "badSurfaceWithSubsets.fms" << endl;
            surfacePtr_->writeSurface("badSurfaceWithSubsets.fms");

            FatalErrorIn
            (
                "void rollingMillMesh::crossSectionMeshGenerator"
                "::createOctree()"
            ) << "Cannot generate a 2D mesh. Geometry is not a 2D surface"
              << exit(FatalError);
        }
    }

    if( surfacePtr_->featureEdges().size() != 0 )
    {
        //- get rid of duplicate triangles as they cause strange problems
//        triSurfaceCleanupDuplicateTriangles
//        (
//            const_cast<triSurf&>(*surfacePtr_)
//        );

        //- create surface patches based on the feature edges
        //- and update the meshDict based on the given data
        triSurfacePatchManipulator manipulator(*surfacePtr_);

        const triSurf* surfaceWithPatches =
            manipulator.surfaceWithPatches(&meshDict_);

        //- delete the old surface and assign the new one
        deleteDemandDrivenData(surfacePtr_);
        surfacePtr_ = surfaceWithPatches;
    }

    if( meshDict_.found("anisotropicSources") )
    {
        surfaceMeshGeometryModification surfMod(*surfacePtr_, meshDict_);

        modSurfacePtr_ = surfMod.modifyGeometry();

        meshOctreePtr_ = new meshOctree(*modSurfacePtr_, true);
    }
    else
    {
        meshOctreePtr_ = new meshOctree(*surfacePtr_, true);
    }

    meshOctreeCreator(*meshOctreePtr_, meshDict_).createOctreeBoxes();

    if( userSettingsDict_.found("contactCellSize") )
    {
        const scalar cSize =
            readScalar(userSettingsDict_.lookup("contactCellSize"));

        const boundBox& rootBox = meshOctreePtr_->rootBox();
        const scalar s = rootBox.max().x() - rootBox.min().x();

        scalar minDist(VGREAT);
        contactRefLevel_ = 0;

        label currLevel(0);

        bool improve;
        do
        {
            improve = false;

            const scalar d = mag(cSize - s / pow(2.0, currLevel));
            if( d < minDist )
            {
                minDist = d;
                contactRefLevel_ = currLevel;

                improve = true;
            }

            ++currLevel;
        } while( improve );

        maxContactRefLevel_ = contactRefLevel_;

        Info << "Contact ref level " << label(contactRefLevel_) << endl;

        enforceContactCellSize_ = true;
    }

    if( userSettingsDict_.found("minContactCellSize") )
    {
        const scalar cSize =
            readScalar(userSettingsDict_.lookup("minContactCellSize"));

        const boundBox& rootBox = meshOctreePtr_->rootBox();
        const scalar s = rootBox.max().x() - rootBox.min().x();

        scalar minDist(VGREAT);
        maxContactRefLevel_ = 0;

        label currLevel(0);

        bool improve;
        do
        {
            improve = false;

            const scalar d = mag(cSize - s / pow(2.0, currLevel));
            if( d < minDist )
            {
                minDist = d;
                maxContactRefLevel_ = currLevel;

                improve = true;
            }

            ++currLevel;
        } while( improve );

        maxContactRefLevel_ = max(maxContactRefLevel_, contactRefLevel_);

        Info << "Max contact ref level " << label(maxContactRefLevel_) << endl;
    }

    //- detect if there exist any empty surface patches
    HashSet<const word> emptySurfacePatches;

    const triSurf& surf = meshOctreePtr_->surface();
    labelList nTriasInPatch(surf.patches().size(), 0);
    forAll(surf, tI)
        ++nTriasInPatch[surf[tI].region()];

    forAll(nTriasInPatch, patchI)
    {
        if( nTriasInPatch[patchI] == 0 )
            emptySurfacePatches.insert(surf.patches()[patchI].name());
    }

    if( emptySurfacePatches.size() != 0 )
    {
        surf.writeSurface("surfWithEmptyPatches.stl");

        FatalError << "Geometry definition error!! "
            << "Patches with zero faces detected. This causes infinite loops."
            << " Patch names are " << emptySurfacePatches << exit(FatalError);
    }
}

void rollingMillMesh::crossSectionMeshGenerator::parseUserDictionary()
{
    if( userSettingsDict_.found("minNumFacesBetweenFeatures") )
    {
        minNumFacesInPatch_ =
            readLabel(userSettingsDict_.lookup("minNumFacesBetweenFeatures"));
    }

    if( userSettingsDict_.found("nAdditionalLayers") )
    {
        nAdditionalLayers_ =
            readLabel(userSettingsDict_.lookup("nAdditionalLayers"));
    }

    if( userSettingsDict_.found("maxCellSize") )
    {
        const scalar cs = readScalar(userSettingsDict_.lookup("maxCellSize"));

        meshDict_.add("maxCellSize", cs, true);
    }
    else
    {
        FatalError << "Cannot find a mandatory setting"
            << " maxCellSize in dictionary "
            << userSettingsDict_ << exit(FatalError);
    }
}

void rollingMillMesh::crossSectionMeshGenerator::createCartesianMesh()
{
    octreePtr_ = meshOctreePtr_;

    deleteDemandDrivenData(meshPtr_);
    meshPtr_ = new polyMeshGen(runTime_, timeStep_, regionName_/"polyMesh");

    //- create polyMesh from octree boxes
    try
    {
        cartesianMeshExtractor cme(*meshOctreePtr_, meshDict_, *meshPtr_);

        cme.createMesh();
    }
    catch(...)
    {}
}

bool rollingMillMesh::crossSectionMeshGenerator::surfacePreparation()
{
    //- removes unnecessary cells and morph the boundary
    //- such that there is only one boundary face per cell
    //- It also checks topology of cells after morphing is performed
    bool changed;

    do
    {
        changed = false;

        if( meshPtr_->cells().size() == 0 )
            return false;

        checkIrregularSurfaceConnections checkConnections(*meshPtr_);
        if( checkConnections.checkAndFixIrregularConnections() )
            changed = true;

        if( checkNonMappableCellConnections(*meshPtr_).removeCells() )
            changed = true;

        if( checkCellConnectionsOverFaces(*meshPtr_).checkCellGroups() )
            changed = true;
    } while( changed );

    checkBoundaryFacesSharingTwoEdges(*meshPtr_).improveTopology();

    if( meshPtr_->cells().size() == 0 )
        return false;

    return true;
}

void rollingMillMesh::crossSectionMeshGenerator::mapMeshToSurface()
{
    //- calculate mesh surface
    meshSurfaceEngine* msePtr = new meshSurfaceEngine(*meshPtr_);

    //- pre-map mesh surface
    meshSurfaceMapper2D mapper(*msePtr, *octreePtr_);

    mapper.adjustZCoordinates();

    mapper.preMapVertices();

    //- map mesh surface on the geometry surface
    mapper.mapVerticesOntoSurface();

    deleteDemandDrivenData(msePtr);
}

void rollingMillMesh::crossSectionMeshGenerator::extractPatches()
{
    meshSurfaceEdgeExtractor2D
    (
        *meshPtr_,
        *octreePtr_
    ).distributeBoundaryFaces();
}

void rollingMillMesh::crossSectionMeshGenerator::mapEdgesAndCorners()
{
    meshSurfaceEdgeExtractor2D(*meshPtr_, *octreePtr_).remapBoundaryPoints();
}

void rollingMillMesh::crossSectionMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(*meshPtr_);
    meshSurfaceOptimizer optimizer(mse, *octreePtr_);
    optimizer.optimizeSurface2D();
    optimizer.optimizeLowQualitySurface2D();
    optimizer.untangleSurface2D();
}

void rollingMillMesh::crossSectionMeshGenerator::generateBoundaryLayers()
{
    boundaryLayers bl(*meshPtr_);

    bl.activate2DMode();

    bl.addLayerForAllPatches();

    if( modSurfacePtr_ )
    {
        polyMeshGenGeometryModification meshMod(*meshPtr_, meshDict_);

        //- revert the mesh into the original space
        meshMod.revertGeometryModification();

        //- contruct a new octree from the input surface
        octreePtr_ = new meshOctree(*surfacePtr_, true);
        meshOctreeCreator(*octreePtr_).createOctreeWithRefinedBoundary(20);

        mapEdgesAndCorners();

        optimiseMeshSurface();
    }
}

void rollingMillMesh::crossSectionMeshGenerator::refBoundaryLayers()
{
    if( meshDict_.isDict("boundaryLayers") )
    {
        refineBoundaryLayers refLayers(*meshPtr_);

        refineBoundaryLayers::readSettings(meshDict_, refLayers);

        refLayers.activate2DMode();

        refLayers.refineLayers();

        meshSurfaceEngine mse(*meshPtr_);
        meshSurfaceOptimizer optimizer(mse, *octreePtr_);

        optimizer.untangleSurface2D();
    }
}

void rollingMillMesh::crossSectionMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(*meshPtr_, meshDict_, true);
}

void rollingMillMesh::crossSectionMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(*meshPtr_).renumberMesh();
}

bool rollingMillMesh::crossSectionMeshGenerator::checkGeometryDeviation()
{
    Info << "Checking geometry deviation" << endl;

    //- calculate deviation of face centres from the original geometry
    meshSurfaceEngine mse(*meshPtr_);
    const faceList::subList& bFaces = mse.boundaryFaces();
    const vectorLongList& fCentres = mse.faceCentres();
    const labelLongList& facePatch = mse.boundaryFacePatches();

    const pointFieldPMG& points = mse.mesh().points();
    const PtrList<boundaryPatch>& patches = mse.mesh().boundaries();

    const triSurf& surf = octreePtr_->surface();
    const geometricSurfacePatchList& sPatches = surf.patches();

    //- check if the contact patch exists in the surface mesh
    bool hasContactInGeometry(false);

    //- indices of contact patches
    labelHashSet contactIds;
    forAllConstIter(HashSet<const word>, contactPatchNames_, cIt)
    {
        forAll(sPatches, patchI)
        {
            if( sPatches[patchI].name().find(cIt.key()) != word::npos )
                hasContactInGeometry = true;;
        }

        forAll(patches, patchI)
        {
            if( patches[patchI].patchName().find(cIt.key()) != word::npos )
                contactIds.insert(patchI);
        }
    }

    if( !hasContactInGeometry )
    {
        FatalError << "No contact patches available in the input geometry."
            << " Available patches are " << sPatches << abort(FatalError);
    }

    //- find patches with insufficient number of faces
    std::set<label> coarsePatches;
    forAll(patches, patchI)
    {
        const boundaryPatch& patch = patches[patchI];

        if( patch.patchSize() < minNumFacesInPatch_ )
        {
            Info << "Patch " << patch.patchName() << " has less than "
                 << minNumFacesInPatch_ << " faces" << endl;

            coarsePatches.insert(patchI);
        }
    }

    //- calculate distances between the geometry and the mesh
    bool modified(false);

    //- find the octree leaves than need to be refined
    labelList refineBox(meshOctreePtr_->numberOfLeaves());

    boolList hasMeshInVicinity(surf.size()), detected(surf.size());

    scalar maxDist(0.0);

    DynList<label> leavesInRange, facets;

    # ifdef USE_OMP
    # pragma omp parallel private(leavesInRange, facets)
    # endif
    {
        scalar dMax(0.0);

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1) nowait
        # endif
        forAll(refineBox, leafI)
            refineBox[leafI] = 0;

        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(hasMeshInVicinity, tI)
        {
            const label patchI = surf[tI].region();

            if
            (
                !contactIds.found(patchI) &&
                (coarsePatches.find(patchI) != coarsePatches.end())
            )
            {
                //- refine coarse patches that are not in the contact region
                hasMeshInVicinity[tI] = false;
            }
            else
            {
                hasMeshInVicinity[tI] = true;
            }

            //- initialise detected to false
            detected[tI] = false;
        }

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 50)
        # endif
        forAll(bFaces, bfI)
        {
            const face& bf = bFaces[bfI];

            //- sometimes the mesh has zero cells and all surface
            //- patches do not exist in the mesh
            if( patches.size() <= sPatches.size() )
                continue;

            const word& pName = patches[facePatch[bfI]].patchName();

            //- find surface facets near the current boundary face
            const scalar r = mag(points[bf[0]] - fCentres[bfI]);

            leavesInRange.clear();
            octreePtr_->findLeavesInSphere(fCentres[bfI], r, leavesInRange);

            std::set<label> nearFacets;
            forAll(leavesInRange, i)
            {
                facets.clear();
                octreePtr_->containedTriangles(leavesInRange[i], facets);

                forAll(facets, j)
                    nearFacets.insert(facets[j]);
            }

            //- find the nearest surface face to the face centre
            label nearest(-1);
            scalar dSq(VGREAT);

            forAllConstIter(std::set<label>, nearFacets, it)
            {
                const label tI = *it;

                //- filter out triangles that are out-of-range
                const point np =
                    help::nearestPointOnTheTriangle(tI, surf, fCentres[bfI]);

                const scalar distSq = magSqr(np - fCentres[bfI]);

                if( distSq > sqr(r) )
                    continue;

                //- mark all surface facets in that patch as marked
                if( sPatches[surf[tI].region()].name() == pName )
                {
                    detected[tI] = true;

                    if( !contactIds.found(surf[tI].region()) )
                        continue;

                    if( distSq < dSq )
                    {
                        nearest = tI;
                        dSq = distSq;
                    }
                }
            }

            if( nearest != -1 )
            {
                dMax = max(dMax, sqrt(dSq));

                if
                (
                    (dSq > sqr(geometryTolerance_)) ||
                    (
                        coarsePatches.find(surf[nearest].region()) !=
                        coarsePatches.end()
                    )
                )
                {
                    hasMeshInVicinity[nearest] = false;
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        for(label leafI=0;leafI<meshOctreePtr_->numberOfLeaves();++leafI)
        {
            const meshOctreeCubeBasic& oc = meshOctreePtr_->returnLeaf(leafI);

            //- check if the box contains the contact patch
            bool hasContact(false);
            facets.clear();
            meshOctreePtr_->containedTriangles(leafI, facets);

            forAll(facets, i)
            {
                const label tI = facets[i];

                if( contactIds.found(surf[tI].region()) )
                {
                    hasContact = true;
                    break;
                }
            }

            //- check if the cell size is larger than requested
            if( hasContact && (oc.level() < contactRefLevel_) )
            {
                modified = true;
                refineBox[leafI] = 1;
                continue;
            }

            //- do not allow the cell size smaller than the minContactCellSize
            if( hasContact && (oc.level() >= maxContactRefLevel_)  )
                continue;

            forAll(facets, i)
            {
                if( !hasMeshInVicinity[facets[i]] || !detected[facets[i]] )
                {
                    modified = true;
                    refineBox[leafI] = 1;
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp critical(maxDist)
        # endif
        {
            maxDist = max(maxDist, dMax);
        }
    }

    Info << "Max deviation from geometry " << maxDist << endl;

    if( modified )
    {
        Info << "Mesh deviates more than the prescribed tolerance" << endl;

        if( enforceUniformRefinement_ )
        {
            for(label leafI=0;leafI<octreePtr_->numberOfLeaves();++leafI)
            {
                if( octreePtr_->hasContainedTriangles(leafI) )
                {
                    DynList<label> ct;
                    octreePtr_->containedTriangles(leafI, ct);

                    bool inContact(false);
                    forAll(ct, i)
                    {
                        if( contactIds.found(surf[ct[i]].region()) )
                        {
                            inContact = true;
                            break;
                        }
                    }

                    if( inContact )
                    {
                        //- leaf has at least one triangle in the contact patch
                        refineBox[leafI] = 1;
                    }
                }
            }
        }

        //- refine the octree in the locations where the mesh is too coarse
        meshOctreeModifier octreeModifier(*meshOctreePtr_);
        octreeModifier.markAdditionalLayers(refineBox, nAdditionalLayers_);
        octreeModifier.refineSelectedBoxes(refineBox);

        //- update inside/outside information
        meshOctreeInsideOutside insideOutsideCalculator(*meshOctreePtr_);
    }
    else
    {
        Info << "Mesh is within the tolerance" << endl;
    }

    if( octreePtr_ != meshOctreePtr_ )
    {
        deleteDemandDrivenData(octreePtr_);
    }
    else
    {
        octreePtr_ = NULL;
    }

    # ifdef DEBUG
    Info << "Checking octree ref levels " << endl;

    for(label leafI=0;leafI<meshOctreePtr_->numberOfLeaves();++leafI)
    {
        if( meshOctreePtr_->returnLeaf(leafI).level() > maxContactRefLevel_ )
        {
            Info << "Contact ref level " << label(maxContactRefLevel_) << endl;
            FatalError << "Too fine leaf " << leafI << " coordinates "
                 << meshOctreePtr_->returnLeaf(leafI) << abort(FatalError);
        }
    }
    # endif

    return modified;
}

void rollingMillMesh::crossSectionMeshGenerator::generate2DMesh()
{
    try
    {
        createOctree();

        do
        {
            createCartesianMesh();

            //- check the surface of the mesh
            //- do not execute the rest of the workflow if there are no cells
            if( !surfacePreparation() )
                continue;

            mapMeshToSurface();

            extractPatches();

            mapEdgesAndCorners();

            optimiseMeshSurface();

            generateBoundaryLayers();

            optimiseMeshSurface();

            refBoundaryLayers();

            mapEdgesAndCorners();

        } while( checkGeometryDeviation() );

        replaceBoundaries();

        renumberMesh();
    }
    catch(...)
    {
        Warning << "Something unpredicted occured during meshing" << endl;
    }
}

rollingMillMesh::crossSectionMeshGenerator::crossSectionMeshGenerator
(
    const triSurf* surfPtr,
    const IOdictionary& meshDict,
    const dictionary& userDict,
    const Time& runTime,
    const fileName regionName,
    const fileName timeStep,
    const HashSet<const word>& contactPatchNames,
    const bool enforceUniformRefinement
)
:
    runTime_(runTime),
    regionName_(regionName),
    timeStep_(timeStep),
    contactPatchNames_(contactPatchNames),
    surfacePtr_(new triSurf(*surfPtr)),
    modSurfacePtr_(NULL),
    meshDict_(meshDict),
    userSettingsDict_(userDict),
    geometryTolerance_(1e-4),
    minNumFacesInPatch_(2),
    nAdditionalLayers_(0),
    meshPtr_(NULL),
    meshOctreePtr_(NULL),
    octreePtr_(NULL),
    contactRefLevel_(0),
    maxContactRefLevel_(255),
    enforceUniformRefinement_(enforceUniformRefinement),
    enforceContactCellSize_(false)
{
    if( meshDict.found("geometryTolerance") )
        geometryTolerance_ = readScalar(meshDict.lookup("geometryTolerance"));

    cleanupSmallPatches();

    parseUserDictionary();

    generate2DMesh();
}

rollingMillMesh::crossSectionMeshGenerator::~crossSectionMeshGenerator()
{
    deleteDemandDrivenData(meshOctreePtr_);
    deleteDemandDrivenData(octreePtr_);
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(modSurfacePtr_);
}

polyMeshGen* rollingMillMesh::crossSectionMeshGenerator::meshPtr() const
{
    if( !meshPtr_ )
        FatalError << "Volume mesh is not generated" << exit(FatalError);

    return meshPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
