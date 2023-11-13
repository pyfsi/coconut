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

#include "cartesianMeshGenerator.H"
#include "triSurf.H"
#include "triSurfacePatchManipulator.H"
#include "demandDrivenData.H"
#include "meshOctreeCreator.H"
#include "cartesianMeshExtractor.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceMapper.H"
#include "edgeExtractor.H"
#include "meshSurfaceEdgeExtractorNonTopo.H"
#include "meshOptimizer.H"
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

#include "meshSurfaceCheckInvertedVertices.H"
#include "polyMeshGenChecks.H"

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void cartesianMeshGenerator::createCartesianMesh()
{
    //- create polyMesh from octree boxes
    cartesianMeshExtractor cme(*octreePtr_, meshDict_, mesh_);

    if( meshDict_.found("decomposePolyhedraIntoTetsAndPyrs") )
    {
        if( readBool(meshDict_.lookup("decomposePolyhedraIntoTetsAndPyrs")) )
            cme.decomposeSplitHexes();
    }

    cme.createMesh();
}

void cartesianMeshGenerator::surfacePreparation()
{
    //- removes unnecessary cells and morph the boundary
    //- such that there is only one boundary face per cell
    //- It also checks topology of cells after morphing is performed
    bool changed;

    do
    {
        changed = false;

        //- check and repair non-manifold connections at the surface of the mesh
        checkIrregularSurfaceConnections checkConnections(mesh_);
        if( checkConnections.checkAndFixIrregularConnections() )
            changed = true;

        //- check if all cell cells can be visited over their neighbours
        //- removes cells until there is only one group remaining
        if( checkCellConnectionsOverFaces(mesh_).checkCellGroups() )
            changed = true;
    } while( changed );

    //- sometimes boundary faces share two common edges
    //- splits such faces into triangles
    checkBoundaryFacesSharingTwoEdges(mesh_).improveTopology();

    if
    (
        (returnReduce(mesh_.cells().size(), sumOp<label>()) == 0) ||
        (returnReduce(mesh_.faces().size(), sumOp<label>()) == 0)
    )
    {
        throw std::string("The mesh has no cells and faces!!");
    }
}

void cartesianMeshGenerator::mapMeshToSurface()
{
    //- calculate mesh surface
    meshSurfaceEngine mse(mesh_);

    //- backup surface vertices
    {
        polyMeshGenModifier meshModifier(mesh_);
        const labelLongList& bPoints = mse.boundaryPoints();

        forAll(bPoints, bpI)
            meshModifier.backupPoint(bPoints[bpI]);
    }

    //- pre-map mesh surface
    meshSurfaceMapper mapper(mse, *octreePtr_);
    mapper.preMapVertices();

    //- untangle surface faces
    meshSurfaceOptimizer(mse, *octreePtr_).optimizeSurface();

    //- map mesh surface on the geometry surface
    mapper.mapVerticesOntoSurface();

    //- untangle surface faces
    meshSurfaceOptimizer(mse, *octreePtr_).optimizeSurface();
}

void cartesianMeshGenerator::extractPatches()
{
    edgeExtractor extractor(mesh_, *octreePtr_);

    Info << "Extracting edges" << endl;
    extractor.extractEdges();

    extractor.updateMeshPatches();
}

void cartesianMeshGenerator::mapEdgesAndCorners()
{
    meshSurfaceEdgeExtractorNonTopo(mesh_, *octreePtr_);

    optimiseMeshSurface();
}

void cartesianMeshGenerator::optimiseMeshSurface()
{
    meshSurfaceEngine mse(mesh_);
    meshSurfaceOptimizer surfOpt(mse, *octreePtr_);
    surfOpt.optimizeSurface();
}

void cartesianMeshGenerator::generateBoundaryLayers()
{
    //- add boundary layers
    boundaryLayers bl(mesh_);
    bl.addLayerForAllPatches();
}

void cartesianMeshGenerator::refBoundaryLayers()
{
    if( meshDict_.isDict("boundaryLayers") )
    {
        refineBoundaryLayers refLayers(mesh_);

        refineBoundaryLayers::readSettings(meshDict_, refLayers);

        refLayers.refineLayers();

        labelLongList pointsInLayer;
        refLayers.pointsInBndLayer(pointsInLayer);

        meshOptimizer mOpt(mesh_);
        mOpt.lockPoints(pointsInLayer);
        mOpt.untangleBoundaryLayer();
    }
}

void cartesianMeshGenerator::optimiseFinalMesh()
{
    //- untangle the surface if needed
    bool enforceConstraints(false);
    if( meshDict_.found("enforceGeometryConstraints") )
    {
        enforceConstraints =
            readBool(meshDict_.lookup("enforceGeometryConstraints"));
    }

    if( true )
    {
        meshSurfaceEngine mse(mesh_);
        meshSurfaceOptimizer surfOpt(mse, *octreePtr_);

        if( enforceConstraints )
            surfOpt.enforceConstraints();

        surfOpt.optimizeSurface();
    }

    deleteDemandDrivenData(octreePtr_);

    //- final optimisation
    meshOptimizer optimizer(mesh_);
    if( enforceConstraints )
        optimizer.enforceConstraints();

    optimizer.optimizeMeshFVWithTopoChanges();

    //- check if the mesh has been untangled
    labelHashSet badFaces;
    if
    (
        returnReduce
        (
            polyMeshGenChecks::findBadFaces(mesh_, badFaces),
            sumOp<label>()
        )
    )
    {
        Warning << "The mesh cannot be untangled. Please improve your settings"
          << " and try again. Exitting..." << endl;

        const label badId = mesh_.addFaceSubset("invalidFaces");
        const label badCellsId = mesh_.addCellSubset("invalidCells");
        const labelLongList& owner = mesh_.owner();
        const labelLongList& neighbour = mesh_.neighbour();

        forAllConstIter(labelHashSet, badFaces, it)
        {
            mesh_.addFaceToSubset(badId, it.key());

            mesh_.addCellToSubset(badCellsId, owner[it.key()]);

            if( neighbour[it.key()] < 0 )
                continue;

            mesh_.addCellToSubset(badCellsId, neighbour[it.key()]);
        }

        //- write the mesh and stop the meshing process
        mesh_.write();

        const std::string message("Mesh cannot be untangled. IT IS NOT VALID!");
        throw message;
    }
    badFaces.clear();

    //- get rid of low quality cells
    optimizer.optimizeLowQualityFaces();

    if( modSurfacePtr_ )
    {
        //- revert the mesh to the original coordinates
        //- when anisotropic meshing is used
        polyMeshGenGeometryModification meshMod(mesh_, meshDict_);

        //- revert the mesh into the original space
        meshMod.revertGeometryModification();

        //- delete modified surface mesh
        deleteDemandDrivenData(modSurfacePtr_);

        //- project the mesh back to the surface mesh
        deleteDemandDrivenData(octreePtr_);
        octreePtr_ = new meshOctree(*surfacePtr_);

        meshOctreeCreator
        (
            *octreePtr_,
            meshDict_
        ).createOctreeWithRefinedBoundary(20, 30);

        //- calculate mesh surface
        meshSurfaceEngine mse(mesh_);

        //- pre-map mesh surface
        meshSurfaceMapper mapper(mse, *octreePtr_);

        //- map mesh surface on the geometry surface
        mapper.mapVerticesOntoSurface();

        //- clear out addressing
        mesh_.clearAddressingData();
    }

    optimizer.untangleMeshFVWithTopoChanges();
    optimizer.optimizeBoundaryLayer(true);
    optimizer.untangleMeshFVWithTopoChanges();
}

void cartesianMeshGenerator::replaceBoundaries()
{
    renameBoundaryPatches rbp(mesh_, meshDict_);
}

void cartesianMeshGenerator::renumberMesh()
{
    polyMeshGenModifier(mesh_).renumberMesh();
}

void cartesianMeshGenerator::finalCleanup()
{
    polyMeshGenModifier meshModifier(mesh_);

    meshModifier.clearPointBackup();
    meshModifier.unlockPoints();

    meshModifier.clearAll();
}

void cartesianMeshGenerator::generateMesh()
{
    if( controller_.runCurrentStep("templateGeneration") )
    {
        createCartesianMesh();
    }

    if( controller_.runCurrentStep("surfaceTopology") )
    {
        surfacePreparation();
    }

    if( controller_.runCurrentStep("surfaceProjection") )
    {
        mapMeshToSurface();
    }

    if( controller_.runCurrentStep("patchAssignment") )
    {
        extractPatches();
    }

    if( controller_.runCurrentStep("edgeExtraction") )
    {
        mapEdgesAndCorners();
    }

    if( controller_.runCurrentStep("boundaryLayerGeneration") )
    {
        generateBoundaryLayers();
    }

    if( controller_.runCurrentStep("meshOptimisation") )
    {
        optimiseFinalMesh();
    }

    if( controller_.runCurrentStep("boundaryLayerRefinement") )
    {
        refBoundaryLayers();
    }

    renumberMesh();

    replaceBoundaries();

    finalCleanup();

    controller_.workflowCompleted();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cartesianMeshGenerator::cartesianMeshGenerator(const Time& time)
:
    db_(time),
    surfacePtr_(NULL),
    modSurfacePtr_(NULL),
    meshDict_
    (
        IOobject
        (
            "meshDict",
            db_.system(),
            db_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    octreePtr_(NULL),
    mesh_(time),
    controller_(mesh_)
{
    try
    {
        if( true )
        {
            checkMeshDict cmd(meshDict_);
        }

        fileName surfaceFile = meshDict_.lookup("surfaceFile");
        if( Pstream::parRun() )
            surfaceFile = ".."/surfaceFile;

        surfacePtr_ = new triSurf(db_.path()/surfaceFile);

        if( true )
        {
            //- save meta data with the mesh (surface mesh + its topology info)
            triSurfaceMetaData sMetaData(*surfacePtr_);
            const dictionary& surfMetaDict = sMetaData.metaData();

            mesh_.metaData().remove("surfaceFile");
            mesh_.metaData().remove("surfaceMeta");
            mesh_.metaData().add("surfaceFile", surfaceFile, true);
            mesh_.metaData().add("surfaceMeta", surfMetaDict, true);
        }

        if( surfacePtr_->featureEdges().size() != 0 )
        {
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

            octreePtr_ = new meshOctree(*modSurfacePtr_);
        }
        else
        {
            octreePtr_ = new meshOctree(*surfacePtr_);
        }

        meshOctreeCreator(*octreePtr_, meshDict_).createOctreeBoxes();

        generateMesh();
    }
    catch(const std::string& message)
    {
        Info << message << endl;
    }
    catch(...)
    {
        Warning << "Meshing process terminated!" << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cartesianMeshGenerator::~cartesianMeshGenerator()
{
    deleteDemandDrivenData(surfacePtr_);
    deleteDemandDrivenData(modSurfacePtr_);
    deleteDemandDrivenData(octreePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void cartesianMeshGenerator::writeMesh() const
{
    mesh_.write();
}

const polyMeshGen& cartesianMeshGenerator::mesh() const
{
    return mesh_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
