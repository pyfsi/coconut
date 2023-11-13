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
#include "triSurfaceCleanupDuplicateTriangles.H"
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
#include "triSurfModifier.H"
#include "meshSurfaceDistanceFromGeometry.H"
#include "meshOctreeModifier.H"
#include "meshOctreeInsideOutside.H"
#include "revolve2DMesh.H"
#include "extrude2DMesh.H"
#include "wireBlockMeshGenerator.H"

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rollingMillMesh::rollingMillMesh
(
    const Time& time,
    const fileName regionName,
    const fileName timeStep
)
:
    db_(time),
    regionName_(regionName),
    timeStep_(timeStep.empty()?"constant":timeStep),
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
    globalMeshPtr_(NULL),
    patchHandler_(meshDict_),
    geomHandler_(time, meshDict_, patchHandler_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rollingMillMesh::~rollingMillMesh()
{
    deleteDemandDrivenData(globalMeshPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rollingMillMesh::generateWire()
{
    try
    {
        //- read the wire properties from disk
        geomHandler_.parseWireDictionary(meshDict_);

        //- performs meshing of wire
        generateWireMesh();

        //- performs meshing of wire coating
        generateWireCoatingMesh();

        //- renumber the mesh
        polyMeshGenModifier(*globalMeshPtr_).renumberMesh();
    }
    catch(const std::string& message)
    {
        Info << message << endl;
    }
    catch(...)
    {
        WarningIn
        (
            "void rollingMillMesh::generateWire()"
        ) << "Meshing process terminated!" << endl;
    }
}

void rollingMillMesh::generateRollers()
{
    try
    {
        if( !geomHandler_.isRollSetupValid() )
        {
            FatalErrorIn
            (
                "void rollingMillMesh::generateRollers()"
            ) << "rollSetup is set to none. Cannot generate roller meshes."
              << exit(FatalError);
        }

        //- performs meshing of rolls
        generateRollerMeshes();

        //- renumber the mesh
        polyMeshGenModifier(*globalMeshPtr_).renumberMesh();
    }
    catch(const std::string& message)
    {
        Info << message << endl;
    }
    catch(...)
    {
        WarningIn
        (
            "void rollingMillMesh::generateRollers()"
        ) << "Meshing process terminated!" << endl;
    }
}

void rollingMillMesh::generateDie()
{
    try
    {
        if( geomHandler_.isRollSetupValid() )
        {
            FatalErrorIn
            (
                "void rollingMillMesh::generateRollers()"
            ) << "rollSetup is not none. Cannot generate a die mesh."
              << exit(FatalError);
        }

        //- read the die properties from disk
        geomHandler_.parseDieDictionary(meshDict_);

        //- performs meshing of a die
        generateDieMesh();

        //- performs meshing of a casing
        generateCasingMesh();

        //- renumber the mesh
        polyMeshGenModifier(*globalMeshPtr_).renumberMesh();
    }
    catch(const std::string& message)
    {
        Info << message << endl;
    }
    catch(...)
    {
        WarningIn
        (
            "void rollingMillMesh::generateDie()"
        ) << "Meshing process terminated!" << endl;
    }
}

void rollingMillMesh::writeMesh() const
{
    if( !globalMeshPtr_ )
    {
        Warning << "Mesh was not generated!!" << endl;
        return;
    }

    globalMeshPtr_->write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
