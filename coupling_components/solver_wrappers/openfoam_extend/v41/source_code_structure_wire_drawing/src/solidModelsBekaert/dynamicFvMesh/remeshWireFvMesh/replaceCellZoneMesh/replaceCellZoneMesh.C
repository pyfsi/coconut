/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "replaceCellZoneMesh.H"
#include "newFvMeshSubset.H"
#include "volFields.H"
#include "meshTriangulation.H"
//#include "cartesianMeshGenerator.H"
#include "triSurf.H"
#include "triSurfaceDetectFeatureEdges.H"
#include "linearNormal.H"
#include "faceMesh.H"
#include "extrudedMesh.H"
#include "mergePolyMesh.H"
#include "perfectInterface.H"
#include "polyTopoChanger.H"
#include "patchManipulationFunctions.H"
#include "removeCells.H"
#include "distributedTriSurfaceMesh.H"
#include "runSystemCommandWithPOpen.H"
#include "nagataPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(replaceCellZoneMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::replaceCellZoneMesh::createOldGlobalMesh() const
{
    if (oldGlobalMeshPtr_.valid())
    {
        FatalErrorIn(type() + "::createOldGlobalMesh()")
            << "Pointer already set!" << abort(FatalError);
    }

    if (Pstream::master())
    {
        // Reconstruct the mesh from the latest time
        oldMesh_.write();
        chDir(globalTime().path());
        runSystemCommandWithPOpen
        (
            "reconstructParMeshZones -latestTime",
            "log.reconstructParMeshZones_remesh_oldGlobalMesh_time"
          + oldMesh_.time().timeName(),
            type(),
            true
        );
        chDir(oldMesh_.time().path());
    }

    // Force slave procs to wait for master to finish reconstructing
    // the mesh
    syncMPI();

    oldGlobalMeshPtr_.set
    (
        new fvMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                globalTime().timeName(),
                globalTime(),
                IOobject::MUST_READ
            )
        )
    );
}


void Foam::replaceCellZoneMesh::createNewCellZoneMesh() const
{
    if (newCellZoneMeshPtr_.valid())
    {
        FatalErrorIn(type() + "::createNewCellZoneMesh()")
            << "Pointer already set!" << abort(FatalError);
    }

    // Take references for convenience
    const Time& runTime = oldMesh_.time();


    // Approach:
    // Note: currently the procedure is specifically for a wire, but this can
    // be generalised without too much difficulty
    // 4. Create region-of-interest mesh (sub-set of the cellZone)
    // 5. Flatten the oldInternalFaces (downstream) patch
    // 6. Write out STL of entire region-of-interest boundary
    // 7. Create new mesh with cartesianMesh
    // 8. Extrude downstream patch to create extruded section of wire
    // 9. Merge extruded mesh with new wire mesh
    // 10. Stitch interface between extruded and remeshed wire regions


    // 4. Create region-of-interest mesh (sub-set of the wire)
    if (debug)
    {
        InfoIn(type() + "::createWireMesh() const")
            << nl << "4. Creating region-of-interest mesh (sub-set of the wire)"
            << endl;
    }

    // Create a sub mesh of part of the wire
    newFvMeshSubset subsetMesh
    (
        IOobject
        (
            "regionOfInterest",
            runTime.timeName(),
            oldMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        oldMesh_
    );

    // Define a list with 1s for all cells in the region-of-interest in the
    // cellZone and 0s for all the rest of cells
    // region-of-interest is a sub-set of the wire
    // We will use a volScalarField for parallel exchange
    labelList region(oldMesh_.nCells(), 0);
    volScalarField regionField
    (
        IOobject
        (
            "regionOfInterestField",
            runTime.timeName(),
            oldMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        oldMesh_,
        dimensionedScalar("zero", dimless, 0.0)
    );
    scalarField& regionFieldI = regionField.internalField();

    const label cellZoneID = oldMesh_.cellZones().findZoneID(cellZoneName_);
    if (cellZoneID == -1)
    {
        FatalErrorIn(type() + "::createNewCellZoneMesh()")
            << "cell zone " << cellZoneName_ << " not found!"
            << abort(FatalError);
    }

    // List of cells in the cellZone of interest
    const labelList& cellIDs = oldMesh_.cellZones()[cellZoneID];

    // Find the maximum X coordinate of any point in the cellZone
    const labelListList& cellPoints = oldMesh_.cellPoints();
    const pointField& points = oldMesh_.points();
    scalar maxXInOldCellZone = -GREAT;
    forAll(cellIDs, cI)
    {
        const label cellID = cellIDs[cI];
        const labelList& curCellPoints = cellPoints[cellID];
        forAll(curCellPoints, cpI)
        {
            const label pointID = curCellPoints[cpI];

            maxXInOldCellZone = max(maxXInOldCellZone, points[pointID].x());
        }
    }

    // Find global maximum
    reduce(maxXInOldCellZone, maxOp<scalar>());

    // Set cutting plane X coordinate
    scalar choppingPlaneXCoord = GREAT;
    if (dict_.found("choppingPlaneXCoordinate"))
    {
        if (debug)
        {
            Info<< "    Looking up choppingPlaneXCoordinate from the remesh "
                << "dict" << endl;
        }

        // Lookup chopping plane x coordinate
        choppingPlaneXCoord =
            readScalar(dict_.lookup("choppingPlaneXCoordinate"));

        // Place reasonable bounds on the chopping plane
        if (choppingPlaneXCoord > 100 || choppingPlaneXCoord < 0)
        {
            FatalErrorIn(type() + "::createNewCellZoneMesh() const")
                << "choppingPlaneXCoordinate should be greater than 0.0 and "
                << "less than 100.0" << abort(FatalError);
        }
    }
    else
    {
        // Define the default cutting plane X coordinate as half the max X
        // coordinate

        choppingPlaneXCoord = 0.5*maxXInOldCellZone;

        if (debug)
        {
            Info<< "    Setting choppingPlaneXCoordinate to half of the max X "
                << "coordinate in the " << cellZoneName_ << " cellZone: "
                << choppingPlaneXCoord << endl;
        }
    }

    const vectorField& CI = oldMesh_.C().internalField();

    forAll(cellIDs, cI)
    {
        const label cellID = cellIDs[cI];

        // Mark this cell as it is in the wire and region-of-interest
        // This chopped cells might end up castellated i.e. there might be a
        // step in the interface, so we will correct for this below
        if (CI[cellID].x() < choppingPlaneXCoord)
        {
            region[cellID] = 1;
        }
    }

    // Pass region list to region field and sync parallel boundaries
    forAll(regionFieldI, cellI)
    {
        regionFieldI[cellI] = region[cellI];
    }
    regionField.correctBoundaryConditions();

    // Correct interface cells so that we have one continuous flat layer of
    // cells downstream at the region of interest

    // Find cells at the downstream interface

    labelHashSet interfaceCellsSet;

    const labelListList& cellCells = oldMesh_.cellCells();
    const cellList& cellFaces = oldMesh_.cells();

    forAll(region, cellI)
    {
        // Check only cells in the region of interest
        if (region[cellI] == 1)
        {
            // Check cell-neighbours on this proc
            const labelList& curCellCells = cellCells[cellI];
            forAll(curCellCells, ccI)
            {
                const label neiCellID = curCellCells[ccI];
                if (region[neiCellID] == 0)
                {
                    // cellI is at the interface
                    interfaceCellsSet.insert(cellI);
                    break;
                }
            }
        }
    }

    // Check processor boundaries
    forAll(regionField.boundaryField(), patchI)
    {
        if (regionField.boundaryField()[patchI].coupled())
        {
            const unallocLabelList& faceCells =
                oldMesh_.boundaryMesh()[patchI].faceCells();

            // The neighbour cell-centre field is stored on the patch
            const scalarField& neiRegion = regionField.boundaryField()[patchI];

            forAll(faceCells, faceI)
            {
                const label cellID = faceCells[faceI];

                if (region[cellID] == 1)
                {
                    // Check neighbour across proc boundary
                    if (neiRegion[faceI] == 0)
                    {
                        if (!interfaceCellsSet.found(cellID))
                        {
                            interfaceCellsSet.insert(cellID);
                        }
                    }
                }
            }
        }
    }

    const labelList interfaceCells = interfaceCellsSet.toc();

    // Recursively correct cells at interface

    if (debug)
    {
        Info<< "    Recursively removing steps in the interface layer" << endl;
    }

    int cellHasBeenRemovedFromInterface = 0;
    const vector xVector = vector(1, 0, 0);
    label i = 0;
    do
    {
        cellHasBeenRemovedFromInterface = 0;

        forAll(interfaceCells, cI)
        {
            const label cellID = interfaceCells[cI];

            // Note: we will not update the interfaceCells list, so we need to
            // check if this cell is still at the kept-side of the interface
            if (region[cellID] ==  1)
            {
                // Check neighbour cells on the current processor
                const labelList& curCellCells = cellCells[cellID];
                forAll(curCellCells, ccI)
                {
                    const label neiCellID = curCellCells[ccI];
                    if (region[neiCellID] == 0)
                    {
                        vector d = CI[neiCellID] - CI[cellID];
                        d /= mag(d);

                        // Compare d vector to the x vector
                        if ((d & xVector) < 0.7071)
                        {
                            // Remove cell
                            region[cellID] = 0;
                            regionFieldI[cellID] = 0.0;
                            cellHasBeenRemovedFromInterface++;
                            break;
                        }
                    }
                }

                // Check neighbour cells on the neighbour processor
                const labelList& curCellFaces = cellFaces[cellID];
                forAll(curCellFaces, cfI)
                {
                    const label faceID = curCellFaces[cfI];

                    // Check if it is a boundary face
                    if (!oldMesh_.isInternalFace(faceID))
                    {
                        // Find patchID and check if it is a coupled patch
                        const label patchID =
                            oldMesh_.boundaryMesh().whichPatch(faceID);

                        if (oldMesh_.boundaryMesh()[patchID].coupled())
                        {
                            const processorPolyPatch& pp =
                                refCast<const processorPolyPatch>
                                (
                                    oldMesh_.boundaryMesh()[patchID]
                                );

                            // Get local faceID
                            const label localFaceID = faceID - pp.start();

                            // Check if the neighbour cell has region equal to
                            // zero
                            if
                            (
                                regionField.boundaryField()
                                [
                                    patchID
                                ][localFaceID] < SMALL
                            )
                            {
                                // Get delta vector for current cell to
                                // neighbour cell centre
                                vector d =
                                    pp.neighbFaceCellCentres()[localFaceID]
                                  - CI[cellID];
                                d /= mag(d);

                                if ((d & xVector) < 0.7071)
                                {
                                    // Remove cell
                                    region[cellID] = 0;
                                    regionFieldI[cellID] = 0.0;
                                    cellHasBeenRemovedFromInterface++;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Sync field at parallel boundaries
        regionField.correctBoundaryConditions();

        i++;

        if (debug)
        {
            Info<< "    Iter = " << i
                << ", removing " << cellHasBeenRemovedFromInterface
                << " cells from the interface" << endl;
        }
    }
    while (cellHasBeenRemovedFromInterface > 0);

    // Take a sub-set of the wire, which contains cells in the region-of-
    // interest
    subsetMesh.setLargeCellSubset(region, 1);

    // Take a reference to the regionOfInterestMesh
    fvMesh& regionOfInterestMesh = subsetMesh.subMesh();


    // 5. Flatten the oldInternalFaces (downstream) patch

    if (debug)
    {
        Info<< nl << "5. flattening the oldInternalFaces (downstream) patch"
            << endl;
    }

    const label oldInternalFacesPatchID =
        regionOfInterestMesh.boundaryMesh().findPatchID("oldInternalFaces");

    if (oldInternalFacesPatchID == -1)
    {
        FatalErrorIn(type() + "::createNewCellZoneMesh()")
            << "Cannot find oldInternalFaces patches in the region-of-interest"
            << " mesh" << abort(FatalError);
    }

    // Find maximum x coordinate on the patch
    const pointField& localPoints =
        regionOfInterestMesh.boundaryMesh()
        [
            oldInternalFacesPatchID
        ].localPoints();

    const scalar maxPointX = gMax(localPoints.component(vector::X));

    // Project all points on the patch to a plane with a normal in the x
    // direction and goes through the point with the maximum x coordinate
    pointField newPoints = regionOfInterestMesh.points();
    const labelList& meshPoints =
        regionOfInterestMesh.boundaryMesh()
        [
            oldInternalFacesPatchID
        ].meshPoints();

    forAll(meshPoints, pI)
    {
        const label pointID = meshPoints[pI];

        newPoints[pointID].x() = maxPointX;
    }

    // Update the regionOfInterestMesh
    // Move the mesh points
    regionOfInterestMesh.movePoints(newPoints);
    regionOfInterestMesh.moving(false);
    regionOfInterestMesh.changing(false);
    regionOfInterestMesh.setPhi().writeOpt() = IOobject::NO_WRITE;

    if (debug)
    {
        Info<< "    Writing " << regionOfInterestMesh.name() << " to time = "
            << runTime.value() << endl;
        regionOfInterestMesh.write();
    }


    // 6. Write out STL of entire region-of-interest boundary
    const word boundarySurfaceNameFms = writeBoundarySTL(regionOfInterestMesh);


    // 7. Create new mesh with cartesianMesh

    if (debug)
    {
        Info<< nl << "7. Creating new mesh with cartesianMesh" << endl;
    }

    // Make a copy of the meshDict file currently in the case as we will
    // overwrite it

    if (debug)
    {
        Info<< "    Copying meshDict to meshDict.beforeRemesh" << endl;
    }

    // Only the master proc copies and moves files
    if (Pstream::master())
    {
        cp
        (
            runTime.caseSystem()/"meshDict",
            runTime.caseSystem()/"meshDict.beforeRemesh"
        );
        cp
        (
            runTime.caseSystem()/"meshDict.gz",
            runTime.caseSystem()/"meshDict.beforeRemesh.gz"
        );
    }

    // Force slave procs to wait
    syncMPI();

    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            runTime.caseSystem(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    // Add settings to mesh dict
    meshDict.add("maxCellSize", maxCellSize());
    meshDict.add
    (
        "surfaceFile",
        fileName("triSurfaces"/fileName(boundarySurfaceNameFms).name())
    );

    // Only the master proc writes the meshDict
    if (Pstream::master())
    {
        if (debug)
        {
            Info<< "    Writing system/meshDict" << endl;
        }

        meshDict.regIOobject::write();
    }


    // Create Cartesian mesh
    // cfMesh writes the mesh to constant so we will make a copy of the
    // orginal mesh and move it back after moving the new mesh

    // Make a copy of the original mesh
    rmDir(runTime.path()/runTime.constant()/"polyMesh.org");
    mv
    (
        runTime.path()/runTime.constant()/"polyMesh",
        runTime.path()/runTime.constant()/"polyMesh.org"
    );

    if (debug)
    {
        Info<< "    Running cartesianMesh" << endl;
    }

    // For now, we will call cartesianMesh from the command line, although it is
    // possible to run it from memory
    // We will trick cartesianMeshGenerator into running in serial
    //const bool parRun = Pstream::parRun();
    // Disable parRun switch
    // if (Pstream::parRun())
    // {
    //     Pstream::parRun() = false;
    // }
    //cartesianMeshGenerator cmg(runTime);
    //cmg.writeMesh();
    // Re-enable parRun switch
    // if (parRun)
    // {
    //     Pstream::parRun() = true;
    // }

    // Each processor will make the same mesh in their own proc directory
    runSystemCommandWithPOpen
    (
        "cartesianMesh",
        "log.cartesianMesh_remesh_time" + runTime.timeName(),
        type(),
        true
    );

    syncMPI();


    // Create time directory
    mkDir(runTime.path()/runTime.timeName());

    // Move polyMesh directory from constant to the time directory
    rmDir(runTime.path()/runTime.timeName()/"polyMesh");
    rmDir(runTime.path()/runTime.timeName()/"remeshedRegionOfInterestSubset");
    mkDir(runTime.path()/runTime.timeName()/"remeshedRegionOfInterestSubset");
    mv
    (
        runTime.path()/runTime.constant()/"polyMesh",
        runTime.path()/runTime.timeName()
       /"remeshedRegionOfInterestSubset/polyMesh"
    );

    // Move original mesh back
    mv
    (
        runTime.path()/runTime.constant()/"polyMesh.org",
        runTime.path()/runTime.constant()/"polyMesh"
    );


    // Copy the mesh dict
    if (Pstream::parRun())
    {
        if (debug)
        {
            Info<< "Copying meshDict to meshDict.remesh" << endl;
        }
        cp
        (
            runTime.caseSystem()/"meshDict",
            runTime.caseSystem()/"meshDict.remesh"
        );
        cp
        (
            runTime.caseSystem()/"meshDict.gz",
            runTime.caseSystem()/"meshDict.remesh.gz"
        );

        if (debug)
        {
            Info<< "Copying meshDict.beforeRemesh to meshDict" << endl;
        }

        rm(fileName(runTime.caseSystem()/"meshDict"));
        cp
        (
            runTime.caseSystem()/"meshDict.beforeRemesh",
            runTime.caseSystem()/"meshDict"
        );
        cp
        (
            runTime.caseSystem()/"meshDict.beforeRemesh.gz",
            runTime.caseSystem()/"meshDict.gz"
        );
    }

    // 8. Extrude downstream patch to create extruded section of wire
    if (debug)
    {
        Info<< nl << "8. Creating extruded section of wire mesh" << endl;
    }

    // Read remeshed region of wire
    polyMesh remeshedRegionOfInterestWire
    (
        IOobject
        (
            "remeshedRegionOfInterestSubset",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // Create extrudeDict
    dictionary extrudeDict;

    dictionary linearNormalCoeffs;

    // Extrude cellZone mesh such that it will have the same final maximum X
    // coordinate
    const scalar extrudeThickness =
        maxXInOldCellZone
      - gMax(remeshedRegionOfInterestWire.points().component(vector::X));

    if (debug)
    {
        Info<< "    Extruding the wire section by " << extrudeThickness << endl;
    }
    linearNormalCoeffs.add("thickness", extrudeThickness);
    extrudeDict.add("linearNormalCoeffs", linearNormalCoeffs);

    // Calculate the number of extruded layers
    if (dict_.found("nLayers"))
    {
        if (debug)
        {
            Info<< "    Setting nLayers from remeshing dict" << endl;
        }

        extrudeDict.add("nLayers", readInt(dict_.lookup("nLayers")));
    }
    else
    {
        const int nLayers = int(extrudeThickness/maxCellSize()) + 1;

        if (debug)
        {
            Info<< "    Setting nLayers to " << nLayers << endl;
        }

        extrudeDict.add("nLayers", nLayers);
    }

    // Create a linearNormal extrude model
    extrudeModels::linearNormal extrudeModel(extrudeDict);


    const label remeshedOldInternalFacesPatchID =
        remeshedRegionOfInterestWire.boundaryMesh().findPatchID
        (
            "oldInternalFaces"
        );

    if (remeshedOldInternalFacesPatchID == -1)
    {
        FatalErrorIn(type() + "::createNewCellZoneMesh()")
            << "Cannot find oldInternalFaces patches in the remeshed "
            "region-of-interest mesh" << abort(FatalError);
    }

    autoPtr<faceMesh> fMesh;
    {
        const polyPatch& pp =
            remeshedRegionOfInterestWire.boundaryMesh()
            [
                remeshedOldInternalFacesPatchID
            ];
        fMesh.reset(new faceMesh(pp.localFaces(), pp.localPoints()));
    }

    extrudedMesh extrudMesh
    (
        IOobject
        (
            "extrudedWire",
            runTime.timeName(),
            runTime
        ),
        fMesh(),
        extrudeModel
    );

    if (debug)
    {
        Info<< "    Writing " << extrudMesh.name() << " mesh" << endl;
        extrudMesh.write();
    }


    // 9. Merge extruded mesh with new wire mesh
    if (debug)
    {
        Info<< nl << "9. Merging extruded wire and remeshed wire" << endl;
    }
    newCellZoneMeshPtr_.set
    (
        new mergePolyMesh(remeshedRegionOfInterestWire)
    );
    mergePolyMesh& mergedWireMesh = newCellZoneMeshPtr_();

    mergedWireMesh.addMesh(extrudMesh);
    mergedWireMesh.merge();

    mergedWireMesh.rename("newWireMerged");

    if (debug)
    {
        Info<< "    Writing " << mergedWireMesh.name() << " to time = "
            << runTime.value() << endl;
        mergedWireMesh.write();
    }


    // 10. Stitch interface between extruded and remeshed wire regions
    if (debug)
    {
        Info<< nl << "10. Stitch interface between extruded and remeshed wire"
            << endl;
    }

    // Create sticher
    polyTopoChanger stitcher(mergedWireMesh);
    stitcher.setSize(1);

    const polyBoundaryMesh& patches = mergedWireMesh.boundaryMesh();

    // Lookup patch IDs of two patches to stich
    const label stitchPatch1ID = patches.findPatchID("oldInternalFaces");
    if (stitchPatch1ID == -1)
    {
        FatalErrorIn(type() + "::createNewCellZoneMesh()")
            << "Cannot find 'oldInternalFaces' patch in the newWire"
            << abort(FatalError);
    }

    const label stitchPatch2ID = patches.findPatchID("originalPatch");
    if (stitchPatch2ID == -1)
    {
        FatalErrorIn(type() + "::createNewCellZoneMesh()")
            << "Cannot find 'originalPatch' patch in the newWire"
            << abort(FatalError);
    }

    // Make list of masterPatch (stitchPatch1) faces
    labelList isf(patches[stitchPatch1ID].size());
    forAll(isf, i)
    {
        isf[i] = patches[stitchPatch1ID].start() + i;
    }

    const word cutZoneName("originalCutFaceZone");
    List<faceZone*> fz
    (
        1,
        new faceZone
        (
            cutZoneName,
            isf,
            boolList(isf.size(), false),
            0,
            mergedWireMesh.faceZones()
        )
    );

    mergedWireMesh.addZones(List<pointZone*>(0), fz, List<cellZone*>(0));

    // Add the perfect interface mesh modifier
    stitcher.set
    (
        0,
        new perfectInterface
        (
            "stitchedInterface",
            0,
            stitcher,
            cutZoneName,
            patches[stitchPatch1ID].name(),
            patches[stitchPatch2ID].name()
        )
    );

    // Execute all polyMeshModifiers
    autoPtr<mapPolyMesh> morphMap = stitcher.changeMesh();

    mergedWireMesh.movePoints(morphMap->preMotionPoints());

    mergedWireMesh.rename("newWireMergedAndStitched");
    polyMesh& newWireMergedAndStitched = mergedWireMesh;


    // 10b. Merge 'sides' and wireContact patch
    if (debug)
    {
        Info<< nl << "10b. Move all faces in the 'sides' patch to the "
            << "wireContact patch" << endl;
    }

    patchManipulationFunctions::mergeTwoPatches
    (
        newWireMergedAndStitched,
        dict_.lookupOrDefault<word>
        (
            "wireContactPatchName", "wireContact1"
        ),
        "sides"
    );

    // 10c. Remove size zero patches

    if (debug)
    {
        Info<< nl << "10c. Remove all size zero patches" << endl;
    }

    patchManipulationFunctions::removePatchesWithNoFaces
    (
        newWireMergedAndStitched
    );

    // 10d. Rename otherSide to downstream

    if (debug)
    {
        Info<< nl << "10d. Renaming downstream patch" << endl;
    }

    const label otherSidePatchID =
        newWireMergedAndStitched.boundaryMesh().findPatchID("otherSide");

    if (otherSidePatchID == -1)
    {
        FatalErrorIn(type() + "::createNewCellZoneMesh()")
            << "Cannot find 'otherSide' patch in the newWireMergedAndStitched"
            << abort(FatalError);
    }

    const_cast<word&>
    (
        newWireMergedAndStitched.boundaryMesh()[otherSidePatchID].name()
    ) = dict_.lookupOrDefault<word>
    (
        "wireDownstreamPatchName",
        "wireDownStream1"
    );

    // Clear zones
    newWireMergedAndStitched.removeZones();

    // 10e. Add wire cell zone

    if (debug)
    {
        Info<< nl << "10e. Add all cells in newWireMergedAndStitched to a "
            << "cellZone" << endl;
    }

    newWireMergedAndStitched.cellZones().setSize(1);
    {
        labelList cellZoneIDs(newWireMergedAndStitched.nCells(), 0);
        forAll(cellZoneIDs, cellI)
        {
            cellZoneIDs[cellI] = cellI;
        }
        newWireMergedAndStitched.cellZones().set
        (
            0,
            new cellZone
            (
                cellZoneName_,
                cellZoneIDs,
                0,
                newWireMergedAndStitched.cellZones()
            )
        );
    }

    // 10f. Add upstream and downstream faceZones

    if (debug)
    {
        Info<< nl << "10f. Create upstream and downstream patch faceZones"
            << endl;
    }

    patchManipulationFunctions::addFaceZoneFromPatch
    (
        newWireMergedAndStitched,
        upstreamPatchName_,
        upstreamFaceZoneName_,
        true
    );

    patchManipulationFunctions::addFaceZoneFromPatch
    (
        newWireMergedAndStitched,
        downstreamPatchName_,
        downstreamFaceZoneName_,
        true
    );

    // Set all zones to write
    newWireMergedAndStitched.faceZones().writeOpt() = IOobject::AUTO_WRITE;
    newWireMergedAndStitched.cellZones().writeOpt() = IOobject::AUTO_WRITE;

    if (debug)
    {
        Info<< "    Writing " << newWireMergedAndStitched.name()
            << " to time = " << runTime.value() << endl;
        newWireMergedAndStitched.write();
    }
}


void Foam::replaceCellZoneMesh::createNewGlobalMesh() const
{
    if (newGlobalMeshPtr_.valid())
    {
        FatalErrorIn(type() + "::createNewGlobalMesh()")
            << "Pointer already set!" << abort(FatalError);
    }

    // 11. Make a copy of original mesh and remove old wire
    // 12. Add new wire mesh to originalMeshWithoutWire
    // 13. reorder newGlobalMesh patches to be in the same order as the oldMesh

    // Take references for convenience
    const Time& runTime = oldMesh_.time();


    // 11. Make a copy of original mesh and remove old wire
    if (debug)
    {
        Info<< nl << "11. Making a copy of old mesh and removing the wire"
            << endl;
    }

    // Create a copy of the old mesh
    polyMesh oldMeshCopy
    (
        IOobject
        (
            "oldMeshCopy",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        xferCopy(oldMesh_.points()),
        xferCopy(oldMesh_.faces()),
        xferCopy(oldMesh_.faceOwner()),
        xferCopy(oldMesh_.faceNeighbour()),
        true
    );

    // Add the boundary
    List<polyPatch*> boundaryMeshCopy(oldMesh_.boundaryMesh().size());
    forAll(boundaryMeshCopy, patchI)
    {
        boundaryMeshCopy[patchI] =
            oldMesh_.boundaryMesh()[patchI].clone
            (
                oldMeshCopy.boundaryMesh()
            ).ptr();
    }
    oldMeshCopy.addPatches(boundaryMeshCopy);

    // Add the zones
    List<pointZone*> pz(oldMesh_.pointZones().size());
    forAll(pz, pzI)
    {
        pz[pzI] =
            new pointZone
            (
                oldMesh_.pointZones()[pzI].name(),
                oldMesh_.pointZones()[pzI],
                pzI,
                oldMeshCopy.pointZones()
            );
    }
    List<faceZone*> fz(oldMesh_.faceZones().size());
    forAll(fz, fzI)
    {
        fz[fzI] =
            new faceZone
            (
                oldMesh_.faceZones()[fzI].name(),
                oldMesh_.faceZones()[fzI],
                oldMesh_.faceZones()[fzI].flipMap(),
                fzI,
                oldMeshCopy.faceZones()
            );
    }
    List<cellZone*> cz(oldMesh_.cellZones().size());
    forAll(cz, czI)
    {
        cz[czI] =
            new cellZone
            (
                oldMesh_.cellZones()[czI].name(),
                oldMesh_.cellZones()[czI],
                czI,
                oldMeshCopy.cellZones()
            );
    }
    oldMeshCopy.addZones(pz, fz, cz);


    // Topo changer that will actually change the mesh
    directTopoChange meshChanger(oldMeshCopy);

    // Create cell remover
    removeCells cellRemover(oldMeshCopy, true);

    const label cellZoneID = oldMesh_.cellZones().findZoneID(cellZoneName_);

    if (cellZoneID == -1)
    {
        FatalErrorIn(type() + "::createNewGlobalMesh() const")
            << "cell zone " << cellZoneName_ << " not found!"
            << abort(FatalError);
    }

    // Set old wire cells to be removed from oldMeshCopy
    cellRemover.setRefinement
    (
        oldMesh_.cellZones()[cellZoneID], // list of cellIDs to remove
        labelList(0),     // facesToExpose: none
        labelList(0),     // patchIDs for facesToExpose: none
        meshChanger
    );

    // Remove old cellZone cells from oldMeshCopy

    if (debug)
    {
        Info<< "    Removing the old cellZone cells from the old mesh" << endl;
    }

    autoPtr<mapPolyMesh> map = meshChanger.changeMesh(oldMeshCopy, false);

    oldMeshCopy.rename("oldMeshWithoutWire");
    polyMesh& oldMeshWithoutWire = oldMeshCopy;

    // Note: for now, we must write this mesh as mergePolyMesh only has a read
    // from disk constructor
    WarningIn(type() + "::createNewGlobalMesh()")
        << "For now we must write " << oldMeshWithoutWire.name() << " to disk"
         << " because the mergePolyMesh class only reads from disk!" << endl;
    Info<< "    Writing " << oldMeshWithoutWire.name() << " to time = "
        << runTime.value() << endl;
    oldMeshWithoutWire.write();


    // 12. Add new wire mesh to originalMeshWithoutWire

    if (debug)
    {
        Info<< nl << "12. Adding new wire to oldMeshWithoutWire"
            << endl;

    }

    //autoPtr<mergePolyMesh> globalMeshPtr;
    if (Pstream::parRun())
    {
        // Reconstruct the oldMeshWithoutWire

        // Note: be careful of the distinction between globalTime and
        // runTime (local time)

        // Reconstruct the region

        if (debug)
        {
            Info<< "    Running reconstructParMeshZones" << endl;
        }

        if (Pstream::master())
        {
            WarningIn(type() + "::createNewGlobalMesh(...)")
                << "Removing V0 and meshPhi from " << globalTime().timePath()
                << " due to a bug in fvMesh.C constructor for sub-regions"
                << endl;
            rm(globalTime().timePath()/"V0");
            rm(globalTime().timePath()/"V0.gz");
            rm(globalTime().timePath()/"meshPhi");
            rm(globalTime().timePath()/"meshPhi.gz");

            chDir(globalTime().path());
            runSystemCommandWithPOpen
            (
                "reconstructParMeshZones -region " + oldMeshWithoutWire.name()
              + " -latestTime",
                "log.reconstructParMeshZones_remesh_time"
              + runTime.timeName(),
                type(),
                true
            );
            chDir(runTime.path());
        }

        // Force slave procs to wait for master to finish reconstructing
        // the oldMeshWithoutWire
        syncMPI();

        // Read mesh from serial case

        chDir(globalTime().constant());
        runSystemCommandWithPOpen
        (
            "\\ln -s ../" + Foam::name(runTime.value()) + "/"
          + oldMeshWithoutWire.name() + " .",
            "log.ln_remesh_time" + runTime.timeName(),
            type(),
            false // it is OK if this fails e.g. if the link already exists
        );
        chDir(runTime.path());

        if (debug)
        {
            Info<< "    Reading global mesh" << endl;
        }

        suspendMPI();
        newGlobalMeshPtr_.set
        (
            new mergePolyMesh
            (
                IOobject
                (
                    oldMeshWithoutWire.name(),
                    globalTime().timeName(),
                    globalTime(),
                    IOobject::MUST_READ
                )
            )
        );
        resumeMPI();
    }
    else
    {
        newGlobalMeshPtr_.set
        (
            new mergePolyMesh
            (
                oldMeshWithoutWire
            )
        );
    }


    // Create new cellZone mesh
    const polyMesh& newWireMesh = newCellZoneMesh();

    newGlobalMeshPtr_().addMesh(newWireMesh);
    newGlobalMeshPtr_().merge();

    newGlobalMeshPtr_().rename("newGlobalMesh");

    if (debug)
    {
        Info<< "    Writing " << newGlobalMeshPtr_().name() << " to time = "
            << runTime.value() << endl;

        if (Pstream::master())
        {
            newGlobalMeshPtr_().write();
        }
    }

    // 13. reorder newGlobalMesh patches to be in the same order as the oldMesh
    // patches

    if (debug)
    {
        Info<< nl << "13. Reorder patches consistently" << endl;
    }

    patchManipulationFunctions::reorderPatchesConsistently
    (
        oldMesh_, newGlobalMeshPtr_()
    );

    if (debug)
    {
        Info<< "    Writing " << newGlobalMeshPtr_().name() << " to time = "
            << runTime.value() << endl;

        if (Pstream::master())
        {
            newGlobalMeshPtr_().write();
        }
    }
}


void Foam::replaceCellZoneMesh::calcMaxCellSize() const
{
    if (maxCellSize_ > 0.0)
    {
        FatalErrorIn(type() + "calcMaxCellSize() const")
            << "maxCellSize has already been set!" << abort(FatalError);
    }

    if (dict_.found("maxCellSize"))
    {
        if (debug)
        {
            Info<< "    Reading maxCellSize from remesh dict" << endl;
        }

        maxCellSize_ = readScalar(dict_.lookup("maxCellSize"));

    }
    else
    {
        // Estimate max cell size from the oldMesh
        maxCellSize_ = Foam::cbrt(gMax(oldMesh_.V()));

        if (debug)
        {
            Info<< "    Estimating maxCellSize from the old mesh: "
                << maxCellSize_ << endl;
        }
    }
}


const Foam::word Foam::replaceCellZoneMesh::writeBoundarySTL
(
    const fvMesh& mesh
) const
{
    if (debug)
    {
        Info<< nl << "Writing boundary STL and FMS surface files" << endl;
    }

    // Take references for convenience
    const Time& runTime = oldMesh_.time();


    // Triangulate all faces of the mesh
    meshTriangulation triangulatedFaces
    (
        mesh,
        -1,    // internal faces patch
        boolList(mesh.nCells(), true), // included cells
        false   // face centre decomposition
    );

    // Create triangulated surface of just the boundary

    // Mark all boundary faces
    const label nFaces = triangulatedFaces.size();
    const label nIntFaces = triangulatedFaces.nInternalFaces();

    boolList includedFaces(nFaces, false);
    labelList boundaryPointMap(triangulatedFaces.nPoints(), -1);
    labelList boundaryFaceMap(nFaces, -1);

    const labelList& faceMap = triangulatedFaces.faceMap();

    for (int triFaceI = nIntFaces; triFaceI < nFaces; triFaceI++)
    {
        const label faceID = faceMap[triFaceI];

        // Skip faces on empty and processor boundaries
        const label patchID =
            mesh.boundaryMesh().whichPatch(faceID);

        if
        (
            mesh.boundaryMesh()[patchID].type()
         != emptyPolyPatch::typeName
         && mesh.boundaryMesh()[patchID].type()
         != processorPolyPatch::typeName
        )
        {
            includedFaces[triFaceI] = true;
        }
    }

    triSurface triangulatedBoundaryFaces =
        triangulatedFaces.subsetMesh
        (
            includedFaces,
            boundaryPointMap,
            boundaryFaceMap
        );

    mkDir(runTime.path()/"triSurfaces");

    fileName boundarySurfaceName =
        runTime.path()/"triSurfaces"/"meshBoundary_"
      + Foam::name(runTime.timeIndex());


    // Use Nagata patches
    if (dict_.lookupOrDefault<Switch>("useNagataPatches", true))
    {
        nagataPatch nagataSurf(triangulatedBoundaryFaces);

        // Remove empty patches from the STL
        nagataSurf.removeEmptyPacthes(triangulatedBoundaryFaces);

        // Refine and project surface
        triSurface refinedAndProjectedSurf =
            nagataSurf.refineAndProjectPatches
            (
                wordList
                (
                    1,
                    dict_.lookupOrDefault<word>
                    (
                        "wireContactPatchName", "wireContact1"
                    )
                ),
                2 // nRefinementLevels
            );

        triangulatedBoundaryFaces = refinedAndProjectedSurf;
    }

    // Write STL
    if (Pstream::parRun())
    {
        // Reconstruct STL before writing
        // We will reconstruct the entire STL on all processors

        // Calculate global bounding box
        List<treeBoundBox> globalBoundBox
        (
            1, treeBoundBox(boundBox(mesh.points(), true))
        );

        // Create dictionary for distribution
        dictionary triSurfaceDict;
        triSurfaceDict.add("bounds", globalBoundBox);
        //meshBb[Pstream::myProcNo()]);
        triSurfaceDict.add("distributionType", "follow");
        triSurfaceDict.add("mergeDistance", SMALL);

        if (debug)
        {
            Info<< "    Creating the distributedTriSurfaceMesh" << endl;
        }

        distributedTriSurfaceMesh parallelTriangulatedBoundaryFaces
        (
            IOobject
            (
                "parallelTriangulatedBoundaryFaces.stl",
                runTime.path(),
                "triSurfaces",
                oldMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            triangulatedBoundaryFaces,
            triSurfaceDict
        );

        autoPtr<mapDistribute> faceMap;
        autoPtr<mapDistribute> pointMap;

        if (debug)
        {
            Info<< "    Reconstructing the STL surface" << endl;
        }

        parallelTriangulatedBoundaryFaces.distribute
        (
            globalBoundBox, //meshBb[Pstream::myProcNo()],
            true, // keepNonMapped,
            faceMap,
            pointMap
        );
        faceMap.clear();
        pointMap.clear();

        // Write the reconstructed STL file
        // It is OK for every processor to write the file, as they write to
        // their own local processor directory

        if (debug)
        {
            Info<< "    Writing boundary surface mesh to "
                << boundarySurfaceName << endl;
        }

        parallelTriangulatedBoundaryFaces.triSurface::write
        (
            boundarySurfaceName + ".stl"
        );

        if (debug)
        {
            Info<< "    Writing serial boundary surface mesh to "
                << word("serial_" + boundarySurfaceName) << endl;

            triangulatedBoundaryFaces.triSurface::write
            (
                fileName
                (
                    runTime.path()/"triSurfaces"/"serial_meshBoundary_"
                  + Foam::name(runTime.timeIndex()) + ".stl"
                )
            );
        }
    }
    else
    {
        // Serial run
        if (debug)
        {
            Info<< "Writing boundary surface mesh to " << boundarySurfaceName
                << ".stl" << endl;
        }
        triangulatedBoundaryFaces.write(boundarySurfaceName + ".stl");
    }


    // Detect edges on the STL file and convert to FMS file

    if (debug)
    {
        Info<< "    Detecting surface feature edges" << endl;
    }

    // triSurf does not run in parallel; however, in our case we have
    // reconstructed the STL so it is the same on each processor. So we will
    // trick triSurf into thinking it is running in serial on each processor by
    // turning off the Pstream::parRun() switch before calling it and then
    // turning the switch back on again

    triSurf stlSurface(boundarySurfaceName + ".stl");

    const bool parRun = Pstream::parRun();
    if (Pstream::parRun())
    {
        Pstream::parRun() = false;
    }

    triSurfaceDetectFeatureEdges edgeDetector(stlSurface);
    edgeDetector.detectFeatureEdges();

    if (parRun)
    {
        Pstream::parRun() = true;
    }

    const word boundarySurfaceNameFms = boundarySurfaceName.name() + ".fms";

    if (debug)
    {
        Info<< "    Writing : " << boundarySurfaceNameFms << endl;
    }

    stlSurface.writeSurface
    (
        fileName(runTime.path()/"triSurfaces"/boundarySurfaceNameFms)
    );

    return boundarySurfaceNameFms;
}


void Foam::replaceCellZoneMesh::suspendMPI() const
{
    if (!parRun_)
    {
        FatalErrorIn(type() + "::suspendMPI() const")
            << "This function can only be called in parallel!"
            << abort(FatalError);
    }

    Pstream::parRun() = false;
}


void Foam::replaceCellZoneMesh::resumeMPI() const
{
    if (!parRun_)
    {
        FatalErrorIn(type() + "::resumeMPI() const")
            << "This function can only be called in parallel!"
            << abort(FatalError);
    }

    Pstream::parRun() = true;
}


void Foam::replaceCellZoneMesh::syncMPI() const
{
    // Perform a global reduction to force the synchronisation of all procs
    returnReduce(bool(true), orOp<bool>());
}


void Foam::replaceCellZoneMesh::makeGlobalTime() const
{
    if (!parRun_)
    {
        FatalErrorIn(type() + "::makeGlobalTime() const")
            << "This function can only be called in parallel!"
            << abort(FatalError);
    }

    if (globalTimePtr_.valid())
    {
        FatalErrorIn(type() + "::makeGlobalTime() const")
            << "Poiner already set!"
            << abort(FatalError);
    }

    const Time& localTime = oldMesh_.time();

    suspendMPI();

    globalTimePtr_.set
    (
        new Time
        (
            localTime.rootPath(),
            localTime.caseName().components()[0],
            "system",
            "constant",
            false // enable function objects
        )
    );

    globalTimePtr_().setTime(localTime);

    resumeMPI();
}


const Foam::Time& Foam::replaceCellZoneMesh::globalTime() const
{
    if (globalTimePtr_.empty())
    {
        makeGlobalTime();
    }

    return globalTimePtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::replaceCellZoneMesh::replaceCellZoneMesh
(
    fvMesh& mesh,
    const word& cellZoneName,
    const word& upstreamPatchName,
    const word& upstreamFaceZoneName,
    const word& downstreamPatchName,
    const word& downstreamFaceZoneName,
    const dictionary& dict
)
:
    oldMesh_(mesh),
    oldGlobalMeshPtr_(),
    cellZoneName_(cellZoneName),
    upstreamPatchName_(upstreamPatchName),
    upstreamFaceZoneName_(upstreamFaceZoneName),
    downstreamPatchName_(downstreamPatchName),
    downstreamFaceZoneName_(downstreamFaceZoneName),
    dict_(dict),
    newCellZoneMeshPtr_(),
    newGlobalMeshPtr_(),
    maxCellSize_(-1.0),
    parRun_(Pstream::parRun()),
    globalTimePtr_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::replaceCellZoneMesh::~replaceCellZoneMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::replaceCellZoneMesh::oldGlobalMesh() const
{
    if (!parRun_)
    {
        return oldMesh_;
    }

    if (oldGlobalMeshPtr_.empty())
    {
        createOldGlobalMesh();
    }

    return oldGlobalMeshPtr_();
}


Foam::fvMesh& Foam::replaceCellZoneMesh::oldGlobalMesh()
{
    if (!parRun_)
    {
        return oldMesh_;
    }

    if (oldGlobalMeshPtr_.empty())
    {
        createOldGlobalMesh();
    }

    return oldGlobalMeshPtr_();
}


const Foam::polyMesh& Foam::replaceCellZoneMesh::newCellZoneMesh() const
{
    if (newCellZoneMeshPtr_.empty())
    {
        createNewCellZoneMesh();
    }

    return newCellZoneMeshPtr_();
}


const Foam::polyMesh& Foam::replaceCellZoneMesh::newGlobalMesh() const
{
    if (newGlobalMeshPtr_.empty())
    {
        createNewGlobalMesh();
    }

    return newGlobalMeshPtr_();
}


Foam::scalar Foam::replaceCellZoneMesh::maxCellSize() const
{
    if (maxCellSize_ < 0.0)
    {
        calcMaxCellSize();
    }

    return maxCellSize_;
}


void Foam::replaceCellZoneMesh::clearOut()
{
    oldGlobalMeshPtr_.clear();
    newCellZoneMeshPtr_.clear();
    newGlobalMeshPtr_.clear();
    globalTimePtr_.clear();
}

// ************************************************************************* //
