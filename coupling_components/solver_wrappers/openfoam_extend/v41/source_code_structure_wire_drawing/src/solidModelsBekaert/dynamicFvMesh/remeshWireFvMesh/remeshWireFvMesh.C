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

#include "remeshWireFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "replaceCellZoneMesh.H"
#include "mergePolyMesh.H"
#include "directTopoChange.H"
#include "MapVolFields.H"
#include "MapConsistentVolFields.H"
#include "UnMapped.H"
#include "meshToMesh.H"
#include "polyModifyFace.H"
#include "geometricFieldContainer.H"
#include "runSystemCommandWithPOpen.H"
#include "pointMeshMapper.H"
#include "nonLinearGeometry.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(remeshWireFvMesh, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        remeshWireFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::remeshWireFvMesh::checkRemeshWire() const
{
    // Check remeshing is active and check remeshing frequency
    if (remeshWire_ && (time().timeIndex() % remeshFrequency_ == 0))
    {
        return true;
    }

    return false;
}


void Foam::remeshWireFvMesh::writeMeshFields(const fvMesh& mesh) const
{
    // First, write old mesh to disk and clear out mesh geometry beforehand

    if (debug)
    {
        Info<< "    Writing mesh to time = " << mesh.time().value()
            << endl;
    }

    mesh.write();

    // Write all of the volFields in the objectRegistry
    WriteAllObjects<scalar, fvPatchField, volMesh>(mesh);
    WriteAllObjects<vector, fvPatchField, volMesh>(mesh);
    WriteAllObjects<tensor, fvPatchField, volMesh>(mesh);
    WriteAllObjects<symmTensor, fvPatchField, volMesh>(mesh);
    WriteAllObjects<diagTensor, fvPatchField, volMesh>(mesh);
    WriteAllObjects<sphericalTensor, fvPatchField, volMesh>(mesh);

    WarningIn(type() + "::writeMeshFields(...)")
        << "Removing V0 from time = " << mesh.time().value()
        << " due to a bug in fvMesh.C constructor when reading V0 and meshPhi "
        << "from sub-regions" << endl;
    rm(mesh.time().timePath()/"V0");
    rm(mesh.time().timePath()/"V0.gz");
    rm(mesh.time().timePath()/"meshPhi");
    rm(mesh.time().timePath()/"meshPhi.gz");

    // If we are running in parallel, then we will reconstruct the fields onto
    // the global mesh
    if (Pstream::parRun() && Pstream::master())
    {
        if (debug)
        {
            Info<< "    Reconstructing fields onto global mesh" << endl;
        }

        chDir(mesh.time().path()/"..");
        runSystemCommandWithPOpen
        (
            "reconstructParMeshZones -latestTime",
            "log.reconstructParMeshZones_writeFields_"
           + Foam::name(mesh.time().timeIndex()),
            type(),
            true
        );
        chDir(mesh.time().path());
    }
}


void Foam::remeshWireFvMesh::mapFieldsAndOverwriteTargetMesh
(
    const polyMesh& newMesh,
    fvMesh& oldMesh
)
{
    // 14. Map fields from complete old mesh to complete new mesh
    // 15. Overwrite the old mesh with the new mesh
    // 16. Clear out mesh info
    // 17. Overwrite the old fields with the new fields


    // Record the oldMesh info, before we change oldMesh
    //const meshInfo oldMeshInfo(targetMesh);


    // 14. Map fields from complete old mesh to complete new mesh
    if (debug)
    {
        Info<< nl << "14. Mapping fields from old mesh to new mesh"
            << endl;
    }

    // Transfer the polyMesh to a fvMesh

    WarningIn(type() + "::mapFieldsAndOverwriteTargetMesh()")
        << "Currently newMesh must be written to disk to transfer polyMesh to "
        << "fvMesh"<< endl;

    if (Pstream::master())
    {
        newMesh.write();

        WarningIn(type() + "::mapFieldsAndOverwriteTargetMesh(...)")
            << "Removing meshPhi from time = " << newMesh.time().value()
            << " due to a bug in fvMesh.C constructor when reading meshPhi"
            << " from sub-regions" << endl;
        rm(newMesh.time().timePath()/"meshPhi");
        rm(newMesh.time().timePath()/"meshPhi.gz");
    }

    // Force procs to synchronise
    returnReduce(bool(true), orOp<bool>());

    fvMesh newFMesh(newMesh);

    // Force procs to synchronise
    returnReduce(bool(true), orOp<bool>());

    // Map all fields from old mesh to the new mesh
    // The master writes the fields after mapping
    chDir(newMesh.time().path());
    mapConsistentMesh(oldMesh, newFMesh);
    chDir(oldMesh.time().path());

    WarningIn(type() + "::mapFieldsAndOverwriteTargetMesh()")
        << "TO-DO: perform direct mapping on non-wire field!" << endl;

    // Force procs to synchronise
    returnReduce(bool(true), orOp<bool>());

    if (debug)
    {
        Info<< "    Writing " << newFMesh.name() << " mesh to time = "
            << time().value() << endl;
    }

    if (Pstream::master())
    {
        newFMesh.write();
    }

    // Force procs to synchronise
    returnReduce(bool(true), orOp<bool>());

    // Force procs to synchronise
    returnReduce(bool(true), orOp<bool>());

    // 15. Overwrite the old mesh with the new mesh
    if (debug)
    {
        Info<< nl << "15. Overwrite oldMesh with newMesh" << endl;
    }

    // Update oldMesh

    labelList patchSizes(newFMesh.boundaryMesh().size(), 0);
    labelList patchStarts(newFMesh.boundaryMesh().size(), 0);
    forAll(patchSizes, patchI)
    {
        patchSizes[patchI] = newFMesh.boundaryMesh()[patchI].size();
        patchStarts[patchI] = newFMesh.boundaryMesh()[patchI].start();
    }

    oldMesh.resetPrimitives
    (
        xferCopy(newFMesh.points()),
        xferCopy(newFMesh.faces()),
        xferCopy(newFMesh.faceOwner()),
        xferCopy(newFMesh.faceNeighbour()),
        patchSizes,
        patchStarts,
        true // valid boundary
    );

    // Copy the zones from newFMesh to oldMesh
    copyZones(newFMesh, oldMesh);


    // 15b. Create a map from the old mesh to the new mesh
    // TO-DO: work-in-progress: we should be able to use a mapPolyMesh rather
    // than using the ReplaceOldMeshFieldsWithNewMeshFields functions
    // if (debug)
    // {
    //    Info<< nl << "15b. Creating a map from the old mesh to the new mesh"
    //        << endl;
    // }
    // // Maps: from current mesh to old mesh
    // labelList pointMap(oldMesh.nPoints(), -1);
    // labelList faceMap(oldMesh.nFaces(), -1);
    // labelList cellMap(oldMesh.nCells(), -1);
    // // Reverse maps: from old mesh to current mesh
    // labelList reversePointMap(oldMeshInfo.nPoints_, -1);
    // labelList reverseFaceMap(oldMeshInfo.nFaces_, -1);
    // labelList reverseCellMap(oldMeshInfo.nCells_, -1);
    // // Set the reversePointMap etc?
    // mapPolyMesh mpm
    // (
    //     oldMesh,
    //     oldMeshInfo.nPoints_,
    //     oldMeshInfo.nFaces_,
    //     oldMeshInfo.nCells_,
    //     pointMap,
    //     List<objectMap>(), // pointsFromPoints,
    //     faceMap,
    //     List<objectMap>(), // facesFromPoints,
    //     List<objectMap>(), // facesFromEdges,
    //     List<objectMap>(), // facesFromFaces,
    //     cellMap,
    //     List<objectMap>(), // cellsFromPoints,
    //     List<objectMap>(), // cellsFromEdges,
    //     List<objectMap>(), // cellsFromFaces,
    //     List<objectMap>(), // cellsFromCells,
    //     reversePointMap,
    //     reverseFaceMap,
    //     reverseCellMap,
    //     labelHashSet(), // flipFaceFlux,
    //     labelListList(0), // patchPointMap,
    //     labelListList(0), // pointZoneMap,
    //     labelListList(0), // faceZonePointMap,
    //     labelListList(0), // faceZoneFaceMap,
    //     labelListList(0), // cellZoneMap,
    //     pointField(oldMesh.nPoints(), vector::zero), // preMotionPoints,
    //     oldMeshInfo.patchStarts_, // oldPatchStarts,
    //     oldMeshInfo.patchNMeshPoints_ // oldPatchNMeshPoints
    // );


    if (debug)
    {
        Info<< "    Writing " << this->name() << " mesh and fields to time = "
             << time().value() << endl;
    }

    if (Pstream::master())
    {
        oldMesh.write();
    }

    // Force sync after mesh write
    returnReduce(bool(true), orOp<bool>());

    // 16. Clear out mesh info e.g. delta vectors
    if (debug)
    {
       Info<< nl << "16. Clear-out mesh geometry and addressing" << endl;
    }

    // Clear out mesh data before performing mapping, as we do not want to map
    // mesh fields: this is important, otherwise non-mesh fields can get messed
    // up due to memory issues
    oldMesh.clearOut();
    newFMesh.clearOut();

    // Force parallel sync
    returnReduce(bool(true), orOp<bool>());

    // 17. Overwrite the old fields with the new fields
    if (debug)
    {
        Info<< nl << "17. Overwrite old fields with new fields" << endl;
    }

    replaceMeshFields(oldMesh, newFMesh, true);

    // Force sync after writing fields
    returnReduce(bool(true), orOp<bool>());
}


void Foam::remeshWireFvMesh::mapConsistentMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget
)
{
    // Create the interpolation scheme
    meshToMesh meshToMeshInterp(meshSource, meshTarget);

    if (debug)
    {
        Info<< nl
            << "Consistently creating and mapping fields for time "
            << meshSource.time().timeName() << nl << endl;
    }

    // In parallel, we must read the global fields from disk, whereas in serial
    // we look fields up from the registry
    if (Pstream::parRun())
    {
        // Search for list of objects for this time
        IOobjectList objects(meshSource, meshSource.time().timeName());

        // Map volFields form disk
        // ~~~~~~~~~~~~~
        newMapConsistentVolFields<scalar>(objects, meshToMeshInterp);
        newMapConsistentVolFields<vector>(objects, meshToMeshInterp);
        newMapConsistentVolFields<sphericalTensor>(objects, meshToMeshInterp);
        newMapConsistentVolFields<symmTensor>(objects, meshToMeshInterp);
        newMapConsistentVolFields<tensor>(objects, meshToMeshInterp);
    }
    else
    {
        // Map volFields from memory
        // ~~~~~~~~~~~~~
        newMapConsistentVolFields<scalar>(meshToMeshInterp);
        newMapConsistentVolFields<vector>(meshToMeshInterp);
        newMapConsistentVolFields<sphericalTensor>(meshToMeshInterp);
        newMapConsistentVolFields<symmTensor>(meshToMeshInterp);
        newMapConsistentVolFields<tensor>(meshToMeshInterp);
    }

    // This doesn't do anything
    // {
    //     // Search for list of target objects for this time
    //     IOobjectList objects(meshTarget, meshTarget.time().timeName());

    //     // Mark surfaceFields as unmapped
    //     // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //     UnMapped<surfaceScalarField>(objects);
    //     UnMapped<surfaceVectorField>(objects);
    //     UnMapped<surfaceSphericalTensorField>(objects);
    //     UnMapped<surfaceSymmTensorField>(objects);
    //     UnMapped<surfaceTensorField>(objects);

    //     // Mark pointFields as unmapped
    //     // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //     UnMapped<pointScalarField>(objects);
    //     UnMapped<pointVectorField>(objects);
    //     UnMapped<pointSphericalTensorField>(objects);
    //     UnMapped<pointSymmTensorField>(objects);
    //     UnMapped<pointTensorField>(objects);
    // }

    // We have no Lagrangian particles
    //mapLagrangian(meshToMeshInterp);
}


void Foam::remeshWireFvMesh::correctMaterialIndexFields
(
    const fvMesh& mesh
) const
{
    // 18. Update material indices fields
    if (debug)
    {
        Info<< nl << "18. Correcting material indices fields from cellZones"
            << endl;
    }

    // Lookup materials index field
    volScalarField& materials = const_cast<volScalarField&>
    (
        mesh.lookupObject<volScalarField>("materials")
    );
    scalarField& materialsI = materials.internalField();
    materialsI = -1;

    // Reset materials field from cellZones
    forAll(mesh.cellZones(), cellZoneI)
    {
        const labelList& curCellZoneIDs = mesh.cellZones()[cellZoneI];

        // Lookup "curMaterial" index field: this is 1 in cells of the "current"
        // material and 0 elsewhere
        volScalarField& curMat = const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>
            (
                "curMaterial" + Foam::name(cellZoneI)
            )
        );
        scalarField& curMatI = curMat.internalField();
        curMatI = 0;

        forAll(curCellZoneIDs, cI)
        {
            const label cellID = curCellZoneIDs[cI];
            materialsI[cellID] = cellZoneI;
            curMatI[cellID] = 1.0;
        }

        if (Pstream::master())
        {
            if (debug)
            {
                Info<< "    Writing " << curMat.name()
                    << " to time = " << time().value() << endl;
            }
            curMat.write();
        }
    }

    if (gMin(materialsI) < 0)
    {
        FatalErrorIn(type() + "::update()")
            << "There are some cells not in any cellZone!"
            << abort(FatalError);
    }

    // Correct material boundaries (zero gradient)
    materials.correctBoundaryConditions();

    if (debug)
    {
        Info<< "    Writing " << materials.name() << " to time = "
            << time().value() << endl;
    }

    if (Pstream::master())
    {
        materials.write();
    }

    // Force sync of procs after write
    returnReduce(bool(true), orOp<bool>());
}


void Foam::remeshWireFvMesh::updateBoundaryConditionsAndInterp
(
    const fvMesh& mesh,
    const remeshWireFvMesh::meshInfo& oldMeshInfo,
    const bool onlyMasterWritesFields
) const
{
    if (debug)
    {
        Info<< nl << "19. Updating boundary conditions" << endl;
    }

    // Create a map that will be just used for mapping boundaries
    mapPolyMesh mpm
    (
        mesh,
        oldMeshInfo.nPoints_,
        oldMeshInfo.nFaces_,
        oldMeshInfo.nCells_,
        labelList(mesh.nPoints(), -1), // pointMap,
        List<objectMap>(), // pointsFromPoints,
        labelList(mesh.nFaces(), -1), // faceMap,
        List<objectMap>(), // facesFromPoints,
        List<objectMap>(), // facesFromEdges,
        List<objectMap>(), // facesFromFaces,
        labelList(mesh.nCells(), -1), // cellMap,
        List<objectMap>(), // cellsFromPoints,
        List<objectMap>(), // cellsFromEdges,
        List<objectMap>(), // cellsFromFaces,
        List<objectMap>(), // cellsFromCells,
        labelList(oldMeshInfo.nPoints_, -1), // reversePointMap,
        labelList(oldMeshInfo.nFaces_, -1), // reverseFaceMap,
        labelList(oldMeshInfo.nCells_, -1), // reverseCellMap,
        labelHashSet(), // flipFaceFlux,
        labelListList(0), // patchPointMap,
        labelListList(0), // pointZoneMap,
        labelListList(0), // faceZonePointMap,
        labelListList(0), // faceZoneFaceMap,
        labelListList(0), // cellZoneMap,
#if FOAMEXTEND > 40
        boolList(mesh.boundaryMesh().size(), false), // resetPatchFlag
#endif
        pointField(mesh.nPoints(), vector::zero), // preMotionPoints,
        oldMeshInfo.patchStarts_, // oldPatchStarts,
        oldMeshInfo.patchNMeshPoints_ // oldPatchNMeshPoints
    );

    // Create face mapper
    faceMapper fMapper(mpm);

    // Lookup DU field
    volVectorField& DU = const_cast<volVectorField&>
    (
        mesh.lookupObject<volVectorField>("DU")
    );

    // Map each patch on the DU field boundary
    forAll(mesh.boundary(), patchI)
    {
        // Create a patch mapper
        fvPatchMapper patchMapper(mesh.boundary()[patchI], fMapper);

        DU.boundaryField()[patchI].autoMap(patchMapper);
    }

    // grad(DU) is needed for correcting the boundaries, so we will create it if
    // it does not exit
    autoPtr<volTensorField> gradDUPtr;
    if (!mesh.foundObject<volTensorField>("grad(DU)"))
    {
        gradDUPtr.set(new volTensorField(fvc::grad(DU)));
    }
    DU.correctBoundaryConditions();

    // Let the pointMesh know that the mesh has changed
    const_cast<pointBoundaryMesh&>
    (
        pointMesh::New(mesh).boundary()
#if FOAMEXTEND > 40
    ).updateMesh(*this);
#else
    ).updateMesh();
#endif

    // Force sync of procs after write
    returnReduce(bool(true), orOp<bool>());

    if (Pstream::master() || !onlyMasterWritesFields)
    {
        DU.write();
        //pointDU.write();
    }

    WarningIn(type() + "::updateBoundaryConditions(...)")
        << "Only DU boundary conditions are mapped! 'T' is not yet!" << endl
        << "Accumulated fields in DU boundary conditions e.g. totalDisp"
        << " will be reset to zero!" << endl;

    // 7. Let volToPoint interpolator know that the mesh has been updated
    // Note: it does not actually use the mapPolyMesh
    if (debug)
    {
        Info<< nl << "20. Updating volToPoint interpolator" << endl;
    }

    WarningIn(type() + "::updateBoundaryConditions(...)")
        << "The mesh smoother was not updated!" << endl;

    volToPointInterp().updateMesh(mpm);
}


void Foam::remeshWireFvMesh::copyZones
(
    const fvMesh& meshSource,
    fvMesh& meshTarget
) const
{
    List<pointZone*> newPointZones(meshSource.pointZones().size());
    List<faceZone*> newFaceZones(meshSource.faceZones().size());
    List<cellZone*> newCellZones(meshSource.cellZones().size());

    // Copy pointZones from meshSource and re-order consistently with
    // current (old) mesh
    if (meshTarget.pointZones().size() != meshSource.pointZones().size())
    {
        FatalErrorIn(type() + "::copyZones(...)")
            << "The number of pointZones in the current mesh and meshSource"
            << " are different!" << nl
            << meshSource.name() << ", nZones = "
            << meshSource.pointZones().size() << nl
            << meshTarget.name() << ", nZones = "
            << meshTarget.pointZones().size()
            << abort(FatalError);
    }

    forAll(meshTarget.pointZones(), pointZoneI)
    {
        const word& pzNameI = meshTarget.pointZones()[pointZoneI].name();

        const label newPointZoneID =
            meshSource.pointZones().findZoneID(pzNameI);

        if (newPointZoneID == -1)
        {
            FatalErrorIn(type() + "::copyZones(...)")
                << "pointZone " << pzNameI << " not found in meshSource!"
                << abort(FatalError);
        }

        newPointZones[pointZoneI] =
            new pointZone
            (
                pzNameI,
                meshSource.pointZones()[newPointZoneID],
                pointZoneI, // pointZoneID
                meshTarget.pointZones()
            );
    }

    // Copy faceZones from meshSource and re-order consistently with
    // current (old) mesh
    if (meshTarget.faceZones().size() != meshSource.faceZones().size())
    {
        FatalErrorIn(type() + "::copyZones(...)")
            << "The number of faceZones in the current mesh and meshSource"
            << " are different!" << nl
            << meshSource.name() << ", nZones = "
            << meshSource.faceZones().size() << nl
            << meshTarget.name() << ", nZones = "
            << meshTarget.faceZones().size()
            << abort(FatalError);
    }

    forAll(meshTarget.faceZones(), faceZoneI)
    {
        const word& fzNameI = meshTarget.faceZones()[faceZoneI].name();

        const label newFaceZoneID =
            meshSource.faceZones().findZoneID(fzNameI);

        if (newFaceZoneID == -1)
        {
            FatalErrorIn(type() + "::copyZones(...)")
                << "faceZone " << fzNameI << " not found in meshSource!"
                << abort(FatalError);
        }

        newFaceZones[faceZoneI] =
            new faceZone
            (
                fzNameI,
                meshSource.faceZones()[newFaceZoneID],
                meshSource.faceZones()[newFaceZoneID].flipMap(),
                faceZoneI, // faceZoneID
                meshTarget.faceZones()
            );
    }

    // Copy cellZones from meshSource and re-order consistently with
    // current (old) mesh
    if (meshTarget.cellZones().size() != meshSource.cellZones().size())
    {
        FatalErrorIn(type() + "::copyZones(...)")
            << "The number of cellZones in the current mesh and meshSource"
            << " are different!" << nl
            << meshSource.name() << ", nZones = "
            << meshSource.cellZones().size() << nl
            << meshTarget.name() << ", nZones = "
            << meshTarget.cellZones().size()
            << abort(FatalError);
    }

    forAll(meshTarget.cellZones(), cellZoneI)
    {
        const word& czNameI = meshTarget.cellZones()[cellZoneI].name();

        const label newCellZoneID =
            meshSource.cellZones().findZoneID(czNameI);

        if (newCellZoneID == -1)
        {
            FatalErrorIn(type() + "::copyZones(...)")
                << "cellZone " << czNameI << " not found in meshSource!"
                << abort(FatalError);
        }

        newCellZones[cellZoneI] =
            new cellZone
            (
                czNameI,
                meshSource.cellZones()[newCellZoneID],
                cellZoneI, // cellZoneID
                meshTarget.cellZones()
            );
    }

    // Remove all zones from current mesh
    meshTarget.removeZones();

    // Add zones to current mesh
    meshTarget.addZones(newPointZones, newFaceZones, newCellZones);
}


void Foam::remeshWireFvMesh::decomposeGlobalMesh
(
    const fvMesh& globalMesh,
    fvMesh& localMesh,
    const meshInfo& oldLocalMeshInfo
)
{
    // Approach:
    // 1. Use decomposePar to decompose the global mesh and fields on disk
    //        - issue: this will overwrite proc dirs
    //        - solution: mv proc dirs and mv contents back after
    // 2. Update proc mesh from disk
    // 3. Update proc fields from disk


    // 1. Use decomposePar to decompose the global mesh and fields on disk

    if (Pstream::master())
    {
        for (int procI = 0; procI < Pstream::nProcs(); procI++)
        {
            const fileName procDir =
                fileName
                (
                    globalMesh.time().path()/"processor" + Foam::name(procI)
                );

            // Remove the current time-step from the processor directory
            rmDir(procDir/localMesh.time().timeName());

            // Before decomposing, we will move the current processor
            // directories as we do not want to lose all the previous
            // time-step results
            cp(procDir, procDir + "_tmp");
        }
    }

     // Force sync
    returnReduce(bool(true), orOp<bool>());

    // decomposePar does not have a latestTime option so we will use our
    // modified version which does have this option
    if (Pstream::master())
    {
        runSystemCommandWithPOpen
        (
            "decomposeParWithLatestTimeOption -force -latestTime",
            "log.decomposeParWithLatestTimeOption_remesh_timeIndex"
           + Foam::name(localMesh.time().timeIndex()),
            type(),
            true
        );

        // decomposePar will place the mesh in the zero time directory, so we
        // will move it to the latest time-step
        for (int procI = 0; procI < Pstream::nProcs(); procI++)
        {
            const fileName procDir =
                fileName
                (
                    globalMesh.time().path()/"processor" + Foam::name(procI)
                );
            const fileName procTmpDir = fileName(procDir + "_tmp");

            // Move polyMesh to the latest time directory
            mv
            (
                procDir/"0/polyMesh",
                procDir/localMesh.time().timeName()/"polyMesh"
            );

            // Remove 0 directory
            if (isDir(procDir/"0"))
            {
                rmDir(procDir/"0");
            }

            // Remove constant directory
            if (isDir(procDir/"constant"))
            {
                rmDir(procDir/"constant");
            }

            // Replace constant with constant_tmp in the old proc directory
            rmDir(procTmpDir/word("constant"));
            mv(procTmpDir/word("constant_tmp"), procTmpDir/word("constant"));

            // Move the previous time-steps back into the procDir
            const fileNameList contents =
                readDir(procTmpDir, fileName::DIRECTORY, true);

            forAll(contents, i)
            {
                mv(procTmpDir/contents[i], procDir/contents[i]);
            }

            // Remove procTmpDir
            rmDir(procTmpDir);
        }
    }


    // Force sync
    returnReduce(bool(true), orOp<bool>());


    // 2. Update proc mesh from disk

    // Read mesh from disk
    fvMesh localMeshFromDisk
    (
        IOobject
        (
            polyMesh::defaultRegion,
            localMesh.time().timeName(),
            localMesh.time(),
            IOobject::MUST_READ
        )
    );

    // Copy localMeshFromDisk to localMesh
    labelList patchSizes(localMeshFromDisk.boundaryMesh().size(), 0);
    labelList patchStarts(localMeshFromDisk.boundaryMesh().size(), 0);
    forAll(patchSizes, patchI)
    {
        patchSizes[patchI] = localMeshFromDisk.boundaryMesh()[patchI].size();
        patchStarts[patchI] = localMeshFromDisk.boundaryMesh()[patchI].start();
    }

    localMesh.resetPrimitives
    (
        xferCopy(localMeshFromDisk.points()),
        xferCopy(localMeshFromDisk.faces()),
        xferCopy(localMeshFromDisk.faceOwner()),
        xferCopy(localMeshFromDisk.faceNeighbour()),
        patchSizes,
        patchStarts,
        true // valid boundary
    );

    // Copy the zones from localMeshFromDisk to meshTarget
    copyZones(localMeshFromDisk, localMesh);

    // Clear out demand driven data
    localMeshFromDisk.clearOut();
    localMesh.clearOut();

    // 3. Update proc fields from disk

    // Read fields
    autoPtr<geometricFieldContainer> localGeomFieldsFromDiskPtr;
    localGeomFieldsFromDiskPtr.set
    (
        new geometricFieldContainer(localMeshFromDisk)
    );

    // Replace fields for localMesh
    replaceMeshFields(localMesh, localMeshFromDisk, false);

    // Clear fields from disk as they are no longer needed
    localGeomFieldsFromDiskPtr.clear();

    // Update DU boundary conditions
    updateBoundaryConditionsAndInterp(localMesh, oldLocalMeshInfo, false);

    // Force sync
    returnReduce(bool(true), orOp<bool>());
}


void Foam::remeshWireFvMesh::replaceMeshFields
(
    const fvMesh& meshTarget,
    const fvMesh& meshSource,
    const bool writeFields
)
{
    if (debug)
    {
        Info<< "    Updating and writing volFields" << endl;
    }

    ReplaceOldMeshFieldsWithNewMeshFields
    <
        scalar, fvPatchField, volMesh, fvMesh
    >
    (
        meshTarget, meshSource, meshTarget, meshSource, writeFields
    );
    ReplaceOldMeshFieldsWithNewMeshFields
    <
        vector, fvPatchField, volMesh, fvMesh
    >
    (
        meshTarget, meshSource, meshTarget, meshSource, writeFields
    );
    ReplaceOldMeshFieldsWithNewMeshFields
    <
        tensor, fvPatchField, volMesh, fvMesh
    >
    (
        meshTarget, meshSource, meshTarget, meshSource, writeFields
    );
    ReplaceOldMeshFieldsWithNewMeshFields
    <
        symmTensor, fvPatchField, volMesh, fvMesh
    >
    (
        meshTarget, meshSource, meshTarget, meshSource, writeFields
    );
    ReplaceOldMeshFieldsWithNewMeshFields
    <
        diagTensor, fvPatchField, volMesh, fvMesh
    >
    (
        meshTarget, meshSource, meshTarget, meshSource, writeFields
    );
    ReplaceOldMeshFieldsWithNewMeshFields
    <
        sphericalTensor, fvPatchField, volMesh, fvMesh
    >
    (
        meshTarget, meshSource, meshTarget, meshSource, writeFields
    );

    if (debug)
    {
        Info<< "    Updating  and writing surfaceFields" << endl;
    }

    ReplaceOldMeshFieldsWithNewMeshFields
    <
        scalar, fvsPatchField, surfaceMesh, fvMesh
    >
    (
        meshTarget, meshSource, meshTarget, meshSource, writeFields
    );
    ReplaceOldMeshFieldsWithNewMeshFields
    <
        vector, fvsPatchField, surfaceMesh, fvMesh
    >
    (
        meshTarget, meshSource, meshTarget, meshSource, writeFields
    );
    ReplaceOldMeshFieldsWithNewMeshFields
    <
        tensor, fvsPatchField, surfaceMesh, fvMesh
    >
    (
        meshTarget, meshSource, meshTarget, meshSource, writeFields
    );
    ReplaceOldMeshFieldsWithNewMeshFields
    <
        symmTensor, fvsPatchField, surfaceMesh, fvMesh
    >
    (
        meshTarget, meshSource, meshTarget, meshSource, writeFields
    );
    ReplaceOldMeshFieldsWithNewMeshFields
    <
        diagTensor, fvsPatchField, surfaceMesh, fvMesh
    >
    (
        meshTarget, meshSource, meshTarget, meshSource, writeFields
    );
    ReplaceOldMeshFieldsWithNewMeshFields
    <
        sphericalTensor, fvsPatchField, surfaceMesh, fvMesh
    >
    (
        meshTarget, meshSource, meshTarget, meshSource, writeFields
    );

    if (debug)
    {
        Info<< "    Updating and writing pointFields" << endl;
    }

    const pointMesh& meshTargetPMesh = pointMesh::New(meshTarget);
    const pointMesh& meshSourcePMesh = pointMesh::New(meshSource);
    ReplaceOldMeshPointFieldsWithNewMeshPointFields<scalar>
    (
        meshTarget, meshSource, meshTargetPMesh, meshSourcePMesh, writeFields
    );
    ReplaceOldMeshPointFieldsWithNewMeshPointFields<vector>
    (
        meshTarget, meshSource, meshTargetPMesh, meshSourcePMesh, writeFields
    );
    ReplaceOldMeshPointFieldsWithNewMeshPointFields<tensor>
    (
        meshTarget, meshSource, meshTargetPMesh, meshSourcePMesh, writeFields
    );
    ReplaceOldMeshPointFieldsWithNewMeshPointFields<symmTensor>
    (
        meshTarget, meshSource, meshTargetPMesh, meshSourcePMesh, writeFields
    );
    ReplaceOldMeshPointFieldsWithNewMeshPointFields<diagTensor>
    (
        meshTarget, meshSource, meshTargetPMesh, meshSourcePMesh, writeFields
    );
    ReplaceOldMeshPointFieldsWithNewMeshPointFields<sphericalTensor>
    (
        meshTarget, meshSource, meshTargetPMesh, meshSourcePMesh, writeFields
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::remeshWireFvMesh::remeshWireFvMesh(const IOobject& io)
:
    layerAdditionFvMesh(io),
    remeshWire_(dict().lookup("remeshWire")),
    remeshFrequency_(readInt(remeshDict().lookup("remeshFrequency"))),
    keepIntermediateMeshesAndFields_
    (
        dict().lookupOrDefault<Switch>("keepIntermediateMeshesAndFields", false)
    ),
    geomFieldsPtr_()
{
    if (remeshFrequency_ < 1)
    {
        FatalErrorIn(type() + "::" + type() + "()")
            << "remeshFrequency must be an integer greater than 0"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::remeshWireFvMesh::~remeshWireFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::remeshWireFvMesh::update()
{
    layerAdditionFvMesh::update();

    if (!checkRemeshWire())
    {
        Info<< type() << ": wire is not re-meshed" << nl << endl;
        return false;
    }
    else
    {
        Info<< type() << ": wire will be re-meshed" << nl << endl;
    }

    // Parallel is WIP
    // if (Pstream::parRun())
    // {
    //     FatalErrorIn(type() + "::update()")
    //         << "For now, only serial runs are allowed!" << abort(FatalError);
    // }

    // Method
    // 1. Store info about original mesh before we make changes
    // 2. Create new wire mesh
    // 3. Map fields from newMesh and overwrite the mesh with the newMesh
    // 4. Clear out mesh replacer
    // 5. Update material indices fields
    // 6. Update boundary conditions and volToPoint interpolator

    // 1. Store info about original mesh before we make changes
    const meshInfo oldLocalMeshInfo(*this);

    if (Pstream::parRun())
    {
        // Make a copy of processor constant directory
        cp
        (
            fileName(time().path()/time().constant()),
            fileName(time().path()/time().constant() + "_tmp")
        );
    }

    // 2. Create the new mesh, where the wire cellZone mesh has been replaced
    replaceCellZoneMesh meshReplacer
    (
        *this,
        movingCellsName(),
        wireInletName(),
        wireInletFaceZoneName(),
        wireOutletName(),
        wireOutletFaceZoneName(),
        remeshDict()
    );

    // Store maxCellSize, if it was automatically calculated, to avoid the
    // maxCellSize progressively increasing with each remesh
    if (!remeshDict().found("maxCellSize"))
    {
        remeshDict().set("maxCellSize", meshReplacer.maxCellSize());
    }

    // 3. Map fields from newMesh and overwrite the mesh with the newMesh
    clearOut();
    writeMeshFields(*this);
    autoPtr<meshInfo> oldGlobalMeshInfoPtr;
    if (Pstream::parRun())
    {
        fvMesh& oldGlobalMesh = meshReplacer.oldGlobalMesh();
        oldGlobalMeshInfoPtr.set(new meshInfo(oldGlobalMesh));

        // Read all fields from disk and place them in the global mesh
        // objectRegistry
        if (geomFieldsPtr_.valid())
        {
            geomFieldsPtr_.clear();
        }
        chDir(time().path()/"..");
        oldGlobalMesh.solutionDict().subDict("solidMechanics").set
        (
            "nonLinear",
            nonLinearGeometry::nonLinearNames_
            [
                nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF
            ]
        );
        geomFieldsPtr_.set(new geometricFieldContainer(oldGlobalMesh));
        chDir(time().path());
    }
    else
    {
        oldGlobalMeshInfoPtr.set(new meshInfo(*this));
    }
    mapFieldsAndOverwriteTargetMesh
    (
        meshReplacer.newGlobalMesh(),
        meshReplacer.oldGlobalMesh() // old global mesh
    );

    // 4. Correct material index fields
    correctMaterialIndexFields(meshReplacer.oldGlobalMesh());

    // 5. Update boundary conditions e.g. the traction field in
    // solidTraction on the DU field will have the wrong size because we
    // did not map them
    updateBoundaryConditionsAndInterp
    (
        meshReplacer.oldGlobalMesh(), oldGlobalMeshInfoPtr(), true
    );

    if (Pstream::parRun())
    {
        // Decompose the global mesh and fields them and distribute them
        decomposeGlobalMesh
        (
            meshReplacer.oldGlobalMesh(), *this, oldLocalMeshInfo
        );
    }

    // 6. Clear global fields and global meshes
    geomFieldsPtr_.clear();
    meshReplacer.clearOut();

    // 7. Remove the intermediate meshing fields, unless asked otherwise
    if (!keepIntermediateMeshesAndFields_)
    {
        Info<< "Removing intermediate meshes and fields" << endl;
        rmDir(fileName(time().path()/time().timeName()));
    }

    return true;
}


// ************************************************************************* //
