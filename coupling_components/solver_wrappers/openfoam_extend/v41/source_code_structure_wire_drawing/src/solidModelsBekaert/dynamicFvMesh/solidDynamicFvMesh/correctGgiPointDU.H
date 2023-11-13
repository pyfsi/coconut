{
    // Correct ggi points
    // Lets final value be equal to average of master and slave
    forAll(DU.boundaryField(), patchI)
    {
        if
        (
            isA<materialGgiFvPatchVectorField>(DU.boundaryField()[patchI])
        )
        {
            const ggiPolyPatch& ggiPatch =
                refCast<const ggiPolyPatch>(this->boundaryMesh()[patchI]);

            if (ggiPatch.master())
            {
                // Slave point field
                label slavePatchID = ggiPatch.shadowIndex();
                label slaveZoneID = ggiPatch.shadow().zoneIndex();
                label masterZoneID = ggiPatch.zoneIndex();

                const labelList& slaveMeshPoints =
                    this->boundaryMesh()[slavePatchID].meshPoints();

                label globalSlaveZoneIndex =
                    findIndex(globalFaceZones(), slaveZoneID);

                Field<vector> slaveZonePointDU
                (
                    this->faceZones()[slaveZoneID]().nPoints(),
                    pTraits<vector>::zero
                );

                if (globalSlaveZoneIndex != -1)
                {
                    // global face zone

                    const labelList& curPointMap =
                        globalToLocalFaceZonePointMap()[globalSlaveZoneIndex];

                    const labelList& slaveZoneMeshPoints =
                        this->faceZones()[slaveZoneID]().meshPoints();

                    Field<vector> slaveZonePointDUGlobal
                    (
                        slaveZoneMeshPoints.size(),
                        pTraits<vector>::zero
                    );

                    //- Inter-proc points are shared by multiple procs
                    //  pointNumProc is the number of procs
                    //  which a point lies on
                    scalarField pointNumProcs(slaveZoneMeshPoints.size(), 0);

                    forAll(slaveZonePointDUGlobal, globalPointI)
                    {
                        label localPoint = curPointMap[globalPointI];

                        if (slaveZoneMeshPoints[localPoint] < this->nPoints())
                        {
                            label procPoint = slaveZoneMeshPoints[localPoint];

                            slaveZonePointDUGlobal[globalPointI] =
                                pointDUI[procPoint];

                            pointNumProcs[globalPointI] = 1;
                        }
                    }

                    if (Pstream::parRun())
                    {
                        reduce
                        (
                            slaveZonePointDUGlobal,
                            sumOp<Field<vector> >()
                        );
                        reduce(pointNumProcs, sumOp<scalarField>());

                        //- now average the displacement between all procs
                        slaveZonePointDUGlobal /= pointNumProcs;
                    }

                    forAll(slaveZonePointDU, globalPointI)
                    {
                        label localPoint = curPointMap[globalPointI];

                        slaveZonePointDU[localPoint] =
                            slaveZonePointDUGlobal[globalPointI];
                    }
                }
                else
                {
                    slaveZonePointDU =
                        Field<vector>
                        (
                            pointDUI,
                            this->faceZones()[slaveZoneID]().meshPoints()
                        );
                }

                Field<vector> slaveZonePointDUAtMaster =
                    ggiPatch.patchToPatch()
                   .slaveToMasterPointInterpolate(slaveZonePointDU);

                // Master point field

                const labelList& masterMeshPoints =
                    this->boundaryMesh()[patchI].meshPoints();

                Field<vector> slavePatchPointDUAtMaster
                (
                    this->boundaryMesh()[patchI].nPoints(),
                    pTraits<vector>::zero
                );

                forAll(slavePatchPointDUAtMaster, pointI)
                {
                    label curMeshPointID = masterMeshPoints[pointI];

                    label curMasterZonePointID =
                        this->faceZones()[masterZoneID]()
                       .whichPoint(curMeshPointID);

                    slavePatchPointDUAtMaster[pointI] =
                        slaveZonePointDUAtMaster[curMasterZonePointID];
                }

                Field<vector> masterPatchPointDU(pointDUI, masterMeshPoints);

                // New master point field
                Field<vector> newMasterPatchPointDU =
                    0.5*(masterPatchPointDU + slavePatchPointDUAtMaster);

                forAll(masterMeshPoints, pI)
                {
                    pointDUI[masterMeshPoints[pI]] = newMasterPatchPointDU[pI];
                }

                // Now correct slave side

                Field<vector> newMasterZonePointDU
                (
                    this->faceZones()[masterZoneID]().nPoints(),
                    pTraits<vector>::zero
                );

                label globalMasterZoneIndex =
                    findIndex(globalFaceZones(), masterZoneID);

                if (globalMasterZoneIndex != -1)
                {
                    // global face zone

                    const labelList& curPointMap =
                        globalToLocalFaceZonePointMap()[globalMasterZoneIndex];

                    const labelList& masterZoneMeshPoints =
                        this->faceZones()[masterZoneID]().meshPoints();

                    Field<vector> masterZonePointDUGlobal
                    (
                        masterZoneMeshPoints.size(),
                        pTraits<vector>::zero
                    );

                    //- Inter-proc points are shared by multiple procs
                    //  pointNumProc is the number of procs
                    //  which a point lies on
                    scalarField pointNumProcs(masterZoneMeshPoints.size(), 0);

                    forAll(masterZonePointDUGlobal, globalPointI)
                    {
                        label localPoint = curPointMap[globalPointI];

                        if (masterZoneMeshPoints[localPoint] < this->nPoints())
                        {
                            label procPoint = masterZoneMeshPoints[localPoint];

                            masterZonePointDUGlobal[globalPointI] =
                                pointDUI[procPoint];

                            pointNumProcs[globalPointI] = 1;
                        }
                    }

                    if (Pstream::parRun())
                    {
                        reduce
                        (
                            masterZonePointDUGlobal,
                            sumOp<Field<vector> >()
                        );
                        reduce(pointNumProcs, sumOp<scalarField>());

                        //- now average the displacement between all procs
                        masterZonePointDUGlobal /= pointNumProcs;
                    }

                    forAll(newMasterZonePointDU, globalPointI)
                    {
                        label localPoint = curPointMap[globalPointI];

                        newMasterZonePointDU[localPoint] =
                            masterZonePointDUGlobal[globalPointI];
                    }
                }
                else
                {
                    newMasterZonePointDU =
                        Field<vector>
                        (
                            pointDUI,
                            this->faceZones()[masterZoneID]().meshPoints()
                        );
                }

                // New slave point field
                Field<vector> newSlaveZonePointDU =
                    ggiPatch.patchToPatch()
                   .masterToSlavePointInterpolate(newMasterZonePointDU);

                Field<vector> newSlavePatchPointDU
                (
                    this->boundaryMesh()[slavePatchID].nPoints(),
                    pTraits<vector>:: zero
                );

                forAll(newSlavePatchPointDU, pointI)
                {
                    label curMeshPointID = slaveMeshPoints[pointI];

                    label curSlaveZonePointID =
                        this->faceZones()[slaveZoneID]()
                       .whichPoint(curMeshPointID);

                    newSlavePatchPointDU[pointI] =
                        newSlaveZonePointDU[curSlaveZonePointID];
                }

                forAll(slaveMeshPoints, pI)
                {
                    pointDUI[slaveMeshPoints[pI]] = newSlavePatchPointDU[pI];
                }
            }
        }
    }
}
