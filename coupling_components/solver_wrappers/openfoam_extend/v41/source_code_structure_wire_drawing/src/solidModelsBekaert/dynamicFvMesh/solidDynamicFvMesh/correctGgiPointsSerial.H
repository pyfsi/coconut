{
    // Correct ggi points
    // Lets final value be equal to average of master and slave
    forAll(mesh().boundaryMesh(), patchI)
    {
        if ( isA<ggiPolyPatch>(mesh().boundaryMesh()[patchI]) )
        {
            const ggiPolyPatch& ggiPatch =
                refCast<const ggiPolyPatch>(mesh().boundaryMesh()[patchI]);

            if (ggiPatch.master())
            {
                // Slave point field
                label slavePatchID = ggiPatch.shadowIndex();

                const labelList& slaveMeshPoints = 
                    mesh().boundaryMesh()[slavePatchID].meshPoints();

                Field<Type> slavePatchPf(pfI, slaveMeshPoints);

                Field<Type> slavePatchPfAtMaster =
                    ggiPatch.patchToPatch()
                   .slaveToMasterPointInterpolate(slavePatchPf);

                // Master point field

                const labelList& masterMeshPoints = 
                    mesh().boundaryMesh()[patchI].meshPoints();
                  
                Field<Type> masterPatchPf(pfI, masterMeshPoints);

                // New master point field
                Field<Type> newMasterPatchPf =
                    0.5*(masterPatchPf + slavePatchPfAtMaster);

                forAll(masterMeshPoints, pI)
                {
                    pfI[masterMeshPoints[pI]] = newMasterPatchPf[pI];
                }

                // New slave point field
                Field<Type> newSlavePatchPf =
                    ggiPatch.patchToPatch()
                   .masterToSlavePointInterpolate(newMasterPatchPf);

                forAll(slaveMeshPoints, pI)
                {
                    pfI[slaveMeshPoints[pI]] = newSlavePatchPf[pI];
                }
            } 
        }
    }
}
