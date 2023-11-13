/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

Description
    Merge two meshes and their fields.

    This is just the mergeMeshes utility with the added functionality of merging
    the fields of the meshes.

    The mesh is overwritten.

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "mergePolyMesh.H"
#include "MergeVolFields.H"
#include "MergeSurfaceFields.H"

using namespace Foam;


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("noFields", "");
    argList::validOptions.insert("latestTime", "");
    argList::validOptions.insert("overwrite", "");
    argList::validOptions.insert("writeToDefaultRegion", "");

#   include "setRoots.H"
#   include "createTimes.H"

    const bool noFields = args.optionFound("noFields");
    const bool overwrite = args.optionFound("overwrite");
    const bool writeToDefaultRegion = args.optionFound("writeToDefaultRegion");

    if (args.optionFound("latestTime"))
    {
        const instantList masterTimes = runTimeMaster.times();
        runTimeMaster.setTime
        (
            masterTimes[masterTimes.size() - 1],
            masterTimes.size() - 1
        );
    }

    Info<< "Reading master mesh from time = " << runTimeMaster.timeName()
        << endl;

    Info<< "Create mesh\n" << endl;
    mergePolyMesh masterMesh
    (
        IOobject
        (
            masterRegion,
            runTimeMaster.timeName(),
            runTimeMaster
        )
    );

    const word oldInstance = masterMesh.pointsInstance();

    // Set time to last time step for meshToAdd
    instantList Times = runTimeToAdd.times();
    label endTime = Times.size() - 1;
    runTimeToAdd.setTime(Times[endTime], endTime);

    Info<< "Reading mesh to add for time = " << runTimeToAdd.timeName()
        << endl;

    Info<< "Create mesh\n" << endl;
    polyMesh meshToAdd
    (
        IOobject
        (
            addRegion,
            runTimeToAdd.timeName(),
            runTimeToAdd
        )
    );

    // Merge the meshes
    masterMesh.addMesh(meshToAdd);
    masterMesh.merge();

    //runTimeMaster++;

    if (overwrite)
    {
        masterMesh.setInstance(oldInstance);
    }
    else
    {
        runTimeMaster++;
    }

    if (writeToDefaultRegion)
    {
        masterMesh.rename(polyMesh::defaultRegion);
    }

    Info<< "Writing combined mesh to " << runTimeMaster.timeName() << endl;
    masterMesh.polyMesh::write();


    // Now we lookup the fields in the masterMesh and the meshToAdd

    fvMesh mesh(masterMesh);
    fvMesh slaveMesh(meshToAdd);

    if (!noFields)
    {
        // Search for list of source objects in both meshes
        IOobjectList objects(mesh, masterMesh.time().timeName());
        IOobjectList slaveObjects
        (
            slaveMesh, slaveMesh.time().timeName()
        );

        Info<< nl << "Merging volfields " << endl;

        // Merge volFields of slave mesh into the new mesh
        MergeVolFields<scalar>
        (
            objects, slaveObjects, mesh, slaveMesh
        );
        MergeVolFields<vector>
        (
            objects, slaveObjects, mesh, slaveMesh
        );
        MergeVolFields<tensor>
        (
            objects, slaveObjects, mesh, slaveMesh
        );
        MergeVolFields<symmTensor>
        (
            objects, slaveObjects, mesh, slaveMesh
        );

        Info<< "\nMerging surfacefields " << endl;

        // Merge surfaceFields of slave mesh into the new mesh
        MergeSurfaceFields<scalar>
        (
            objects, slaveObjects, mesh, slaveMesh
        );
        MergeSurfaceFields<vector>
        (
            objects, slaveObjects, mesh, slaveMesh
        );
        MergeSurfaceFields<tensor>
        (
            objects, slaveObjects, mesh, slaveMesh
        );
        MergeSurfaceFields<symmTensor>
        (
            objects, slaveObjects, mesh, slaveMesh
        );

        // Fix F/J/etc fields for rolling/drawing restart
        // Now done in solver and mechanicalLaws classes
        //#       include "correctFields.H"
    }

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
