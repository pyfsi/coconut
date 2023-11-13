/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    scalePointsRadially

Description
    Scale points radially that intersect the two given patches.

    This utility was intends to be used to generate an initial contraction/notch
    in a cylindrical tensile test.

    It is assumed that the x axis is the cylinder axis.

    The user must supply two patch names and also the percentage decrease in
    radius, e.g. for 1%

        scalePointsRadially patch1 patch2 0.01

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("patchName1");
    argList::validArgs.append("patchName1");
    argList::validArgs.append("reductionFactor");
    argList::validOptions.insert("latestTime", "");

#   include "setRootCase.H"
#   include "createTime.H"

    runTime.functionObjects().off();

    if (args.optionFound("latestTime"))
    {
        instantList Times = runTime.times();
        runTime.setTime(Times[Times.size() - 1], Times.size() - 1);
    }

#   include "createMesh.H"

    // Store old point instance to overwrite mesh
    const word oldInstance = mesh.pointsInstance();

    // Read patch names

    const word patchName1(args.additionalArgs()[0]);
    const label patchID1 = mesh.boundaryMesh().findPatchID(patchName1);
    if (patchID1 == -1)
    {
        FatalError
            << "Patch " << patchName1 << " not found!" << abort(FatalError);
    }

    const word patchName2(args.additionalArgs()[1]);
    const label patchID2 = mesh.boundaryMesh().findPatchID(patchName2);
    if (patchID2 == -1)
    {
        FatalError
            << "Patch " << patchName2 << " not found!" << abort(FatalError);
    }

    // Read percentage change in radius (decrease is assummed positive)
    const scalar reductionFactor
    (
        readScalar(IStringStream(args.additionalArgs()[2])())
    );
    if (reductionFactor >= 1.0)
    {
        FatalError
            << "The reductionFactor should be less than 1.0!"
            << abort(FatalError);
    }

    // Find points and scale them

    pointField newPoints = mesh.points();

    const labelList& meshPoints = mesh.boundaryMesh()[patchID1].meshPoints();
    const labelListList& pointFaces = mesh.pointFaces();

    forAll(meshPoints, meshPointI)
    {
        const label pointID = meshPoints[meshPointI];

        const labelList& curPointFaces = pointFaces[pointID];

        // Find faces on the boundary and check if they are on patch 2
        forAll(curPointFaces, pfI)
        {
            const label faceID = curPointFaces[pfI];

            if (!mesh.isInternalFace(faceID))
            {
                if (mesh.boundaryMesh().whichPatch(faceID) == patchID2)
                {
                    // This point should be scaled
                    const vector curPoint = newPoints[pointID];

                    // Vector from axis to the point
                    const vector d = vector(0, curPoint.y(), curPoint.z());

                    // Current radius from x axis
                    const scalar r = mag(d);

                    // New radius
                    const scalar newR = (1.0 - reductionFactor)*r;

                    // New d vector
                    const vector newD = newR*d/r;

                    // New point
                    newPoints[pointID].y() = newD.y();
                    newPoints[pointID].z() = newD.z();

                    break;
                }
            }
        }
    }

    // Move the mesh
    mesh.movePoints(newPoints);

    // Write the mesh
    mesh.setInstance(oldInstance);
    Info<< "Writing the mesh to " << oldInstance << endl;
    mesh.write();


    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
