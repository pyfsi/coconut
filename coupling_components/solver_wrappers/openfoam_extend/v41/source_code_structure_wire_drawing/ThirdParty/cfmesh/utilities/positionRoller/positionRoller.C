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
    positionRoller

Description
    Moving roller to touch wire surface.

Author
    Philip Cardiff UCD
    Zeljko Tukovic FSB/UCD

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "boolList.H"
#include "emptyPolyPatch.H"
#include "twoDPointCorrector.H"
#include "pointFields.H"
#include "volFields.H"
#include "standAlonePatch.H"
#include "patchToPatchInterpolation.H"
#include "helperFunctionsAux.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    argList::validArgs.append("roller");
    argList::validArgs.append("wire");
    argList::validArgs.append("rollerContactPatchName");
    argList::validArgs.append("wireContactPatchName");
    argList args(argc, argv);

    if( args.additionalArgs().size() < 4 )
    {
        FatalError
            << "At least two regions should be specified"
            << abort(FatalError);
    }

    const word rollerName = args.additionalArgs()[0];
    const word wireName = args.additionalArgs()[1];
    const word rollerContactPatchName = args.additionalArgs()[2];
    const word wireContactPatchName = args.additionalArgs()[3];

    if
    (
        !Switch(rollerName == "topRoll") &&
        !Switch(rollerName == "bottomRoll") &&
        !Switch(rollerName == "rightRoll") &&
        !Switch(rollerName == "leftRoll")
    )
    {
        FatalErrorIn
        (
            "positionRoller"
        )   << "roller " << rollerName << " is unknown. It should be "
            "topRoll, bottomRoll, leftRoll or rightRoll" << abort(FatalError);
    }


    Info << nl << "Positioning roller " << rollerName << " to touch "
         << wireName << nl << endl;

#   include "createTime.H"

    Foam::fvMesh wireMesh
    (
        Foam::IOobject
        (
            wireName,
            "constant"/wireName,
            runTime
        )
    );

    const label wirePatchID =
        wireMesh.boundaryMesh().findPatchID(wireContactPatchName);

    if (wirePatchID == -1)
    {
        FatalErrorIn
        (
            "positionRoller"
        )   << "wire contact patch not found!" << abort(FatalError);
    }

    // Create a standAlonePatch for the wire surface
    standAlonePatch wirePatch
    (
        wireMesh.boundaryMesh()[wirePatchID].localFaces(),
        wireMesh.boundaryMesh()[wirePatchID].localPoints()
    );

    // Define the positioning tolerance for moving the roller
    // We will define it as a fraction of the smallest face on the
    // wire patch
    const scalar positioningTol =
        1e-6*Foam::sqrt
        (
            min(mag(wireMesh.boundaryMesh()[wirePatchID].faceAreas()))
        );

    // We will iteratively move the roller until it is "just touching"
    // the wire surface.
    // We will do this by calculting the distances from the roller surface
    // to the wire surface and then moving the roller away from the wire
    // (positive y direction for the topRoll, negative y for the bottom
    // roll, etc.)

    scalar rollerMinY = GREAT;
    scalar rollerMaxY = -GREAT;
    scalar rollerMinZ = GREAT;
    scalar rollerMaxZ = -GREAT;

    vector translation = vector::zero;
    const int maxPositioningCorrectors = 10;
    int i = 0;

    do
    {
        Info<< "            Iteration " << (i + 1) << endl;

        // Reset translation to zero
        translation = vector::zero;

        // Read roller mesh from dict
        // This will re-read the mesh from disk if it has changed
        Foam::fvMesh rollerMesh
        (
            Foam::IOobject
            (
                rollerName,
                "constant"/rollerName,
                runTime
            )
        );

        // Lookup roller contact patch index
        const label rollerPatchID =
            rollerMesh.boundaryMesh().findPatchID
            (
                word(rollerContactPatchName)
            );

        if (rollerPatchID == -1)
        {
            FatalErrorIn
            (
                "positionRoller"
            )   << "roller contact patch not found!" << abort(FatalError);
        }

        // Create a standAlonePatch for the roller surface
        standAlonePatch rollerPatch
        (
            rollerMesh.boundaryMesh()[rollerPatchID].localFaces(),
            rollerMesh.boundaryMesh()[rollerPatchID].localPoints()
        );

        // Create patchToPatch object to calculate point distances between
        // the wire and roller patches
        //IMPORTANT: VISIBLE algorithm option disabled. CHECKS STILL NEEDED ON PROFILED ROLLS!!!
        PatchToPatchInterpolation<standAlonePatch, standAlonePatch>
            patchToPatch
            (
                rollerPatch,
                wirePatch
                //intersection::VISIBLE
            );

        // Store all rollerPatch boundary points for calculating the roller
        // position. This is needed to determine the displacement needed to
        // reach the required roll gap.
        const vectorField boundaryPoints = rollerPatch.localPoints();

        // Calculation of point distances from the roller patch to the wire
        // patch and convert to vector distances by multiply by the point
        // normals
        const vectorField wirePointDist =
            patchToPatch.pointDistanceToIntersection()*
            wirePatch.pointNormals();


        scalar translationY = 0;
        scalar translationZ = 0;

        vector rollerPosition(vector::zero);

        // Find the most negative distance (when there is no crossover this
        // is the least positive distance)
        if (rollerName == "topRoll")
        {
            scalar minTranslationY = GREAT;

            // Find minimum distance between wire and roller
            forAll(wirePointDist, pI)
            {
                if (mag(wirePointDist[pI]) < 1e14)
                {
                    minTranslationY =
                        min(minTranslationY, wirePointDist[pI].y());
                }
            }

            translationY = -minTranslationY;

            // Get the absolute position of the roller
            forAll(boundaryPoints, pI)
            {
                if (boundaryPoints[pI].y() < rollerMinY)
                {
                    rollerMinY = boundaryPoints[pI].y();
                }
            }
            rollerPosition = vector(0, rollerMinY, 0);

        }
        else if (rollerName == "bottomRoll")
        {
            scalar maxTranslationY = -GREAT;

            forAll(wirePointDist, pI)
            {
                if (mag(wirePointDist[pI]) < 1e14)
                {
                    maxTranslationY =
                        max(maxTranslationY, wirePointDist[pI].y());
                }
            }

            translationY = -maxTranslationY;

            forAll(boundaryPoints, pI)
            {
                if (boundaryPoints[pI].y() > rollerMaxY)
                {
                    rollerMaxY = boundaryPoints[pI].y();
                }
            }
            rollerPosition = vector(0, rollerMaxY, 0);
        }
        else if (rollerName == "rightRoll")
        {
            scalar minTranslationZ = GREAT;

            forAll(wirePointDist, pI)
            {
                if (mag(wirePointDist[pI]) < 1e14)
                {
                    minTranslationZ =
                        min(minTranslationZ, wirePointDist[pI].z());
                }
            }

            translationZ = -minTranslationZ;

            forAll(boundaryPoints, pI)
            {
                if (boundaryPoints[pI].z() < rollerMinZ)
                {
                    rollerMinZ = boundaryPoints[pI].z();
                }
            }
            rollerPosition = vector(0, 0, rollerMinZ);
        }
        else if (rollerName == "leftRoll")
        {
            scalar maxTranslationZ = -GREAT;

            forAll(wirePointDist, pI)
            {
                if (mag(wirePointDist[pI]) < 1e14)
                {
                    maxTranslationZ =
                        max(maxTranslationZ, wirePointDist[pI].z());
                }
            }

            translationZ = -maxTranslationZ;

            forAll(boundaryPoints, pI)
            {
                if (boundaryPoints[pI].z() > rollerMaxZ)
                {
                    rollerMaxZ = boundaryPoints[pI].z();
                }
            }
            rollerPosition = vector(0, 0, rollerMaxZ);
        }
        else
        {
            FatalErrorIn
            (
                "positionRoller"
            )   << "Positioning not implemented for roller type "
                << rollerName
                << abort(FatalError);
        }

        // Calculate translation vector
        translation = vector(0, translationY, translationZ);

        Info<< "            Translating by: " << translation << endl;

        // Move the roller mesh
        int result = help::runProcess
        (
            std::string
            (
                "transformPoints -region " + rollerName
              + " -translate "
              + "\"("
              + Foam::name(translation.x()) + " "
              + Foam::name(translation.y()) + " "
              + Foam::name(translation.z())
              + ")\""
            ).c_str(),
            std::string("log.transformPoints_" + rollerName).c_str()
        );
        if (result > 0)
        {
            FatalErrorIn
                (
                    "positionRoller"
                )
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }

    }
    while
    (
        i++ < maxPositioningCorrectors &&
        mag(translation) > positioningTol
    );

    Info<< "            Roller positioned to a tolerance of "
        << translation << endl;

    return 0;
}



// ************************************************************************* //
