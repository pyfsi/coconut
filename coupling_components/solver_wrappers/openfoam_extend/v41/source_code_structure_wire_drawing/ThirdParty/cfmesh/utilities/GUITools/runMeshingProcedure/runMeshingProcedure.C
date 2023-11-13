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
    runMeshingProcedure

Description
    Automated meshing procedure for given Bekaert mesh setup.

Author
    Mario Baburic, Creative Fields

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "clock.H"
#include "fileName.H"
#include "IOdictionary.H"
#include "foamTime.H"
#include "fvMeshSubset.H"
#include "polyMesh.H"
#include "dictionaryEntry.H"
#include "IFstream.H"
#include "IOobjectList.H"
#include "fileStat.H"
#include "helperFunctionsAux.H"

#include <string.h>
#include <cstdlib>
#include <unistd.h>
#include <fstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Optional argument
    // The utility cannot be run in parallel
    argList::noParallel();
    argList::validOptions.insert("overwrite", "");
    argList::validOptions.insert("positionRollers", "");
    argList::validOptions.insert("clean", "");
    argList::validOptions.insert("smooth", "word");
    argList args(argc, argv);

    // Overwrite wire / position rollers / delete submeshes
    bool overwriteWire = args.optionFound("overwrite");
    bool positionRollers = args.optionFound("positionRollers");
    bool deleteSubmeshes = args.optionFound("clean");
    bool smooth = args.optionFound("smooth");

    // Smooth the wire mesh and fields in the flow direction
    scalar newWireLength = -1;
    if (smooth)
    {
        newWireLength = strtod(args.option("smooth").c_str(),NULL);
        if (newWireLength <= 0)
        {
            FatalError
                << "Inappropriate -smooth argument. Try with positive number. "
                << abort(FatalError);
        }
    }
    else
    {
        FatalError
            << "Missing -smooth argument."
            << abort(FatalError);
    }

    // Checking case folder
    fileName passDir = cwd();
    help::checkFile(passDir/"0","Missing 0 directory");
    help::checkFile(passDir/"system","Missing system directory");
    help::checkFile(passDir/"system/controlDict","Missing controlDict file");
    help::checkFile(passDir/"constant","Missing constant directory");

    // Creating time object
    Info << nl << "Creating time object" << endl;
    Time runTime
    (
        Time::controlDictName,
        passDir.path(),
        passDir.name(),
        "system",
        "constant",
        false
    );

    // Checking pass
    IOdictionary controlDict
    (
        IOobject
        (
            "controlDict",
            "system",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    word passType = controlDict.lookup("application");
    if (passType != "drawingPass" && passType != "rollingPass")
    {
        FatalErrorIn("runMeshingProcedure")
            << nl << "    >>>> Invalid pass type. Aborting. "
            << nl << "    Pass type: " << passType
            << abort(FatalError);
    }
    else
    {
        Info << nl << "Pass type: " << passType << endl;
    }

    // Determine wire directory name
    fileName wireDir ="wire";
    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            "system",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    if( meshDict.found("wireDir") )
    {
        wireDir = word(meshDict.lookup("wireDir"));
        Info <<  nl << "Wire directory (meshDict): " << wireDir << endl;
    }
    else
    {
        Info << nl << "Wire directory (default): " << wireDir << endl;
    }

    // Determine roll directory name
    fileName rollDir ="roller";
    if( meshDict.found("rollDir") )
    {
        rollDir = word(meshDict.lookup("rollDir"));
        if (passType == "rollingPass")
            Info <<  nl << "Roller directory (meshDict): " << rollDir << endl;
    }
    else
    {
        if (passType == "rollingPass")
            Info << nl << "Roller directory (default): " << rollDir << endl;
    }

    // Determine die/casing directory name
    fileName dieDir ="die";
    fileName casingDir ="casing";
    if( meshDict.found("dieDir") )
    {
        dieDir = word(meshDict.lookup("dieDir"));

        if (passType == "drawingPass")
            Info <<  "Die directory (meshDict): " << dieDir << endl;
    }
    else
    {
        if (passType == "drawingPass")
            Info << "Die directory (default): " << dieDir << endl;
    }

    // Check for roll setup
    word rollSetup = "-";
    word singleRollerContactPatch = "-";
    fileName rollerName = "-";
    if (passType == "rollingPass")
    {
        if( meshDict.found("rollSetup") )
        {
            rollSetup = word(meshDict.lookup("rollSetup"));
            Info << nl <<  "Roll setup: " << rollSetup << endl;
        }
        else
        {
            FatalErrorIn("runMeshingProcedure")
                << nl << "    >>>> System command error. Aborting. "
                << nl << "    Roll setup not found."
                << abort(FatalError);
        }

        if( meshDict.found("singleRollerContactPatch") )
        {
            singleRollerContactPatch = word(meshDict.lookup("singleRollerContactPatch"));
            Info << nl <<  "Single roller contact patch: " << singleRollerContactPatch << endl;
        }
        else
        {
            FatalErrorIn("runMeshingProcedure")
                << nl << "    >>>> System command error. Aborting. "
                << nl << "    Single roller contact patch not found."
                << abort(FatalError);
        }
    }

    // Capture result from calling system commands. If non-zero value is
    // returned, give an error.
    int result(0);

    // Automated meshing procedure, depending on the pass type
    Info << nl << "Mesh preparation procedure:" << endl;

    // Cleaning submeshes
    Info << "   Cleaning sub-meshes ...   (Running application deleteSubMeshes)" << endl;
    result =
        help::runProcess
        (
            "deleteSubMeshes -all", "log.deleteSubMeshes_initial"
        );
    if (result > 0)
    {
        FatalErrorIn("runMeshingProcedure --> deleteSubMeshes")
            << nl << "    >>>> System command error. Aborting. "
            << abort(FatalError);
    }

    // Check for previous pass existence
    Info << "   Check for wire mesh in previous pass ...   (Running application prepareWireFromPrevious)" << endl;
    result =
        help::runProcess
        (
            std::string("prepareWireFromPrevious -remesh -wire " + wireDir + " -smooth " + Foam::name(newWireLength).c_str()).c_str(), "log.prepareWireFromPrevious"
        );
    if (result > 0)
    {
        FatalErrorIn("runMeshingProcedure --> prepareWireFromPrevious")
            << nl << "    >>>> System command error. Aborting. "
            << abort(FatalError);
    }

    // Case hygiene prior to meshing procedure
    result = help::runProcess("rm -rf 0/materials*");
    result = help::runProcess("rm -rf 0/cellToRegion*");

    // Check if wire mesh exists
    overwriteWire = true;
    if (fileStat(passDir/"constant"/wireDir/"polyMesh").isValid() != 0)
    {
        overwriteWire = false;
    }

    // Check for wire mesh existance
    if (!overwriteWire)
    {
        help::checkFile(passDir/"constant"/wireDir,"Missing wire directory. Remedy ==> CREATE/OVERWRITE WIRE MESH");
        help::checkFile(passDir/"constant"/wireDir/"polyMesh","Missing wire mesh. Remedy ==> CREATE/OVERWRITE WIRE MESH");
        Info << "   !!! Wire mesh found and copied from previous pass !!!" << endl;
    }

    // Generate wire mesh
    if (overwriteWire)
    {
        //help::runProcess("rm -rf constant"/wireDir);
        Info << "   Generating wire mesh ...   (Running application generateWireMesh)" << endl;
        result =
            help::runProcess
            (
                std::string("generateWireMesh -region " + wireDir).c_str(), "log.generateWireMesh"
            );
        if (result > 0)
        {
            FatalErrorIn("runMeshingProcedure --> generateWireMesh")
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }
    }

    if (passType == "drawingPass")
    {

        // Generate die mesh
        Info << "   Generating die mesh ...   (Running application generateDieMesh)" << endl;
        result =
            help::runProcess
            (
                std::string("generateDieMesh -region " + dieDir).c_str(), "log.generateDieMesh"
            );
        if (result > 0)
        {
            FatalErrorIn("runMeshingProcedure --> generateDieMesh")
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }

        // Merge sub-meshes
        Info << "   Merging wire and die meshes ...   (Running application mergeSubMeshes)" << endl;
        result =
            help::runProcess
            (
                std::string("mergeSubMeshes " + dieDir + " " + wireDir).c_str(), std::string("log.mergeSubMeshes_" + dieDir + "_" + wireDir).c_str()
            );
        if (result > 0)
        {
            FatalErrorIn("runMeshingProcedure --> mergeSubMeshes")
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }

    }
    else
    {

        // Generate roller mesh
        Info << "   Generating roller mesh ...   (Running application generate3DRollerMesh)" << endl;
        result =
            help::runProcess
            (
                std::string("generate3DRollerMesh -region " + rollDir).c_str(), "log.generate3DRollerMesh"
            );
        if (result > 0)
        {
            FatalErrorIn("runMeshingProcedure --> generate3DRollerMesh")
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }

        // Merge sub-meshes
        Info << "   Merging wire and roller meshes ...   (Running application mergeSubMeshes)" << endl;
        result =
            help::runProcess
            (
                std::string("mergeSubMeshes " + rollDir + " " + wireDir).c_str(), std::string("log.mergeSubMeshes_" + rollDir + "_" + wireDir).c_str()
            );
        if (result > 0)
        {
            FatalErrorIn("runMeshingProcedure --> mergeSubMeshes")
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }

        // Position rollers
        word wireContactPatch = "-";
        word rollerContactPatch = "-";
        if (positionRollers && rollSetup != "threeRolls")
        {

            // Extract sub-meshes
            Info << "   Extract sub-meshes ...   (Running application getSubMeshes)" << endl;
            result =
                help::runProcess
                (
                    std::string("getSubMeshes").c_str(),  "log.getSubMeshes"
                );
            if (result > 0)
            {
                FatalErrorIn("runMeshingProcedure --> getSubMeshes")
                    << nl << "    >>>> System command error. Aborting. "
                    << abort(FatalError);
            }

            // Read contact patches
            if (meshDict.subDict("rollingMillPatchNames").found("wireContactPatchName"))
            {
                wireContactPatch = word(meshDict.subDict("rollingMillPatchNames").lookup("wireContactPatchName"));

                //Info << "     Wire contact patch: " << wireContactPatch << endl;
            }
            else
            {
                FatalErrorIn("runMeshingProcedure")
                    << nl << "    >>>> System command error. Aborting. "
                    << nl << "    Wire contact patch not found."
                    << abort(FatalError);
            }

            if (meshDict.subDict("rollingMillPatchNames").found("rollerContactPatchName"))
            {
                rollerContactPatch = word(meshDict.subDict("rollingMillPatchNames").lookup("rollerContactPatchName"));

                //Info << "     Roller contact patch: " << rollerContactPatch << endl;
            }
            else
            {
                FatalErrorIn("runMeshingProcedure")
                    << nl << "    >>>> System command error. Aborting. "
                    << nl << "    Roller contact patch not found."
                    << abort(FatalError);
            }

            // Positon meshes depending on the roll setup
            if (rollSetup == "singleRoll")
            {
                if (meshDict.subDict("singleRollDict").found("topRoll"))
                {
                    rollerName = "topRoll";

                }
                else if (meshDict.subDict("singleRollDict").found("bottomRoll"))
                {
                    rollerName = "bottomRoll";

                }
                else if (meshDict.subDict("singleRollDict").found("rightRoll"))
                {
                    rollerName = "rightRoll";

                }
                else if (meshDict.subDict("singleRollDict").found("leftRoll"))
                {
                    rollerName = "leftRoll";

                }
                else
                {
                    FatalErrorIn("runMeshingProcedure")
                        << nl << "    >>>> System command error. Aborting. "
                        << nl << "    Single roll not found."
                        << abort(FatalError);
                }

                // Position sub-mesh
                Info << "   Position sub-mesh ...   (Running application positionRoller)" << endl;
                if(singleRollerContactPatch == "no")
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerName + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                else
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }

                if (result > 0)
                {
                    FatalErrorIn("runMeshingProcedure --> positionRoller")
                        << nl << "    >>>> System command error. Aborting. "
                        << abort(FatalError);
                }
            }
            else if (rollSetup == "edgeRolls")
            {
                // Position sub-mesh
                rollerName = "rightRoll";
                Info << "   Position rightRoll sub-mesh ...   (Running application positionRoller)" << endl;
                if(singleRollerContactPatch == "no")
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerName + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                else
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                if (result > 0)
                {
                    FatalErrorIn("runMeshingProcedure --> positionRoller")
                        << nl << "    >>>> System command error. Aborting. "
                        << abort(FatalError);
                }

                // Position sub-mesh
                rollerName = "leftRoll";
                Info << "   Position leftRoll sub-mesh ...   (Running application positionRoller)" << endl;
                if(singleRollerContactPatch == "no")
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerName + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                else
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                if (result > 0)
                {
                    FatalErrorIn("runMeshingProcedure --> positionRoller")
                        << nl << "    >>>> System command error. Aborting. "
                        << abort(FatalError);
                }
            }
            else if (rollSetup == "sideRolls")
            {
                // Position sub-mesh
                rollerName = "topRoll";
                Info << "   Position topRoll sub-mesh ...   (Running application positionRoller)" << endl;
                if(singleRollerContactPatch == "no")
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerName + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                else
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                if (result > 0)
                {
                    FatalErrorIn("runMeshingProcedure --> positionRoller")
                        << nl << "    >>>> System command error. Aborting. "
                        << abort(FatalError);
                }

                // Position sub-mesh
                rollerName = "bottomRoll";
                Info << "   Position bottomRoll sub-mesh ...   (Running application positionRoller)" << endl;
                if(singleRollerContactPatch == "no")
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerName + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                else
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                if (result > 0)
                {
                    FatalErrorIn("runMeshingProcedure --> positionRoller")
                        << nl << "    >>>> System command error. Aborting. "
                        << abort(FatalError);
                }

                // Position sub-mesh
                rollerName = "rightRoll";
                Info << "   Position rightRoll sub-mesh ...   (Running application positionRoller)" << endl;
                if(singleRollerContactPatch == "no")
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerName + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                else
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                if (result > 0)
                {
                    FatalErrorIn("runMeshingProcedure --> positionRoller")
                        << nl << "    >>>> System command error. Aborting. "
                        << abort(FatalError);
                }

                // Position sub-mesh
                rollerName = "leftRoll";
                Info << "   Position leftRoll sub-mesh ...   (Running application positionRoller)" << endl;
                if(singleRollerContactPatch == "no")
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerName + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                else
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                if (result > 0)
                {
                    FatalErrorIn("runMeshingProcedure --> positionRoller")
                        << nl << "    >>>> System command error. Aborting. "
                        << abort(FatalError);
                }
            }
            else if (rollSetup == "twoHighRolls")
            {
                // Position sub-mesh
                rollerName = "topRoll";
                Info << "   Position topRoll sub-mesh ...   (Running application positionRoller)" << endl;
                if(singleRollerContactPatch == "no")
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerName + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                else
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                if (result > 0)
                {
                    FatalErrorIn("runMeshingProcedure --> positionRoller")
                        << nl << "    >>>> System command error. Aborting. "
                        << abort(FatalError);
                }

                // Position sub-mesh
                rollerName = "bottomRoll";
                Info << "   Position bottomRoll sub-mesh ...   (Running application positionRoller)" << endl;
                if(singleRollerContactPatch == "no")
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerName + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                else
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                if (result > 0)
                {
                    FatalErrorIn("runMeshingProcedure --> positionRoller")
                        << nl << "    >>>> System command error. Aborting. "
                        << abort(FatalError);
                }
            }
            else if (rollSetup == "threeRolls")
            {
                // Position sub-mesh
                rollerName = "topRoll";
                Info << "   Position topRoll sub-mesh ...   (Running application positionRoller)" << endl;
                if(singleRollerContactPatch == "no")
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerName + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                else
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                if (result > 0)
                {
                    FatalErrorIn("runMeshingProcedure --> positionRoller")
                        << nl << "    >>>> System command error. Aborting. "
                        << abort(FatalError);
                }

                // Position sub-mesh
                rollerName = "bottomRightRoll";
                Info << "   Position bottomRightRoll sub-mesh ...   (Running application positionRoller)" << endl;
                if(singleRollerContactPatch == "no")
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerName + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                else
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                if (result > 0)
                {
                    FatalErrorIn("runMeshingProcedure --> positionRoller")
                        << nl << "    >>>> System command error. Aborting. "
                        << abort(FatalError);
                }

                // Position sub-mesh
                rollerName = "bottomLeftRoll";
                Info << "   Position bottomLeftRoll sub-mesh ...   (Running application positionRoller)" << endl;
                if(singleRollerContactPatch == "no")
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerName + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                else
                {
                    result =
                        help::runProcess
                        (
                            std::string("positionRoller " + rollerName + " " + wireDir + " " + rollerContactPatch + " " + wireContactPatch).c_str(), std::string("log.positionRoller_" + rollerName).c_str()
                        );
                }
                if (result > 0)
                {
                    FatalErrorIn("runMeshingProcedure --> positionRoller")
                        << nl << "    >>>> System command error. Aborting. "
                        << abort(FatalError);
                }
            }
            else
            {
                FatalErrorIn("runMeshingProcedure")
                    << nl << "    >>>> System command error. Aborting. "
                    << nl << "    Unrecognised roll setup: " << rollSetup
                    << abort(FatalError);
            }

            // Merge all sub-meshes
            if (overwriteWire)
            {
                Info << "   Merging all sub-meshes into main mesh ...   (Running application mergeSubMeshesAllPass)" << endl;
                result =
                    help::runProcess
                    (
                        std::string("mergeSubMeshesAllPass").c_str(), "log.mergeSubMeshesAllPass"
                    );
            }
            else
            {

                Info << "   Merging roller sub-meshes ...   (Running application mergeSubMeshesAllPass)" << endl;
                result =
                    help::runProcess
                    (
                        std::string("mergeSubMeshesAllPass -exclude").c_str(), "log.mergeSubMeshesAllPass"
                    );
            }
            if (result > 0)
            {
                FatalErrorIn("runMeshingProcedure --> mergeSubMeshesAllPass")
                    << nl << "    >>>> System command error. Aborting. "
                    << abort(FatalError);
            }

        }
        else if (positionRollers && rollSetup == "threeRolls")
        {
            Info << "   !!! Positioning skipped for threeRolls setup !!!" << endl;
        }

    }

    // Case hygiene
    result = help::runProcess("rm -rf 0/materials*");
    result = help::runProcess("rm -rf 0/cellToRegion*");

    // Merge roller and wire sub-meshes and fields
    if (!overwriteWire)
    {
        const fileName source = runTime.path()/"source";
        const fileName target = runTime.path()/"target";

        // Create source directory
        result = help::runProcess("rm -rf target");
        mkDir("target");
        chDir(target);
        mkDir("constant");
        mkDir("0");
        cp(runTime.path()/"constant"/"polyMesh",target/"constant"/"polyMesh");
        cp(runTime.path()/"system",target/"system");
        chDir(passDir);

        // Merging
        help::checkFile(source);
        help::checkFile(target);
        Info << "   Merging roller and wire meshes and fields ...   (Running application mergeMeshesAndFields)" << endl;
        result =
            help::runProcess
            (
                std::string("mergeMeshesAndFields . target . source").c_str(), "log.mergeMeshesAndFields"
            );
        if (result > 0)
        {
            FatalErrorIn("runMeshingProcedure --> mergeMeshesAndFields")
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }

        // Copy merged mesh and fields to constant and 0
        result = help::runProcess("rm -rf constant/polyMesh");
        result = help::runProcess("rm -rf 0");
        cp(target/"constant"/"polyMesh",passDir/"constant"/"polyMesh");
        cp(target/"0",passDir/"0");

        // Clean source and target dirs
        result = help::runProcess("rm -rf source");
        result = help::runProcess("rm -rf target");
    }

    // Renumbering mesh
    Info << "   Renumbering main mesh ...   (Running application renumberMesh)" << endl;
    result =
        help::runProcess
        (
            std::string("renumberMesh -overwrite").c_str(), "log.renumberMesh"
        );
    if (result > 0)
    {
        FatalErrorIn("runMeshingProcedure --> renumberMesh")
            << nl << "    >>>> System command error. Aborting. "
            << abort(FatalError);
    }

    // Writing face zones
    Info << "   Sets to zones ...   (Running application setsToZones)" << endl;
    result =
        help::runProcess
        (
            std::string("setsToZones").c_str(), "log.setsToZones"
        );
    if (result > 0)
    {
        FatalErrorIn("runMeshingProcedure --> setsToZones")
            << nl << "    >>>> System command error. Aborting. "
            << abort(FatalError);
    }

    // Write material fields
    Info << "   Writing material fields ...   (Running application setMatFromCellZones)" << endl;
    result =
        help::runProcess
        (
            std::string("setMatFromCellZones").c_str(), "log.setMatFromCellZones"
        );
    if (result > 0)
    {
        FatalErrorIn("runMeshingProcedure --> setMatFromCellZones")
            << nl << "    >>>> System command error. Aborting. "
            << abort(FatalError);
    }

    Info << "Done" << endl;

    // Cleaning sub-mesh directories
    if (deleteSubmeshes)
    {
        Info << nl << "Cleaning sub-meshes: " << endl;
        Info << nl << "   Running application deleteSubMeshes" << endl;
        result =
            help::runProcess
            (
                std::string("deleteSubMeshes -all").c_str(), "log.deleteSubMeshes_final"
            );
        if (result > 0)
        {
            FatalErrorIn("runMeshingProcedure --> deleteSubMeshes")
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }
        Info << "Done" << endl;
    }

    // Print out date and time
    Info << nl << nl << nl << "Running automated Beakert meshing procedure completed on " << clock::date() << " at "
         << clock::clockTime() << endl;

    Info << "END" << nl << endl;

    return 0;
}


// ************************************************************************* //
