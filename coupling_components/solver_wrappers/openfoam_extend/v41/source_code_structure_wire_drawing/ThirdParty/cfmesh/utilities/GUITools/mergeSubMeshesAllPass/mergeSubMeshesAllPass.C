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
    mergeSubMeshesAllPass

Description
    Merge all pass sub-meshes.

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


#include <cstdlib>

using namespace Foam;

void mergeSubMeshes(const fileName&, const fileName&, const fileName&);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Optional argument
    // The utility cannot be run in parallel
    argList::noParallel();
    argList::validOptions.insert("clean", "");
    argList::validOptions.insert("exclude", "");
    argList args(argc, argv);

    // Cleaning submeshes argument
    bool deleteSubmeshes = args.optionFound("clean");
    bool excludeWire = args.optionFound("exclude");

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
        FatalErrorIn("mergeSubMeshesAllPass")
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

    // Start merging procedure, depending on the pass type
    if (passType == "drawingPass")
    {
        mergeSubMeshes(passDir,dieDir,wireDir);
        if (fileStat(passDir/"constant"/casingDir).isValid() != 0)
        {
            mergeSubMeshes(passDir,casingDir,"polyMeshMerge");
        }
    }
    else
    {
        // Roll setup
        word rollSetup = "-";
        fileName rollerName = "-";
        fileName merged = "-";
        if( meshDict.found("rollSetup") )
        {
            rollSetup = word(meshDict.lookup("rollSetup"));
            Info << nl <<  "Roll setup: " << rollSetup << endl;
        }
        else
        {
            FatalErrorIn("mergeSubMeshesAllPass")
                << nl << "    >>>> System command error. Aborting. "
                << nl << "    Roll setup not found."
                << abort(FatalError);
        }

        // Merging meshes depending on the roll setup
        if (rollSetup == "singleRoll")
        {
            if (meshDict.subDict("singleRollDict").found("topRoll"))
            {
                rollerName = "topRoll";
                if (!excludeWire)
                    mergeSubMeshes(passDir,rollerName,wireDir);

            }
            else if (meshDict.subDict("singleRollDict").found("bottomRoll"))
            {
                rollerName = "bottomRoll";
                if (!excludeWire)
                    mergeSubMeshes(passDir,rollerName,wireDir);
            }
            else if (meshDict.subDict("singleRollDict").found("rightRoll"))
            {
                rollerName = "rightRoll";
                if (!excludeWire)
                    mergeSubMeshes(passDir,rollerName,wireDir);
            }
            else if (meshDict.subDict("singleRollDict").found("leftRoll"))
            {
                rollerName = "leftRoll";
                if (!excludeWire)
                    mergeSubMeshes(passDir,rollerName,wireDir);
            }
            else
            {
                FatalErrorIn("mergeSubMeshesAllPass")
                    << nl << "    >>>> System command error. Aborting. "
                    << nl << "    Single roll not found."
                    << abort(FatalError);
            }
        }
        else if (rollSetup == "edgeRolls")
        {
            mergeSubMeshes(passDir,"rightRoll","leftRoll");
            if (!excludeWire)
                mergeSubMeshes(passDir,wireDir,"polyMeshMerge");

        }
        else if (rollSetup == "sideRolls")
        {
            mergeSubMeshes(passDir,"rightRoll","leftRoll");
            mergeSubMeshes(passDir,"topRoll","polyMeshMerge");
            mergeSubMeshes(passDir,"bottomRoll","polyMeshMerge");
            if (!excludeWire)
                mergeSubMeshes(passDir,wireDir,"polyMeshMerge");
        }
        else if (rollSetup == "twoHighRolls")
        {
            mergeSubMeshes(passDir,"topRoll","bottomRoll");
            if (!excludeWire)
                mergeSubMeshes(passDir,wireDir,"polyMeshMerge");
        }
        else if (rollSetup == "threeRolls")
        {
            mergeSubMeshes(passDir,"topRoll","bottomRightRoll");
            mergeSubMeshes(passDir,"bottomLeftRoll","polyMeshMerge");
            if (!excludeWire)
                mergeSubMeshes(passDir,wireDir,"polyMeshMerge");
        }
        else
        {
            FatalErrorIn("mergeSubMeshesAllPass")
                << nl << "    >>>> System command error. Aborting. "
                << nl << "    Unrecognised roll setup: " << rollSetup
                << abort(FatalError);
        }
    }

    int result(0);

    // Case hygiene
    result = help::runProcess("rm -rf 0/materials*");
    result = help::runProcess("rm -rf 0/cellToRegion*");

    // Additional mesh manipualtions
    if(false)
    {

        // Renumbering mesh
        Info << nl << "Renumbering main mesh ...   (Running application renumberMesh)" << endl;
        result =
            help::runProcess
            (
                "renumberMesh -overwrite", "log.renumberMesh"
            );
        if (result > 0)
        {
            FatalErrorIn("mergeSubMeshesAllPass --> renumberMesh")
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }

        // Writing face zones
        Info << "Sets to zones ...   (Running application setsToZones)" << endl;
        result =
            help::runProcess
            (
                "setsToZones", "log.setsToZones"
            );
        if (result > 0)
        {
            FatalErrorIn("mergeSubMeshesAllPass --> setsToZones")
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }

        // Write material fields
        Info << "Writing material fields ...   (Running application setMatFromCellZones)" << endl;
        result =
            help::runProcess
            (
                "setMatFromCellZones", "log.setMatFromCellZones"
            );
        if (result > 0)
        {
            FatalErrorIn("mergeSubMeshesAllPass --> setMatFromCellZones")
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }

    }

    // Cleaning sub-mesh directories
    if (deleteSubmeshes)
    {
        Info << nl << "Cleaning submeshes (constant directory): " << endl;
        Info << "   Running application deleteSubMeshes" << endl;
        result =
            help::runProcess
            (
                "deleteSubMeshes -constant", "log.deleteSubMeshes"
            );
        if (result > 0)
        {
            FatalErrorIn("mergeSubMeshesAllPass --> deleteSubMeshes")
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }
        Info << "Done" << endl;
    }

    // Print out date and time
    Info << nl << nl << nl << "Merging all rolling pass sub-meshes completed on " << clock::date() << " at "
         << clock::clockTime() << endl;

    Info << "END" << nl << endl;

    return 0;
}

void mergeSubMeshes(const fileName& dir, const fileName& file1, const fileName& file2)
{
    int result(0);

    help::checkFile(dir/"constant"/file1, "Missing submesh. Remedy ==> CREATE/EXTRACT SUBMESHES");
    help::checkFile(dir/"constant"/file2, "Missing submesh. Remedy ==> CREATE/EXTRACT SUBMESHES");

    Info << nl << "Merging " << file1 << " & " << file2 << " ..." << endl;

    // Cleaning constant/polyMesh
    Info << "   Cleaning constant/polyMesh" << endl;
    result = help::runProcess("rm -rf constant/polyMesh");

    // Merge
    Info << "   Running application mergeSubMeshes" << endl;
    result =
        help::runProcess
        (
            std::string("mergeSubMeshes " + file1 + " " + file2).c_str(), std::string("log.mergeSubMeshes_" + file1 + "_" + file2).c_str()
        );
    if (result > 0)
    {
        FatalErrorIn("mergeSubMeshesAllPass --> mergeSubMeshes")
            << nl << "    >>>> System command error. Aborting."
            << abort(FatalError);
    }

    // Copying constant/polyMesh
    if (fileStat(dir/"constant/polyMeshMerge").isValid() == 0)
    {
        result = help::runProcess("rm -rf constant/polyMeshMerge");
        mkDir("constant/polyMeshMerge");
    }
    Info << "   Copying constant/polyMesh to constant/polyMeshMerge" << endl;
    result = help::runProcess("rm -rf constant/polyMeshMerge/polyMesh");
    help::checkFile(dir/"constant/polyMesh");
    cp("constant/polyMesh","constant/polyMeshMerge");

    Info << "Done" << endl;
}

// ************************************************************************* //
