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
    getSubMeshes

Description
    Getting sub-meshes from the main mesh and putting them into constant directory.

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Optional argument
    // The utility cannot be run in parallel
    argList::noParallel();
    argList args(argc, argv);

    // Capture result from calling system commands. If non-zero value is
    // returned, give an error.
    int result(0);

    // Checking case folder
    fileName passDir = cwd();
    help::checkFile(passDir/"0","Missing 0 directory");
    help::checkFile(passDir/"system","Missing system directory");
    help::checkFile(passDir/"system/controlDict","Missing controlDict file");
    help::checkFile(passDir/"constant","Missing constant directory");
    help::checkFile(passDir/"constant/polyMesh","Missing polyMesh directory");

    // Case hygiene prior to procedure
    result = help::runProcess("rm -rf 0/materials*");
    result = help::runProcess("rm -rf 0/cellToRegion*");

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
        FatalErrorIn("getSubMeshes")
            << nl << "    >>>> Invalid pass type. Aborting. "
            << nl << "    Pass type: " << passType
            << abort(FatalError);
    }
    else
    {
        Info << nl << "Pass type: " << passType << endl;
    }

    // Split the sub-meshes from the base mesh
    Info << nl << "Running application splitMeshRegions" << endl;
    result =
        help::runProcess
        (
            "splitMeshRegions -overwrite", "log.splitMeshRegions"
        );
    if (result > 0)
    {
        FatalErrorIn("getSubMeshes --> splitMeshRegions")
            << nl << "    >>>> System command error. Aborting. "
            << abort(FatalError);
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

    // Copy submeshes to constant, depeneding on the pass type
    help::checkFile(passDir/"0"/wireDir,wireDir/" not present in 0 folder");
    if (passType == "drawingPass")
    {
        help::checkFile(passDir/"0"/dieDir,dieDir/" not present in 0 folder");

        Info << nl << "Copying submeshes to constant:" << endl;

        cp(passDir/"0"/dieDir,"constant");
        Info << "   --> " << dieDir << endl;

        if (fileStat(passDir/"0"/casingDir).isValid() != 0)
        {
            cp(passDir/"0"/casingDir,"constant");
            Info << "   --> " << casingDir << endl;
        }
    }
    else
    {
        word rollSetup = "-";
        fileName rollerName = "-";
        if( meshDict.found("rollSetup") )
        {
            rollSetup = word(meshDict.lookup("rollSetup"));
            Info << nl <<  "Roll setup: " << rollSetup << endl;
        }
        else
        {
            FatalErrorIn("getSubMeshes")
                << nl << "    >>>> System command error. Aborting. "
                << nl << "    Roll setup not found."
                << abort(FatalError);
        }

        Info << nl << "Copying submeshes to constant:" << endl;

        if (rollSetup == "singleRoll")
        {

            rollerName = "topRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }

            rollerName = "bottomRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }

            rollerName = "rightRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }

            rollerName = "leftRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }
        }
        else if (rollSetup == "edgeRolls")
        {
            rollerName = "rightRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }

            rollerName = "leftRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }
        }
        else if (rollSetup == "sideRolls")
        {
            rollerName = "topRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }

            rollerName = "bottomRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }

            rollerName = "rightRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }

            rollerName = "leftRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }
        }
        else if (rollSetup == "twoHighRolls")
        {
            rollerName = "topRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }

            rollerName = "bottomRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }
        }
        else if (rollSetup == "threeRolls")
        {
            rollerName = "topRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }

            rollerName = "bottomLeftRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }

            rollerName = "bottomRightRoll";
            if (fileStat(passDir/"0"/rollerName).isValid() != 0)
            {
                cp(passDir/"0"/rollerName,"constant");
                Info << "   --> " << rollerName << endl;
            }
        }
        else
        {
            FatalErrorIn("getSubMeshes")
                << nl << "    >>>> System command error. Aborting. "
                << nl << "    Unrecognised roll setup: " << rollSetup
                << abort(FatalError);
        }

    }
    cp(passDir/"0"/wireDir,"constant");
    Info << "   --> " << wireDir << endl;

    // Cleaning
    Info << nl << "Final cleanup (0 directory)" << endl;
    Info << "   Running application deleteSubMeshes" << endl;
    result =
        help::runProcess
        (
            "deleteSubMeshes -zero", "log.deleteSubMeshes"
        );
    if (result > 0)
    {
        FatalErrorIn("getSubMeshes --> deleteSubMeshes")
            << nl << "    >>>> System command error. Aborting. "
            << abort(FatalError);
    }
    Info << "Done" << endl;


    // Print out date and time
    Info << nl << nl << nl << "Extracting sub-meshes completed on " << clock::date() << " at "
         << clock::clockTime() << endl;

    Info << "END" << nl << endl;

    return 0;
}

// ************************************************************************* //
