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
    deleteSubMeshes

Description
    Delete sub-meshes in 0 and/or constant directory.

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

void deleteSubMeshes(const fileName&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Optional argument
    // The utility cannot be run in parallel
    argList::noParallel();
    argList::validOptions.insert("constant", "");
    argList::validOptions.insert("zero", "");
    argList::validOptions.insert("all", "");
    argList args(argc, argv);

    // Deletion flags
    bool deleteSubmeshesConstant = args.optionFound("constant");
    bool deleteSubmeshesZero = args.optionFound("zero");
    bool deleteSubmeshesAll = args.optionFound("all");

    // Checking case folder
    fileName passDir = cwd();
    //help::checkFile(passDir/"0","Missing 0 directory");
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
        FatalErrorIn("deleteSubMeshes")
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



    // Cleaning sub-mesh directories
    Info << nl << "Cleaning sub-meshes: " << endl;
    word dirClean = "-";
    if (deleteSubmeshesConstant || deleteSubmeshesAll)
    {
        dirClean = "constant";
        //Info << nl << "Cleaning submeshes in constant: " << endl;
        deleteSubMeshes(passDir/dirClean/"topRoll");
        deleteSubMeshes(passDir/dirClean/"bottomRoll");
        deleteSubMeshes(passDir/dirClean/"rightRoll");
        deleteSubMeshes(passDir/dirClean/"leftRoll");
        deleteSubMeshes(passDir/dirClean/"bottomRightRoll");
        deleteSubMeshes(passDir/dirClean/"bottomLeftRoll");
        deleteSubMeshes(passDir/dirClean/"polyMeshMerge");
        deleteSubMeshes(passDir/dirClean/wireDir);
        deleteSubMeshes(passDir/dirClean/dieDir);
        deleteSubMeshes(passDir/dirClean/casingDir);
        deleteSubMeshes(passDir/dirClean/rollDir);
        //Info << "Done" << endl;
    }
    if (deleteSubmeshesZero || deleteSubmeshesAll)
    {
        dirClean = "0";
        //Info << nl << "Cleaning submeshes in 0: " << endl;
        deleteSubMeshes(passDir/dirClean/"topRoll");
        deleteSubMeshes(passDir/dirClean/"bottomRoll");
        deleteSubMeshes(passDir/dirClean/"rightRoll");
        deleteSubMeshes(passDir/dirClean/"leftRoll");
        deleteSubMeshes(passDir/dirClean/"bottomRightRoll");
        deleteSubMeshes(passDir/dirClean/"bottomLeftRoll");
        deleteSubMeshes(passDir/dirClean/"polyMeshMerge");
        deleteSubMeshes(passDir/dirClean/wireDir);
        deleteSubMeshes(passDir/dirClean/dieDir);
        deleteSubMeshes(passDir/dirClean/casingDir);
        deleteSubMeshes(passDir/dirClean/rollDir);
        //Info << "Done" << endl;
    }
    if (deleteSubmeshesAll)
    {
        dirClean = "target";
        //Info << nl << "Cleaning submeshes in target: " << endl;
        deleteSubMeshes(passDir/dirClean/"constant");
        deleteSubMeshes(passDir/dirClean/"0");
        //Info << "Done" << endl;

        dirClean = "source";
        //Info << nl << "Cleaning submeshes in source: " << endl;
        deleteSubMeshes(passDir/dirClean/"constant");
        deleteSubMeshes(passDir/dirClean/"0");
        //Info << "Done" << endl;
    }
    Info << "Done" << endl;

    // Print out date and time
    Info << nl << nl << nl << "Deleting sub-meshes completed on " << clock::date() << " at "
         << clock::clockTime() << endl;

    Info << "END" << nl << endl;

    return 0;
}

void deleteSubMeshes(const fileName& dir)
{
    int result(0);
    if (fileStat(dir).isValid() != 0)
    {
        Info << "   " << dir << endl;
        result = help::runProcess(std::string("rm -rf " + dir).c_str());
    }
}

// ************************************************************************* //
