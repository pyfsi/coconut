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
    prepareWireFromPrevious

Description
    Prepares wire from the end-time of the previous process line pass.

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

void findPassDirs(fileName&, fileName&);
void remeshWire(const Time&, const word&, const fileName&, const word& keyword = "geometryFile");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Optional argument
    // The utility cannot be run in parallel
    argList::noParallel();
    argList::validOptions.insert("wire", "word");
    argList::validOptions.insert("smooth", "word");
    argList::validOptions.insert("from", "word");
    argList::validOptions.insert("remesh", "");
    argList args(argc, argv);

    // Capture result from calling system commands. If non-zero value is
    // returned, give an error.
    int result(0);

    // Determine directories
    fileName passDir, passDirFrom;
    if (args.optionFound("from"))
    {
        passDir=cwd();
        passDirFrom = passDir.path() + "/" + args.option("from");

    }
    else
    {
        findPassDirs(passDir,passDirFrom);
    }

    // Exit in case no previous pass
    if (passDirFrom == "none")
    {
        Info << nl << "No previous pass detected.";

        Info << nl << nl << nl << "Preparing wire from previous pass completed on " << clock::date() << " at "
             << clock::clockTime() << endl;

        Info << "END" << nl << endl;

        return 0;
    }

    // Check for directory existence / case validity
    help::checkFile(passDir);
    help::checkFile(passDirFrom, "Invalid/nonexistent previous pass.");
    help::checkFile(passDirFrom/"system", "Invalid case.");
    help::checkFile(passDirFrom/"constant", "Invalid case.");
    help::checkFile(passDirFrom/"0", "Invalid case.");
    Info << nl << "Pass directories check:   OK";
    Info << nl << "          Current pass:   " << passDir.name();
    Info << nl << "         Previous pass:   " << passDirFrom.name() << endl;

    // Creating time objects
    Info << nl << "Creating time objects" << endl;
    Time runTime
    (
        Time::controlDictName,
        passDir.path(),
        passDir.name(),
        "system",
        "constant",
        false
    );
    Time runTimeFrom
    (
        Time::controlDictName,
        passDirFrom.path(),
        passDirFrom.name(),
        "system",
        "constant",
        false
    );

    // Determine the wire directory
    fileName wireDir;
    if (args.optionFound("wire"))
    {
        wireDir = args.option("wire");
        Info << nl << "Wire directory (argument): " << wireDir << endl;
    }
    else
    {
        // Deafault wire name
        wireDir ="wire";

        // Look if there is wireDir entry in meshDict
        IOdictionary meshDict
        (
            IOobject
            (
                "meshDict",
                "system",
                runTimeFrom,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        if( meshDict.found("wireDir") )
        {
            wireDir = word(meshDict.lookup("wireDir"));
            Info <<  "Wire directory (meshDict): " << wireDir << endl;
        }
        else
        {
            Info << nl << "Wire directory (default): " << wireDir << endl;
        }
    }


    // PREVIOUS PASS OPERATIONS
    Info << nl << nl << "* * * * * * * * * * * *" << endl;
    Info << "PREVIOUS PASS: " << passDirFrom.name() << endl;
    Info << nl << "   Location: " << passDirFrom << endl;
    chDir(passDirFrom);

    // Split the wire sub-mesh from the base mesh
    Info << nl << "   Running application splitMeshRegionsLastTimeStep" << endl;
    result =
        help::runProcess
        (
            "splitMeshRegionsLastTimeStep", "log.splitMeshRegionsLastTimeStep"
        );
    if (result > 0)
    {
        FatalErrorIn("prepareWireFromPrevious --> splitMeshRegionsLastTimeStep")
            << nl << "    >>>> System command error. Aborting. "
            << abort(FatalError);
    }

    // Smooth the wire mesh and fields in the flow direction
    if (args.optionFound("smooth"))
    {
        scalar newWireLength = strtod(args.option("smooth").c_str(),NULL);
        if (newWireLength > 0)
        {
            Info << "   Running application smoothWireMeshFields " << newWireLength << endl;
            result = help::runProcess
                (
                    std::string("smoothWireMeshFields " + Foam::name(newWireLength)).c_str() , "log.smoothWireMeshFields"
                );
            if (result > 0)
            {
                FatalErrorIn("prepareWireFromPrevious --> smoothWireMeshFields")
                    << nl << "    >>>> System command error. Aborting. "
                    << abort(FatalError);
            }
        }
        else
        {
            FatalError
                << "Inappropriate -smooth argument. Try with positive number. "
                << abort(FatalError);
        }
    }

    // Set time to latest time
    Time& runTimeLatest = runTimeFrom;
    const instantList& times = runTimeLatest.times();
    const label latestTimeID = times.size() -1;
    runTimeLatest.setTime(times[latestTimeID], latestTimeID);

    // Create wire mesh object
    help::checkFile(runTimeLatest.path()/runTimeLatest.timeName()/wireDir,"Invalid wire directory.");
    help::checkFile(runTimeLatest.path()/runTimeLatest.timeName()/wireDir/"polyMesh","No mesh in wire directory.");
    Info << nl << "   Reading wire mesh for previous pass from time = " << runTimeLatest.timeName() << endl;
    fvMesh wireMeshFrom
    (
        IOobject
        (
            wireDir,
            runTimeLatest.timeName(),
            runTimeLatest,
            IOobject::MUST_READ
        )
    );

    if (runTimeLatest.timeName() != wireMeshFrom.pointsInstance())
    {
        FatalError
            << "There is no " << wireDir << " directory for time = " << runTimeLatest.timeName()
            << abort(FatalError);
    }

    Info << nl << "   Wire mesh info:" << endl;
    Info << "      Cells:      " << wireMeshFrom.nCells() << endl;
    Info << "      Faces:      " << wireMeshFrom.nFaces() << endl;
    Info << "      Points:     " << wireMeshFrom.nPoints() << endl;
    Info << "      Read from:  " << runTimeLatest.path() + "/" + wireMeshFrom.pointsInstance() + "/" + wireMeshFrom.meshDir() << endl;

    // Remesh wire
    if (args.optionFound("remesh"))
    {
        // Extract wire profile & remesh wire
        Info << nl << "   Running application extractWireProfile" << endl;
        word wireProfile;
        result = help::runProcess
            (
                std::string("extractWireProfile -region " + wireDir).c_str(), "log.extractWireProfile"
            );
        if (result > 0)
        {
            FatalErrorIn("prepareWireFromPrevious --> extractWireProfile")
                << nl << "    >>>> System command error. Aborting. "
                << abort(FatalError);
        }

        wireProfile = "wireProfile_" +  passDirFrom.name() + ".fms";

        remeshWire(runTimeLatest, wireProfile, wireDir, "surfaceFile");
    }

    // CURRENT PASS OPERATIONS
    Info << nl << nl << "* * * * * * * * * * * *" << endl;
    Info << "CURRENT PASS: " << passDir.name() << endl;
    Info << nl << "   Location: " << passDir << endl;
    chDir(passDir);

    // Copy the wire region mesh from the last time step of the
    // previous pass into the constant directory of the current
    // pass
    const fileName targetDirPrev = runTimeLatest.path()/"target";
    const fileName targetDirCur = runTime.path()/"source";
    help::checkFile(targetDirPrev);
    result = help::runProcess("rm -rf source");
    Info << nl <<  "   Copying target case from:" << nl
            << "    " << targetDirPrev << nl
            << "   to:" << nl
            << "    " << targetDirCur << endl;
    cp(targetDirPrev, targetDirCur);

    // Copy the wire region mesh from the last time step of the
    // previous pass into the constant directory of the current
    // pass
    const fileName wireMeshDirPrev = runTimeLatest.path()/runTimeLatest.timeName()/wireDir/"polyMesh";
    const fileName wireMeshDirCur = runTime.path()/"constant"/wireDir/"polyMesh";
    help::checkFile(wireMeshDirPrev);
    result = help::runProcess(std::string("rm -rf constant"/wireDir).c_str());
    Info << nl <<  "   Copying wire mesh from:" << nl
            << "    " << wireMeshDirPrev << nl
            << "   to:" << nl
            << "    " << wireMeshDirCur << endl;
    cp(wireMeshDirPrev, wireMeshDirCur);

    // Copy the wire region fields from target
    // into the 0 directory of the current pass
//    const fileName wireFieldsDirPrev = runTime.path()/"target"/"0";
//    const fileName wireFieldsDirCur = runTime.path()/"0";
//    help::checkFile(wireFieldsDirPrev);
//    result = help::runProcess("rm -rf 0");
//    Info << nl << "   Copying wire fields from:" << nl
//            << "    " << wireFieldsDirPrev << nl
//            << "   to:" << nl
//            << "    " << wireFieldsDirCur << endl;
//    cp(wireFieldsDirPrev, wireFieldsDirCur);

    // Remove new latest time step created in previous pass + remove old fields which cause problem on start of new case
    Info << nl << "   Final cleanup" << endl;
    result = help::runProcess(std::string("rm -rf "/runTimeLatest.path()/runTimeLatest.timeName()).c_str());
    result = help::runProcess(std::string("rm -rf "/targetDirPrev).c_str());
    result = help::runProcess("rm -rf 0/polyMesh");
    chDir("target/0");
    result = help::runProcess("rm -rf DU*");
    result = help::runProcess("rm -rf rho*");
    result = help::runProcess("rm -rf U*");
    result = help::runProcess("rm -rf polyMesh");

    // Print out date and time
    Info << nl << nl << nl << "Preparing wire from previous pass completed on " << clock::date() << " at "
         << clock::clockTime() << endl;

    Info << "END" << nl << endl;

    return 0;
}

void findPassDirs(fileName& dir, fileName& dirFrom)
{
    // Current pass directory
    dir=cwd();
    Foam::string dirName(dir.name());

    // Determine pass number
    int num = -1;
    int underscore = -1;
    for (int i = dirName.length()-1; i >= 0.;i--)
    {
        if(dirName[i] == '_')
        {
            underscore = i;
            break;
        }
    }
    Foam::string passNumString(dirName(underscore + 1, dirName.length() - underscore - 1));
    num = int(strtod(passNumString.c_str(),NULL));
    if(num < 2)
    {
//        FatalError
//            << "The utility cannot be used in current directory. There is no previous pass."
//            << abort(FatalError);
        dirFrom = "none";
        return;
    }

    // Previous pass directory
    string dirFromName(dirName(0,underscore+1) + name(num-1));
    dirFrom = dir.path() + "/" + dirFromName;

}


void remeshWire(const Time& runTime, const word& wireProfile, const fileName& wireDir, const word& keyword)
{
    Info << nl << "   Re-meshing the wire and mapping fields" << endl;

    // Capture result from calling system commands. If non-zero value is
    // returned, give an error.
    int result(0);

    // Prepare directories
    Info << "      Preparing directories ..." << endl;
    const fileName oldWire =
        runTime.path()/runTime.timeName()/wireDir;
    const fileName source = runTime.path()/"source";
    const fileName target = runTime.path()/"target";
    chDir(runTime.path());
    result = help::runProcess("rm -rf source");
    result = help::runProcess("rm -rf target");
    result = help::runProcess("rm -rf system_copy");
    mkDir("source");
    chDir(source);
    mkDir("constant");

    // Copy the fields to be mapped. This needs to be aligned with the
    // setFieldsDict for initialisation in the target case.
    // Also copy the old mesh
    Info << "      Copying fields and data ..." << endl;

    help::checkFile(oldWire);
    cp(oldWire,"0");

    // Remove (old) fields which cause problems
    chDir(source/"0");
    result = help::runProcess("rm -rf DU*");
    result = help::runProcess("rm -rf rho*");
    result = help::runProcess("rm -rf U*");
    result = help::runProcess("rm -rf J*");
    chDir(source);

    help::checkFile(cwd()/"0"/"polyMesh");
    result = help::runProcess("mv 0/polyMesh constant/polyMesh");

    chDir(runTime.path());

    help::checkFile(source);
    cp(source, target);

    help::checkFile(cwd()/wireProfile, "Wire profile file is missing.");
    cp(wireProfile, target);

    help::checkFile(cwd()/"system");
    cp("system", source);

    // Update meshDict entry for re-meshing
    cp("system", "system_copy");
    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            "system",
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );
    if( meshDict.found("wireMeshDict") )
    {
        meshDict.subDict("wireMeshDict").remove("geometryFile");
        meshDict.subDict("wireMeshDict").remove("surfaceFile");
        meshDict.subDict("wireMeshDict").remove("dxfFile");
        meshDict.subDict("wireMeshDict").add(keyword, wireProfile, true);
    }
    meshDict.regIOobject::write();

    cp("system",target);
    result = help::runProcess("rm -r system");

    help::checkFile(cwd()/"system_copy");
    cp("system_copy", "system");

    // Remshing wire & mapping fields
    chDir(target);

    Info << "         Running application setFields" << endl;
    result = help::runProcess("setFields", "log.setFields");
    if (result > 0)
    {
        FatalErrorIn("prepareWireFromPrevious --> setFields")
            << nl << "    >>>> System command error. Aborting. "
            << abort(FatalError);
    }

    Info << "         Running application generateWireMesh" << endl;
    result = help::runProcess("generateWireMesh", "log.generateWireMesh");
    if (result > 0)
    {
        FatalErrorIn("prepareWireFromPrevious --> generateWireMesh")
            << nl << "    >>>> System command error. Aborting. "
            << abort(FatalError);
    }

    Info << "         Running application mapFields" << endl;   
    result = help::runProcess("mapFields ../source", "log.mapFields");
    if (result > 0)
    {
        FatalErrorIn("prepareWireFromPrevious --> mapFields")
            << nl << "    >>>> System command error. Aborting. "
            << abort(FatalError);
    }

    // Directory management
    Info << "      Final directory management ..." << endl;

    help::checkFile(cwd()/"constant"/"polyMesh");
    result = help::runProcess("cp -r constant/polyMesh 0/");

    chDir(runTime.path()/runTime.timeName());
    result = help::runProcess(std::string("rm -rf " + wireDir).c_str());
    chDir(runTime.path());

    cp(target/"0", runTime.path()/runTime.timeName());
    mv(runTime.path()/runTime.timeName()/"0", oldWire);

    // Cleaning unnecessary directories
    result = help::runProcess("rm -rf system_copy");
    result = help::runProcess("rm -rf source");

}



// ************************************************************************* //
