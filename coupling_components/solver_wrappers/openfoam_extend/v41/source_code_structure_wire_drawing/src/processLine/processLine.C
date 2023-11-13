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

\*---------------------------------------------------------------------------*/

#include "processLine.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processLine, 0);
}


// * * * * * * * * * *  Private Member Functions  * * * * * * * * * * * * * * //

void Foam::processLine::makeRunTime() const
{
    if (runTimePtr_)
    {
        FatalErrorIn("void Foam::processLine::makeRunTime() const")
            << "pointer already set" << abort(FatalError);
    }

    // We must add the folowing two items to the controlDict to keep the Time
    // object happy
    controlDict_.set("deltaT", 1.0);
    controlDict_.set("writeInterval", 1.0);

    runTimePtr_ =
        new Time
        (
            controlDict_,               // controlDict
            cwd(),                      // root directory
            "./",                       // case directory
            ".",                        // system directory name
            ".",                        // constant directory name
            false                       // functionObjects disabled
        );
}


const Foam::Time& Foam::processLine::runTime() const
{
    if (!runTimePtr_)
    {
        makeRunTime();
    }

    return *runTimePtr_;
}


void Foam::processLine::makeDict() const
{
    if (dictPtr_)
    {
        FatalErrorIn("void Foam::processLine::makeDict() const")
            << "pointer already set" << abort(FatalError);
    }

    // Read the processLineDict
    dictPtr_ =
        new IOdictionary
        (
            IOobject
            (
                word(processLineName_ + "Dict"),
                fileName("dictFiles"),
                runTime(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
}


Foam::dictionary& Foam::processLine::dict() const
{
    if (!dictPtr_)
    {
        makeDict();
    }

    return *dictPtr_;
}


void Foam::processLine::makePasses()
{
    if (passes_.size() != 0)
    {
        FatalError
            << "pointer list already set!" << abort(FatalError);
    }

    // Read the passes dict
    const dictionary& passesDict(dict().subDict("passes"));

    // Get the list of names of the passes
    const wordList passNames = passesDict.toc();

    // Check all passes are specified as subDicts
    forAll(passNames, passI)
    {
        if (!passesDict.isDict(passNames[passI]))
        {
            FatalError
                << "All entries in the passes dict should be pass subDicts, but"
                << " pass number " << (passI + 1) << " is not!"
                << abort(FatalError);
        }
    }

    Info<< "Found " << passNames.size() << " passes for process line "
        << dataContainer().returnKey<word>("processLineName") << ":" << endl;

    // Create a list processPass objects
    passes_.setSize(passNames.size());

    forAll(passes_, passI)
    {
        const word& passName = passNames[passI];

        // Create each pass object
        passes_.set
        (
            passI,
            new processPass
            (
                *this,
                passName,
                passesDict.subDict(passName),
                passI,
                &dataContainer()
             )
        );
    }
}


Foam::PtrList<Foam::processPass>& Foam::processLine::passes()
{
    if (passes_.size() == 0)
    {
        makePasses();
    }

    return passes_;
}


void Foam::processLine::copyWireMeshAndFields
(
    const label fromPassID, const label toPassID
)
{
    // Find the latest time-step in the previous pass

    // Time object from previous pass
    Time& runTimePrev = passes()[fromPassID].runTime();

    // List of time-steps from previous pass
    const instantList& times = runTimePrev.times();

    // Index of the latest time-step
    const label latestTimeID = times.size() - 1;

    // Set the previous time to be the latest time-step
    runTimePrev.setTime(times[latestTimeID], latestTimeID);


    // Copy the wire region mesh from the last time step of the
    // previous pass into the constant directory of the current
    // pass

    {
        // Wire mesh directory from the previous pass
        const fileName wireMeshDirPrev =
            runTimePrev.path()/runTimePrev.timeName()
            /passes()[fromPassID].wire().name()/"polyMesh";

        // Wire mesh directory for the current pass
        const fileName wireMeshDirCur =
            passes()[toPassID].runTime().path()/"constant"
            /passes()[toPassID].wire().name()/"polyMesh";

        Info<< "    Copying wire mesh from:" << nl
            << "        " << wireMeshDirPrev << nl
            << "    to:" << nl
            << "        " << wireMeshDirCur << endl;

        cp(wireMeshDirPrev, wireMeshDirCur);
    }


    // Copy the wire region fields from the last time step of the
    // previous pass into the 0 directory of the current pass

    {
        // Wire mesh directory from the previous pass
        const fileName wireFieldsDirPrev =
            runTimePrev.path()/runTimePrev.timeName()
            /passes()[fromPassID].wire().name();

        // Wire mesh directory for the current pass
        const fileName wireFieldsDirCur =
            passes()[toPassID].runTime().path()/"0";

        Info<< "    Copying wire fields from:" << nl
            << "        " << wireFieldsDirPrev << nl
            << "    to:" << nl
            << "        " << wireFieldsDirCur << endl;

        cp(wireFieldsDirPrev, wireFieldsDirCur);


        Info<< "    Removing polyMesh, DU, U and rho fields from 0 folder"
            << endl;

        const fileName fieldDir =
            wireFieldsDirCur/passes()[fromPassID].wire().name();

        // Remove wire polyMesh from 0 folder in current pass to ensure the wire
        // mesh from the constant folder is subsequently used
        rmDir(fileName(fieldDir/"polyMesh"));

        // Remove old fields which cause problem on start of new case
        rm(fileName(fieldDir/"DU"));
        rm(fileName(fieldDir/"DU_0"));
        rm(fileName(fieldDir/"DU_0_0"));
        rm(fileName(fieldDir/"DU_0_0_0"));
        rm(fileName(fieldDir/"U"));
        rm(fileName(fieldDir/"U_0"));
        rm(fileName(fieldDir/"U_0_0"));
        rm(fileName(fieldDir/"U_0_0_0"));
        rm(fileName(fieldDir/"rho"));
        rm(fileName(fieldDir/"rho_0"));
        rm(fileName(fieldDir/"rho_0_0"));
        rm(fileName(fieldDir/"rho_0_0_0"));
    }

    // Reset the fromPass time state to the the initial time
    runTimePrev.setTime(times[0], 0);
}


void Foam::processLine::remeshWire(const label passID)
{
    // Take reference to the pass
    processPass& currentPass = passes()[passID];

    // Take reference to the runTime for this pass
    Time& runTime = currentPass.runTime();

    // Take reference to the wire for this pass
    wirePassItem& wire = currentPass.wire();

    if
    (
        dataContainer_.processLineProgressDict().lookupOrDefault<bool>
        (
            currentPass.name() + "_remeshWire", false
        )
    )
    {
        Info<< nl
            << "INTER_PASS: "
            << "Wire remeshing in pass " << currentPass.name()
            << " has already been performed: skipping remeshWire" << nl
            << "--------------------" << endl;

        return;
    }

    if (wire.meshInputDict().found("remeshAfterwards"))
    {
        FatalErrorIn("void Foam::processLine::remeshWire(const Time& runTime)")
            << "'remeshAfterwards' option is no longer used for the wire."
            << " The wire is now remeshed at the start of a pass rather than at"
            << " the end of a pass. Please remove this option from the current "
            << " pass and add the option 'remeshFromPreviousPass' to the wire "
            << " mesh input dict for the next pass."
            << endl;
    }

    if
    (
        wire.meshInputDict().lookupOrDefault<Switch>
        (
            "remeshFromPreviousPass", "no"
        )
    )
    {
        Info<< nl
            << "INTER_PASS: "
            << " re-meshing the wire in pass "
            << passes()[passID].name() << nl
            << "--------------------"
            << endl;

        const fileName oldWire =
            runTime.path()/runTime.timeName()/wire.name();

        const fileName source = runTime.path()/"source";

        const fileName target = runTime.path()/"target";

        // Always remove the source and target if they exists, in case a
        // previous remesh failed
        if (isDir(source))
        {
            dataContainer().runSystemCommand
            (
                "rm -rf " + source,
                "log.rm_source",
                "void Foam::processLine::remeshWire(const label passID)"
            );
        }
        if (isDir(target))
        {
            dataContainer().runSystemCommand
            (
                "rm -rf " + target,
                "log.rm_target",
                "void Foam::processLine::remeshWire(const label passID)"
            );
        }

        cp(oldWire, source);
        cp("system", source);

        dictionary controlDict(runTime.controlDict());

        // Copy the fields and mesh to be mapped
        chDir(source);
        mkDir("0");
        mkDir("constant");

        if (isDir("polyMesh"))
        {
            mv(cwd()/"polyMesh", cwd()/"constant"/"polyMesh");
        }
        else
        {
            FatalErrorIn(type())
                << "There is no wire mesh in the last time-step!"
                << abort(FatalError);
        }

        // Remove fields that cause problems for mapping
        dataContainer().runSystemCommand
        (
            "rm cellToRegion.gz materials.gz frictionCoeff.gz *_0.gz DU.gz",
            "log.rm_fieldsThatCauseProblemsForMapping",
            "void Foam::processPass::run()",
            false
        );

        // Move all other fields to time 0 for mapping
        dataContainer().runSystemCommand
        (
            "mv *.gz 0",
            "log.mv_fieldsForMapping",
            "void Foam::processPass::run()",
            false
        );

        chDir(runTime.path());

        // Copy the source directory to the target directory
        cp(source, target);

        // PC, MC, 13/Mar/18
        // Remove fields in the target 0 directory as they will be created by
        // the mapFields utility
        dataContainer().runSystemCommand
        (
            "rm -rf target/0/*",
            "log.rm_targetFieldsIn0",
            "void Foam::processPass::run()",
            false
        );

        // Copy a setFieldsDict and mapFieldsDict from the dictFiles directory
        const fileName setFieldsDict =
            dataContainer().returnKey<fileName>
            (
                "processLineRootDir"
            )/"dictFiles"/"defaultDicts"/"setFieldsDict";

        const fileName mapFieldsDict =
            dataContainer().returnKey<fileName>
            (
                "processLineRootDir"
            )/"dictFiles"/"defaultDicts"/"mapFieldsDict";

        // fileName profileFileName("wireProfile.prof");
        // cp(profileFileName, target);

        Info<< "            Copying setFieldsDict to ./target/system" << endl;
        cp(setFieldsDict, target/"system");

        Info<< "            Copying mapFieldsDict to ./target/system" << endl;
        cp(mapFieldsDict, target/"system");

        // Copy in default meshDict
        dataContainer().runSystemCommand
        (
            "cp ../../dictFiles/defaultDicts/meshDict system/",
            "log.cp_meshDict_remesh_wire",
            "void Foam::processLine::remeshWire(const Time& runTime)"
        );

        // Make a copy of the original meshDict
        dataContainer().runSystemCommand
        (
            "cp system/meshDict system/meshDict.orig",
            "log.cp_meshDict_meshDictOrig",
            "void Foam::processLine::remeshWire(const Time& runTime)"
        );

        // Write out new wire meshDict
        const fileName profileFileName =
            "wireProfile_" + runTime.caseName() + ".fms";
        wire.writeMeshDict(runTime, profileFileName);

        // Copy meshDict to the target case for remeshing later, but also leave
        // a copy in the main system folder
        cp("system/meshDict", target/"system"/"meshDict");

        // Extract wire profile from downstream patch
        dataContainer().runSystemCommand
        (
            word
            (
                "extractWireProfile -region " +
                wire.name()
            ),
            "log.extractWireProfile",
            "void Foam::processLine::remeshWire(const Time& runTime)"
        );

        // Copy the profile to the target case where it will be used to create
        // a new wire mesh
        cp(profileFileName, target/profileFileName);

        // Replace the original meshDict: this should be done after
        // extractWireProfile above, as extractWireProfile needs the target
        // meshDict
        dataContainer().runSystemCommand
        (
            "mv system/meshDict.orig system/meshDict",
            "log.mv_meshDictOrg_meshDict",
            "void Foam::processLine::remeshWire(const Time& runTime)"
        );

        chDir(target);

        // setFields needed for merging fields after mapping
        dataContainer().runSystemCommand
        (
            "setFields",
            "log.setFields",
            "void Foam::processLine::remeshWire(const Time& runTime)"
        );

        // Generate new wire mesh
        dataContainer().runSystemCommand
        (
            "export OMP_NUM_THREADS=1; "
            "generateWireMesh",
            "log.generateWireMesh",
            "void Foam::processLine::remeshWire(const Time& runTime)"
        );

        // Before using mapFields with the "consistent" option, we must make
        // sure that the patches of the source and target mesh are in the same
        // order
        dataContainer().runSystemCommand
        (
            "reorderPatchesConsistent ../source ",
            "log.reorderPatchesConsistent",
            "void Foam::processLine::remeshWire(const Time& runTime)"
        );

        // Map fields from old to new wire mesh
        // PC, 02/Mar/18: mapConservativeFields some times misses cells near
        // the boundary, which can cause the solver to crash; mapFields is
        // only first order but is robust
        dataContainer().runSystemCommand
        (
            //"mapConservativeFields ../source &> log.mapFields",
            "mapFields ../source -consistent",
            "log.mapFields",
            "void Foam::processLine::remeshWire(const Time& runTime)"
        );

        cp("constant/polyMesh", "0/polyMesh");

        chDir(runTime.path()/runTime.timeName());

        if (isDir(wire.name()))
        {
            rmDir(wire.name());
        }

        chDir(runTime.path());

        cp(target/"0", runTime.path()/runTime.timeName());

        mv(runTime.path()/runTime.timeName()/"0", oldWire);

        // Record that remeshWire finished correctly to allow restart
        dataContainer_.processLineProgressDict().set
        (
            word(currentPass.name() + "_remeshWire"), true
        );
        dataContainer_.processLineProgressDict().regIOobject::write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processLine::processLine
(
   const fileName& processLineName,
   const bool restart
)
:
    processLineName_(processLineName),
    controlDict_(),
    runTimePtr_(NULL),
    dictPtr_(NULL),
    passes_(),
    dataContainer_(runTime())
{
    // PC, 2/Mar/18: calling a system command creates a separate bash session,
    // so this command will not affect the subsequent cfMesh calls; instead we
    // need to call this command at the same time as the cfMesh call
    // Force cfMesh to run in serial as the parallel option is less stable for
    // now.
    // dataContainer_.runSystemCommand
    // (
    //     "export OMP_NUM_THREADS=1",
    //     "Foam::processLine::processLine"
    // );

    // Getting general settings for the process line
    dataContainer().addKey<word>("processLineName", processLineName);

    const word processLineSymmetry =
        word(generalSettings().lookup("processLineSymmetry"));

    dataContainer().addKey<word>("processLineSymmetry", processLineSymmetry);

    if (processLineSymmetry == "axiSymmetric")
    {
        dataContainer().addKey<scalar>
        (
            "wedgeAngle",
            readScalar(generalSettings().lookup("wedgeAngle"))
        );
    }

    dataContainer().addKey<scalar>
    (
        "lineSpeed",
        readScalar(generalSettings().lookup("lineSpeed"))
    );

    dataContainer().addKey<scalar>
    (
        "ambientTemperature",
        readScalar(generalSettings().lookup("ambientTemperature"))
    );

    dataContainer().addKey<dictionary>
    (
        "patchNames",
        generalSettings().subDict("patchNames")
    );

    dataContainer().addKey<dictionary>
    (
        "defaultFields",
        generalSettings().subDict("defaultFields")
    );

    Info<< nl << "Creating processLine: "
        << word(processLineName) << nl << endl;

    // Adding default directories
    dataContainer().addKey<fileName>
    (
        "processLineRootDir",
        cwd()
    );

    dataContainer().addKey<fileName>
    (
       "processLineDir",
       fileName(cwd()/processLineName)
    );

    // Setup the process line directory
    bool restartedProcessLine = false;
    if (isDir(dataContainer().returnKey<fileName>("processLineDir")))
    {
        char response;

        if (restart)
        {
            response = 'R';
        }
        else
        {
            Info<< "Previous processLine directory found: "
                << word(dataContainer().returnKey<fileName>("processLineDir"))
                << nl << "To restart the process line type 'R',"
                << "to start a new process line and remove the previous line "
                << "type 'Y': ";
            cin >> response;
        }

        if (response == 'Y')
        {
            Info<< "Removing directory "
                << word(dataContainer().returnKey<fileName>("processLineDir"))
                << endl;
            rmDir(dataContainer().returnKey<fileName>("processLineDir"));

            // Clear the progress dict
            dataContainer().processLineProgressDict().clear();
        }
        else if (response == 'R')
        {
            Info<< "Restarting the process line" << endl;
            restartedProcessLine = true;
        }
        else
        {
            Info<< "Previous process line directory not removed" << endl;

            FatalErrorIn
            (
                "Foam::processLine::processLine"
                "(const fileName& processLineName)"
            )   << "Process line directory already exists!"
                << abort(FatalError);
        }
    }

    if (!restartedProcessLine)
    {
        // Create the process line directory
        Info<< "Creating empty directory "
            << word(dataContainer().returnKey<fileName>("processLineDir"))
            << nl << endl;

        mkDir(dataContainer().returnKey<fileName>("processLineDir"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::processLine::~processLine()
{}


// * * * * * * * * * * * * * * * Public Member Functions  * * *  * * * * * * //

Foam::dictionary& Foam::processLine::generalSettings()
{
    return dict().subDict("generalSettings");
}

void Foam::processLine::setup()
{
    Info<< nl << "******************************************" << endl;
    Info<< "        SETUP THE PROCESS LINE" << endl;
    Info<< "******************************************" << nl << endl;

    Info<< "Process line symmetry: "
        << word(dataContainer().returnKey<word>("processLineSymmetry")) << endl;

    // Setup the passes
    forAll(passes(), passI)
    {
        if
        (
            dataContainer_.processLineProgressDict().lookupOrDefault<bool>
            (
                passes()[passI].name(), false
            )
        )
        {
            Info<< passes()[passI].name() << " has already been run!" << endl;
        }
        else
        {
            const fileName previousPassDir =
                fileName
                (
                    dataContainer().returnKey<fileName>("processLineDir")
                   /passes()[passI].name()
                );

            if (isDir(previousPassDir))
            {
                Info<< "Removing previous pass directory: "
                    << previousPassDir << endl;
                rmDir(previousPassDir);
            }
            passes()[passI].setup();
        }
    }
}


void Foam::processLine::setupAndRun()
{
    Info<< nl << nl << "****************************************" << endl;
    Info<<  "         SETUP AND RUN THE PROCESS LINE" << endl;
    Info<< "****************************************" << endl;

    Info<< "Process line symmetry: "
        << word(dataContainer().returnKey<word>("processLineSymmetry")) << endl;

    // Run each passes in turn
    forAll(passes(), passI)
    {
        Info<< nl << "PASS: " << passes()[passI].name() << endl;
        Info<< "--------------------" << endl;
        Info<< nl << ">> CASE PREPARATION" << endl;

        if
        (
            dataContainer_.processLineProgressDict().lookupOrDefault<bool>
            (
                passes()[passI].name(), false
            )
        )
        {
            Info<< passes()[passI].name() << " has already been run!" << endl;
        }
        else
        {
            // Delete pass directory if it already exists e.g. from a previous
            // filed run
            const fileName passDir =
                fileName
                (
                    dataContainer().returnKey<fileName>("processLineDir")
                   /passes()[passI].name()
                );

            if (isDir(passDir))
            {
                Info<< "Removing pass directory from previous run: "
                    << passDir << endl;
                rmDir(passDir);
            }

            // Setup the pass directories and dictionaries and meshes
            passes()[passI].setup();

            // Copy wire from previous pass
            if (passI > 0)
            {
                // Copy wire submesh from previous pass in current pass
                if (!passes()[passI].wire().createMesh())
                {
                    // Optionally: remesh wire at the end of the previous pass
                    remeshWire(passI - 1);

                    Info<< "    Copying wire from previous pass" << endl
                        << "    Current and previous pass names: "
                        << passes()[passI].name() << ", "
                        << passes()[passI - 1].name() << endl;

                    // We must copy the wire mesh and fields from the latest
                    // time in the previous pass to the current pass
                    copyWireMeshAndFields(passI - 1, passI);
                }
                else
                {
                    FatalErrorIn("void Foam::processLine::setupAndRun()")
                        << "Wire mesh created for " << passes()[passI].name()
                        << ", with pass number " << (passI + 1)
                        << " greater than 1!"
                        << abort(FatalError);
                }
            }

            Info<< nl
                << "PASS: "<< passes()[passI].name()
                << " : final preparations and running the solver" << nl
                << "--------------------"
                << endl;

            // Final preparations, run the solver and post-processing
            passes()[passI].run();

            // Record that this pass finished correctly to allow restart
            dataContainer_.processLineProgressDict().set
            (
                passes()[passI].name(), true
            );
            dataContainer_.processLineProgressDict().regIOobject::write();
        }
    }
}


void Foam::processLine::postProcess()
{
    Info<< nl << "Process line post processing" << nl << endl;

    Info<< "    Generating PDF report" << endl;

    // Open the report file

    chDir(runTime().path());

    //dataContainer().printData();

    const fileName reportName
    (
        dataContainer().returnKey<fileName>("processLineDir")
       /dataContainer().returnKey<word>("processLineName") + ".pdf"
    );
    OFstream report(reportName);

    report
        << "\\documentclass{article}" << nl
        << "\\usepackage{siunitx}" << nl
        << "\\usepackage{graphicx}" << nl
        << "\\usepackage{natbib}" << nl
        << "\\usepackage{amsmath}" << nl
        << "\\setlength\\parindent{0pt} % "
        << "Removes all indentation from paragraphs" << nl
        << "%\\usepackage{times} % Uncomment to use the Times New Roman font"
        << nl
        << nl
        << "\\title{Bekaert report \\\\ process line \\\\ model results}" << nl
        << "\\author{John \\textsc{Smith}}" << nl
        << "\\date{\\today}" << nl
        << nl
        << "\\begin{document}" << nl
        << "\\maketitle" << nl
        << "\\begin{center}" << nl
        << "\\begin{tabular}{l r}" << nl
        << "Date Performed: & January 1, 2012 \\\\" << nl
        << "Partners: & James Smith \\\\ % Partner names" << nl
        << "& Mary Smith \\\\" << nl
        << "Instructor: & Professor Smith % Instructor/supervisor" << nl
        << "\\end{tabular}" << nl
        << "\\end{center}" << nl
        << nl
        << "% \\begin{abstract}" << nl
        << "Abstract/summary text goes here" << nl
        << "% \\end{abstract}" << nl
        << "\\section{Section heading}" << nl
        << "Lots of really important information will go here!" << nl
        << "\\end{document}" << nl
        << endl;

    // Info<< "    Running pdflatex" << endl;
    // chDir(dataContainer().returnKey<fileName>("processLineDir"));

    // dataContainer_.runSystemCommand
    // (
    //     word("pdflatex " + reportName.name() + " >& log.pdflatex"),
    //     "void Foam::processLine::postProcess()"
    // );
}

// ************************************************************************* //
