/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    coconut_foamRun

Description
    Loads and executes an OpenFOAM solver module either specified by the
    optional \c solver entry in the \c controlDict or as a command-line
    argument.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

Usage
    \b foamRun [OPTION]

      - \par -solver <name>
        Solver name

      - \par -libs '(\"lib1.so\" ... \"libN.so\")'
        Specify the additional libraries loaded

    Example usage:
      - To run a \c rhoPimpleFoam case by specifying the solver on the
        command line:
        \verbatim
            foamRun -solver fluid
        \endverbatim

      - To update and run a \c rhoPimpleFoam case add the following entries to
        the controlDict:
        \verbatim
            application     foamRun;

            solver          fluid;
        \endverbatim
        then execute \c foamRun

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "solver.H"
#include "pimpleSingleRegionControl.H"
#include "setDeltaT.H"

using namespace Foam;

//#include "Time.H"
#include "fsiDisplacement.H"
#include "waitForSync.H"
#include <unistd.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "solver",
        "name",
        "Solver name"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    // Read the solverName from the optional solver entry in controlDict
    word solverName
    (
        runTime.controlDict().lookupOrDefault("solver", word::null)
    );

    // Optionally reset the solver name from the -solver command-line argument
    args.optionReadIfPresent("solver", solverName);

    // Check the solverName has been set
    if (solverName == word::null)
    {
        args.printUsage();

        FatalErrorIn(args.executable())
            << "solver not specified in the controlDict or on the command-line"
            << exit(FatalError);
    }
    else
    {
        // Load the solver library
        solver::load(solverName);
    }

    // Create the default single region mesh
    #include "createMesh.H"

    // Instantiate the selected solver
    autoPtr<solver> solverPtr(solver::New(solverName, mesh));
    solver& solver = solverPtr();

    // Create the outer PIMPLE loop and control structure
    pimpleSingleRegionControl pimple(solver.pimple);

    // Set the initial time-step
    setDeltaT(runTime, solver);

    // Define coupling iteration counter
    unsigned int iteration;
    iteration = 0;

    // Read CoCoNuT controls
    #include "readCoconutControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting loop\n" << endl;

    while (true)
    {
        usleep(1000); // Expressed in microseconds

        if (exists("next.coco"))
        {
            waitForSync("next"); // Keep the sync at the beginning of the block

            runTime.run();

            solver.preSolve();

            // Adjust the time-step according to the solver maxDeltaT -- NOT SUPPORTED BY COCONUT
            adjustDeltaT(runTime, solver);

            runTime++;

            // Reset coupling iteration counter
            iteration = 0;

            Info<< "Time = " << runTime.userTimeName() << nl << endl;
        }

        if (exists("continue.coco"))
        {
            waitForSync("continue");

            // Increment coupling iteration counter
            iteration++;
            Info << "Coupling iteration = " << iteration << nl << endl;

            // Define movement of the coupling interface
            forAll(boundaryNames, s)
            {
                word boundaryName = boundaryNames[s];
                ApplyFSIPointDisplacement(mesh, boundaryName);
            }

            // PIMPLE corrector loop
            while (pimple.loop())
            {
                solver.moveMesh();
                solver.fvModels().correct();
                solver.prePredictor();
                solver.momentumPredictor();
                solver.thermophysicalPredictor();
                solver.pressureCorrector();
                solver.postCorrector();

                // Check coupling convergence
                if (checkCouplingConvergence && pimple.firstIter())
                {
                    #include "checkCouplingConvergence.H"
                }
            }

            solver.postSolve();

            // Return the coupling interface output
            #include "executeCoconutFunctionObjects.H"

            Info << "Coupling iteration " << iteration << " end" << nl << endl;

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        }

        if (exists("save.coco"))
        {
            waitForSync("save");

            runTime.write();
        }

        if (exists("stop.coco"))
        {
            waitForSync("stop");

            break;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
