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
    setupProcessLine

Description
    Setup all ncessary files for the virtual process line.

Author
    Philip Cardiff, UCD
    Peter De Jaeger, Bekaert.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "processLine.H"
#include "dirent.h"
#include "clock.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("restart", "");
    argList::validOptions.insert("setupCasesOnly", "");
    argList args(argc, argv);

    // Check if restart option is specified
    bool restart = false;
    if (args.optionFound("restart"))
    {
        restart = true;
    }

    // Create processLines
    processLine processLine1("processLine", restart);

    if (args.optionFound("setupCasesOnly"))
    {
        // Setup all the passes but do not run them
        processLine1.setup();
    }
    else
    {
        // Setup and runs all the passes
        processLine1.setupAndRun();
    }

    // Perform post-processing steps
    processLine1.postProcess();

    // Print out date and time
    Info<< nl << "Process line completed on " << clock::date() << " at "
        << clock::clockTime() << nl << endl;

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
