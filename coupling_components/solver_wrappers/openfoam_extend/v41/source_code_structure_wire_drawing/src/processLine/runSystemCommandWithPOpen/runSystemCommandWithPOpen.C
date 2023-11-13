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

\*---------------------------------------------------------------------------*/

#include "runSystemCommandWithPOpen.H"
#include <stdio.h>
#include "messageStream.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::runSystemCommandWithPOpen
(
    const word& systemCommand,
    const word& logName,
    const word& functionCalling,
    const bool checkForSystemError
)
{
    // From helperFunctionsAuxI.H in cfMesh code: avoid problems with
    // running system commands on different systems and shells

    Info<< "            Running system command: " << systemCommand << endl;

    // Append stderr to the log file. A problem seems to be that the error
    // code does not get passed so the while loop below gets stuck...
    // const word systemCommandAndRedirectStdErr = systemCommand + " 2>&1";

    // Cast word to char* for command
    //const char* progName = systemCommandAndRedirectStdErr.c_str();
    const char* progName = systemCommand.c_str();

    // Execute the command but do not wait for it to finish
    FILE *handle = popen(progName, "r");

    // If handle is zero, it means the command gave an error
    if (handle == 0)
    {
        if (checkForSystemError)
        {
            FatalError
                << nl
                << "    >>> Failed system command: "
                << word(progName)
                << abort(FatalError);
        }
        else
        {
            Warning
                << nl
                << "    >>> Failed system command: "
                << word(progName)
                << endl;
        }
    }

    // Create and write the log file

    FILE* filePtr = NULL;
    if (logName != "")
    {
        filePtr = fopen(logName.c_str(), "wa");
    }

    char buf[128];
    size_t readn;
    while ( (readn = fread(buf, 1, sizeof(buf), handle)) > 0 )
    {
        fwrite(buf, 1, readn, filePtr?filePtr:stdout);

        // PC, PDJ, 12Apr18: force buffer to be written to file
        fflush(filePtr);
    }

    // Flush all data before trying to close the file
    fflush(filePtr);

    if (pclose(handle) && checkForSystemError)
    {
        FatalError
            << nl
            << "    >>> There was an error running the system command: "
            << word(progName)
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
