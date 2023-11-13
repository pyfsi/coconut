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

#include "dataHandler.H"
#include "objectRegistry.H"
#include "runSystemCommandWithPOpen.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dataHandler, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dataHandler::dataHandler(const Time& runTime)
:
    dataContainerDict_(),
    processLineProgressDict_
    (
        IOobject
        (
            "processLineProgress",
            ".",
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        dictionary("processLineProgressDict")
    )
{
    processLineProgressDict_.regIOobject::write();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::dataHandler::~dataHandler()
{}


// * * * * * * * * * * * * * * * Public Member Functions  * * *  * * * * * * //

void Foam::dataHandler::printData()
{
    Info<< "dataContainerDict_: " << nl << dataContainerDict_ << endl;
}

void Foam::dataHandler::runSystemCommand
(
    const word& systemCommand,
    const word& logName,
    const word& functionCalling,
    const bool checkForSystemError
)
{
    runSystemCommandWithPOpen
    (
        systemCommand,
        logName,
        functionCalling,
        checkForSystemError
    );
}

// ************************************************************************* //
