/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Application
    Generates a mesh of a wire

Description
    Generates a 3D mesh representing a wire

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "rollingMillMesh.H"
#include "helperFunctions.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"

    //- wire mesher cannot be run in parallel
    argList::noParallel();

    //- get times list
    instantList Times = runTime.times();

    label startTime = -1;

    // check -time and -latestTime options
    if (args.optionFound("time"))
    {
        Foam::scalar timeValue = args.optionRead<scalar>("time");

        startTime = Foam::Time::findClosestTimeIndex(Times, timeValue);
    }

    if (args.optionFound("latestTime"))
    {
        startTime = Times.size() - 1;
    }

    //- read the requested time step to write the mesh
    fileName timeStep = "";
    if( startTime != -1 )
    {
        runTime.setTime(Times[startTime], startTime);
        timeStep = runTime.timeName();

        Info << "Write to time step " << timeStep << endl;
    }

    //- read the region of the mesh
    word regionName = "";
    if( args.optionReadIfPresent("region", regionName) )
    {
        Info<< "Wire mesh will be saved in region " << regionName << nl << endl;
    }

    rollingMillMesh rmg(runTime, regionName, timeStep);

    rmg.generateWire();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s\n"
        << "ClockTime = " << runTime.elapsedClockTime() << " s" << endl;

    rmg.writeMesh();

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
