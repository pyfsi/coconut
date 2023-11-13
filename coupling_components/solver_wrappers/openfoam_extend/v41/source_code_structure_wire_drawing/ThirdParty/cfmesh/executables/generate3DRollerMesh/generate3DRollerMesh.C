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
    generate3DRollerMesh

Description
    Generates a 3D mesh of a roller

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "rollingMillMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"

    //- the mesher cannot be run using MPI
    argList::noParallel();

    word regionName = "";
    if( args.optionReadIfPresent("region", regionName) )
    {
        Info<< "Wire mesh will be saved in region " << regionName << nl << endl;
    }

    rollingMillMesh rmg(runTime, regionName);

    rmg.generateRollers();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s\n"
        << "ClockTime = " << runTime.elapsedClockTime() << " s" << endl;

    rmg.writeMesh();

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
