/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2005-2007 Franjo Juretic
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    Test for smoothers

Description
    - reads the mesh and tries to untangle negative volume cells

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMeshGen.H"
#include "revolve2DMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    polyMeshGen pmg(runTime);

    Info << "Starting reading mesh" << endl;
    pmg.read();
    Info << "Finished reading mesh" << endl;

    const point origin(0., -0.1, 0.);
    const vector axis(1., 0., 0.);
    const word pName = "bottomEmptyFaces";

    revolve2DMesh revolver(pmg);
    revolver.setOrigin(origin);
    revolver.setRotationAxis(axis);
    revolver.setRevolvingPatch(pName);

    revolver.setCircumResolution(20);
    revolver.setIntervalResolution(0, M_PI/4, 20);
    revolver.setIntervalResolution(M_PI, 1.5*M_PI, 30);
    revolver.setMaxRatio(1.2);

    revolver.generateRevolvedMesh();

    pmg.write();

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
