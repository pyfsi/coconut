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

Description
    Scales the mesh into other units.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMeshGen.H"
#include "polyMeshGenModifier.H"
#include "helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::validArgs.clear();
    argList::noParallel();
    argList::validArgs.append("regions");
    argList::validOptions.insert("noFunctionObjects", "");

    argList args(argc, argv, false, true);

    if( args.additionalArgs().size() < 2 )
    {
        FatalError
            << "At least two regions should be specified"
            << abort(FatalError);
    }

#   include "createTime.H"

    Info << "Reading sub-meshes" << endl;
    DynList<word> submeshes;
    forAll(args.additionalArgs(), argI)
        submeshes.append(args.additionalArgs()[argI]);

    //- read the mesh from disk
    polyMeshGen pmg(runTime);

    polyMeshGenModifier meshModifier(pmg);

    forAll(submeshes, i)
    {
        Info << "Merging sub-mesh " << submeshes[i] << endl;
        polyMeshGen subMesh
        (
            runTime,
            "constant",
            submeshes[i]+"/polyMesh"
        );

        subMesh.read();

        meshModifier.addMesh(subMesh);
    }

    pmg.write();

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
