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
    Reads the surface mesh, remove the selected facets
    and writes the modified mesh into a new file

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMeshGen.H"
#include "revolve2DMesh.H"
#include "demandDrivenData.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validOptions.insert("rotationAxis", "vector");
    argList::validOptions.insert("origin", "point");
    argList::validOptions.insert("startAngle", "scalar");
    argList::validOptions.insert("endAngle", "scalar");
    argList::validOptions.insert("circumResolution", "integer");
    argList::validOptions.insert("activePatch", "word");

    argList args(argc, argv);

    # include "createTime.H"

    polyMeshGen pmg(runTime);
    pmg.read();

    revolve2DMesh revolver(pmg);

    //- set the rotation origin
    point origin = vector::zero;
    if( args.options().found("origin") )
    {
        origin = vector(IStringStream(args.options()["origin"])());

        Info << "Origin of rotation set to " << origin << endl;
    }
    else
    {
        Info << "Default origin " << origin << endl;
    }

    revolver.setOrigin(origin);

    //- set the rotation axis
    vector rotationAxis(1., 0., 0.);
    if( args.options().found("rotationAxis") )
    {
        rotationAxis = vector(IStringStream(args.options()["rotationAxis"])());

        Info << "Rotation axis set to " << rotationAxis << endl;
    }
    else
    {
        Info << "Default rotation axis " << rotationAxis << endl;
    }
    revolver.setRotationAxis(rotationAxis);

    //- set the resolution
    scalar startAngle(0.0);
    if( args.options().found("startAngle") )
    {
        startAngle = readScalar(IStringStream(args.options()["startAngle"])());
        startAngle *= M_PI / 180.0;

        Info << "Start angle is set to " << startAngle << " radians" << endl;
    }
    else
    {
        Info << "Default starting angle " << startAngle << " radians" << endl;
    }

    scalar endAngle(2.0 * M_PI);
    if( args.options().found("endAngle") )
    {
        endAngle = readScalar(IStringStream(args.options()["endAngle"])());
        endAngle *= M_PI / 180.0;

        Info << "End angle is set to " << endAngle << " radians" << endl;
    }
    else
    {
        Info << "Default end angle " << endAngle << " radians" << endl;
    }

    label circumResolution(80);
    if( args.options().found("circumResolution") )
    {
        circumResolution =
            readLabel(IStringStream(args.options()["circumResolution"])());

        Info << "Circumferential resolution is set to "
             << circumResolution << endl;
    }
    else
    {
        Info << "Default circumferential resolution "
             << circumResolution << endl;
    }

    revolver.setIntervalResolution(startAngle, endAngle, circumResolution);

    //- set the revolving patch
    word revolvingPatchName("bottomEmptyFaces");
    if( args.options().found("activePatch") )
    {
        revolvingPatchName = args.options()["activePatch"];

        Info << "Revolving patch is set to " << revolvingPatchName << endl;
    }
    else
    {
        Info << "Default revolvng patch is " << revolvingPatchName << endl;
    }

    revolver.setRevolvingPatch(revolvingPatchName);

    //- generate the revolved mesh and write the mesh
    revolver.generateRevolvedMesh();

    pmg.write();

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
