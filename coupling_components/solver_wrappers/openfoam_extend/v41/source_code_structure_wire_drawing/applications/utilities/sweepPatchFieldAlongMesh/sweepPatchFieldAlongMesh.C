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
    sweepPatchFieldAlongMesh

Description
    Reads a patch field and then sweeps this field along mesh cell-streamlines,
    then writes out the updated field.

    This only really makes sense for meshes that are generated in a
    swept/extruded manner.


Author
   Philip Cardiff, UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wireStreamlines.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("fieldName");
    argList::validArgs.append("patchName");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    // Read arguments
    const word fieldName(args.additionalArgs()[0]);
    const word patchName(args.additionalArgs()[1]);

    // Check the patch name
    const label patchID = mesh.boundaryMesh().findPatchID(patchName);

    // Create a mesh streamline object
    wireStreamlines streamlines(mesh, patchName);

    // Get the cell streamlines
    const List< SLList<label> >& cellStreamlines =
        streamlines.cellStreamlines();

    Info<< "There are " << cellStreamlines[0].size()
        << " cells along the first mesh streamline" << endl;

    // Read the yield strength file
    Info<< "Reading field " << fieldName << endl;
    volScalarField vf
    (
        IOobject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Set first layer of cells values from the patch: these will be the first
    // cells in the streamlines
    scalarField& vfI = vf.internalField();
    const scalarField& vfP = vf.boundaryField()[patchID];

    // Set values along streamlines
    // Note: there is one streamline for each face on the patch
    forAll(cellStreamlines, faceI)
    {
        // Convert to standard list
        const labelList curStreamline(cellStreamlines[faceI]);

        forAll(curStreamline, cI)
        {
            const label cellID = curStreamline[cI];
            vfI[cellID] = vfP[faceI];
        }
    }

    // Extrapolated (zero-gradient) boundary values
    forAll(vf.boundaryField(), patchI)
    {
        if (vf.boundaryField()[patchI].coupled())
        {
            vf.boundaryField()[patchI].evaluate();
        }
        else if (vf.boundaryField()[patchI].type() != "empty")
        {
            vf.boundaryField()[patchI] =
                vf.boundaryField()[patchI].patchInternalField();
        }
    }

    // Write the field
    Info<< "Writing field " << fieldName << endl;
    vf.write();

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
