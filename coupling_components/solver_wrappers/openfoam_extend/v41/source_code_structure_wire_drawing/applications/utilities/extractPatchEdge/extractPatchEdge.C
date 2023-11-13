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
    extractPatchEdge

Description
    Extract the list of points on the boundary of a patch and write them to a
    text file. The points are in order around the boundary.

    Operates on the latest time by default.

Author
    Philip Cardiff, UCD. All rights reserved.
    David McAuliffe, Resero. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("patchName", "patchName");

#   include "setRootCase.H"
#   include "createTime.H"

    // Read patchName
    word patchName = "wireDownStream1";
    if (args.optionFound("patchName"))
    {
        patchName = word(args.optionLookup("patchName")());
    }

    // Set the time to the latest time
    instantList Times = runTime.times();
    runTime.setTime(Times[Times.size() - 1], Times.size() - 1);

    // Read the mesh
#   include "createMesh.H"

    // FInd the index fo the patch of interest
    const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    if (patchID == -1)
    {
        FatalError
            << "Patch " << patchName << " not found!" << abort(FatalError);
    }

    // Take a reference to the patch
    const polyPatch& ppatch = mesh.boundaryMesh()[patchID];

    // Get a list of boundary edges on the patch
    const labelListList& edgeLoops = ppatch.edgeLoops();

    // Give an error if there are multiple disconnnected regions on the patch
    if (edgeLoops.size() != 1)
    {
        FatalError
            << "There are multiple disconnnected regions on the patch. This "
            << "utility is only implemented for a patch that is continous in"
            << "space" << nl
            << "There are currently " << edgeLoops.size() << " edge loops"
            << abort(FatalError);
    }

    // Take a reference to the one edge loop
    const labelList& curEdgeLoop = edgeLoops[0];

    // Get list of points on patch
    const pointField& localPoints = ppatch.localPoints();

    // Open a file
    word fName("extractPatchEdge_" + word(runTime.caseName()) + ".txt");
    Info<< "Writing patch boundary points to file " << fName << endl;
    OFstream outFile(fName);

    // Write ordered points to file
    forAll(curEdgeLoop, pI)
    {
        const label pointID = curEdgeLoop[pI];

        outFile
            << localPoints[pointID] << endl;
    }

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
