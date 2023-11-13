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
    testNagataPatch

Description
    Simple utility for testing Nagata patch implementation.

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "nagataPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("surfaceFile");
    argList::validArgs.append("nRefinmentLevels");
    argList::validOptions.insert("patches", "patches");
    argList::validOptions.insert("epsilon1", "epsilon1");
    argList::validOptions.insert("epsilon2", "epsilon1");

#   include "setRootCase.H"


    // Read triSurface
    fileName surfName(args.additionalArgs()[0]);
    Info<< nl << "Reading " << surfName << endl;
    triSurface triSurf(surfName);

    // Read list of patches to refine and project
    wordList patchesToProject;
    if (args.optionFound("patches"))
    {
        patchesToProject = wordList(args.optionLookup("patches")());
        Info<< "    refining the patches = " << patchesToProject << endl;
    }
    else
    {
        Info<< "    refining all the patches" << endl;

        // Mark all patches for projection
        patchesToProject.resize(triSurf.patches().size(), "");
        forAll(triSurf.patches(), patchI)
        {
            patchesToProject[patchI] = triSurf.patches()[patchI].name();
        }
    }


    // Read curvature control parameters
    scalar epsilon1 = 0.05;
    if (args.optionFound("epsilon1"))
    {
        epsilon1 = readScalar(args.optionLookup("epsilon1")());
    }
    scalar epsilon2 = 0.025;
    if (args.optionFound("epsilon2"))
    {
        epsilon2 = readScalar(args.optionLookup("epsilon2")());
    }


    // Create nagata patch
    nagataPatch nagataSurf(triSurf, epsilon1, epsilon2);


    // Read number of refinement levels
    //const int nRefinmentLevels(IStringStream(args.additionalArgs()[1])());
    const int nRefinmentLevels
    (
        readInt(IStringStream(args.additionalArgs()[1])())
    );

    if (nRefinmentLevels < 1)
    {
        FatalError
            << "nRefinmentLevels should be greater than or equal to 1" << endl;
    }


    // Refine and project the patches
    triSurface refinedAndProjectedSurf =
        nagataSurf.refineAndProjectPatches(patchesToProject, nRefinmentLevels);


    // Write out the refined and projected surface
    const fileName newSurfName(surfName.lessExt() + "refinedAndProjected.stl");
    Info<< nl << "Writing " << newSurfName << nl << endl;
    refinedAndProjectedSurf.write(newSurfName);

    Info<< nl << "End" << nl << endl;

    return(0);
}


// ************************************************************************* //
