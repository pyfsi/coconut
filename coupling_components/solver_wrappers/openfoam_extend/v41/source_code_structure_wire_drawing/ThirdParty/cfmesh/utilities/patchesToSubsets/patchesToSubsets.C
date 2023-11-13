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
    Converts specified patches into subsets

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurf.H"
#include "triFaceList.H"
#include "labelLongList.H"
#include "IFstream.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");
    argList::validOptions.insert("patchNames", "list of patch names");
    argList::validOptions.insert
    (
        "matchingString",
        "patch nemes matching the string"
    );
    argList::validOptions.insert("subsetName", "name of the resulting subset");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);

    if( outFileName.ext() != "fms" )
        Warning << "The subsets cann only be saved in the .fms format" << endl;

    triSurf origSurf(inFileName);

    DynList<word> patchNames;

    //- check patch names
    if( args.options().found("patchNames") )
    {
        IStringStream is(args.options()["patchNames"]);
        patchNames = wordList(is);

        Info << "Adding patches with names " << patchNames << endl;
    }
    else if( args.options().found("matchingString") )
    {
        const word ms = args.options()["matchingString"];

        Info << "Adding patches matching the string " << ms << endl;

        forAll(origSurf.patches(), patchI)
        {
            const word& pName = origSurf.patches()[patchI].name();
            if( pName.find(ms) != word::npos )
                patchNames.append(pName);
        }
    }
    else
    {
        Info << "Adding all patches" << endl;
        patchNames.setSize(origSurf.patches().size());
        forAll(origSurf.patches(), patchI)
            patchNames[patchI] = origSurf.patches()[patchI].name();

    }

    //- resulting subset name
    label subsetId(-1);
    if( args.options().found("subsetName") )
    {
        const word sName = args.options()["subsetName"];
        subsetId = origSurf.addFacetSubset(sName);
    }

    forAll(patchNames, patchI)
    {
        labelLongList subsetFacets;
        forAll(origSurf, triI)
        {
            if( origSurf[triI].region() == patchI )
                subsetFacets.append(triI);
        }

        label sId = subsetId;
        if( sId < 0 )
            sId = origSurf.addFacetSubset(patchNames[patchI]);

        forAll(subsetFacets, i)
            origSurf.addFacetToSubset(subsetId, subsetFacets[i]);
    }

    origSurf.writeSurface(outFileName);

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
