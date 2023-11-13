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
    Creates a cell subset from a face subset.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMeshGenModifier.H"
#include "meshOptimizer.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.clear();
    argList::validArgs.append("face subset name");

#   include "setRootCase.H"
#   include "createTime.H"

    const word sName = args.args()[1];

    //- load the mesh from disk
    polyMeshGen pmg(runTime);
    pmg.read();

    const labelLongList& owner = pmg.owner();
    const labelLongList& neighbour = pmg.neighbour();

    const label fId = pmg.faceSubsetIndex(sName);

    if( fId < 0 )
    {
        Warning << "Face subset " << sName << " does not exist" << endl;
        Info << "\nEnd" << endl;
        return 0;
    }

    const label cId = pmg.addCellSubset(sName+"_cells");

    labelLongList facesInSubset;
    pmg.facesInSubset(fId, facesInSubset);

    forAll(facesInSubset, i)
    {
        const label faceI = facesInSubset[i];

        pmg.addCellToSubset(cId, owner[faceI]);

        if( neighbour[faceI] < 0 )
            continue;

        pmg.addCellToSubset(cId, neighbour[faceI]);
    }

    Info << "Writing mesh" << endl;
    pmg.write();

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
