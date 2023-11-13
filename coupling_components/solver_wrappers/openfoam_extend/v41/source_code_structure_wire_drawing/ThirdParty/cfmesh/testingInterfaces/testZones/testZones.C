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

    //- find available cell zones
    DynList<word> names;
    pmg.cellZoneNames(names);
    Info << "Cell zones have names " << names << endl;

    forAll(pmg.cells(), cellI)
        Info << "Cell " << cellI << " is in zone "
             << pmg.cellZone(cellI) << endl;

    const label cId = pmg.addCellZone("myFirstCellZone");
    forAll(pmg.cells(), cellI)
    {
        if( cellI % 3 == 0 )
            pmg.addCellToZone(cId, cellI);
    }

    //- find available face zones
    pmg.faceZoneNames(names);
    Info << "Face zones have names " << names << endl;

    forAll(pmg.faces(), faceI)
        Info << "Face " << faceI << " is in zone "
             << pmg.faceZone(faceI) << " with flip "
             << pmg.isFaceFlipped(faceI) << endl;

    const label fId = pmg.addFaceZone("myFirstFaceZone");

    forAll(pmg.faces(), faceI)
    {
        if( faceI % 4 == 0 )
            pmg.addFaceToZone(fId, faceI);
    }

    //- find available point zones
    pmg.pointZoneNames(names);
    Info << "Point zones have names " << names << endl;

    forAll(pmg.points(), pointI)
        Info << "Point " << pointI << " is in zone "
             << pmg.pointZone(pointI) << endl;

    const label pId = pmg.addPointZone("myFirstPointZone");
    Info << "pId " << pId << endl;
    forAll(pmg.points(), pointI)
    {
        if( pointI % 2 )
            continue;

        pmg.addPointToZone(pId, pointI);
    }

    pmg.write();

    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
