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
    Reads a mesh with many zones and translates the vertices in each zone
    by a given translation vector

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMeshGenModifier.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void translateZone
(
    polyMeshGen& mesh,
    const word zoneName,
    const vector& translationVec
)
{
    const label zoneId = mesh.cellZoneIndex(zoneName);

    if( zoneId < 0 )
    {
        Warning << "Zone " << zoneName << " does no exist. Skipping.." << endl;

        return;
    }

    Info << "Translating vertices in zone " << zoneName
         << " in direction " << translationVec << endl;

    const faceListPMG& faces = mesh.faces();
    const cellListPMG& cells = mesh.cells();

    //- find points in the zone
    std::set<label> pointsInZone;
    forAll(cells, cellI)
    {
        if( mesh.cellZone(cellI) == zoneId )
        {
            const cell& c = cells[cellI];

            forAll(c, fI)
            {
                const face& f = faces[c[fI]];

                forAll(f, pI)
                {
                    pointsInZone.insert(f[pI]);
                }
            }
        }
    }

    //- translate points
    pointFieldPMG& points = polyMeshGenModifier(mesh).pointsAccess();

    forAllConstIter(std::set<label>, pointsInZone, it)
    {
        points[*it] += translationVec;
    }
}

int main(int argc, char *argv[])
{
    argList::validOptions.insert("2DLayers", "bool");

#   include "setRootCase.H"
#   include "createTime.H"

    IOdictionary translationDict
    (
        IOobject
        (
            "translationDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if( !translationDict.found("translateCellZones") )
    {
        FatalError << "translationCellZones does no exist in translationDict"
                   << exit(FatalError);
    }

    //- load the mesh from disk
    polyMeshGen pmg(runTime);
    pmg.read();

    PtrList<entry> translateZones(translationDict.lookup("translateCellZones"));

    forAll(translateZones, zI)
    {
        const dictionary& dict = translateZones[zI].dict();

        word zoneName;
        vector translationVec;

        if( dict.found("zoneName") )
        {
            zoneName = word(dict.lookup("zoneName"));
        }
        else
        {
            FatalError << "zoneName keyword not present in " << dict
                       << abort(FatalError);
        }

        if( dict.found("translationVector") )
        {
            translationVec = vector(dict.lookup("translationVector"));
        }
        else
        {
            FatalError << "translationVector keyword not present in " << dict
                       << abort(FatalError);
        }

        translateZone(pmg, zoneName, translationVec);
    }

    Info << "Writing mesh" << endl;
    pmg.write();

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
