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
    A testing interface for mesh quality optimisation

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "cartesianMeshGenerator.H"

#include "Random.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    //- generate a reference mesh
    cartesianMeshGenerator cmgRef(runTime);

    const polyMeshGen& meshRef = cmgRef.mesh();

    for(label i=0;i<2;++i)
    {
        cartesianMeshGenerator cmgComp(runTime);

        const polyMeshGen& meshComp = cmgComp.mesh();

        bool areIdentical(true);

        //- compare vertices
        if( meshRef.points().size() != meshComp.points().size() )
        {
            Pout << "Number of points does not match" << endl;
            areIdentical = false;
        }

        if( areIdentical )
        {
            forAll(meshRef.points(), pointI)
            {
                if
                (
                    mag(meshRef.points()[pointI] - meshComp.points()[pointI]) >
                    SMALL
                )
                {
                    Pout << "Points " << pointI << " are not identical "
                         << meshRef.points()[pointI] << " "
                         << meshComp.points()[pointI] << endl;

                    areIdentical = false;
                }
            }
        }

        if( meshRef.faces().size() != meshComp.faces().size() )
        {
            Pout << "Number of faces does not match" << endl;
            areIdentical = false;
        }

        if( areIdentical )
        {
            forAll(meshRef.faces(), faceI)
                if( meshRef.faces()[faceI] != meshComp.faces()[faceI] )
                    Pout << "Faces " << faceI << " are not identical "
                         << meshRef.faces()[faceI] << " "
                         << meshComp.faces()[faceI] << endl;
        }

        if( meshRef.cells().size() != meshComp.cells().size() )
        {
            Pout << "Number of cells does not match" << endl;
            areIdentical = false;
        }

        if( areIdentical )
        {
            forAll(meshRef.cells(), cellI)
                if( meshRef.cells()[cellI] != meshComp.cells()[cellI] )
                {
                    Pout << "Cells " << cellI << " are not identical "
                         << meshRef.cells()[cellI] << " "
                         << meshComp.cells()[cellI] << endl;

                    areIdentical = false;
                }
        }

        if( meshRef.boundaries().size() != meshRef.boundaries().size() )
        {
            Pout << "Number of patches does not match" << endl;
            areIdentical = false;
        }

        if( areIdentical )
        {
            forAll(meshRef.boundaries(), patchI)
            {
                if
                (
                    meshRef.boundaries()[patchI].patchStart() !=
                    meshComp.boundaries()[patchI].patchStart()
                )
                {
                    Pout << "Patches " << patchI << " are not identical "
                         << meshRef.boundaries()[patchI] << " "
                         << meshComp.boundaries()[patchI] << endl;

                    areIdentical = false;
                }
            }
        }

        reduce(areIdentical, minOp<bool>());

        if( areIdentical )
        {
            Info << "Meshes are identical" << endl;
        }
        else
        {
            Info << "Meshes are not identical" << endl;
        }
    }

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
