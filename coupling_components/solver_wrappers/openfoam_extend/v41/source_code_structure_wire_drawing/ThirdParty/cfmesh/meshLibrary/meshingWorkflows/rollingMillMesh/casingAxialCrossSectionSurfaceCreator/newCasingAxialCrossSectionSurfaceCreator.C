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

\*---------------------------------------------------------------------------*/

#include "casingAxialCrossSectionSurfaceCreator.H"
#include "dictionary.H"
#include "error.H"

//#define DEBUGDie

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

autoPtr<Foam::casingAxialCrossSectionSurfaceCreator>
Foam::casingAxialCrossSectionSurfaceCreator::New
(
    const rollingMillPatchNamesHandler& patchNames,
    const dictionary& dict,
    const scalar tol
)
{
    if( debug )
    {
        Info<< "casingAxialCrossSectionSurfaceCreator::New("
            << "const rollingMillPatchNamesHandler& patchNames,"
            << " const dictionary& dict,"
            << " const scalar tol) : "
            << "constructing casingAxialCrossSectionSurfaceCreator"
            << endl;
    }

    // default type is self
    word cmType(typeName_());
    if( dict.found("type") )
    {
        dict.lookup("type") >> cmType;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(cmType);

    if( cstrIter == dictionaryConstructorTablePtr_->end() )
    {
        FatalIOErrorIn
        (
            "casingAxialCrossSectionSurfaceCreator::New("
            " const rollingMillPatchNamesHandler& patchNames,"
            " const dictionary& dict, const scalar tol)",
            dict
        )   << "Unknown casingAxialCrossSectionSurfaceCreator type "
            << cmType << nl << nl
            << "Valid casingAxialCrossSectionSurfaceCreator types are :" << nl
            << "[default: " << typeName_() << "]"
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<casingAxialCrossSectionSurfaceCreator>
    (
        cstrIter()(patchNames, dict, tol)
    );
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
