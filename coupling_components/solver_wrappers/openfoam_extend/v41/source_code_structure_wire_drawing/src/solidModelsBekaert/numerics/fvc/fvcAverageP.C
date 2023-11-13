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

\*---------------------------------------------------------------------------*/

#include "fvcAverageP.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >average
(
    const GeometricField<Type, pointPatchField, pointMesh>& pf
)
{
    const fvMesh& mesh =
        pf.db().parent().objectRegistry::lookupObject<fvMesh>
        (
            pf.mesh().mesh().name()
        );

    tmp<GeometricField<Type, fvPatchField, volMesh> > tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "average(" + pf.name() + ')',
                pf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("0", pf.dimensions(), pTraits<Type>::zero)
        )
    );
    Field<Type>& vfI = tvf().internalField();

    const labelListList& cellPoints = mesh.cellPoints();

    const Field<Type>& pfI = pf.internalField();

    forAll(vfI, cellI)
    {
        const labelList& curCellPoints = cellPoints[cellI];

        Type& curV = vfI[cellI];

        forAll(curCellPoints, cpI)
        {
            const label pointID = curCellPoints[cpI];

            curV += pfI[pointID];
        }

        curV /= curCellPoints[cpI];
    }

    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();

    forAll(tvf().boundaryField(), patchI)
    {
        Field<Type>& patchV = tvf().boundaryField()[patchI];

        if
        (
            mesh.boundaryMesh()[patchI].type() != "empty"
         && !mesh.boundaryMesh()[patchI].coupled()
        )
        {
            const faceList& localFaces = mesh.boundaryMesh()[patchI];
            const label start = mesh.boundaryMesh()[patchI].start();

            forAll(patchV, fI)
            {
                const label faceID = start + fI;
                const face& curFace = faces[faceID];

                patchV[faceI] = curFace.average(points, pfI);
            }
        }
    }

    tvf().correctBoundaryConditions();

    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
