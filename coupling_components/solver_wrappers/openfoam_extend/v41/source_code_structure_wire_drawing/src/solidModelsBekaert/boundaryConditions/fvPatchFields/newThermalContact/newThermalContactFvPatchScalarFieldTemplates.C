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

#include "newThermalContactFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::newThermalContactFvPatchScalarField::
zoneField
(
    const label zoneIndex,
    const label patchIndex,
    const Field<Type>& patchField
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const faceZone& fZone = mesh.faceZones()[zoneIndex];

    tmp<Field<Type> > tZoneField
    (
        new Field<Type>(fZone.size(), pTraits<Type>::zero)
    );
    scalarField& zField = tZoneField();

    if (!Pstream::parRun())
    {
        zField = patchField;
    }
    else
    {
        const label patchStart
            = mesh.boundaryMesh()[patchIndex].start();

        // put local patchField into global zoneField
        forAll(patchField, i)
        {
            zField[fZone.whichFace(patchStart + i)] = patchField[i];
        }

        reduce(zField, sumOp<Field<Type> >());
    }

    return tZoneField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::newThermalContactFvPatchScalarField::
patchField
(
    const label patchIndex,
    const label zoneIndex,
    const Field<Type>& zoneField
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const faceZone& fZone = mesh.faceZones()[zoneIndex];

    tmp<Field<Type> > tPatchField
    (
        new Field<Type>
        (
            mesh.boundaryMesh()[patchIndex].size(),
            pTraits<Type>::zero
        )
    );
    Field<Type>& pField = tPatchField();

    if (!Pstream::parRun())
    {
        pField = zoneField;
    }
    else
    {
        const label patchStart
            = mesh.boundaryMesh()[patchIndex].start();

        // take local patchField from the global zoneField
        forAll(pField, i)
        {
            pField[i] =
                zoneField[fZone.whichFace(patchStart + i)];
        }
    }

    return tPatchField;
}


// ************************************************************************* //
