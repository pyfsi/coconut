/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mechanicalLaw.H"

// * * * * * * * * * * *  Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::mechanicalLaw::SetFieldPtr
(
    GeometricField<Type, fvPatchField, volMesh>*& fieldPtr,
    const word& fieldName,
    const dimensioned<Type>& defaultValue,
    const word& patchTypes
)
{
    if (fieldPtr)
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "void Foam::mechanicalLaw::SetFieldPtr\n"
            "(\n"
            "    (Type*)& fieldPtr,\n"
            "    const word& fieldName,\n"
            "    const dimensioned<Type>& defaultValue\n"
            ")"
        )   << "pointer already set" << abort(FatalError);
    }

    if
    (
        mesh().foundObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            fieldName
        )
    )
    {
        // Lookup the field from the objectRegistry and set the pointer
        fieldPtr =
           &const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (
                mesh().lookupObject
                <
                    GeometricField<Type, fvPatchField, volMesh>
                >(fieldName)
            );

        Info<< "Looking up field " << fieldName << " from the objectRegistry"
            << " and setting the field pointer" << endl;
    }
    else
    {
        // Create the field and read from the disk if present
        fieldPtr =
            new GeometricField<Type, fvPatchField, volMesh>
            (
                IOobject
                (
                    fieldName,
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                defaultValue,
                patchTypes
            );

        Info<< "Creating field " << fieldName << endl;
    }
}


template<class Type>
void Foam::mechanicalLaw::SetFieldPtr
(
    GeometricField<Type, fvsPatchField, surfaceMesh>*& fieldPtr,
    const word& fieldName,
    const dimensioned<Type>& defaultValue,
    const word& patchTypes
)
{
    if (fieldPtr)
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "void Foam::mechanicalLaw::SetFieldPtr\n"
            "(\n"
            "    (Type*)& fieldPtr,\n"
            "    const word& fieldName,\n"
            "    const dimensioned<Type>& defaultValue\n"
            ")"
        )   << "pointer already set" << abort(FatalError);
    }

    if
    (
        mesh().foundObject<GeometricField<Type, fvsPatchField, surfaceMesh> >
        (
            fieldName
        )
    )
    {
        // Lookup the field from the objectRegistry and set the pointer
        fieldPtr =
           &const_cast<GeometricField<Type, fvsPatchField, surfaceMesh>&>
            (
                mesh().lookupObject
                <
                    GeometricField<Type, fvsPatchField, surfaceMesh>
                >(fieldName)
            );

        Info<< "Looking up field " << fieldName << " from the objectRegistry"
            << " and setting the field pointer" << endl;
    }
    else
    {
        // Create the field and read from the disk if present
        fieldPtr =
            new GeometricField<Type, fvsPatchField, surfaceMesh>
            (
                IOobject
                (
                    fieldName,
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                defaultValue,
                patchTypes
            );

        Info<< "Creating field " << fieldName << endl;
    }
}


template<class Type>
void Foam::mechanicalLaw::CleanFieldPtr
(
    GeometricField<Type, fvPatchField, volMesh>*& fieldPtr,
    const word& fieldName
)
{
    // Check if the field exists in the object registry
    if
    (
        mesh().foundObject<GeometricField<Type, fvPatchField, volMesh> >
        (
            fieldName
        )
    )
    {
        // Lookup the field from the objectRegistry and delete it

        Info<< "Removing the field " << fieldName << " from the objectRegistry"
            << " and setting the field pointer to NULL" << endl;

        fieldPtr =
           &const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (
                mesh().lookupObject
                <
                    GeometricField<Type, fvPatchField, volMesh>
                >(fieldName)
            );

        deleteDemandDrivenData(fieldPtr);
    }
    else
    {
        Info<< "The field " << fieldName << " does not exist in the "
            "objectRegistry: setting the field pointer to NULL" << endl;

        // Set the pointer to NULL
        fieldPtr = NULL;
    }
}


template<class Type>
void Foam::mechanicalLaw::CleanFieldPtr
(
    GeometricField<Type, fvsPatchField, surfaceMesh>*& fieldPtr,
    const word& fieldName
)
{
    // Check if the field exists in the object registry
    if
    (
        mesh().foundObject<GeometricField<Type, fvsPatchField, surfaceMesh> >
        (
            fieldName
        )
    )
    {
        // Lookup the field from the objectRegistry and delete it

        Info<< "Removing the field " << fieldName << " from the objectRegistry"
            << " and setting the field pointer to NULL" << endl;

        fieldPtr =
           &const_cast<GeometricField<Type, fvsPatchField, surfaceMesh>&>
            (
                mesh().lookupObject
                <
                    GeometricField<Type, fvsPatchField, surfaceMesh>
                >(fieldName)
            );

        deleteDemandDrivenData(fieldPtr);
    }
    else
    {
        Info<< "The field " << fieldName << " does not exist in the "
            "objectRegistry: setting the field pointer to NULL" << endl;

        // Set the pointer to NULL
        fieldPtr = NULL;
    }
}


// ************************************************************************* //
