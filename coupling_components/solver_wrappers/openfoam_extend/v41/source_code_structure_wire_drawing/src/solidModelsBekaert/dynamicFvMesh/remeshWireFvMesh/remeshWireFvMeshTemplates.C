/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "remeshWireFvMesh.H"
#include "calculatedFvPatchFields.H"
#include "GeometricField.H"
#include "meshToMesh.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class Type>
void Foam::remeshWireFvMesh::newMapConsistentVolFields
(
    const IOobjectList& objects,
    const meshToMesh& meshToMeshInterp
)
{
    const fvMesh& meshSource = meshToMeshInterp.fromMesh();
    const fvMesh& meshTarget = meshToMeshInterp.toMesh();

    const word fieldClassName
    (
        GeometricField<Type, fvPatchField, volMesh>::typeName
    );

    IOobjectList fields = objects.lookupClass(fieldClassName);

    // Remove oldTime fields as they will be read with the primary field
    // for
    // (
    //     IOobjectList::iterator fieldIter = fields.begin();
    //     fieldIter != fields.end();
    //     ++fieldIter
    // )
    // {
    //     const word& fieldName = fieldIter()->name();
    //     const label fieldNameSize = fieldName.size();

    //     if (fieldNameSize > 2)
    //     {
    //         if
    //         (
    //             fieldName[fieldNameSize - 2] == '_'
    //          && fieldName[fieldNameSize - 1] == '0'
    //         )
    //         {
    //             Info<< "Ignoring " << fieldName << endl;
    //             fields.erase(fieldIter);
    //         }
    //     }
    // }

    for
    (
        IOobjectList::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        if (debug)
        {
            Info<< "    interpolating " << fieldIter()->name()
                << endl;
        }

        // Read field
        GeometricField<Type, fvPatchField, volMesh> fieldSource
        (
            *fieldIter(),
            meshSource
        );

        IOobject fieldTargetIOobject
        (
            fieldIter()->name(),
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (fieldTargetIOobject.headerOk())
        {
            // Read fieldTarget
            GeometricField<Type, fvPatchField, volMesh> fieldTarget
            (
                fieldTargetIOobject,
                meshTarget
            );

            // Interpolate field
            meshToMeshInterp.interpolate
            (
                fieldTarget,
                fieldSource,
                meshToMesh::INTERPOLATE
            );

            // // Write field
            // fieldTarget.write();

            // Only the master writes the field
            returnReduce(bool(true), orOp<bool>());
            if (Pstream::master())
            {
                if (fieldTarget.name()[0] != '(')
                {
                    fieldTarget.write();
                }
            }
            returnReduce(bool(true), orOp<bool>());
        }
        else
        {
            fieldTargetIOobject.readOpt() = IOobject::NO_READ;

            // Interpolate field
            GeometricField<Type, fvPatchField, volMesh> fieldTarget
            (
                fieldTargetIOobject,
                meshToMeshInterp.interpolate
                (
                    fieldSource,
                    meshToMesh::INTERPOLATE
                )
            );

            // // Write field
            // fieldTarget.write();

            // Only the master writes the field
            returnReduce(bool(true), orOp<bool>());
            if (Pstream::master())
            {
                if (fieldTarget.name()[0] != '(')
                {
                    fieldTarget.write();
                }
            }
            returnReduce(bool(true), orOp<bool>());
        }
    }
}


template<class Type>
void Foam::remeshWireFvMesh::newMapConsistentVolFields
(
    const meshToMesh& meshToMeshInterp
)
{
    const fvMesh& meshSource = meshToMeshInterp.fromMesh();
    const fvMesh& meshTarget = meshToMeshInterp.toMesh();

    // Read volField objects from object registry
    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields
    (
        meshSource.thisDb().objectRegistry::template lookupClass
        <GeometricField<Type, fvPatchField, volMesh> >()
    );

    for
    (
        typename HashTable
        <
            const GeometricField<Type, fvPatchField, volMesh>*
        >::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        if (debug)
        {
            Info<< "    interpolating " << fieldIter()->name()
                << endl;
        }

        // Read field
        GeometricField<Type, fvPatchField, volMesh>& fieldSource =
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (*fieldIter());

        IOobject fieldTargetIOobject
        (
            fieldIter()->name(),
            meshTarget.time().timeName(),
            meshTarget,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (fieldTargetIOobject.headerOk())
        {
            // Read fieldTarget
            GeometricField<Type, fvPatchField, volMesh> fieldTarget
            (
                fieldTargetIOobject,
                meshTarget
            );

            // Interpolate field
            meshToMeshInterp.interpolate
            (
                fieldTarget,
                fieldSource,
                meshToMesh::INTERPOLATE
            );

            // // Write field
            // fieldTarget.write();

            // Only the master writes the field
            returnReduce(bool(true), orOp<bool>());
            if (Pstream::master())
            {
                if (fieldTarget.name()[0] != '(')
                {
                    fieldTarget.write();
                }
            }
            returnReduce(bool(true), orOp<bool>());
        }
        else
        {
            fieldTargetIOobject.readOpt() = IOobject::NO_READ;

            // Interpolate field
            GeometricField<Type, fvPatchField, volMesh> fieldTarget
            (
                fieldTargetIOobject,
                meshToMeshInterp.interpolate
                (
                    fieldSource,
                    meshToMesh::INTERPOLATE
                )
            );

            // // Write field
            // fieldTarget.write();

            // Only the master writes the field
            returnReduce(bool(true), orOp<bool>());
            if (Pstream::master())
            {
                if (fieldTarget.name()[0] != '(')
                {
                    fieldTarget.write();
                }
            }
            returnReduce(bool(true), orOp<bool>());
        }
    }
}


template
<
    class Type,
    template<class> class PatchFieldType,
    class MeshType,
    class FieldMeshType
>
void Foam::remeshWireFvMesh::ReplaceOldMeshFieldsWithNewMeshFields
(
    const fvMesh& oldMesh,
    const fvMesh& newMesh,
    const FieldMeshType& oldFieldMesh,
    const FieldMeshType& newFieldMesh,
    const bool writeFields
)
{
    // Read volField objects from object registry
    HashTable<const GeometricField<Type, PatchFieldType, MeshType>*> fields
    (
        oldMesh.thisDb().objectRegistry::template lookupClass
        <GeometricField<Type, PatchFieldType, MeshType> >()
    );

    // Remove oldTime fields as hey will be read with the primary field
    for
    (
        typename HashTable
        <
            const GeometricField<Type, PatchFieldType, MeshType>*
        >::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        // GeometricField<Type, PatchFieldType, MeshType>& field =
        //     const_cast<GeometricField<Type, PatchFieldType, MeshType>&>
        //     (*fieldIter());

        // field.storeOldTimes();

        // const word& fieldName = fieldIter()->name();
        // const label fieldNameSize = fieldName.size();

        // if (fieldNameSize > 2)
        // {
        //     if
        //     (
        //         fieldName[fieldNameSize - 2] == '_'
        //      && fieldName[fieldNameSize - 1] == '0'
        //     )
        //     {
        //         if (debug)
        //         {
        //             Info<< "Ignoring " << fieldName << endl;
        //         }

        //         fields.erase(fieldIter);
        //     }
        // }
    }

    // Loop through all of the volFields
    for
    (
        typename HashTable
        <
            const GeometricField<Type, PatchFieldType, MeshType>*
        >::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        // Read field
        GeometricField<Type, PatchFieldType, MeshType>& field =
            const_cast<GeometricField<Type, PatchFieldType, MeshType>&>
            (*fieldIter());

        // Should we skip oldTime and prevIter fields and update them with the
        // fields themselves? For now, we will leave it.

        if (debug)
        {
            Info<< "        Updating " << field.name();
        }

        // Set field to zero and correct length
        // Note: multiplying by 1.0 fixes some strange memory freeing error
        // with surface fields
        field =
            1.0*GeometricField<Type, PatchFieldType, MeshType>
            (
                IOobject
                (
                    field.name(),
                    oldMesh.time().timeName(),
                    oldMesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                oldFieldMesh,
                dimensioned<Type>
                (
                    "zero",
                    field.dimensions(),
                    pTraits<Type>::zero
                )
            );

        // Attempt to read newMesh field from disk
        IOobject newMeshFieldIOobject
        (
            field.name(),
            newMesh.time().timeName(),
            newMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        // Check if the field exists on the newMesh
        if (newMeshFieldIOobject.headerOk())
        {
            // Read new mesh field from disk
            const GeometricField<Type, PatchFieldType, MeshType>
                newMeshField
                (
                    newMeshFieldIOobject,
                    newFieldMesh
                );

            // Overwrite internal field
            field.internalField() = newMeshField.internalField();

            // Overwrite boundary field
            forAll(oldMesh.boundaryMesh(), patchI)
            {
                field.boundaryField()[patchI] =
                    Field<Type>(newMeshField.boundaryField()[patchI]);
            }
        }
        else
        {
            if (debug)
            {
                Info<< "    field resized and initialised to zero" << endl;
            }
        }

        // No need to write
        if (writeFields)
        {
            if (debug)
            {
                Info<< ":    writing to time = " << oldMesh.time().value()
                    << endl;
            }

            returnReduce(bool(true), orOp<bool>());
            if (Pstream::master())
            {
                field.write();
            }
            returnReduce(bool(true), orOp<bool>());
        }
    }
}


template<class Type>
void Foam::remeshWireFvMesh::ReplaceOldMeshPointFieldsWithNewMeshPointFields
(
    const fvMesh& oldMesh,
    const fvMesh& newMesh,
    const pointMesh& oldPMesh,
    const pointMesh& newPMesh,
    const bool writeFields
)
{
    // Read objects from object registry
    HashTable<const GeometricField<Type, pointPatchField, pointMesh>*> fields
    (
        oldMesh.thisDb().objectRegistry::template lookupClass
        <GeometricField<Type, pointPatchField, pointMesh> >()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.
    for
    (
        typename HashTable
        <
            const GeometricField<Type, pointPatchField, pointMesh>*
        >::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, pointPatchField, pointMesh>& field =
            const_cast<GeometricField<Type, pointPatchField, pointMesh>&>
            (*fieldIter());

        field.storeOldTimes();
    }

    // Loop through all of the fields
    for
    (
        typename HashTable
        <
            const GeometricField<Type, pointPatchField, pointMesh>*
        >::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        // Read field
        GeometricField<Type, pointPatchField, pointMesh>& field =
            const_cast<GeometricField<Type, pointPatchField, pointMesh>&>
            (*fieldIter());

        if (debug)
        {
            Info<< "        Updating " << field.name();
        }

        // Set field to zero and correct length
        // Note: multiplying by 1.0 fixes some strange memory freeing error with
        // surface fields: maybe this is resolved...
        field =
            1.0*GeometricField<Type, pointPatchField, pointMesh>
            (
                IOobject
                (
                    field.name(),
                    oldMesh.time().timeName(),
                    oldMesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                oldPMesh,
                dimensioned<Type>
                (
                    "zero",
                    field.dimensions(),
                    pTraits<Type>::zero
                )
            );

        // Attempt to read newMesh field from disk
        IOobject newMeshFieldIOobject
        (
            field.name(),
            newMesh.time().timeName(),
            newMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        // Check if the field exists on the newMesh
        if (newMeshFieldIOobject.headerOk())
        {
            // Read new mesh field from disk
            const GeometricField<Type, pointPatchField, pointMesh> newMeshField
            (
                newMeshFieldIOobject,
                newPMesh
            );

            // Overwrite internal field
            field.internalField() = newMeshField.internalField();

            // No boundary update for point fields
        }
        else
        {
            if (debug)
            {
                Info<< "    field resized and initialised to zero" << endl;
            }
        }

        if (writeFields)
        {
            if (debug)
            {
                Info<< ":    writing to time = " << oldMesh.time().value()
                    << endl;
            }

            returnReduce(bool(true), orOp<bool>());
            if (Pstream::master())
            {
                field.write();
            }
            returnReduce(bool(true), orOp<bool>());
        }
    }
}


template<class Type, template<class> class PatchFieldType, class MeshType>
void Foam::remeshWireFvMesh::WriteAllObjects(const fvMesh& mesh) const
{
    // Read objects from object registry
    HashTable<const GeometricField<Type, PatchFieldType, MeshType>*> fields
    (
        mesh.thisDb().objectRegistry::template lookupClass
        <GeometricField<Type, PatchFieldType, MeshType> >()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.
    for
    (
        typename HashTable
        <
            const GeometricField<Type, PatchFieldType, MeshType>*
        >::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        // GeometricField<Type, PatchFieldType, MeshType>& field =
        //     const_cast<GeometricField<Type, PatchFieldType, MeshType>&>
        //     (*fieldIter());
        //field.storeOldTimes();

        const word& fieldName = fieldIter()->name();
        const label fieldNameSize = fieldName.size();

        if (fieldNameSize > 2)
        {
            if
            (
                fieldName[fieldNameSize - 2] == '_'
             && fieldName[fieldNameSize - 1] == '0'
            )
            {
                if (debug)
                {
                    Info<< "Ignoring " << fieldName << endl;
                }

                fields.erase(fieldIter);
            }
        }
    }

    // Loop through all of the volFields
    for
    (
        typename HashTable
        <
            const GeometricField<Type, PatchFieldType, MeshType>*
        >::iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        // Read field
        GeometricField<Type, PatchFieldType, MeshType>& field =
            const_cast<GeometricField<Type, PatchFieldType, MeshType>&>
            (*fieldIter());

        // Write the field to disk, if its name does not start with a
        // bracket
        if (field.name().size() > 0)
        {
            if (field.name()[0] != '(')
            {
                if (debug)
                {
                    Info<< "    Writing " << field.name() << endl;
                }
                field.write();
            }
        }
    }
}


// ************************************************************************* //
