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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template
<
    class Type,
    template<class> class PatchFieldType,
    class MeshType,
    class FieldMeshType
>
void Foam::geometricFieldContainer::readFields
(
    PtrList< GeometricField<Type, PatchFieldType, MeshType> >& ptrList,
    const IOobjectList& objects,
    const FieldMeshType& mesh
)
{
    // Lookup class name e.g. volScalarField, pointVectorField
    const word fieldClassName
    (
        GeometricField<Type, PatchFieldType, MeshType>::typeName
    );

    // Find the objects that are of this field type
    IOobjectList fieldNames = objects.lookupClass(fieldClassName);

    // Remove oldTime fields as they will be read with the primary field
    for
    (
        IOobjectList::iterator fieldIter = fieldNames.begin();
        fieldIter != fieldNames.end();
        ++fieldIter
    )
    {
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
                Info<< "Ignoring " << fieldName << endl;
                fieldNames.erase(fieldIter);
            }
        }
    }

    // Set size of storage list
    ptrList.setSize(fieldNames.size());

    // Read in fields and store them
    label i = 0;
    for
    (
        IOobjectList::iterator fieldIter = fieldNames.begin();
        fieldIter != fieldNames.end();
        ++fieldIter
    )
    {
        //if (debug)
        {
            Info<< "    Reading " << fieldIter()->name() << endl;
        }

        // Read field and store in list
        ptrList.set
        (
            i,
            new GeometricField<Type, PatchFieldType, MeshType>
            (
                *fieldIter(),
                mesh
            )
        );

        i++;
    }
}


// ************************************************************************* //
