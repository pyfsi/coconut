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

#ifndef MergeSurfaceFields_H
#define MergeSurfaceFields_H

#include "GeometricField.H"
#include "IOobjectList.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void MergeSurfaceFields
(
    const IOobjectList& masterObjects,
    const IOobjectList& slaveObjects,
    const fvMesh& masterMesh,
    const fvMesh& slaveMesh
)
{
    word fieldClassName
    (
        GeometricField<Type, fvsPatchField, surfaceMesh>::typeName
    );

    IOobjectList masterFields = masterObjects.lookupClass(fieldClassName);
    IOobjectList slaveFields = slaveObjects.lookupClass(fieldClassName);

    // Method
    // Loop through all slave fields and check if the master has the same
    // fields, if so then we merge the fields, if not then we create a new field
    // in the master mesh with the slave values.

    for
    (
        IOobjectList::iterator fieldIter = slaveFields.begin();
        fieldIter != slaveFields.end();
        ++fieldIter
    )
    {
        IOobject slaveFieldIOobject
        (
            fieldIter()->name(),
            slaveMesh.time().timeName(),
            slaveMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        IOobject masterFieldIOobject
        (
            fieldIter()->name(),
            masterMesh.time().timeName(),
            // (
            //    masterMesh.time().value() - masterMesh.time().deltaT().value()
            // ),
            masterMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (slaveFieldIOobject.headerOk())
        {
            // Read slave field
            GeometricField<Type, fvsPatchField, surfaceMesh> slaveField
            (
                *fieldIter(),
                slaveMesh
            );

            // Read master field if it exists otherwise create it

            GeometricField<Type, fvsPatchField, surfaceMesh>* masterFieldPtr
                = NULL;

            if (masterFieldIOobject.headerOk())
            {
                Info<< "    Mapping " << fieldIter()->name() << endl;

                masterFieldPtr =
                    new GeometricField<Type, fvsPatchField, surfaceMesh>
                    (
                        masterFieldIOobject,
                        masterMesh
                    );
            }
            else
            {
                Info<< "    Creating " << fieldIter()->name() << endl;

                masterFieldPtr =
                    new GeometricField<Type, fvsPatchField, surfaceMesh>
                    (
                        IOobject
                        (
                            fieldIter()->name(),
                            masterMesh.time().timeName(),
                            masterMesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE
                        ),
                        masterMesh,
                        dimensioned<Type>
                        (
                            "zero", slaveField.dimensions(), pTraits<Type>::zero
                        )
                    );
            }

            GeometricField<Type, fvsPatchField, surfaceMesh>& masterField =
                *masterFieldPtr;

            // Merge fields

            // The slave cell values are at the end of the list
            const label nOldIntFaces =
                masterMesh.nInternalFaces() - slaveMesh.nInternalFaces();
            forAll(slaveField.internalField(), faceI)
            {
                masterField.internalField()[faceI + nOldIntFaces] =
                    slaveField.internalField()[faceI];
            }

            // Map patches with the same name and same size
            forAll(slaveField.boundaryField(), patchI)
            {
                if (!slaveField.boundaryField()[patchI].coupled())
                {
                    // Check of the master has a patch with the same name
                    label masterPatchID = masterMesh.boundaryMesh().findPatchID
                    (
                        slaveMesh.boundaryMesh()[patchI].name()
                    );

                    if (masterPatchID != -1)
                    {
                        // Check if patch sizes are the same
                        if
                        (
                            slaveMesh.boundaryMesh()[patchI].size()
                            == masterMesh.boundaryMesh()[masterPatchID].size()
                        )
                        {
                            forAll(slaveField.boundaryField()[patchI], faceI)
                            {
                                masterField.boundaryField()
                                [
                                    masterPatchID
                                ][faceI] =
                                    slaveField.boundaryField()[patchI][faceI];
                            }
                        }
                        // else
                        // {
                        //     Warning
                        //         << "Master and slave both have patch "
                        //         << slaveMesh.boundaryMesh()[patchI].name()
                        //         << " but they have a different size: "
                        //         << " boundary patch not merged!" << endl;
                        // }
                    }
                }
            }

            // Write master field
            masterField.write();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
