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

#ifndef MergeVolFields_H
#define MergeVolFields_H

#include "GeometricField.H"
//#include "meshToMesh.H"
#include "IOobjectList.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void MergeVolFields
(
    const IOobjectList& masterObjects,
    const IOobjectList& slaveObjects,
    const fvMesh& masterMesh,
    const fvMesh& slaveMesh
)
{
    word fieldClassName
    (
        GeometricField<Type, fvPatchField, volMesh>::typeName
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
        if
        (
            fieldIter()->name() != "U_0"
         && fieldIter()->name() != "U_0_0"
         && fieldIter()->name() != "U_0_0_0"
         && fieldIter()->name() != "rho_0"
         && fieldIter()->name() != "rho_0_0"
         && fieldIter()->name() != "DU_0"
         && fieldIter()->name() != "DU_0_0"
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
                masterMesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            );

            if (slaveFieldIOobject.headerOk())
            {
                // Read slave field
                GeometricField<Type, fvPatchField, volMesh> slaveField
                    (
                        *fieldIter(),
                        slaveMesh
                    );

                // Read master field if it exists otherwise create it

                GeometricField<Type, fvPatchField, volMesh>* masterFieldPtr =
                    NULL;

                if (masterFieldIOobject.headerOk())
                {
                    Info<< "    Mapping " << fieldIter()->name() << endl;

                    masterFieldPtr =
                        new GeometricField<Type, fvPatchField, volMesh>
                        (
                            masterFieldIOobject,
                            masterMesh
                        );
                }
                else
                {
                    Info<< "    Creating " << fieldIter()->name() << endl;

                    masterFieldPtr =
                        new GeometricField<Type, fvPatchField, volMesh>
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
                                "zero",
                                slaveField.dimensions(),
                                pTraits<Type>::zero
                            )
                        );
                }

                GeometricField<Type, fvPatchField, volMesh>& masterField =
                    *masterFieldPtr;

                // Merge fields
                // The slave cell values are at the end of the list
                const label nOldCells =
                    masterMesh.nCells() - slaveMesh.nCells();
                forAll(slaveField.internalField(), cellI)
                {
                    masterField.internalField()[cellI + nOldCells] =
                        slaveField.internalField()[cellI];
                }

                // Map patches with the same name and same size
                forAll(slaveField.boundaryField(), patchI)
                {
                    if (!slaveField.boundaryField()[patchI].coupled())
                    {
                        // Check of the master has a patch with the same name
                        label masterPatchID =
                            masterMesh.boundaryMesh().findPatchID
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
                                Info<< "        Mapping patch "
                                    << slaveMesh.boundaryMesh()[patchI].name()
                                    << endl;

                                forAll
                                (
                                    slaveField.boundaryField()[patchI],
                                    faceI
                                )
                                {
                                    masterField.boundaryField()
                                        [
                                            masterPatchID
                                        ][faceI] =
                                        slaveField.boundaryField()[patchI]
                                        [
                                            faceI
                                        ];
                                }
                            }
                        }
                    }
                }

                // Do I need to update the BCs?
                //masterField.correctBoundaryConditions();

                // Write master field
                masterField.write();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
