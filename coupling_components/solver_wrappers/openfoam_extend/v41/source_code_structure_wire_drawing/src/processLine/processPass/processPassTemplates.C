/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "processPass.H"
#include "ZoneID.H"
#include "GeometricField.H"
#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class Type>
void Foam::processPass::MergeVolFields
(
    const IOobjectList& baseObjects,
    const IOobjectList& passItemObjects,
    const fvMesh& baseMesh,
    const fvMesh& passItemMesh
)
{
    const word fieldClassName
    (
        GeometricField<Type, fvPatchField, volMesh>::typeName
    );

    // Create a map from the passItem mesh into the baseMesh
    // We assume here that the baseMesh contains a cellZone with the same name
    // as the name of the passItemMesh: this cellZone will be the map from the
    // passItemMesh to the baseMesh

    // Look for a cellZone in the baseMesh with the same name as the
    // passItemMesh
    ZoneID<cellZone> cellZoneID(passItemMesh.name(), baseMesh.cellZones());

    if (!cellZoneID.active())
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "    void Foam::processPass::MergeVolFields\n"
            "(\n"
            "    const IOobjectList& baseObjects,\n"
            "    const IOobjectList& passItemObjects,\n"
            "    const fvMesh& baseMesh,\n"
            "    const fvMesh& passItemMesh\n"
            ")"
        )   << "There is no cellZone in the baseMesh with the same name as the "
            << "passItemMesh (" << passItemMesh.name() << ")"
            << abort(FatalError);
    }

    // Cell map from the passItemMesh to the baseMesh
    const labelList& cellMap = baseMesh.cellZones()[cellZoneID.index()];


    // Method
    // Loop through all passItem fields and check if the base has the same
    // fields, if so then we merge the fields, if not then we create a new field
    // in the base mesh with the passItem values.

    IOobjectList baseFields = baseObjects.lookupClass(fieldClassName);
    IOobjectList passItemFields = passItemObjects.lookupClass(fieldClassName);

    for
    (
        IOobjectList::iterator fieldIter = passItemFields.begin();
        fieldIter != passItemFields.end();
        ++fieldIter
    )
    {
        // IOobject passItemFieldIOobject
        // (
        //     fieldIter()->name(),
        //     passItemMesh.time().timeName(),
        //     passItemMesh,
        //     IOobject::MUST_READ,
        //     IOobject::AUTO_WRITE
        // );

        // Read passItem field
        GeometricField<Type, fvPatchField, volMesh> passItemField
        (
            *fieldIter(),
            passItemMesh
        );

        // Check if this field exists in the base mesh
        IOobject baseFieldIOobject
        (
            fieldIter()->name(),
            baseMesh.time().timeName(),
            baseMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        GeometricField<Type, fvPatchField, volMesh>* baseFieldPtr = NULL;

        if (baseFieldIOobject.headerOk())
        {
            if (debug)
            {
                Info<< "    Mapping " << fieldIter()->name() << endl;
            }

            baseFieldPtr =
                new GeometricField<Type, fvPatchField, volMesh>
                (
                    baseFieldIOobject,
                    baseMesh
                );
        }
        else
        {
            Info<< "        Creating " << fieldIter()->name() << endl;

            baseFieldPtr =
                new GeometricField<Type, fvPatchField, volMesh>
                (
                    IOobject
                    (
                        fieldIter()->name(),
                        baseMesh.time().timeName(),
                        baseMesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    baseMesh,
                    dimensioned<Type>
                    (
                        "zero",
                        passItemField.dimensions(),
                        pTraits<Type>::zero
                    )
                );
        }

        GeometricField<Type, fvPatchField, volMesh>& baseField =
            *baseFieldPtr;


        // Merge fields

        // Merge internalField

        const Field<Type>& passItemFieldI = passItemField.internalField();
        Field<Type>& baseFieldI = baseField.internalField();

        forAll(passItemFieldI, cellI)
        {
            baseFieldI[cellMap[cellI]] = passItemFieldI[cellI];
        }

        // Merge the boundaryField

        // We will only map patches with the same name and same size

        forAll(passItemField.boundaryField(), patchI)
        {
            if (!passItemField.boundaryField()[patchI].coupled())
            {
                // Check if the baseMesh has a patch with the same name
                const label basePatchID =
                    baseMesh.boundaryMesh().findPatchID
                    (
                        passItemMesh.boundaryMesh()[patchI].name()
                    );

                if (basePatchID != -1)
                {
                    // Check if patch sizes are the same
                    if
                    (
                        passItemMesh.boundaryMesh()[patchI].size()
                     == baseMesh.boundaryMesh()[basePatchID].size()
                    )
                    {
                        if(debug)
                        {
                            Info<< "        Mapping patch "
                                << passItemMesh.boundaryMesh()[patchI].name()
                                << endl;
                        }

                        baseField.boundaryField()[basePatchID] =
                            Field<Type>(passItemField.boundaryField()[patchI]);
                    }
                    else
                    {
                        FatalErrorIn
                        (
                            "template<class Type>\n"
                            "    void Foam::processPass::MergeVolFields\n"
                            "(\n"
                            "    const IOobjectList& baseObjects,\n"
                            "    const IOobjectList& passItemObjects,\n"
                            "    const fvMesh& baseMesh,\n"
                            "    const fvMesh& passItemMesh\n"
                            ")"
                        )   << "The patch "
                            << baseMesh.boundaryMesh()[patchI].name()
                            << " in the baseMesh has a different size to "
                            << "the equivalent patch in the passItemMesh!"
                            << abort(FatalError);
                    }
                }
                else
                {
                    FatalErrorIn
                    (
                        "template<class Type>\n"
                        "    void Foam::processPass::MergeVolFields\n"
                        "(\n"
                        "    const IOobjectList& baseObjects,\n"
                        "    const IOobjectList& passItemObjects,\n"
                        "    const fvMesh& baseMesh,\n"
                        "    const fvMesh& passItemMesh\n"
                        ")"
                    )   << "Cannot find patch "
                        << baseMesh.boundaryMesh()[patchI].name()
                        << " in the baseMesh!"
                        << abort(FatalError);
                }
            }
        }

        // No need to update boundary conditions: this will cause an error with
        // DU for example
        //baseField.correctBoundaryConditions();

        // Write base field
        baseField.write();
    }
}


// ************************************************************************* //
