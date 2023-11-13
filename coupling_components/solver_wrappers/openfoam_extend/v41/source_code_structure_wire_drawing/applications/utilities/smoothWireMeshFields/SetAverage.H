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
void SetAverage
(
    const fvMesh& mesh,
    const labelList& cellIDs,
    const labelList& cellIDsInRange
)
{
    // Read volField objects from object registry
    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields
    (
        mesh.thisDb().objectRegistry::template lookupClass
       <GeometricField<Type, fvPatchField, volMesh> >()
    );

    for
    (
        typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        // Read field
        GeometricField<Type, fvPatchField, volMesh>& field =
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (*fieldIter());

        //Info<< "    Setting " << field.name() << endl;

        // Calculate average value along streamline

        Type av = pTraits<Type>::zero;

        forAll(cellIDsInRange, cI)
        {
            av += field[cellIDsInRange[cI]];
        }
        av /= cellIDsInRange.size();

        // Set all values to be average value
        forAll(cellIDs, cI)
        {
            field[cellIDs[cI]] = av;
        }

        // Extrapolated boundary values
        forAll(field.boundaryField(), patchI)
        {
            if (field.boundaryField()[patchI].coupled())
            {
                field.boundaryField()[patchI].evaluate();
            }
            else if (field.boundaryField()[patchI].type() != "empty")
            // (
            //     !field.boundaryField()[patchI].coupled()
            //  && field.boundaryField()[patchI].type() != "empty"
            // )
            {
                field.boundaryField()[patchI] =
                    field.boundaryField()[patchI].patchInternalField();
            }
        }

        //field.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
