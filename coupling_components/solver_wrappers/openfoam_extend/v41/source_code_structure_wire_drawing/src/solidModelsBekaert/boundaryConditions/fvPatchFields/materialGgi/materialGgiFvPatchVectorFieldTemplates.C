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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "materialGgiFvPatchVectorField.H"
#include "symmTransformField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > materialGgiFvPatchVectorField::patchNeighbourField
(
    const Field<Type>& iField
) const
{
    const ggiFvPatch& ggiPatch =
        refCast<const ggiFvPatch>(this->patch());

    // Get shadow face-cells and assemble shadow field
    // This is a patchInternalField of neighbour but access is inconvenient.
    // Assemble by hand.
    // HJ, 27/Sep/2011
    const unallocLabelList& sfc = ggiPatch.shadow().faceCells();

    Field<Type> sField(sfc.size());

    forAll (sField, i)
    {
        sField[i] = iField[sfc[i]];
    }

    tmp<Field<Type> > tpnf(ggiPatch.interpolate(sField));

    return tpnf;
}


template<class Type>
tmp<Field<Type> > materialGgiFvPatchVectorField::patchNeighbourField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
) const
{
    const ggiFvPatch& ggiPatch =
        refCast<const ggiFvPatch>(this->patch());

    // Get shadow field
    Field<Type> sField = sf.boundaryField()[ggiPatch.shadowIndex()];

    tmp<Field<Type> > tpnf(ggiPatch.interpolate(sField));

    return tpnf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
