/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Typedef
    Foam::symmTensor4thOrder

Description
    DevSymmTensor4thOrder of scalars.

SourceFiles
    symmTensor4thOrder.C

\*---------------------------------------------------------------------------*/

#ifndef devSymmTensor4thOrder_H
#define devSymmTensor4thOrder_H

#include "DevSymmTensor4thOrderTemplate.H"
#include "vector.H"
#include "contiguous.H"


#include "objectRegistry.H"
#include "GeometricFields.H"
#include "volMesh.H"
#include "fvMesh.H"
#include "fvPatchFields.H"
#include "volFieldsFwd.H"
#include "calculatedFvPatchFields.H"
#include "fvMatrices.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef DevSymmTensor4thOrder<scalar> devSymmTensor4thOrder;

//- Data associated with symmTensor4thOrder type are contiguous
template<>
inline bool contiguous<devSymmTensor4thOrder>() {return true;}

    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// this is to include everything else that needs to be defined in order to create a field of these things

#endif

// ************************************************************************* //
