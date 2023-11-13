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

Description
    Interpolation class dealing with transfer of data between two
    primitivePatches

\*---------------------------------------------------------------------------*/

#include "PatchToPatchInterpolationTemplateNew.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class FromPatch, class ToPatch>
const Foam::debug::tolerancesSwitch
PatchToPatchInterpolationNew<FromPatch, ToPatch>::directHitTol_
(
    "patchToPatchDirectHit",
    1e-5
);

template<class FromPatch, class ToPatch>
const Foam::debug::tolerancesSwitch
PatchToPatchInterpolationNew<FromPatch, ToPatch>::projectionTol_
(
    "patchToPatchProjectionTol",
    0.05
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FromPatch, class ToPatch>
void PatchToPatchInterpolationNew<FromPatch, ToPatch>::clearOut()
{
    deleteDemandDrivenData(pointAddressingPtr_);
    deleteDemandDrivenData(pointWeightsPtr_);
    deleteDemandDrivenData(pointDistancePtr_);
    deleteDemandDrivenData(faceAddressingPtr_);
    deleteDemandDrivenData(faceWeightsPtr_);
    deleteDemandDrivenData(faceDistancePtr_);
}


template<class FromPatch, class ToPatch>
void PatchToPatchInterpolationNew<FromPatch, ToPatch>::setWeights
(
    labelList* paPtr,
    FieldField<Field, scalar>* pwPtr,
    scalarField* pdPtr,
    labelList* faPtr,
    FieldField<Field, scalar>* fwPtr,
    scalarField* fdPtr
)
{
    clearOut();

    pointAddressingPtr_ = paPtr;
    pointWeightsPtr_ = pwPtr;
    pointDistancePtr_ = pdPtr;
    faceAddressingPtr_ = faPtr;;
    faceWeightsPtr_ = fwPtr;
    faceDistancePtr_ = fdPtr;;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class FromPatch, class ToPatch>
PatchToPatchInterpolationNew<FromPatch, ToPatch>::PatchToPatchInterpolationNew
(
    const FromPatch& fromPatch,
    const ToPatch& toPatch,
    intersection::algorithm alg,
    const intersection::direction dir
)
:
    fromPatch_(fromPatch),
    toPatch_(toPatch),
    alg_(alg),
    dir_(dir),
    pointAddressingPtr_(NULL),
    pointWeightsPtr_(NULL),
    pointDistancePtr_(NULL),
    faceAddressingPtr_(NULL),
    faceWeightsPtr_(NULL),
    faceDistancePtr_(NULL)
{
}


// Construct as copy
template<class FromPatch, class ToPatch>
PatchToPatchInterpolationNew<FromPatch, ToPatch>::PatchToPatchInterpolationNew
(
    const PatchToPatchInterpolationNew<FromPatch, ToPatch>& ppi
)
:
    PatchToPatchInterpolationNewName(),
    fromPatch_(ppi.fromPatch_),
    toPatch_(ppi.toPatch_),
    alg_(ppi.alg_),
    dir_(ppi.dir_),
    pointAddressingPtr_(NULL),
    pointWeightsPtr_(NULL),
    pointDistancePtr_(NULL),
    faceAddressingPtr_(NULL),
    faceWeightsPtr_(NULL),
    faceDistancePtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class FromPatch, class ToPatch>
PatchToPatchInterpolationNew<FromPatch, ToPatch>::~PatchToPatchInterpolationNew()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FromPatch, class ToPatch>
const labelList&
PatchToPatchInterpolationNew<FromPatch, ToPatch>::pointAddr() const
{
    if (!pointAddressingPtr_)
    {
        calcPointAddressing();
    }

    return *pointAddressingPtr_;
}


template<class FromPatch, class ToPatch>
const FieldField<Field, scalar>&
PatchToPatchInterpolationNew<FromPatch, ToPatch>::pointWeights() const
{
    if (!pointWeightsPtr_)
    {
        calcPointAddressing();
    }

    return *pointWeightsPtr_;
}


template<class FromPatch, class ToPatch>
const labelList&
PatchToPatchInterpolationNew<FromPatch, ToPatch>::faceAddr() const
{
    if (!faceAddressingPtr_)
    {
        calcFaceAddressing();
    }

    return *faceAddressingPtr_;
}


template<class FromPatch, class ToPatch>
const FieldField<Field, scalar>&
PatchToPatchInterpolationNew<FromPatch, ToPatch>::faceWeights() const
{
    if (!faceWeightsPtr_)
    {
        calcFaceAddressing();
    }

    return *faceWeightsPtr_;
}


template<class FromPatch, class ToPatch>
const scalarField&
PatchToPatchInterpolationNew<FromPatch, ToPatch>
::pointDistanceToIntersection() const
{
    if (!pointDistancePtr_)
    {
        calcPointAddressing();
    }

    return *pointDistancePtr_;
}


template<class FromPatch, class ToPatch>
const scalarField&
PatchToPatchInterpolationNew<FromPatch, ToPatch>
::faceDistanceToIntersection() const
{
    if (!faceDistancePtr_)
    {
        calcFaceAddressing();
    }

    return *faceDistancePtr_;
}


template<class FromPatch, class ToPatch>
bool PatchToPatchInterpolationNew<FromPatch, ToPatch>::movePoints()
{
    clearOut();

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "CalcPatchToPatchWeightsNew.C"
#   include "PatchToPatchInterpolateNew.C"

// ************************************************************************* //
