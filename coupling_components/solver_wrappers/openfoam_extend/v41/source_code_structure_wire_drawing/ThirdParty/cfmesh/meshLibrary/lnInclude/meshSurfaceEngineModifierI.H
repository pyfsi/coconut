/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "meshSurfaceEngineModifier.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline labelMapType& meshSurfaceEngineModifier::labelMap()
{
    if( !exchangeLabelsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "labelMapType& meshSurfaceEngineModifier::labelMap()"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        labelMapType* mapPtr = new labelMapType();
        labelMapType& m = *mapPtr;

        const DynList<label>& neiProcs =
            surfaceEngine_.bpNeiProcs();
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();

        forAll(neiProcs, i)
        {
            labelListType& l = m[neiProcs[i]];
            l.setSize(globalToLocal.size());
            l.clear();
        }

        exchangeLabelsPtr_ = mapPtr;
    }
    else
    {
        //- clear remaining data
        forAllIter(labelMapType, *exchangeLabelsPtr_, it)
            it->second.clear();
    }

    return *exchangeLabelsPtr_;
}

inline labelledPointMapType& meshSurfaceEngineModifier::labelledPointMap()
{
    if( !exchangeLabelledPointsPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "labelledPointMapType& "
                "meshSurfaceEngineModifier::labelledPointMap()"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        labelledPointMapType* mapPtr = new labelledPointMapType();
        labelledPointMapType& m = *mapPtr;

        const DynList<label>& neiProcs =
            surfaceEngine_.bpNeiProcs();
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();

        forAll(neiProcs, i)
        {
            labelledPointListType& l = m[neiProcs[i]];
            l.setSize(globalToLocal.size());
            l.clear();
        }

        exchangeLabelledPointsPtr_ = mapPtr;
    }
    else
    {
        forAllIter(labelledPointMapType, *exchangeLabelledPointsPtr_, m)
            m->second.clear();
    }

    return *exchangeLabelledPointsPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
