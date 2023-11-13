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

InNamespace
    Foam::flowFactors

Description
    Functions used for calculating flow factors by Patir and Cheng

    Patir, Cheng - An Average Flow Model for Determining Effects of
                   Three-Dimensional Roughness on Partial Hydrodynamic
                   Lubrication, 1978

    Patir, Cheng - Application of Average Flow Model to Lubrication Between
                   Rough Sliding Surfaces, 1979

SourceFiles
    flowFactorsPatirCheng.C

\*---------------------------------------------------------------------------*/

#ifndef flowFactorsPatirCheng_H
#define flowFactorsPatirCheng_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Namespace flowFactors functions Declaration
\*---------------------------------------------------------------------------*/

namespace flowFactors
{
    scalarField pressureFlowFactorsPatirCheng
    (
        const scalar& corrLengthRatio,
        const scalarField& normSurfaceSeparation
    );

    scalarField shearFlowFactorsPatirCheng
    (
        const scalar& corrLengthRatio,
        const scalarField& normSurfaceSeparation
    );

    scalarField shearStressFactorsPatirCheng
    (
        const scalar& corrLengthRatio,
        const scalarField& normSurfaceSeparation
    );

//    scalar shearStressSlidingFactorPatirCheng
//    (
//    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "flowFactorsPatirCheng.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
