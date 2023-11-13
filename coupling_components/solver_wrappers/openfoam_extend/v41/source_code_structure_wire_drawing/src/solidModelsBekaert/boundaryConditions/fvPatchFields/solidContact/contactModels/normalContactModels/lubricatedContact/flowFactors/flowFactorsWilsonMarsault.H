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
    Functions used for calculating flow factors for surfaces with Christensen
    height distribution and large fractional contact areas by Wilson and
    Marsault

    Wilson, Marsault - Partial Hydrodynamic Lubrication With Large Fractional
                       Contact Areas, 1998

SourceFiles
    flowFactorsWilsonMarsault.C

\*---------------------------------------------------------------------------*/

#ifndef flowFactorsWilsonMarsault_H
#define flowFactorsWilsonMarsault_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Namespace flowFactors functions Declaration
\*---------------------------------------------------------------------------*/

namespace flowFactors
{
    scalarField pressureFlowFactorWilsonMarsault
    (
        const scalar& corrLengthRatio,
        const scalarField& normFilmMeanThickness
    );

    scalarField shearFlowFactorsWilsonMarsault
    (
        const scalar& corrLengthRatio,
        const scalarField& normFilmMeanThickness
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "flowFactorsWilsonMarsault.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
