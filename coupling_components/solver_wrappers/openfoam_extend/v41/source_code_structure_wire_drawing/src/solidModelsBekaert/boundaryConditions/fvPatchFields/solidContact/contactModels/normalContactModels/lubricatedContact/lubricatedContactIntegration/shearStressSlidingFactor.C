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

\*---------------------------------------------------------------------------*/

#include "shearStressSlidingFactor.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace lubricatedContactIntegration
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


double shearStressSlidingFactorFunction
(
    double x,
    void* params
)
{
    // Read function parameters
    functionParameters &fParams =
        *reinterpret_cast<functionParameters *>(params);

    const double& h = fParams.h;
    const double& sigma = fParams.sigma;

    return exp(-0.5*sqr(x/sigma))/(h - x);
}

scalar shearStressSlidingFactor
(
    const scalar& surfaceSeparation,
    const scalar& surfaceRoughness,
    gsl_integration_glfixed_table* gaussLegendreTablePtr
)
{
    // Shear stress sliding factor (phi_f) is being calculated using
    // Gauss-Legendre numerical integration from the GSL library

    // Set function parameters
    functionParameters fParams;
    fParams.h = surfaceSeparation;
    fParams.sigma = surfaceRoughness;

    // Set integration functions and parameters
    gsl_function fShearStressSlidingFactor;
    fShearStressSlidingFactor.function =
        &Foam::lubricatedContactIntegration::shearStressSlidingFactorFunction;
    fShearStressSlidingFactor.params =
        reinterpret_cast<void *>(&fParams);

    // Calculate shear stress sliding factor (phi_f)
    scalar shearStressSlidingFactor =
        surfaceSeparation/(surfaceRoughness*Foam::sqrt(2*pi))
       *gsl_integration_glfixed
        (
            &fShearStressSlidingFactor,
            -100*surfaceRoughness,
            surfaceSeparation - (surfaceRoughness/100),
            gaussLegendreTablePtr
        );

    return shearStressSlidingFactor;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace asperityContact

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
