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

#include "flowFactorsPatirCheng.H"
#include "interpolateXY.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace flowFactors
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalarField pressureFlowFactorsPatirCheng
(
    const scalar& corrLengthRatio,
    const scalarField& normSurfaceSeparation
)
{
    // Calculate list of coefficients needed for calculating pressure flow
    // factors

    // Coefficients for calculating pressure flow factors are given for a
    // discrete number of correlation lengths (gamma). For that reason, we need
    // to interpolate pressure flow factors from specified correlation lengths
    // to our correlation length

    // Create array of coefficients
    // NOTE: Added extra coefficient for gamma = (1.0 + SMALL)
    const scalar gammaArray[8] =
        {1.0/9.0, 1.0/6.0, 1.0/3.0, 1.0, (1.0 + SMALL), 3.0, 6.0, 9.0};
    const scalar cArray[8] =
        {1.48, 1.38, 1.18, 0.9, SMALL, 0.225, 0.520, 0.870};
    const scalar rArray[8] = {0.42, 0.42, 0.42, 0.56, 1.5, 1.5, 1.5, 1.5};

    // Find exact, or ceal and floor values of gamma with its appropriate
    // coefficients
    scalarList coeffs(3, 0.0);
    DynamicList<scalarList> coeffsList;

    if (corrLengthRatio <= gammaArray[0])
    {
        coeffs[0] = gammaArray[0];
        coeffs[1] = cArray[0];
        coeffs[2] = rArray[0];

        coeffsList.append(coeffs);
    }
    else if (corrLengthRatio >= gammaArray[7])
    {
        coeffs[0] = gammaArray[7];
        coeffs[1] = cArray[7];
        coeffs[2] = rArray[7];

        coeffsList.append(coeffs);
    }
    else
    {
        for (int i = 1; i < 8; ++i)
        {
            if
            (
                corrLengthRatio > gammaArray[i-1]
             && corrLengthRatio < gammaArray[i]
            )
            {
                // First, floor coefficients
                coeffs[0] = gammaArray[i-1];
                coeffs[1] = cArray[i-1];
                coeffs[2] = rArray[i-1];
                coeffsList.append(coeffs);
                // Second, ceil coefficients
                coeffs[0] = gammaArray[i];
                coeffs[1] = cArray[i];
                coeffs[2] = rArray[i];
                coeffsList.append(coeffs);
            }
            else if (corrLengthRatio == gammaArray[i])
            {
                coeffs[0] = gammaArray[i];
                coeffs[1] = cArray[i];
                coeffs[2] = rArray[i];
                coeffsList.append(coeffs);
            }
        }
    }


    // Calculate and return pressure flow factors
    scalarField pressureFlowFactors(normSurfaceSeparation.size(), 0.0);

    forAll(pressureFlowFactors,faceI)
    {
        // For gamma >= 1, pressure flow factors are valid for H > 0.5
        scalar H = max(normSurfaceSeparation[faceI], 0.5);

        // For gamma < 1, pressure flow factors are valid for H > 1
        if (corrLengthRatio < 1.0)
        {
            H = max(H, 1.0);
        }

        scalarField gammaField(coeffsList.size(), 0.0);
        scalarField pressureFlowFactorField(coeffsList.size(), 0.0);

        forAll(gammaField, gammaI)
        {
            const scalar gamma = coeffsList[gammaI][0];
            const scalar& C = coeffsList[gammaI][1];
            const scalar& r = coeffsList[gammaI][2];

            gammaField[gammaI] = gamma;

            if (corrLengthRatio <= 1.0)
            {
                pressureFlowFactorField[gammaI] = 1.0 - C*exp(-r*H);
            }
            else
            {
                pressureFlowFactorField[gammaI] = 1.0 + C*pow(H, -r);
            }
        }

        pressureFlowFactors[faceI] =
            interpolateXY(corrLengthRatio, gammaField, pressureFlowFactorField);
    }

    return pressureFlowFactors;
}


scalarField shearFlowFactorsPatirCheng
(
    const scalar& corrLengthRatio,
    const scalarField& normSurfaceSeparation
)
{
    // Calculate list of coefficients needed for calculating shear flow factors

    // Coefficients for calculating shear flow factors are given for a discrete
    // number of correlation lengths (gamma). For that reason, we need to
    // interpolate pressure flow factors from specified correlation lengths to
    // our correlation length

    // Create array of coefficients
    const scalar gammaArray[7] =
        {1.0/9.0, 1.0/6.0, 1.0/3.0, 1.0, 3.0, 6.0, 9.0};
    const scalar a1Array[7] = {2.046, 1.962, 1.858, 1.899, 1.560, 1.290, 1.011};
    const scalar alpha1Array[7] = {1.12, 1.08, 1.01, 0.98, 0.85, 0.62, 0.54};
    const scalar alpha2Array[7] = {0.78, 0.77, 0.76, 0.92, 1.13, 1.09, 1.07};
    const scalar alpha3Array[7] = {0.03, 0.03, 0.03, 0.05, 0.08, 0.08, 0.08};
    const scalar a2Array[7] = {1.856, 1.754, 1.561, 1.126, 0.556, 0.388, 0.295};

    // Find exact, or ceal and floor values of gamma with its appropriate
    // coefficients
    scalarList coeffs(6, 0.0);
    DynamicList<scalarList> coeffsList;

    if (corrLengthRatio <= gammaArray[0])
    {
        coeffs[0] = gammaArray[0];
        coeffs[1] = a1Array[0];
        coeffs[2] = alpha1Array[0];
        coeffs[3] = alpha2Array[0];
        coeffs[4] = alpha3Array[0];
        coeffs[5] = a2Array[0];

        coeffsList.append(coeffs);
    }
    else if (corrLengthRatio >= gammaArray[6])
    {
        coeffs[0] = gammaArray[6];
        coeffs[1] = a1Array[6];
        coeffs[2] = alpha1Array[6];
        coeffs[3] = alpha2Array[6];
        coeffs[4] = alpha3Array[6];
        coeffs[5] = a2Array[6];

        coeffsList.append(coeffs);
    }
    else
    {
        for (int i = 1; i < 7; ++i)
        {
            if
            (
                corrLengthRatio > gammaArray[i-1]
             && corrLengthRatio < gammaArray[i]
            )
            {
                // First, floor coefficients
                coeffs[0] = gammaArray[i-1];
                coeffs[1] = a1Array[i-1];
                coeffs[2] = alpha1Array[i-1];
                coeffs[3] = alpha2Array[i-1];
                coeffs[4] = alpha3Array[i-1];
                coeffs[5] = a2Array[i-1];
                coeffsList.append(coeffs);
                // Second, ceil coefficients
                coeffs[0] = gammaArray[i];
                coeffs[1] = a1Array[i];
                coeffs[2] = alpha1Array[i];
                coeffs[3] = alpha2Array[i];
                coeffs[4] = alpha3Array[i];
                coeffs[5] = a2Array[i];
                coeffsList.append(coeffs);
            }
            else if (corrLengthRatio == gammaArray[i])
            {
                coeffs[0] = gammaArray[i];
                coeffs[1] = a1Array[i];
                coeffs[2] = alpha1Array[i];
                coeffs[3] = alpha2Array[i];
                coeffs[4] = alpha3Array[i];
                coeffs[5] = a2Array[i];
                coeffsList.append(coeffs);
            }
        }
    }


    // Calculate and return shear flow factor
    scalarField shearFlowFactors(normSurfaceSeparation.size(), 0.0);

    forAll(shearFlowFactors, faceI)
    {
        // Shear flow factors are valid for H > 0.5
        scalar H = max(normSurfaceSeparation[faceI], 0.5);

        scalarField gammaField(coeffsList.size(), 0.0);
        scalarField shearFlowFactorField(coeffsList.size(), 0.0);

        forAll(gammaField, gammaI)
        {
            const scalar gamma = coeffsList[gammaI][0];
            const scalar& A1 = coeffsList[gammaI][1];
            const scalar& alpha1 = coeffsList[gammaI][2];
            const scalar& alpha2 = coeffsList[gammaI][3];
            const scalar& alpha3 = coeffsList[gammaI][4];
            const scalar& A2 = coeffsList[gammaI][5];

            gammaField[gammaI] = gamma;

            if (H <= 5.0)
            {
                shearFlowFactorField[gammaI] =
                    A1*pow(H, alpha1)*exp(-alpha2*H + alpha3*sqr(H));
            }
            else
            {
                shearFlowFactorField[gammaI] = A2*exp(-0.25*H);
            }
        }

        shearFlowFactors[faceI] =
            interpolateXY(corrLengthRatio, gammaField, shearFlowFactorField);
    }

    return shearFlowFactors;
}


scalarField shearStressFactorsPatirCheng
(
    const scalar& corrLengthRatio,
    const scalarField& normSurfaceSeparation
)
{
    // Calculate list of coefficients needed for calculating shear stress
    // factors

    // Coefficients for calculating shear stress factors are given for a
    // discrete number of correlation lengths (gamma). For that reason, we need
    // to interpolate pressure flow factors from specified correlation lengths
    // to our correlation length

    // Create array of coefficients
    const scalar gammaArray[7] =
        {1.0/9.0, 1.0/6.0, 1.0/3.0, 1.0, 3.0, 6.0, 9.0};
    const scalar a3Array[7] = {14.1, 13.4, 12.3, 11.1, 9.8, 10.1, 8.7};
    const scalar alpha4Array[7] = {2.45, 2.42, 2.32, 2.31, 2.25, 2.25, 2.15};
    const scalar alpha5Array[7] = {2.3, 2.3, 2.3, 2.38, 2.80, 2.90, 2.97};
    const scalar alpha6Array[7] = {0.1, 0.1, 0.1, 0.11, 0.18, 0.18, 0.18};

    // Find exact, or ceal and floor values of gamma with its appropriate
    // coefficients
    scalarList coeffs(5, 0.0);
    DynamicList<scalarList> coeffsList;

    if (corrLengthRatio <= gammaArray[0])
    {
        coeffs[0] = gammaArray[0];
        coeffs[1] = a3Array[0];
        coeffs[2] = alpha4Array[0];
        coeffs[3] = alpha5Array[0];
        coeffs[4] = alpha6Array[0];

        coeffsList.append(coeffs);
    }
    else if (corrLengthRatio >= gammaArray[6])
    {
        coeffs[0] = gammaArray[6];
        coeffs[1] = a3Array[6];
        coeffs[2] = alpha4Array[6];
        coeffs[3] = alpha5Array[6];
        coeffs[4] = alpha6Array[6];

        coeffsList.append(coeffs);
    }
    else
    {
        for (int i = 1; i < 7; ++i)
        {
            if
            (
                corrLengthRatio > gammaArray[i-1]
             && corrLengthRatio < gammaArray[i]
            )
            {
                // First, floor coefficients
                coeffs[0] = gammaArray[i-1];
                coeffs[1] = a3Array[i-1];
                coeffs[2] = alpha4Array[i-1];
                coeffs[3] = alpha5Array[i-1];
                coeffs[4] = alpha6Array[i-1];
                coeffsList.append(coeffs);
                // Second, ceil coefficients
                coeffs[0] = gammaArray[i];
                coeffs[1] = a3Array[i];
                coeffs[2] = alpha4Array[i];
                coeffs[3] = alpha5Array[i];
                coeffs[4] = alpha6Array[i];
                coeffsList.append(coeffs);
            }
            else if (corrLengthRatio == gammaArray[i])
            {
                coeffs[0] = gammaArray[i];
                coeffs[1] = a3Array[i];
                coeffs[2] = alpha4Array[i];
                coeffs[3] = alpha5Array[i];
                coeffs[4] = alpha6Array[i];
                coeffsList.append(coeffs);
            }
        }
    }


    // Calculate and return shear stress factor
    scalarField shearStressFactors(normSurfaceSeparation.size(), 0.0);

    forAll(shearStressFactors, faceI)
    {
        // Shear stress factors are valid for H > 0.5
        scalar H = max(normSurfaceSeparation[faceI], 0.5);

        scalarField gammaField(coeffsList.size(), 0.0);
        scalarField shearStressFactorField(coeffsList.size(), 0.0);

        forAll(gammaField, gammaI)
        {
            const scalar gamma = coeffsList[gammaI][0];
            const scalar& A3 = coeffsList[gammaI][1];
            const scalar& alpha4 = coeffsList[gammaI][2];
            const scalar& alpha5 = coeffsList[gammaI][3];
            const scalar& alpha6 = coeffsList[gammaI][4];

            gammaField[gammaI] = gamma;

            if (H < 7.0)
            {
                shearStressFactorField[gammaI] =
                    A3*pow(H, alpha4)*exp(-alpha5*H + alpha6*sqr(H));
            }
            else
            {
                shearStressFactorField[gammaI] = 0.0;
            }
        }

        shearStressFactors[faceI] =
            interpolateXY(corrLengthRatio, gammaField, shearStressFactorField);
    }

    return shearStressFactors;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowFactors

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
