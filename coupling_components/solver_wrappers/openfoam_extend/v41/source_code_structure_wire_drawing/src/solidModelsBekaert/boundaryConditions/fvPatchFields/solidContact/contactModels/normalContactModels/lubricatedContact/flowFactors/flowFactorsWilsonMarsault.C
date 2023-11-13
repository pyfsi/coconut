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

#include "flowFactorsWilsonMarsault.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace flowFactors
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalarField pressureFlowFactorWilsonMarsault
(
    const scalar& corrLengthRatio,
    const scalarField& normFilmMeanThickness
)
{
    // Non-dimensional pressure flow is a pressure flow factor multiplied by
    // the third potention of normalized film flow thickness, i.e. = phiX*Ht^3

    scalar logGamma = log(corrLengthRatio);
    scalar log9Gamma = log(9*corrLengthRatio);

    scalar Htc = 3.0*(1.0 - pow(0.47476/corrLengthRatio + 1.0, -0.25007));
    scalar a2 = 0.051375*pow3(log9Gamma) - 0.0071901*pow4(log9Gamma);
    scalar a3 =
        1.0019 - 0.17927*logGamma + 0.047583*sqr(logGamma)
      - 0.016417*pow3(logGamma);

    scalarField pressureFlow(normFilmMeanThickness.size(), 0.0);

    forAll(pressureFlow, faceI)
    {
        if (normFilmMeanThickness[faceI] >= Htc)
        {
            pressureFlow[faceI] =
                a2*sqr(normFilmMeanThickness[faceI] - Htc)
              + a3*pow3(normFilmMeanThickness[faceI] - Htc);
        }

    }

    return pressureFlow/pow3(normFilmMeanThickness);
}


scalarField shearFlowFactorsWilsonMarsault
(
    const scalar& corrLengthRatio,
    const scalarField& normFilmMeanThickness
)
{
    // Calculated using spline surface interpolation and polynomial fitting of
    // data in Fig. 4 (Wilson and Marsault). The calculation was done by Matlab
    // R2016b using Curve Fitting Toolbox.

    scalar p00 =    1.036000000;
    scalar p10 =    0.378100000;    scalar p01 =   -4.050000000;
    scalar p20 =   -0.129800000;    scalar p02 =    9.711000000;
    scalar p30 =   -0.014770000;    scalar p03 =  -13.490000000;
    scalar p40 =    0.007365000;    scalar p04 =   10.210000000;
    scalar p50 =   -0.000566800;    scalar p05 =   -3.194000000;
    scalar p11 =    0.681700000;
    scalar p21 =    0.003060000;    scalar p12 =   -1.204000000;
    scalar p31 =   -0.012890000;    scalar p13 =    0.663200000;
    scalar p41 =    0.000511600;    scalar p14 =   -0.109800000;
    scalar p22 =    0.090080000;
    scalar p32 =    0.003367000;    scalar p23 =   -0.046410000;

    if (corrLengthRatio >= 1.0)
    {
        p00 =   0.140300000;
        p10 =   0.744800000;    p01 =   0.109800000;
        p20 =  -0.274400000;    p02 =  -0.104200000;
        p30 =   0.031410000;    p03 =   0.028280000;
        p40 =  -0.000593200;    p04 =  -0.003105000;
        p50 =  -0.000070390;    p05 =   0.000120300;
        p11 =  -0.203000000;
        p21 =   0.059160000;    p12 =   0.026890000;
        p31 =  -0.004833000;    p13 =  -0.001413000;
        p41 =   0.000100600;    p14 =   0.000021830;
        p22 =  -0.005968000;
        p32 =   0.000265600;    p23 =   0.000191600;
    }

    // Calculating shear flow factors

    const scalarField& Ht = normFilmMeanThickness;
    const scalar& gamma = corrLengthRatio;

    // For Ht > 5, equation (22) from Patir & Cheng (1979) is used for
    // calculating shear flow factor: A2*exp(-0.25*Ht)
    // Since function continuity is required, we need to calculate A2
    scalar A2 =
        (
            p00
          + 5.0*(p10 + 5.0*(p20 + 5.0*(p30 + 5.0*(p40 + 5.0*p50))))
          + 5.0*gamma
           *(
                p11
              + 5.0*(p21 + 5.0*(p31 + 5.0*p41))
              + gamma*(p12 + gamma*(p13 + gamma*p14))
              + 5.0*gamma*(p22 + 5.0*p32 + gamma*p23)
            )
          + gamma*(p01 + gamma*(p02 + gamma*(p03 + gamma*(p04 + gamma*p05))))
        )*exp(0.25*5.0);

    // Calculating polynomial
    scalarField shearFlowFactors(normFilmMeanThickness.size(), 0.0);

    forAll(shearFlowFactors, faceI)
    {
        const scalar& fHt = Ht[faceI];

        if (fHt > 5.0)
        {
            shearFlowFactors[faceI] = A2*exp(-0.25*fHt);
        }
        else
        {
            shearFlowFactors[faceI] =
                p00
              + fHt*(p10 + fHt*(p20 + fHt*(p30 + fHt*(p40 + fHt*p50))))
              + fHt*gamma
               *(
                    p11
                  + fHt*(p21 + fHt*(p31 + fHt*p41))
                  + gamma*(p12 + gamma*(p13 + gamma*p14))
                  + fHt*gamma*(p22 + fHt*p32 + gamma*p23)
                )
              + gamma
               *(p01 + gamma*(p02 + gamma*(p03 + gamma*(p04 + gamma*p05))));
        }

        if (shearFlowFactors[faceI] > fHt)
        {
            shearFlowFactors[faceI] = fHt;
        }
    }

    return shearFlowFactors;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flowFactors

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
