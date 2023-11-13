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

#include "LemaitreDamage.H"
#include "addToRunTimeSelectionTable.H"
// #include "zeroGradientFvPatchFields.H"
// #include "transformGeometricField.H"
#include "logVolFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LemaitreDamage, 0);
    addToRunTimeSelectionTable
    (
        damageLaw, LemaitreDamage, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::LemaitreDamage::LemaitreDamage
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const volScalarField& curMaterial
)
:
    damageLaw(name, mesh, dict, curMaterial),
    h_(dict.lookup("h")),
    S0_(dict.lookup("S0")),
    b_(dict.lookup("b")),
    epsilonD_(dict.lookup("epsilonD"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::LemaitreDamage::~LemaitreDamage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LemaitreDamage::updateDamage()
{
    // Lookup Kirchhoff stress field
    const volSymmTensorField& sigma =
        mesh().lookupObject<volSymmTensorField>("sigmaCauchy");

    // Calculate equivalent Cauchy stress
        const volScalarField sigmaEq =
            sqrt((3.0/2.0)*magSqr(dev(sigma)))
          + dimensionedScalar("SMALL", dimPressure, SMALL);

        // Calculate hydrostatic Cauchy stress
        const volScalarField sigmaHyd = (1.0/3.0)*tr(sigma);

        // Calculate constraint/triality
        const volScalarField triaxiality = sigmaHyd/sigmaEq;

    // Lookup plastic strain increment field
    const volScalarField& DEpsilonPEq =
        mesh().lookupObject<volScalarField>("DEpsilonPEq");

     // Lookup plastic strain field
    const volScalarField& epsilonPEq =
        mesh().lookupObject<volScalarField>("epsilonPEq");
  
     const volScalarField lambda =
            mesh().lookupObject<volScalarField>("lambda");

     // Lame's second parameter, aka the shear modulus
        const volScalarField& mu = mesh().lookupObject<volScalarField>("mu");

        // Young's modulus
        const volScalarField E = mu*(3.0*lambda + 2.0*mu)/(lambda + mu);

        // Poisson's ration
        const volScalarField nu = lambda/(2*(lambda + mu));

         // Calculate Lemaitre damage

       // Calculate the Y value
    volScalarField Y
    (
        IOobject
        (
            "Y",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

  volScalarField Ystar
    (
        IOobject
        (
            "Ystar",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );
        scalarField& YI = Y.internalField();
         scalarField& YstarI = Ystar.internalField();
        const scalarField& sigmaEqI = sigmaEq.internalField();
        const scalarField& epsilonPEqI = epsilonPEq.internalField();
     
        const scalarField& triaxialityI = triaxiality.internalField();
        const scalarField& EI = E.internalField();
        const scalarField& nuI = nu.internalField();
        const scalarField& damageLemaitreI = damage().oldTime().internalField();
        const scalar epsilonDValue = epsilonD_.value();

        forAll(YI, cellI)
        {
            // Calculate damage
            if
            (
                triaxialityI[cellI] > (-1.0/3.0)
             && epsilonPEqI[cellI] > epsilonDValue
            )
            {
                // Calculate Lemaitre Y parameter
                const scalar denom =
                    (2.0/3.0)*(1.0 + nuI[cellI])
                  + 3.0*(1.0 - 2.0*nuI[cellI])*pow(triaxialityI[cellI], 2.0);

                if (triaxialityI[cellI] > 0.0)
                {
                    // Damage in tension
                    YI[cellI]=
                       -(
                            pow(sigmaEqI[cellI], 2.0)
                           /(
                               2.0*EI[cellI]
                              *pow(1.0 - damageLemaitreI[cellI], 2.0)
                           )
                        )*denom;
                }
                else
                {
                    // Damage in compression is less than in tension; this
                    // depends on the h_ parameter
                    YI[cellI] =
                        -(
                            h_.value()*pow(sigmaEqI[cellI], 2.0)
                           /(
                                2.0*EI[cellI]
                               *pow
                                (
                                    1.0 - h_.value()*damageLemaitreI[cellI],
                                    2.0
                                )
                            )
                        )*denom;
                }

                YstarI[cellI]=pow(-YI[cellI]/S0_.value(), b_.value());
            }
        }

     Y.correctBoundaryConditions();
     Ystar.correctBoundaryConditions();
     const volScalarField newDamage =
        damage().oldTime()
      + Ystar*DEpsilonPEq;

    // Be careful with multi-material as the damage field is shared by all the
    // damage laws
    damage() = curMaterial()*newDamage + (1.0 - curMaterial())*damage();

}


// ************************************************************************* //
