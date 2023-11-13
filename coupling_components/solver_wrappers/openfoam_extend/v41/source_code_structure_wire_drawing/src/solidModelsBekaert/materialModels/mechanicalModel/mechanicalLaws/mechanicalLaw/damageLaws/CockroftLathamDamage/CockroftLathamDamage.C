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

#include "CockroftLathamDamage.H"
#include "addToRunTimeSelectionTable.H"
// #include "zeroGradientFvPatchFields.H"
// #include "transformGeometricField.H"
#include "logVolFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CockroftLathamDamage, 0);
    addToRunTimeSelectionTable
    (
        damageLaw, CockroftLathamDamage, dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::CockroftLathamDamage::CockroftLathamDamage
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const volScalarField& curMaterial
)
:
    damageLaw(name, mesh, dict, curMaterial),
    Vd_(dict.lookup("criticalDamageVd"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CockroftLathamDamage::~CockroftLathamDamage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CockroftLathamDamage::updateDamage(const volScalarField DEpsilonPEq)
{
    // Lookup Kirchhoff stress field
  //  const volSymmTensorField& tau =
  //      mesh().lookupObject<volSymmTensorField>("tauKirchhoff");
const volSymmTensorField& tau =
        mesh().lookupObject<volSymmTensorField>("sigmaCauchy");


/*const volSymmTensorField& DEpsilonP =
        mesh().lookupObject<volSymmTensorField>("Dp");
    // Lookup Kirchhoff stress field
     volScalarField DEpsilonPEq = sqrt((2.0/3.0)*magSqr(dev(DEpsilonP)));*/
// Read plastic strain tensor
       
//const volScalarField& DEpsilonPEq =
 //       mesh().lookupObject<volScalarField>("DEpsilonPEq");

    // Calculate the maximum (most positive) Kirchhoff principal stresses
    volVectorField eigenVal
    (
        IOobject
        (
            "eigenVal(tau)",
            tau.time().timeName(),
            tau.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", tau.dimensions(), vector::zero)
    );

    volTensorField eigenVec
    (
        IOobject
        (
            "eigenVec(tau)",
            tau.time().timeName(),
            tau.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    // Calculate eigen values and eigen vectors of tau
    eig3Field(tau, eigenVec, eigenVal);

    // Calculate the maximum principal stress
    volScalarField maxPrincipalStress
    (
        IOobject
        (
            "maxPrincipalStress",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    scalarField& maxPrincipalStressI = maxPrincipalStress.internalField();
    vectorField& eigenValI = eigenVal.internalField();

    forAll(maxPrincipalStressI, cellI)
    {
        maxPrincipalStressI[cellI] =
            max
            (
                eigenValI[cellI].x(),
                max(eigenValI[cellI].y(), eigenValI[cellI].z())
            );
    }
    maxPrincipalStress.correctBoundaryConditions();

    // Calculate the equivalent Kirchhoff stress
    const volScalarField sigmaEq =
        sqrt((3.0/2.0)*magSqr(dev(tau)))
        + dimensionedScalar("SMALL", dimPressure, SMALL);

    // Update damage field
    const volScalarField newDamage =
        damage().oldTime()
      + (
            max
            (
                maxPrincipalStress,
                dimensionedScalar("zero", dimPressure, 0.0)
            )*DEpsilonPEq///sigmaEq
        )/Vd_;

    // Be careful with multi-material as the damage field is shared by all the
    // damage laws
    damage() = curMaterial()*newDamage + (1.0 - curMaterial())*damage();
}


// ************************************************************************* //
