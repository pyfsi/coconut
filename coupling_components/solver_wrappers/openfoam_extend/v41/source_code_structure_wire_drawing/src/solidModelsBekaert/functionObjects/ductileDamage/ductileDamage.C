/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "ductileDamage.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "mechanicalModel.H"
#include "logVolFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ductileDamage, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        ductileDamage,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::ductileDamage::writeData()
{
    writeDamageLemaitre();
    writeDamageCockcroftLatham();

    Info<<endl;

    return true;
}


bool Foam::ductileDamage::writeDamageLemaitre()
{
    // Total Lemaitre damage
    volScalarField damageLemaitre
    (
        IOobject
        (
            "damageLemaitre",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    // Add on old accumulated damage
    // Note: the reason we do not store the field is because the function
    // object gets destructed when the time class is destructed and it
    // causes a seg fault when destructing geometric fields members of a
    // function object; I am not sure of an elegant solution but this will
    // do.

    // Careful: we should check if the number of cells changed.
    if (damageLemaitre.internalField().size() != oldDamageLemaitreI_.size())
    {
        WarningIn("bool Foam::ductileDamage::writeDamageLemaitre()")
            << "The number of cells changed: damageLemaitre not written!"
            << endl;
    }
    else
    {
        damageLemaitre.internalField() = oldDamageLemaitreI_;
        damageLemaitre.correctBoundaryConditions();

        // Lookup Cauchy stress field
        const volSymmTensorField& sigma =
            mesh_.lookupObject<volSymmTensorField>("sigmaCauchy");

        // Calculate equivalent Cauchy stress
        const volScalarField sigmaEq =
            sqrt((3.0/2.0)*magSqr(dev(sigma)))
          + dimensionedScalar("SMALL", dimPressure, SMALL);

        // Calculate hydrostatic Cauchy stress
        const volScalarField sigmaHyd = (1.0/3.0)*tr(sigma);

        // Calculate constraint/triality
        const volScalarField triaxiality = sigmaHyd/sigmaEq;

        // Lookup the increment of equivalent plastic strain
        const volScalarField epsilonPEq = this->epsilonPEq();
        const volScalarField DEpsilonPEq = this->DEpsilonPEq();


        // Lookup elastic parameters

        // Lame's first parameter
        const volScalarField lambda =
            mesh_.lookupObject<mechanicalModel>
            (
                "mechanicalProperties"
            ).lambda();

        // Lame's second parameter, aka the shear modulus
        const volScalarField& mu = mesh_.lookupObject<volScalarField>("mu");

        // Young's modulus
        const volScalarField E = mu*(3.0*lambda + 2.0*mu)/(lambda + mu);

        // Poisson's ration
        const volScalarField nu = lambda/(2*(lambda + mu));


        // Calculate Lemaitre damage

        scalarField& damageLemaitreI = damageLemaitre.internalField();
        const scalarField& sigmaEqI = sigmaEq.internalField();
        const scalarField& epsilonPEqI = epsilonPEq.internalField();
        const scalarField& DEpsilonPEqI = DEpsilonPEq.internalField();
        const scalarField& triaxialityI = triaxiality.internalField();
        const scalarField& EI = E.internalField();
        const scalarField& nuI = nu.internalField();

        const scalar epsilonDValue = epsilonD_.value();

        forAll(damageLemaitreI, cellI)
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

                dimensionedScalar Y("zero", dimPressure, 0.0);

                if (triaxialityI[cellI] > 0.0)
                {
                    // Damage in tension
                    Y.value() =
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
                    Y.value() =
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

                damageLemaitreI[cellI] +=
                    (pow(-Y/S0_, b_)*DEpsilonPEqI[cellI]).value();
            }
        }

        damageLemaitre.correctBoundaryConditions();

        // Write field
        if (runTime_.outputTime())
        {
            damageLemaitre.write();
        }

        // Update old field
        oldDamageLemaitreI_ = damageLemaitre.internalField();

        Info<< "Max damage (Lemaitre): " << gMax(damageLemaitre) << endl;
    }

    return true;
}


bool Foam::ductileDamage::writeDamageCockcroftLatham()
{
    // Total Cockcroft-Latham damage
    volScalarField damageCockcroftLatham
    (
        IOobject
        (
            "damageCockcroftLatham",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    // Add on old accumulated damage
    // Note: the reason we do not store the field is because the function
    // object gets destructed when the time class is destructed and it
    // causes a seg fault when destructing geometric fields members of a
    // function object; I am not sure of an elegant solution but this will
    // do.

    // Careful: we should check if the number of cells changed.
    if
    (
        damageCockcroftLatham.internalField().size()
     != oldDamageCockcroftLathamI_.size()
    )
    {
        WarningIn("bool Foam::ductileDamage::writeDamageCockcroftLatham()")
            << "The number of cells changed: damageCockcroftLatham not written!"
            << endl;
    }
    else
    {
        damageCockcroftLatham.internalField() = oldDamageCockcroftLathamI_;
        damageCockcroftLatham.correctBoundaryConditions();

        // Lookup Kirchhoff stress field
        const volSymmTensorField& tau =
            mesh_.lookupObject<volSymmTensorField>("tauKirchhoff");

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
            mesh_,
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
            mesh_,
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
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
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

        // Add on increment of damageLemaitre
        damageCockcroftLatham +=
            (
                max
                (
                    maxPrincipalStress,
                    dimensionedScalar("zero", dimPressure, 0.0)
                )*DEpsilonPEq()/sigmaEq
            )/cockcroftLathamVd_;

        // Write field
        if (runTime_.outputTime())
        {
            damageCockcroftLatham.write();
        }

        // Update old field
        oldDamageCockcroftLathamI_ = damageCockcroftLatham.internalField();

        Info<< "Max damage (Cockcroft-Latham): " << gMax(damageCockcroftLatham)
            << endl;
    }

    return true;
}



Foam::tmp<Foam::volScalarField> Foam::ductileDamage::epsilonPEq() const
{
    tmp<volScalarField> tepsilonPEq
    (
        new volScalarField
        (
            IOobject
            (
                "ductileDamageEpsilonPEq",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (!mesh_.foundObject<volScalarField>("materials"))
    {
        tepsilonPEq() = mesh_.lookupObject<volScalarField>("epsilonPEq");
    }
    else
    {
        const volScalarField& materials =
            mesh_.lookupObject<volScalarField>("materials");

        const int nMat = gMax(materials);

        for (int matI = 0; matI < nMat; matI++)
        {
            if
            (
                mesh_.foundObject<volScalarField>
                (
                    "epsilonPEq" + Foam::name(matI)
                )
            )
            {
                const volScalarField& epsilonPEqMatI =
                    mesh_.lookupObject<volScalarField>
                    (
                        "epsilonPEq" + Foam::name(matI)
                    );

                tepsilonPEq() = max(tepsilonPEq(), epsilonPEqMatI);
            }
        }
    }

    return tepsilonPEq;
}


Foam::tmp<Foam::volScalarField> Foam::ductileDamage::DEpsilonPEq() const
{
    tmp<volScalarField> tDEpsilonPEq
    (
        new volScalarField
        (
            IOobject
            (
                "ductileDamageDEpsilonPEq",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (!mesh_.foundObject<volScalarField>("materials"))
    {
        tDEpsilonPEq() = mesh_.lookupObject<volScalarField>("DEpsilonPEq");
    }
    else
    {
        const volScalarField& materials =
            mesh_.lookupObject<volScalarField>("materials");

        const int nMat = gMax(materials);

        for (int matI = 0; matI < nMat; matI++)
        {
            if
            (
                mesh_.foundObject<volScalarField>
                (
                    "DEpsilonPEq" + Foam::name(matI)
                )
            )
            {
                const volScalarField& DEpsilonPEqMatI =
                    mesh_.lookupObject<volScalarField>
                    (
                        "DEpsilonPEq" + Foam::name(matI)
                    );

                tDEpsilonPEq() = max(tDEpsilonPEq(), DEpsilonPEqMatI);
            }
        }
    }

    return tDEpsilonPEq;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ductileDamage::ductileDamage
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    runTime_(t),
    mesh_
    (
        runTime_.lookupObject<fvMesh>
        (
            dict.lookupOrDefault<word>("region", "region0")
        )
    ),
    oldDamageLemaitreI_(mesh_.nCells(), 0.0),
    oldDamageCockcroftLathamI_(mesh_.nCells(), 0.0),
    // Default Lemaitre values taken from Masse thesis (2010)
    h_
    (
        dict.subDict("Lemaitre").lookupOrDefault<dimensionedScalar>
        (
            "h",
            dimensionedScalar("h", dimPressure, 0.2)
        )
    ),
    S0_
    (
        dict.subDict("Lemaitre").lookupOrDefault<dimensionedScalar>
        (
            "S0",
            dimensionedScalar("S0", dimPressure, 5.23e6)
        )
    ),
    b_
    (
        dict.subDict("Lemaitre").lookupOrDefault<dimensionedScalar>
        (
            "b",
            dimensionedScalar("b", dimless, 3.14)
        )
    ),
    epsilonD_
    (
        dict.subDict("Lemaitre").lookupOrDefault<dimensionedScalar>
        (
            "epsilonD",
            dimensionedScalar("epsilonD", dimless, 0.1)
        )
    ),
    cockcroftLathamVd_
    (
        dict.subDict("Cockcroft-Latham").lookupOrDefault<dimensionedScalar>
        (
            "Vd",
            dimensionedScalar("cockcroftLathamVd", dimless, 1.0)
        )
    )
{
    Info<< "Creating " << this->name() << " function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ductileDamage::start()
{
    return writeData();
}


#if FOAMEXTEND > 40
bool Foam::ductileDamage::execute(const bool forceWrite)
#else
bool Foam::ductileDamage::execute()
#endif
{
    return writeData();
}


bool Foam::ductileDamage::read(const dictionary& dict)
{
    return true;
}

// ************************************************************************* //
