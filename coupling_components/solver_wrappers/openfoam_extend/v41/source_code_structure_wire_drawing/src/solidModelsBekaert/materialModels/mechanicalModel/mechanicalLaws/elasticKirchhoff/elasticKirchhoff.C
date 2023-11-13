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

#include "elasticKirchhoff.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "transformGeometricField.H"
#include "fvc.H"
#include "mechanicalModel.H"
#include "thermalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(elasticKirchhoff, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, elasticKirchhoff, dictionary
    );


    // Store sqrt(2/3) as we use it often
    scalar elasticKirchhoff::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

Foam::volSymmTensorField& Foam::elasticKirchhoff::bEbar()
{
    if (!bEbarPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& Foam::elasticKirchhoff::bEbar()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *bEbarPtr_;
}


const Foam::volSymmTensorField& Foam::elasticKirchhoff::bEbar() const
{
    if (!bEbarPtr_)
    {
        FatalErrorIn
        (
            "Foam::volSymmTensorField& Foam::elasticKirchhoff::bEbar()"
        )   << "pointer not set" << abort(FatalError);
    }

    return *bEbarPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::elasticKirchhoff::elasticKirchhoff
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const label lawIndex
)
:
    mechanicalLaw(name, mesh, dict, lawIndex),
    restarted_
    (
        IOobject
        (
            "sigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ
        ).headerOk()
    ),
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    K_
    (
        mesh.lookupObject<mechanicalModel>("mechanicalProperties").planeStress()
      ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
      : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
    ),
    bEbarPtr_(NULL)
{
    // Create fields
    // We will check if the field already exists in the object registry (i.e. if
    // it has been created by another mechanical law); if it is not found then
    // the field is created; if it is found, then we will set a pointer to it
    // using const_cast. This may not be the most elegant or safe solution but
    // it is OK for now!

    SetFieldPtr<symmTensor>
    (
        bEbarPtr_,
        "bEbar",
        dimensionedSymmTensor("one", dimless, symmTensor(I))
    );

    // On restart, some fields may be defaulted to zero when they should default
    // to I
    if (gMin(mag(bEbar().internalField())) < SMALL)
    {
        WarningIn("elasticKirchhoff::elasticKirchhoff()")
            << "Reseting zero in bEbar fields to I" << endl;

        symmTensorField& bEbarI = bEbar().internalField();

        forAll(bEbarI, cellI)
        {
            if (mag(bEbarI[cellI]) < SMALL)
            {
                bEbarI[cellI] = 1.0*I;
            }
        }

        forAll(bEbar().boundaryField(), patchI)
        {
            symmTensorField& bEbarP = bEbar().boundaryField()[patchI];

            forAll(bEbarP, faceI)
            {
                if (mag(bEbarP[faceI]) < SMALL)
                {
                    bEbarP[faceI] = 1.0*I;
                }
            }
        }

        bEbar().correctBoundaryConditions();
    }

    // Force sotrage of old time bEbar
    bEbar().oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::elasticKirchhoff::~elasticKirchhoff()
{
    CleanFieldPtr<symmTensor>(bEbarPtr_, "bEbar");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::elasticKirchhoff::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            calculatedFvPatchScalarField::typeName
            //zeroGradientFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::elasticKirchhoff::E() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "E",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            E_,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::elasticKirchhoff::nu() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            nu_,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::elasticKirchhoff::Ep() const
{
    notImplemented
    (
        "Foam::tmp<Foam::volScalarField> "
        "Foam::elasticKirchhoff::Ep() const"
    );

    // Keep compiler happy
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "undefined",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            0.0,
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


void Foam::elasticKirchhoff::correct(volSymmTensorField& tau)
{
    const fvMesh& mesh = this->mesh();

    // Lookup relative deformation gradient from the solver
    const volTensorField& relFbar =
        mesh.lookupObject<volTensorField>("relFbar");

    // Update left Cauchy Green strain tensor with volumetric term removed
    const volSymmTensorField newBEbar = transform(relFbar, bEbar().oldTime());
    bEbar() = curMaterial()*newBEbar + (1.0 - curMaterial())*bEbar();

    // Lookup the Jacobian of the deformation gradient from the solver
    const volScalarField& J = mesh.lookupObject<volScalarField>("J");

    // Calculate deviatoric stress
    volSymmTensorField s = mu_*dev(bEbar());

    // Calculate new Kirchhoff stress

    // Add deviatoric stress term
    volSymmTensorField newTau("newTau", s);

    // Add hydrostatic pressure term
    if (mesh.foundObject<volScalarField>("p"))
    {
        // If the pressure field exists in the solver (e.g. from the hybrid
        // approach) then we will use it
        const volScalarField& p = mesh.lookupObject<volScalarField>("p");

        newTau -= p*symmTensor(I);
    }
    else
    {
        // Otherwise calculate pressure as a function of displacement
        newTau += 0.5*K_*(pow(J, 2) - 1.0)*symmTensor(I);
    }

    // Add thermal component, if active
    if (mesh.foundObject<thermalModel>("thermalProperties"))
    {
        const volScalarField& threeKalpha =
            mesh.lookupObject<volScalarField>("(threeK*alpha)");

        const volScalarField& DT = mesh.lookupObject<volScalarField>("DT");

        newTau += threeKalpha*DT*symmTensor(I);
    }

    // Assign Kirchhoff stress
    // For now, to deal with multi-materials, we will multiply by curMaterial
    // index field so only cells in the current material are calculated:
    // we should be able to do this in a better way

    tau = curMaterial()*newTau + (1.0 - curMaterial())*tau;
}


void Foam::elasticKirchhoff::correct(volSymmTensorField& tau, const int flag)
{
    correct(tau);
}


void Foam::elasticKirchhoff::correct(surfaceSymmTensorField& tau)
{
    notImplemented("elasticKirchhoff::correct(surfaceSymmTensorField& tau)");
}


Foam::tmp<Foam::volScalarField>
Foam::elasticKirchhoff::plasticDissipationRate() const
{
    // Elastic material dissipates zero plastic energy
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "plasticDissipationRate",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimForce/(dimArea*dimTime), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


// ************************************************************************* //
