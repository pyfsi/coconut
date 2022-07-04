#include "readThermodynamicProperties.H"

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

mesh.setFluxRequired(p.name());

// Mass flux (corrected by rhoEqn.H)
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

autoPtr<EOSModel> EOSModel =
    EOSModel::New
    (
        "EOSModel",
        thermodynamicProperties,
        p,
        T,
        "rhol",
        "psil"
    );

const volScalarField& rhol(EOSModel->rho());
const volScalarField& psil(EOSModel->psi());

volScalarField rholSat(EOSModel->calcRho(pSat, T));

Info<< "Reading transportProperties\n" << endl;

incompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alphav(mixture.alpha1());
alphav.oldTime();

volScalarField& alphal(mixture.alpha2());

Info<< "Creating compressibilityModel\n" << endl;

autoPtr<mixtureCompressibilityModel> psiModel =
    mixtureCompressibilityModel::New
    (
        "mixComprModel",
        thermodynamicProperties,
        alphav,
        psil
    );

const volScalarField& psi = psiModel->psi();
Info<< "psi max-min : " << max(psi).value()
    << " " << min(psi).value() << endl;

rho == max(alphav*psiv*p + alphal*rhol, rhoMin);

mesh.setFluxRequired(p.name());
mesh.setFluxRequired(rho.name());

// Create incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);

autoPtr<thermalConductivityModel> thermalConductivityModel=
    thermalConductivityModel::New
    (
        "kappal",
        fluidThermalProperties,
        U,
        p,
        rhol,
        T,
        phi
    );

const volScalarField& kappal(thermalConductivityModel->kappa());

autoPtr<heatCapacityModel> heatCapacityModel=
  heatCapacityModel::New
    (
        "Cpl",
        fluidThermalProperties,
        U,
        p,
        rhol,
        T,
        phi
    );

const volScalarField& Cpl(heatCapacityModel->Cp());

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alphal*kappal + alphav*kappav
);

volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alphal*Cpl + alphav*Cpv
);

volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    (hL0 + hEvap + Cpv*(T - T0))*alphav + (hL0 + Cpl*(T - T0) /*+ EOSModel->IntRhoInv()*/)*alphal,
    // CpE * (T - T0),
    T.boundaryField().types()
);
