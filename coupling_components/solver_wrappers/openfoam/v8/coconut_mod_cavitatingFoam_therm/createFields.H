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

volScalarField rhov = psiv*p;

Info<< "" << endl;
Info<< "Reading transportProperties\n" << endl;

compressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alphav(mixture.alpha1());

volScalarField& alphal(mixture.alpha2());

Info<< "" << endl;
Info<< "Creating compressibilityModel\n" << endl;

autoPtr<mixtureCompressibilityModel> psiModel =
    mixtureCompressibilityModel::New
    (
        "mixComprModel",
        thermodynamicProperties,
        alphav,
        psil,
        rholSat
    );

const volScalarField& psi = psiModel->psi();
Info<< "max-min psi: " << max(psi).value()
    << " " << min(psi).value() << endl;

rho == alphav*psiv*p + alphal*rhol;
Info<< "max-min EOS rho: " << max(rho).value()
    << " " << min(rho).value() << endl;

Info<< "" << endl;
Info<< "Creating momentumTransportModel\n" << endl;

// Create incompressible turbulence model
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho*turbulence->nuEff()
);

Info<< "" << endl;
Info<< "Creating thermal models\n" << endl;

autoPtr<thermalConductivityModel> thermalConductivityModel=
    thermalConductivityModel::New
    (
        "kl",
        liquidThermalProperties,
        U,
        p,
        rhol,
        T,
        phi
    );

const volScalarField& kl(thermalConductivityModel->k());

autoPtr<heatCapacityModel> heatCapacityModel=
  heatCapacityModel::New
    (
        "cpl",
        liquidThermalProperties,
        U,
        p,
        rhol,
        T,
        phi
    );

const volScalarField& cpl(heatCapacityModel->cp());

volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alphal*kl + alphav*kv
);

volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (alphal*rhol*cpl + alphav*rhov*cpv)
);
rhoCp.oldTime(); // sets initial value as 'old value' for ddt

Info<< "" << endl;
Info<< "Creating field hl\n" << endl;
volScalarField hl
(
    IOobject
    (
        "hl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    hl0 + cpl*(T - T0) + EOSModel->intDhDp()
);

Info<< "Creating field hv\n" << endl;
volScalarField hv
(
    IOobject
    (
        "hv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    hl0 + hEvap + cpv*(T - T0)
);

// p and T0 contribution to enthalpy
volScalarField rhoH =  alphal*rhol*hl + alphav*rhov*hv;
rhoH.oldTime();
