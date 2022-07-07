Info<< "Reading thermodynamicProperties\n" << endl;

IOdictionary thermodynamicProperties
(
    IOobject
    (
        "thermodynamicProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar psiv
(
    "psiv",
    dimCompressibility,
    thermodynamicProperties
);

dimensionedScalar pSat
(
    "pSat",
    dimPressure,
    thermodynamicProperties
);

dimensionedScalar rhovSat("rhovSat", psiv*pSat);

dimensionedScalar rhoMin
(
    "rhoMin",
    dimDensity,
    thermodynamicProperties
);

dimensionedScalar Cpv
(
    "Cpv",
    dimSpecificHeatCapacity,
    thermodynamicProperties
);

dimensionedScalar kappav
(
    "kappav",
    dimPower/dimLength/dimTemperature,
    thermodynamicProperties
);

dimensionedScalar T0(
    "T0",
    dimTemperature,
    thermodynamicProperties
);

dimensionedScalar hEvap(
    "hEvap",
    dimEnergy/dimMass,
    thermodynamicProperties
);

dimensionedScalar hl0(
    "hl0",
    dimEnergy/dimMass,
    thermodynamicProperties
);

IOdictionary fluidThermalProperties
(
    IOobject
    (
        "fluidThermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Switch thermalSimulation(fluidThermalProperties.lookupOrDefault<Switch>("thermalSimulation", true));
