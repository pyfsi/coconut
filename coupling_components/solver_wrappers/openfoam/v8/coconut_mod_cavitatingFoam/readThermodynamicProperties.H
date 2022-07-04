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

// added

dimensionedScalar CpLiquid
(
    "CpLiquid",
    dimSpecificHeatCapacity,
    thermodynamicProperties
);

dimensionedScalar CpGas
(
    "CpGas",
    dimSpecificHeatCapacity,
    thermodynamicProperties
);

dimensionedScalar KLiquid
(
    "KLiquid",
    dimPower/dimLength/dimTemperature,
    thermodynamicProperties
);

dimensionedScalar KGas
(
    "KGas",
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

dimensionedScalar hL0(
    "hL0",
    dimEnergy/dimMass,
    thermodynamicProperties
);

dimensionedScalar pConst(
    "pConst",
    dimPressure,
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
Switch pEqnBarotropicBase(fluidThermalProperties.lookupOrDefault<Switch>("pEqnBarotropicBase", true));
Switch UEqnNewVersion;(fluidThermalProperties.lookupOrDefault<Switch>("UEqnNewVersion;", true));
// end Peyman