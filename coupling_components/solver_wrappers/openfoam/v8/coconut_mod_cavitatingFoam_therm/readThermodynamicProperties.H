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

Switch thermalSimulation(thermodynamicProperties.lookupOrDefault<Switch>("thermalSimulation", true));

dictionary thermalProperties = thermodynamicProperties.subDict("thermalProperties");

dimensionedScalar T0
(
    "T0",
    dimTemperature,
    thermalProperties
);

dimensionedScalar hEvap
(
    "hEvap",
    dimEnergy/dimMass,
    thermalProperties
);

dimensionedScalar hl0
(
    "hl0",
    dimEnergy/dimMass,
    0
);

if (thermalProperties.found("hl0"))
{
    thermalProperties.lookup("hl0") >> hl0;
}

dictionary vapourThermalProperties = thermalProperties.subDict("vapourThermalProperties");

dimensionedScalar cpv
(
    "cpv",
    dimSpecificHeatCapacity,
    vapourThermalProperties
);

dimensionedScalar kv
(
    "kv",
    dimPower/dimLength/dimTemperature,
    vapourThermalProperties
);

dictionary liquidThermalProperties = thermalProperties.subDict("liquidThermalProperties");
