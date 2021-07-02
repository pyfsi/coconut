/*---------------------------------------------------------------------------*\
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

Global
    readSolidMechanicsControls

Description
    Read solid mechanics solver control parameters

\*---------------------------------------------------------------------------*/

dictionary& solidControl = mesh.solutionDict().subDict("solidMechanics");

const int nCorr(solidControl.lookupOrDefault<int>("nCorrectors", 1000));

const int nCorrT(solidControl.lookupOrDefault<int>("nTCorrectors", nCorr));

const scalar convergenceToleranceU
(
    solidControl.lookupOrDefault<scalar>("DU", 1e-6)
);

const scalar convergenceToleranceT
(
    solidControl.lookupOrDefault<scalar>("T", convergenceToleranceU)
);

const scalar convergenceTolRelU
(
    solidControl.lookupOrDefault<scalar>("relDU", convergenceToleranceU)
);

const scalar convergenceTolRelT
(
    solidControl.lookupOrDefault<scalar>("relT", convergenceToleranceT)
);

scalar convergenceTolMaterial
(
    solidControl.lookupOrDefault<scalar>
    (
        "material", 100*convergenceToleranceU
    )
);

// For backwards compatibility, check if DEpsilonP tolerance is specified
// instead
if (solidControl.found("DEpsilonP"))
{
    Warning
        << "Reading deprecated 'DEpsilonP' tolerance;" << nl
        << "'material' tolerance should now be specified instead of 'DEpsilon'"
        << endl;

    convergenceTolMaterial = readScalar(solidControl.lookup("DEpsilonP"));
}

// Frequency of residual print out
const int infoFrequency(solidControl.lookupOrDefault<int>("infoFrequency", 10));

// Use diffusion equation to initialise the solution field as a predictor step
const Switch predictor
(
    solidControl.lookupOrDefault<Switch>
    (
        "predictor",
        false
    )
);

const Switch allowRestart
(
    solidControl.lookupOrDefault<Switch>
    (
        "allowRestart",
        true
    )
);

if (!allowRestart)
{
    // These fields are only required for restarts
    gradDU.writeOpt() = IOobject::NO_WRITE;
    F.writeOpt() = IOobject::NO_WRITE;
    relF.writeOpt() = IOobject::NO_WRITE;
    relFinv.writeOpt() = IOobject::NO_WRITE;
    J.writeOpt() = IOobject::NO_WRITE;
    relJ.writeOpt() = IOobject::NO_WRITE;
    tau.writeOpt() = IOobject::NO_WRITE;
    if (thermalStress)
    {
        gradTPtr->writeOpt() = IOobject::NO_WRITE;
    }
}

Switch updatePointFields
(
    solidControl.lookupOrDefault<Switch>("updatePointFields", false)
);

// For backwards compatibility, check if DEpsilonP tolerance is specified
// instead
// if (solidControl.found("updatePointDU"))
// {
//     Warning
//         << "Reading deprecated 'updatePointDU' switch;" << nl
//         << "'updatePointFields' should now be specified instead of "
//         << "'updatePointDU'"
//         << endl;

//     updatePointFields = Switch(solidControl.lookup("updatePointDU"));
// }

const Switch stabilisePressure
(
    solidControl.lookupOrDefault<Switch>("stabilisePressure", false)
);

if (stabilisePressure)
{
    //updatePointDU = true; // now not needed
    Info<< nl << "stabilisePressure is on" << endl;
}

// Adjustable time-step
const Switch adjustTimeStep
(
    runTime.controlDict().lookupOrDefault<Switch>("adjustTimeStep", false)
);

// Maximum allowed time-step
const scalar maxDeltaT
(
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT)
);

// Minimum allowed time-step
const scalar minDeltaT
(
    runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", 0.0)
);

if (adjustTimeStep)
{
    Info<< nl << "Adjustable time-step" << nl
        << "    Max deltaT " << maxDeltaT << nl
        << "    Min deltaT " << minDeltaT << nl << endl;
}


// The energy equation time-interval during Eulerian mode: this is how much
// physical time passes in the energy equation per momentum equation
// time-step
const scalar accelerationInterval =
    solidControl.lookupOrDefault<scalar>("accelerationInterval", 0.1);


// Thermal adaptive time-step settings

// Maximum Courant number
const scalar thermalMaxCo =
    runTime.controlDict().lookupOrDefault<scalar>("thermalMaxCo", 1.0);

// Maximum deltaT
scalar thermalMaxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>
    (
        "thermalMaxDeltaT",
        GREAT
    );

// Factor to control the magnitude of the oscillation smoothing/diffusion term
// in the momentum equation. When this factor is 1.0, it can excessively smooth
// regions of local plasticity
if (solidControl.found("smoothScale"))
{
    FatalError
        << "smoothScale should no longer be set in fvSolution; instead it should"
        << " be set individually for each mechanical law in "
        << "mechanicalProperties as 'RhieChowScaleFactor'" << abort(FatalError);
}

// Optional: check if it is an RVE case; if so, read specified average
// total deformation gradient time-series
Switch representativeVolumeElementCase = false;
interpolationTable<tensor> avgFSeries;
if (solidControl.found("representativeVolumeElementCase"))
{
    representativeVolumeElementCase = true;
    avgFSeries = interpolationTable<tensor>
    (
        solidControl.subDict("representativeVolumeElementCase")
    );
}

// ************************************************************************* //
