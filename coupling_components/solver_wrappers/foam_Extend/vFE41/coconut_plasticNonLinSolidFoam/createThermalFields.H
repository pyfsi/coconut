// Thermal fields
// Thermal stress is optional

volScalarField* TPtr = NULL;
volScalarField* ToldPtr = NULL;
//volScalarField* TprevIterPtr = NULL;
thermalModel* thermalPtr = NULL;
volScalarField* CPtr = NULL;
volScalarField* kPtr = NULL;
volScalarField* threeKalphaPtr = NULL;
surfaceScalarField* threeKalphafPtr = NULL;
volScalarField* rhoCPtr = NULL;
surfaceScalarField* rhoCfPtr = NULL;
volScalarField* DTPtr = NULL;
//surfaceScalarField* DTfPtr = NULL;
volScalarField* QTPtr = NULL;
volVectorField* gradTPtr = NULL;
//pointScalarField* pointTPtr = NULL;
//volScalarField* contactPtr = NULL;
surfaceScalarField* phiTPtr = NULL;

// Switch to turn on energy equation
Switch thermalStress
(
    mesh.solutionDict().subDict("solidMechanics").lookupOrDefault<Switch>
    (
        "thermalStress",
        false
    )
);

const scalar thermalStressStartTime
(
    mesh.solutionDict().subDict("solidMechanics").lookupOrDefault<scalar>
    (
        "thermalStressStartTime",
        0.0
    )
);

if (thermalStress)
{
    Info<< "Reading thermal options" << endl;

    if (thermalStressStartTime < GREAT)
    {
        Info<< "    thermalStressStartTime: " << thermalStressStartTime << endl;
    }

    Info<< "Reading field T\n" << endl;

    TPtr =
        new volScalarField
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

    Info<< "Reading thermal model\n" << endl;
    thermalPtr = new thermalModel(*TPtr);

    ToldPtr =
        new volScalarField
        (
            IOobject
            (
                "Told",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimTemperature, 0.0)
        );

    // TprevIterPtr =
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             "TprevIter",
    //             runTime.timeName(),
    //             mesh,
    //             IOobject::NO_READ,
    //             IOobject::NO_WRITE
    //         ),
    //         mesh,
    //         dimensionedScalar("zero", dimTemperature, 0.0)
    //     );

    CPtr = new volScalarField(thermalPtr->C());
    kPtr = new volScalarField(thermalPtr->k());
    threeKalphaPtr =
        new volScalarField
        (
            mechanical.threeK()*thermalPtr->alpha()
        );
    threeKalphafPtr =
        new surfaceScalarField
        (
            "threeKalphaf",
            fvc::interpolate(*threeKalphaPtr, "threeKalpha")
        );
    volScalarField T0 = thermalPtr->T0();
    rhoCPtr = new volScalarField(rho*(*CPtr));
    rhoCfPtr =
        new surfaceScalarField
        (
            "rhoCf",
            fvc::interpolate(*rhoCPtr, "rhoC")
        );

    //DTPtr = new volScalarField("DT",(*TPtr) - (*TPtr).oldTime());
    // At the beginning DT is the difference between the starting
    // temperature and the stress-free temperature. For other time-steps
    // DT is T - T.old
    DTPtr = new volScalarField("DT",(*TPtr) - T0);
    // DTfPtr =
    //     new surfaceScalarField
    //     (
    //         "DTf",
    //         fvc::interpolate(*DTPtr, "DT")
    //     );

    QTPtr =
        new volScalarField
        (
            mechanical.plasticDissipationRate()
        );

    gradTPtr =
        new volVectorField
        (
            IOobject
            (
                "grad(T)",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            ("zero", dimTemperature/dimLength, vector::zero)
        );

    // contactPtr =
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             "contact",
    //             runTime.timeName(),
    //             mesh,
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         mesh,
    //         dimensionedScalar("zero", dimless, -1)
    //     );

    // divThreeKalphaDTPtr =
    //     new volVectorField
    //     (
    //         IOobject
    //         (
    //             "divThreeKalphaDT",
    //             runTime.timeName(),
    //             mesh,
    //             IOobject::NO_READ,
    //             IOobject::NO_WRITE
    //         ),
    //         fvc::div
    //         (
    //         //mesh.Sf()*(*threeKalphafPtr)*fvc::interpolate(*DTPtr, "deltaT")
    //             (*threeKalphaPtr)*(*DTPtr)*I
    //         )
    //     );

    phiTPtr =
        new surfaceScalarField
        (
            IOobject
            (
                "phiT",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh.Sf() & fvc::interpolate(V, "V")*(*rhoCfPtr)
        );
}
