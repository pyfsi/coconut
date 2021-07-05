{
    volScalarField& T = *TPtr;
    const volScalarField& Told = *ToldPtr;
    const volScalarField& k = *kPtr;
    (*rhoCPtr) = rho*(*CPtr);
    (*rhoCfPtr) = fvc::interpolate(*rhoCPtr, "rhoC");
    const volScalarField& rhoC = *rhoCPtr;
    const surfaceScalarField& rhoCf = *rhoCfPtr;
    volScalarField& DT = *DTPtr;
    //surfaceScalarField& DTf = *DTfPtr;
    volScalarField& QT = *QTPtr;
    //volVectorField& gradT = *gradTPtr;
    //pointScalarField& pointT = *pointTPtr;
    //volScalarField& contact = *contactPtr;;
    //const volScalarField& threeKalpha = *threeKalphaPtr;
    surfaceScalarField& phiT = *phiTPtr;
    //const volScalarField Told("Told", 1.0*T);
    //volScalarField& TprevIter = *TprevIterPtr;


    // Plastic dissipated energy
    QT = mechanical.plasticDissipationRate();


    // Solve diffusion equation i.e. pure Lagrangian
    {
        int iCorr = 0;
        lduSolverPerformance solverPerfT;
        scalar relativeResidualT = 1.0;
        bool converged = false;
        lduMatrix::debug = 0;

        Info<< nl << "Solving diffusion energy equation for T" << endl;

        do
        {
            T.storePrevIter();

            // Heat conduction energy equation
            fvScalarMatrix TEqn
            (
                rhoC*fvm::ddt(T)
              - fvm::laplacian(k, T, "laplacian(k,T)")
             == QT
            );

            // Under-relax equation
            TEqn.relax();

            // Solve equation
            solverPerfT = TEqn.solve(mesh.solutionDict().solver("diffusionT"));

            // Under-relax the temperature field
            T.relax();

            // Update gradient of temperature
            // Only needed if newSkewCorrected Laplacian scheme is used
            //gradT = fvc::grad(T);

#           include "checkConvergenceT.H"
        }
        while (!converged && ++iCorr < nCorrT);
    }


    // Solve convection diffusion equation i.e. pure Eulerian
    // By taking an Eulerian perspective, we accelerate the energy equation in
    // time
    // Enable acceleration after 4 thermal diffusion time-steps)
    if (runTime.value() >= (thermalStressStartTime + 4.0*runTime.deltaTValue()))
    {
        Info<< nl << "Solving convection-diffusion energy equation for T"
            << endl;

        // Energy equation sub time loop

        // Update the velocity field: this is the average velocity over the
        // current time-step
        V = DU/runTime.deltaT();

        // Calculate the velocity at the start of the current time-step
        // This is to ensure the velocity flux vectors are consistent with the
        // current mesh configuration (previously we moved the mesh to the
        // mid-point configuration of the time-step)
        const volVectorField VThermal = 0.5*(V.oldTime() + V);

        // Calculate the point velocity field by interpolating from the cells to
        // points
        // Then we will interpolate to faces from the points to calculate the
        // flux
        const newLeastSquaresVolPointInterpolation& volToPointInterp =
            mesh.lookupObject<newLeastSquaresVolPointInterpolation>
            (
                "newLeastSquaresVolPointInterpolation"
            );

        const pointMesh& pMesh = mesh.lookupObject<pointMesh>("pointMesh");

        // Point velocity field
        pointVectorField pointVThermal
        (
            IOobject
            (
                "pointVelocity",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector("zero", dimVelocity, vector::zero),
            calculatedFvPatchVectorField::typeName
        );

        // Interpolate the point velocity field to the faces
        volToPointInterp.interpolate(VThermal, pointVThermal);

        // Calculate the flux of velocity through the faces
        phiT = mesh.Sf() & fvc::interpolate(pointVThermal)*rhoCf;
        // This way may suffer from skewness error
        //phiT = mesh.Sf() & fvc::interpolate(VThermal, "V")*rhoCf;

        // Magnitude of phiT: used to calculate the Courant number
        const surfaceScalarField magPhiT = mag(phiT)/rhoCf;

        // Store quantity that is used to calculate the Courant number
        surfaceScalarField SfUfbyDelta =
            mesh.surfaceInterpolation::deltaCoeffs()*magPhiT;

        // Store the current momentum time-step
        const dimensionedScalar momentumDeltaT = runTime.deltaT();

        // Current thermal time (within this momentum time-step)
        scalar currentThermalTime = 0;

        // Pseudo time loop
        do
        {
            // Caculate the max Courant number
            const scalar deltaT = runTime.deltaT().value();
            const scalar CoNum = max(SfUfbyDelta/mesh.magSf()).value()*deltaT;
            const scalar meanCoNum =
                (sum(SfUfbyDelta)/sum(mesh.magSf())).value()*deltaT;
            const scalar velMag = max(magPhiT/mesh.magSf()).value();

            Info<< nl << "Courant Number mean: " << meanCoNum
                << " max: " << CoNum
                << " velocity magnitude: " << velMag
                << nl << endl;

            // Set the time-step
            const scalar maxDeltaTFact = thermalMaxCo/(CoNum + SMALL);
            const scalar deltaTFact =
                min
                (
                    min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact),
                    4.0
                );

            runTime.setDeltaT
            (
                min
                (
                    deltaTFact*runTime.deltaT().value(),
                    thermalMaxDeltaT
                )
            );

            Info<< "Thermal deltaT = " <<  runTime.deltaT().value()
                << nl << endl;

            // Increase the thermal time
            currentThermalTime += runTime.deltaTValue();

            Info<< "Thermal interval time: " << currentThermalTime
                << nl << endl;

            // Update old time T values
            // We are mimicking runTime++
            T.oldTime() == T;
            T.oldTime().oldTime() == T.oldTime();

            int iCorr = 0;
            lduSolverPerformance solverPerfT;
            scalar relativeResidualT = 1.0;
            bool converged = false;
            lduMatrix::debug = 0;

            do
            {

                // Store T for under-relaxation
                T.storePrevIter();

                // Heat conduction energy equation
                fvScalarMatrix TEqn
                (
                    rhoC*fvm::ddt(T)
                  + fvm::div(phiT, T)
                  - fvm::laplacian(k, T, "laplacian(k,T)")
                 == QT
                );

                // Under-relax equation
                TEqn.relax();

                // Solve equation
                solverPerfT =
                    TEqn.solve(mesh.solutionDict().solver("convectionT"));

                // Under-relax the temperature field
                T.relax();

                // Update gradient of temperature
                // Only needed if newSkewCorrected Laplacian scheme is used
                //gradT = fvc::grad(T);

#               include "checkConvergenceT.H"
            }
            while (!converged && ++iCorr < nCorrT);

            Info<< nl << "Max T: "<< gMax(T) << ", Min T: " << gMin(T)
                << nl << endl;
        }
        while (currentThermalTime < (accelerationInterval - SMALL));

        // Reset the time-step back to the momentum time-step
        runTime.setDeltaT(momentumDeltaT);
    }


    // Update increment of temperature as it may be used bu the mechanical
    // when calculating thermal expansion
    DT = T - Told;
}
