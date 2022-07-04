if(thermalSimulation)
{
    volScalarField mu = mixture.mu() + rho*turbulence->nut();

    thermalConductivityModel->correct();
    heatCapacityModel->correct();

    kappa = alphal*kappal + alphav*kappav ;
    Cp = (alphal*rhol*Cpl + alphav*rhovSat*Cpv)/rho;

    volScalarField hl =
    (
        hl0 + Cpl*(T - T0) + EOSModel->intDhDp()
    );

    volScalarField hv =
    (
        hl0 + hEvap + Cpv*(T - T0)
    );

    volScalarField K = 0.5*magSqr(U);

 	volTensorField gradU = fvc::grad(U);
	volTensorField tau = -1.0*mu*((gradU + gradU.T()) + (2.0/3.0*fvc::div(U))*I);

    fvScalarMatrix hEqn
    (
        (
          rho*Cp*fvm::ddt(T)
        + rho*Cp*fvm::div(phi, T, "div(phi,T)")
        // - fvm::SuSp(fvc::div(phi + fvc::meshPhi(U)), T)
        )  //T* fvc::div(phi )   // the problem related to the ddt scheme
      + (rhovSat*hv
      - rholSat*hl)*(fvc::ddt(alphav) + fvc::div(phi, alphav) - fvc::div(phi/*+ fvc::meshPhi(U)*/)*alphav)
      ==
        fvm::laplacian(kappa, T)
      - (tau && gradU)
      + alphav * fvc::DDt(phi, p)  //third DpDt form*/
    );

    hEqn.relax();
    hEqn.solve();

    T.relax();
    
    Info<< "max-min T: " << max(T).value()
        << " " << min(T).value() << endl;

    T= min(max(T, T0), dimensionedScalar(dimTemperature, 150.0) + T0); // bounded is essential for steady state simulation

    EOSModel->correctRho();
    EOSModel->correctPsi();
}
