{
    thermalConductivityModel->correct();
    heatCapacityModel->correct();

    rhov = psiv*p;

    Info<< "max-min rhov.oldTime(): " << max(rhov.oldTime()).value()
        << " " << min(rhov.oldTime()).value() << endl;
    Info<< "max-min rhov: " << max(rhov).value()
        << " " << min(rhov).value() << endl;

    kappa = alphal*kappal + alphav*kappav ;
    Cp = (alphal*rhol*Cpl + alphav*rhov*Cpv)/rho;

    volScalarField hl =
    (
        hl0 + Cpl*(T - T0) + EOSModel->intDhDp()
    );

    volScalarField hv =
    (
        hl0 + hEvap + Cpv*(T - T0)
    );

 	volTensorField gradU = fvc::grad(U);
 	volScalarField mu = mixture.mu() + rho*turbulence->nut();
	volSymmTensorField tau = 2*mu*dev(symm(gradU));

    Info<< "max-min phi: " << max(tau && gradU).value()
        << " " << min(tau && gradU).value() << endl;

    Info<< "max-min alphal.oldTime(): " << max(alphal.oldTime()).value()
        << " " << min(alphal.oldTime()).value() << endl;
    Info<< "max-min alphal: " << max(alphal).value()
        << " " << min(alphal).value() << endl;

    Info<< "max-min fvc::ddt(alphal): " << max(fvc::ddt(alphal)).value()
        << " " << min(fvc::ddt(alphal)).value() << endl;

    Info<< "max-min fvc::ddt(alphav): " << max(fvc::ddt(alphav)).value()
        << " " << min(fvc::ddt(alphav)).value() << endl;

    Info<< "max-min fvc::ddt(p): " << max(fvc::ddt(p)).value()
        << " " << min(fvc::ddt(p)).value() << endl;

    Info<< "max-min fvc::ddt(psiv, p): " << max(fvc::ddt(psiv, p)).value()
        << " " << min(fvc::ddt(psiv, p)).value() << endl;

    Info<< "max-min fvc::ddt(rhov): " << max(fvc::ddt(rhov)).value()
        << " " << min(fvc::ddt(rhov)).value() << endl;

    Info<< "max-min fvc::ddt(rhol): " << max(fvc::ddt(rhol)).value()
        << " " << min(fvc::ddt(rhol)).value() << endl;

    fvScalarMatrix TEqn
    (
        (
          rho*Cp*fvm::ddt(T)
        + rho*Cp*fvm::div(phi, T, "div(phi,T)")
        - rho*Cp*fvm::Sp(fvc::div(phi), T)
        // - fvm::SuSp(fvc::div(phi + fvc::meshPhi(U)), T)
        )  //T* fvc::div(phi )   // the problem related to the ddt scheme

      //+ (rhovSat*hv - rholSat*hl)*(fvc::ddt(alphav) + fvc::div(phi, alphav) - fvc::div(phi/*+ fvc::meshPhi(U)*/)*alphav)

      + hv*(psiv*fvc::ddt(alphav, p) + fvc::div(phi, alphav*rhov) - alphav*rhov*fvc::div(phi))
      + hl*(fvc::ddt(alphal, rhol) + fvc::div(phi, alphal*rhol) - alphal*rhol*fvc::div(phi))
      + (alphav*rhov*hv + alphal*rhol*hl)*fvc::div(phi)
      ==
        fvm::laplacian(kappa, T)
      + (tau && gradU)
      + (alphav + alphal*EOSModel->betaVolExp()*T)*fvc::DDt(fvc::absolute(phi, U), p) // correct? DDt
    );

    volScalarField temp3 = fvc::ddt(alphal, rhol);
    Info<< "max-min ddt(alphal, rhol): " << max(temp3).value()
        << " " << min(temp3).value() << endl;

    volScalarField temp4 = alphal*fvc::ddt(rhol)-rhol*fvc::ddt(alphav);
    Info<< "max-min ddt(alphal, rhol): " << max(temp4).value()
        << " " << min(temp4).value() << endl;

    volScalarField temp5 = alphal*fvc::ddt(rhol)+rhol*fvc::ddt(alphal);
    Info<< "max-min ddt(alphal, rhol): " << max(temp5).value()
        << " " << min(temp5).value() << endl;

    volScalarField temp7 = fvc::DDt(phi, p);
    Info<< "max-min DDt(p): " << max(temp7).value()
        << " " << min(temp7).value() << endl;

    volScalarField temp6 = fvc::DDt(fvc::absolute(phi, U), p);
    Info<< "max-min DDt(p): " << max(temp6).value()
        << " " << min(temp6).value() << endl;

    volScalarField temp = fvc::ddt(p) + fvc::div(phi, p) - p*fvc::div(phi);
    Info<< "max-min DDt(p): " << max(temp).value()
        << " " << min(temp).value() << endl;

    volScalarField temp0 = fvc::ddt(p) + fvc::div(phi, p) - p*fvc::div(fvc::absolute(phi, U));
    Info<< "max-min DDt(p): " << max(temp0).value()
        << " " << min(temp0).value() << endl;

    volScalarField temp2 = fvc::DDt(phi, p) - (fvc::ddt(p) + fvc::div(phi, p) - p*fvc::div(fvc::absolute(phi, U)));
    Info<< "max-min DDt(p): " << max(temp2).value()
        << " " << min(temp2).value() << endl;
    Info<< temp2.weightedAverage(mesh.V());
    Info<< temp2.weightedAverage(mesh.V()).value();

//    hEqn.relax();
    TEqn.solve();

//    T.relax();
    
    Info<< "max-min T: " << max(T).value()
        << " " << min(T).value() << endl;

    //T= min(max(T, T0), dimensionedScalar(dimTemperature, 150.0) + T0); // bounded is essential for steady state simulation

    EOSModel->correctRho();
    EOSModel->correctPsi();
}
