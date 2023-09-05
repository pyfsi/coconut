{
    thermalConductivityModel->correct();
    heatCapacityModel->correct();

    rhov = psiv*p;

    k = alphal*kl + alphav*kv;
    rhoCp = alphal*rhol*cpl + alphav*rhov*cpv;

    hl = hl0 + cpl*(T - T0) + EOSModel->intDhDp();
    hv = hl0 + hEvap + cpv*(T - T0);
    rhoH =  alphal*rhol*hl + alphav*rhov*hv;

 	volTensorField gradU = fvc::grad(U);
	volSymmTensorField tau = 2*rho*turbulence->nuEff()*dev(symm(gradU));

    fvScalarMatrix TEqn
    (
        fvc::ddt(rhoH)
      + correction(fvm::ddt(rhoCp, T))
      + fvc::div(phi, rhoH)
      + correction(fvm::div(fvc::interpolate(rhoCp)*phi, T))
      ==
        fvm::laplacian(k, T)
      + (tau && gradU)
      + fvc::DDt(phi, p)
    );

    TEqn.solve();
    
    Info<< "max-min T: " << max(T).value()
        << " " << min(T).value() << endl;

    rholSat = EOSModel->calcRho(pSat, T);
    EOSModel->correctRho();
    EOSModel->correctPsi();
}
