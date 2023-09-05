{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvm::div(phi, rho)
    );

    rhoEqn.solve();

    Info<< "max-min rho: " << max(rho).value()
        << " " << min(rho).value() << endl;
}
