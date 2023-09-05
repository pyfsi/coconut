correctUphiBCs(U, phi, true);

{
    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(p.dimensions(), 0),
        pcorrTypes
    );

    surfaceScalarField rhof(fvc::interpolate(rho));
    dimensionedScalar rAUf("rAUf", dimTime, 1.0);

    mesh.setFluxRequired(pcorr.name());

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
          // fvm::ddt(psi, pcorr) +
             fvc::div(phi, rho)
           - fvm::laplacian(rAUf, pcorr)
          ==
             divRhoU()
        );

        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pcorrEqn.flux()/rhof;
        }
    }
}
