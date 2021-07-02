scalar thermalTol = 1e-6;
if (accelerateEnergyEqn)
{
    Info<< nl << "    Thermal Corr = " << iThermalCorr
        << ", Thermal residual = " << thermalResidual << nl << endl;

    if (++iThermalCorr > nThermalCorr)
    {
        Info<< "    Max thermal stress correctors reached" << endl;
        thermalConverged = true;
    }
    else if (thermalResidual < thermalTol)
    {
        Info<< "    Thermal stress converged" << endl;
        thermalConverged = true;
    }
    else
    {
        thermalConverged = false;
    }
}
else
{
    thermalConverged = true;
}
