#include "calculateRelativeResidualT.H"

// If one of the residuals has converged to an order of magnitude
// less than the tolerance then consider the solution converged
if ( iCorr > 0 )// force at leaast 1 outer iteration
{
    if
    (
        solverPerfT.initialResidual() < convergenceToleranceT
     && relativeResidualT < convergenceTolRelT
    )
    {
        Info<< "    Both residuals have converged" << endl;
        converged = true;
    }
    else if (relativeResidualT < 0.1*convergenceTolRelT)
    {
        Info<< "    The relative residual has converged" << endl;
        converged = true;
    }
    else if (solverPerfT.initialResidual() < 0.1*convergenceToleranceT)
    {
        Info<< "    The solver residual has converged" << endl;
        converged = true;
    }
    else
    {
        converged = false;
    }
}
else
{
    converged = false;
}


if (iCorr % infoFrequency == 0 || converged || iCorr == nCorrT)
{
    Info<< "    Time " << runTime.value()
        << ", Corr " << iCorr
        << ", res = " << solverPerfT.initialResidual()
        << ", rel res = " << relativeResidualT
        << ", iters = " << solverPerf.nIterations() << endl;
}
else if (iCorr == nCorrT - 1)
{
    maxTIterReached++;
    Warning
        << "Max iterations reached within energy equation" << endl;
}
