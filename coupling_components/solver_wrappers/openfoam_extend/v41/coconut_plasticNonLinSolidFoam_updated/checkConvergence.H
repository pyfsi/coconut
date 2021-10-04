#include "calculateRelativeResidual.H"
materialResidual = mechanical.residual();

// If one of the residuals has converged to an order of magnitude
// less than the tolerance then consider the solution converged
if
(
    materialResidual < convergenceTolMaterial
 && iCorr > 1 // force at least 2 outer iteration
)
{
    if
    (
        solverPerf.initialResidual() < convergenceToleranceU
        && relativeResidual < convergenceTolRelU
    )
    {
        Info<< "    Both residuals have converged" << endl;
        converged = true;
    }
    else if (relativeResidual < 0.1*convergenceTolRelU)
    {
        Info<< "\tThe relative residual has converged" << endl;
        converged = true;
    }
    else if (solverPerf.initialResidual() < 0.1*convergenceToleranceU)
    {
        Info<< "\tThe solver residual has converged" << endl;
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

// Print out residuals
if (iCorr % infoFrequency == 0 || converged)
{
    Info<< "    Time " << runTime.value()
        << ", Corr " << iCorr
        << ", res = " << solverPerf.initialResidual()
        << ", rel res = " << relativeResidual
        << ", material res = " << materialResidual
        << ", iters = " << solverPerf.nIterations() << endl;
}

if (iCorr == nCorr - 1)
{
    maxUIterReached++;
    Warning
        << "Max iterations reached within momentum equation" << endl;

    // Add flag to the solidControl dict so that the contact boundary conditions
    // can check if the momentum looop did not converge
    solidControl.set("maxUIterReached", true);
}
else if (converged)
{
    solidControl.set("maxUIterReached", false);
}

// If J becomes very small or large, it suggests that the solver is diverging
if (gMin(J) < 0.01 || gMax(J) > 100.0)
{
    int invertedCells = 0;
    forAll(J, cellI)
    {
        if (J[cellI] < 0.01 || J[cellI] > 100.0)
        {
            invertedCells++;
        }
    }

    // Write fields for debugging
    runTime.write();

    FatalErrorIn("checkConvergence.H")
        << "min(J) is less than 0.01 or max(J) greater than 100 in "
        << invertedCells << " cells: this means the solver is diverging!" << nl
        << "The latest fields have been written to time = " << runTime.value()
        << " for debugging."
        << abort(FatalError);
}
