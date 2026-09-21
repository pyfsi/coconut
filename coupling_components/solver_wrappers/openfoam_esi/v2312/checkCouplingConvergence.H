{
    // check converged_ variable
    if (pimple.finalIter())
    {
        Info << "from checkCouplingConvergence.h: PIMPLE Converged in 3 iterations" << endl;
        sendMessage("solver_converged.coco");
    }
    else
    {
        Info<< "from checkCouplingConvergence.h: PIMPLE Did not converge in 3 iterations" << endl;
    }
    sendMessage("check_ready.coco");
}
