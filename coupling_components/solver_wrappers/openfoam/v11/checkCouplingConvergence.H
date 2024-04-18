{
    if (solver.pimple.corrCriteriaSatisfied())
    {
        Info<< solver.pimple.algorithmName() << ": Converged after one iteration" << endl;
        sendMessage("solver_converged.coco");
    }
    else
    {
        Info<< solver.pimple.algorithmName() << ": Not converged after one iteration" << endl;
    }
    sendMessage("check_ready.coco");
}