{
    if (pimple.corrCriteriaSatisfied())
    {
        Info<< pimple.algorithmName() << ": Converged after one iteration" << endl;
        sendMessage("solver_converged.coco");
    }
    else
    {
        Info<< pimple.algorithmName() << ": Not converged after one iteration" << endl;
    }
    sendMessage("check_ready.coco");
}