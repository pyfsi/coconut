{
    if (pimple.finalInnerIter()) // should be pimple.criteriaSatisfied but its protected TODO make inhrtd class ??
    {
        Info<< "PIMPLE: Converged after one iteration" << endl;
        sendMessage("solver_converged.coco");
    }
    else
    {
        Info<< "PIMPLE: Not converged after one iteration" << endl;
    }
    sendMessage("check_ready.coco");
}
