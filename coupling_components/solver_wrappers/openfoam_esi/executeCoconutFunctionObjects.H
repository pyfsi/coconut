forAll(coconutFunctionObjects, i)
{
    word functionObject(coconutFunctionObjects[i]);
    label oi(runTime.functionObjects().findObjectID(functionObject));
    runTime.functionObjects()[oi].execute();
    runTime.functionObjects()[oi].write();
}
