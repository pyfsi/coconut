IOdictionary controlDict
(
    IOobject
    (
        "controlDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    )
);

wordList boundaryNames(controlDict.lookup("boundaryNames"));

wordList coconutFunctionObjects(controlDict.lookup("coconutFunctionObjects"));

bool checkCouplingConvergence(controlDict.lookupOrDefault("checkCouplingConvergence", false));
