// Add this to source code to force a write of the face centres when the solver is started


using namespace Foam;
Info<< "Writing face centres of wall patches\n" << endl;
//change the patch name to your boundary name

//const objectRegistry& db = mesh.thisDb();
//const IOdictionary& controlDict = db.lookupObject<IOdictionary>("controlDict");

//
//const objectRegistry& db = mesh.thisDb();
//const IOdictionary& controlDict = db.lookupObject<IOdictionary>("controlDict");

//IOobject controlDict_IO = IOobject("controlDict", runTime.system(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);
//IOdictionary controlDict(controlDict_IO);

wordList interface(runTime.controlDict().lookup("interface"));

Info << interface << endl;

forAll(interface,i)
{
    word patchLabel = interface[i];
    Info << patchLabel << endl;
    forAll(mesh.boundaryMesh()[patchLabel].faceCentres(), faceI)
    {
        scalar x = mesh.boundaryMesh()[patchLabel].faceCentres()[faceI].x();
        scalar y = mesh.boundaryMesh()[patchLabel].faceCentres()[faceI].y();
        scalar z = mesh.boundaryMesh()[patchLabel].faceCentres()[faceI].z();
        Info<<faceI<<" "<<x<<" "<<y<<" "<<z<<" "<<endl;
    }
}

//
//IOobject controlDict_IO = IOobject("controlDict", runTime.system(),mesh,IOobject::MUST_READ,IOobject::AUTO_WRITE);
//IOdictionary controlDict(controlDict_IO);

//word patch = controlDict.lookupOrDefault("interface", "empty");


/*forAll(mesh.boundaryMesh()[patchI].faceCentres(), faceI)
{
    scalar x = mesh.boundaryMesh()[patchI].faceCentres()[faceI].x();
    scalar y = mesh.boundaryMesh()[patchI].faceCentres()[faceI].y();
    scalar z = mesh.boundaryMesh()[patchI].faceCentres()[faceI].z();
    Info<<faceI<<" "<<x<<" "<<y<<" "<<z<<" "<<endl;
}
*/