volVectorField DU_oldTime
(
    IOobject
    (
        "DU_oldTime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("DU_oldTime", dimLength, vector::zero)

);

// Displacement increment
Info<< "Reading field DU\n" << endl;
volVectorField DU
(
    IOobject
    (
        "DU",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Store old time: needed for some boundary conditions
DU.oldTime().oldTime();

// Gradient of displacement increment
volTensorField gradDU
(
    IOobject
    (
        "grad(DU)",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("zero", dimless, tensor::zero)
);

// Total displacement
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimLength, vector::zero)
);

// Store old time: needed for some boundary conditions
U.oldTime().oldTime();

volTensorField gradU
(
    IOobject
    (
        "grad(U)",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("zero", dimless, tensor::zero)
);

// Velocity
volVectorField V
(
    IOobject
    (
        "Velocity",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
//    fvc::ddt(U)
    DU/runTime.deltaT()
);

// Total deformation gradient
volTensorField F
(
    IOobject
    (
        "F",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("zero", dimless, I)
);

// Jacobian of total deformation gradient
volScalarField J
(
    IOobject
    (
        "J",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    det(F)
);

// Relative deformation gradient
volTensorField relF
(
    IOobject
    (
        "relF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    I + gradDU.T()
);

// Inverse of relative deformation gradient
volTensorField relFinv
(
    IOobject
    (
        "relFinv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    hinv(relF)
);

// Relative Jacobian - the Jacobian of the relative deformation gradient
volScalarField relJ
(
    IOobject
    (
        "relJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    det(relF)
);

// Volume preserving part of relative deformation gradient
volTensorField relFbar
(
    IOobject
    (
        "relFbar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    pow(relJ, -1.0/3.0)*relF
);

// Cauchy stress tensor vol field
volSymmTensorField sigmaCauchy
(
    IOobject
    (
        "sigmaCauchy",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
);

// Kirchhoff stress tensor vol field
volSymmTensorField tau
(
    IOobject
    (
        "tauKirchhoff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
);

// True (Hencky) strain
// volSymmTensorField epsilon
// (
//     IOobject
//     (
//         "epsilon",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedSymmTensor("zero", dimless, symmTensor::zero)
// );

// Create mechanical model
mechanicalModel mechanical(mesh);

// Density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mechanical.rho()
);

volScalarField rho_oldFSIiteration
(
    IOobject
    (
        "rho_oldFSIiteration",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mechanical.rho()
//    dimensionedScalar("rho_oldFSIiteration", dimLength, 0.0)
);
#include "correctFieldsOnRestart.H"

volScalarField mu = mechanical.mu();
volScalarField lambda = mechanical.lambda();
volScalarField RhieChowScaleFactor = mechanical.RhieChowScaleFactor();


bool regionCoupled = false;
wordList twoMuLambdaTypes
(
    mesh.boundary().size(),
    calculatedFvPatchVectorField::typeName
);
forAll(mesh.boundaryMesh(), patchI)
{
    if (isA<regionCouplePolyPatch>(mesh.boundaryMesh()[patchI]))
    {
        twoMuLambdaTypes[patchI] =
            regionCoupledFvPatchScalarField::typeName;
        regionCoupled = true;
    }
}

if (regionCoupled)
{
    // Attach region coupled patches
    #include "attachPatches.H"
}

volScalarField twoMuLambda
(
    IOobject
    (
        "twoMuLambda",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    2.0*mu + lambda,
    twoMuLambdaTypes
);

surfaceScalarField twoMuLambdaf
(
    fvc::interpolate(twoMuLambda, "twoMuLambda")
);


// Detach region coupled patches
if (regionCoupled)
{
    #include "detachPatches.H"
}

// Correct material property fields on bi-material GGI interfaces
bool materialGgi = false;
forAll(DU.boundaryField(), patchI)
{
    if
    (
        isA<materialGgiFvPatchVectorField>
        (
            DU.boundaryField()[patchI]
        )
    )
    {
        mu.boundaryField()[patchI] =
            mu.boundaryField()[patchI].patchInternalField();
        //muf.boundaryField()[patchI] = mu.boundaryField()[patchI];

        lambda.boundaryField()[patchI] =
            lambda.boundaryField()[patchI].patchInternalField();
        //lambdaf.boundaryField()[patchI] = lambda.boundaryField()[patchI];

        // K.boundaryField()[patchI] =
        //     K.boundaryField()[patchI].patchInternalField();
        // Kf.boundaryField()[patchI] = K.boundaryField()[patchI];

        materialGgi = true;
    }
}

if (materialGgi)
{
    if
    (
        word(mesh.schemesDict().gradScheme(gradDU.name()))
     != fv::pointGaussLeastSquaresGrad<vector>::typeName
    )
    {
        FatalError
            << "For cases with internal material interfaces, the "
            << gradDU.name() << " gradScheme must be set to "
            << fv::pointGaussLeastSquaresGrad<vector>::typeName
            << abort(FatalError);
    }
}

#include "createThermalFields.H"

// Lookup pointDU from the dynamicFvMesh
pointVectorField& pointDU =
    const_cast<pointVectorField&>
    (
        mesh.lookupObject<pointVectorField>
        (
            "pointDU"
         )
    );

// Lookup volToPoint interpolator from the dynamicMesh
const newLeastSquaresVolPointInterpolation& volToPointInterp =
    mesh.lookupObject<newLeastSquaresVolPointInterpolation>
    (
        "newLeastSquaresVolPointInterpolation"
    );

// Current area vectors surface field
surfaceVectorField curS = (mesh.Sf() & fvc::interpolate(relFinv*relJ));

// Current Cauchy traction vectors surface field
surfaceVectorField cauchyTraction
(
    "cauchyTraction",
    (curS & fvc::interpolate(tau/J))/mag(curS)
);

// Gradient of displacement increment surface field
surfaceTensorField gradDUf = fvc::interpolate(gradDU);
