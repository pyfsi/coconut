/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    coconut_interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.
	- Adapted for coupling with the pyKratos/CoCoNuT-framework
	  Developed at Ghent University by Laurent De Moerloose (2019)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "Time.H"
#include "fsiDisplacement.H"
#include "waitForSync.H"

#include <unistd.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createDyMControls.H"
    #include "createRDeltaT.H"
    #include "createFields.H"
    #include "createFvOptions.H"

    volScalarField rAU
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1.0)
    );

    #include "CorrectPhi.H"
    #include "createUf.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    unsigned int iteration;
    iteration = 0;

    IOdictionary controlDict(IOobject("controlDict", runTime.system(), mesh, IOobject::MUST_READ,IOobject::NO_WRITE));
    wordList boundaryNames (controlDict.lookup("boundaryNames"));
    wordList coconutFunctionObjects (controlDict.lookup("coconutFunctionObjects"));

    while (true)
    {
        usleep(1000); // Expressed in microseconds

        if (exists("next.coco"))
        {
            waitForSync("next"); // Keep the sync at the beginning of the block

            runTime.run();

            #include "readControls.H"

            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"

            runTime++;
            iteration = 0;

            Info << "Time = " << runTime.timeName() << nl << endl;
        }

        if (exists("continue.coco"))
        {
            waitForSync("continue");

            iteration++;
            Info << "Coupling iteration = " << iteration << nl << endl;

            // Define movement of the coupling interface
            forAll(boundaryNames, s)
            {
                word boundaryName = boundaryNames[s];
                ApplyFSIPointDisplacement(mesh, boundaryName);
            }

            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
               if (pimple.firstIter() || moveMeshOuterCorrectors)
                {

                    scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();
                    mesh.update();

                    Info<< "Execution time for mesh.update() = "
                     << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                     << " s" << endl;

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    // Calculate absolute flux from the mapped surface velocity
                    phi = mesh.Sf() & Uf;

                    if (mesh.changing() && correctPhi)
                    {

                         #include "correctPhi.H"
                    }

                    // Make the flux relative to the mesh motion
                    fvc::makeRelative(phi, U);

                    mixture.correct();

                    if (mesh.changing() && checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }

                #include "alphaControls.H"
                #include "alphaEqnSubCycle.H"

                mixture.correct();

                #include "UEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            // Return the coupling interface output
            #include "executeCoconutFunctionObjects.H"

            Info << "Coupling iteration " << iteration << " end" << nl << endl;
        }

        if (exists("save.coco"))
        {
            waitForSync("save");

            runTime.write(); // OF-command: loops over all objects and requests writing
            // writing is done based on the specific settings of each variable (AUTO_WRITE, NO_WRITE)
        }

        if (exists("stop.coco"))
        {
            waitForSync("stop");
            break;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
