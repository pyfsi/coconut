/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    coconut_pimpleFoam

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids,
    with optional mesh motion and mesh topology changes.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
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

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    unsigned int iteration;
    iteration = 0;

    #include "readCoconutControls.H"

    while (true)
    {
        usleep(1000); // Expressed in microseconds

        if (exists("next.coco"))
        {
            waitForSync("next"); // Keep the sync at the beginning of the block

            runTime.run();

            #include "readDyMControls.H"
            #include "CourantNo.H"
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
                if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
                {
                    // Calculate the mesh motion and update the mesh
                    mesh.update();

                    if (mesh.changing())
                    {
                        MRF.update();

                        if (correctPhi)
                        {
                            // Calculate absolute flux
                            // from the mapped surface velocity
                            phi = mesh.Sf() & Uf();

                            #include "correctPhi.H"

                            // Make the flux relative to the mesh motion
                            fvc::makeRelative(phi, U);
                        }

                        if (checkMeshCourantNo)
                        {
                            #include "meshCourantNo.H"
                        }
                    }
                }

                #include "UEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    laminarTransport.correct();
                    turbulence->correct();
                }

                // Check coupling convergence
                if (checkCouplingConvergence && pimple.firstPimpleIter())
                {
                    #include "checkCouplingConvergence.H"
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
