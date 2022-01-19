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
    coconut_pimpleFoam

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
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
    #include "createControls.H"
    #include "createFields.H"
    #include "createUf.H"
    #include "createFvOptions.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    word prev_runTime;

    unsigned int iteration;
    iteration = 0;

    IOdictionary controlDict(IOobject("controlDict", runTime.system(), mesh, IOobject::MUST_READ,IOobject::NO_WRITE));
    wordList boundary_names (controlDict.lookup("boundary_names"));

    while (true) // NOT runTime.run()
    {
        usleep(1000); // Expressed in microseconds

        if (exists("next.coco"))
        {
            #include "readControls.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"

            prev_runTime = runTime.timeName();
            // For adjustable time steps insert a while (controlDict.deltaT > runTime.deltaTValue()):  pimple loop

            runTime++;
            iteration = 0;

            waitForSync("next");

            Info << "Time = " << runTime.timeName() << nl << endl;
        }

        if (exists("continue.coco"))
        {
            iteration++;
            Info << "Coupling iteration = " << iteration << nl << endl;

            // Define movement of the coupling interface
            forAll(boundary_names, s)
            {
                word boundary_name = boundary_names[s];
                ApplyFSIPointDisplacement(mesh, boundary_name);
            }

            // Calculate the mesh motion and update the mesh
            mesh.update();

            // Calculate absolute flux from the mapped surface velocity
            phi = mesh.Sf() & Uf;

            if (mesh.changing() && correctPhi)
            {
                #include "correctPhi.H"
            }

            // Make the flux relative to the mesh motion
            fvc::makeRelative(phi, U);

            if (mesh.changing() && checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }

            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
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
            }

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            // Return the coupling interface output
            runTime.run();

            waitForSync("continue");

            Info << "Coupling iteration " << iteration << " end" << nl << endl;
        }

        if (exists("save.coco"))
        {
            runTime.write(); // OF-command: loops over all objects and requests writing
            // writing is done based on the specific settings of each variable (AUTO_WRITE, NO_WRITE)

            waitForSync("save");
        }

        if (exists("stop.coco"))
        {
            // remove("stop.coco"); should not be uncommented (has to be seen by all processors)
            runTime.stopAt(Time::saNoWriteNow);
            OFstream outfile ("stop_ready.coco");
            break;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
