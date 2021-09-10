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
    pimpleFoam

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
#include "fixedValuePointPatchField.H"
#include "IOstream.H"
#include "Ostream.H"
#include "forces.H"
#include "fsiDisplacement.H"
#include "Time.H"

#include <stdlib.h>
#include <assert.h>
#include <set>
#include <string>
#include <map>
#include <sstream>
#include <unistd.h>
#include <fstream>
#include <iostream>

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

	runTime.run();
    word prev_runTime;

    unsigned int iteration;
    iteration = 0;

    IOdictionary controlDict(IOobject("controlDict",runTime.system(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE));
    wordList boundary_names (controlDict.lookup("boundary_names"));

    while (true) // NOT (pimple.run(runTime))
    {
        usleep(10000); // Expressed in microseconds

        if (exists("next.coco"))
        {
            #include "readDyMControls.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"

            prev_runTime = runTime.timeName();
            runTime++;
            remove("next.coco");
            OFstream outfile ("next_ready.coco");
            outfile << "next.coco" << endl;
            Info << "Time = " << runTime.timeName() << nl << endl; // Might be deleted when linked to CoCoNuT (which already outputs current time step)
            iteration = 0;
        }

        if (exists("continue.coco"))
        {
            iteration++;
            Info << "Coupling iteration = " << iteration << nl << endl;

            // Define movement of the coupling interface
            forAll(boundary_names, s)
            {
                    word boundary_name = boundary_names[s];
                    ApplyFSIPointDisplacement(mesh, boundary_name, prev_runTime);
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

                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        if (correctPhi)
                        {
                            #include "correctPhi.H"
                        }

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

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
            }
            remove("continue.coco");
            // Return the coupling interface output

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            runTime.run();
            Info << "Coupling iteration " << iteration << " end" << nl << endl;
            OFstream outfile ("continue_ready.coco");
        }

        if (exists("save.coco"))
        {
            runTime.write(); // OF-command: loops over all objects and requests writing - writing is done based on the specific settings of each variable (AUTO_WRITE, NO_WRITE)
            remove("save.coco");
            OFstream outfile ("save_ready.coco");
        }

        if (exists("stop.coco"))
        {
            // remove("stop.coco"); // should not be uncommented
            runTime.stopAt(Time::stopAtControl::noWriteNow);
            OFstream outfile ("stop_ready.coco");
            break;
        }

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
