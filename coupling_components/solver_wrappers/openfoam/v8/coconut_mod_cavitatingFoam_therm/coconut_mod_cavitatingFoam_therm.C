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
    coconut_cavitatingFoam

Description
    Transient cavitation code based on the homogeneous equilibrium model
    from which the compressibility of the liquid/vapour "mixture" is obtained,
    with optional mesh motion and mesh topology changes.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"

#include "mixtureCompressibilityModel.H"
#include "EOSModel.H"
#include "heatCapacityModel.H"
#include "thermalConductivityModel.H"
#include "compressibleTwoPhaseMixture.H"

#include "kinematicMomentumTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"

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
    #include "createControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "createPcorrTypes.H"

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    unsigned int iteration;
    iteration = 0;

    IOdictionary controlDict(IOobject("controlDict", runTime.system(), mesh, IOobject::MUST_READ,IOobject::NO_WRITE));
    wordList boundary_names (controlDict.lookup("boundary_names"));

    autoPtr<volScalarField> divRhoU;

    while (true) // NOT (pimple.run(runTime))
    {
        usleep(1000); // Expressed in microseconds

        if (exists("next.coco"))
        {
            waitForSync("next"); // Keep the sync at the beginning of the block

            #include "readControls.H"

            #include "CourantNo.H"
            #include "setDeltaT.H"

            // Store divRhoU from the previous mesh so that it can be mapped
            // and used in correctPhi to ensure the corrected phi has the
            // same divergence
            if (correctPhi)
            {
                divRhoU = new volScalarField
                (
                    "divRhoU",
                    fvc::div(fvc::absolute(phi, U), rho)
                );
            }

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
            forAll(boundary_names, s)
            {
                word boundary_name = boundary_names[s];
                ApplyFSIPointDisplacement(mesh, boundary_name);
            }

            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
                {
                    // Do any mesh changes
                    mesh.update();

                    if (mesh.changing())
                    {
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

                if
                (
                    !mesh.steady()  // steady operation has not been tested
                 && pimple.firstPimpleIter()
                )
                {
                     #include "rhoEqn.H"
                     #include "alphavPsi.H"
                }

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

                if (thermalSimulation)
                {
                    #include "TEqn.H"
                }

                mixture.correct();
                mu = rho*turbulence->nuEff();
            }

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            // Return the coupling interface output
            runTime.run();

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
