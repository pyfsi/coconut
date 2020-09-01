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
    pimpleDyMFoam.C

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
#include "fixedValuePointPatchField.H"
#include "IOstream.H"
#include "Ostream.H"
#include "forces.H"

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
	
	runTime.run();
    word prev_runTime;

    while (true) // NOT runTime.run()
    {
        usleep(1000); // Expressed in microseconds 


    	if (exists("next.coco"))
		{

    		
			#include "readControls.H"
			#include "CourantNo.H"
			#include "setDeltaT.H"

			 // Get the patch ID: patch name must be "mantle"
            // This should be generalised and a loop should be provided for the different patches of interest

    		prev_runTime = runTime.timeName();
    		Info << prev_runTime << nl << endl;

    		runTime++;
    		remove("next.coco");
    		OFstream outfile ("next_ready.coco");
    		outfile << "Joris says: good job on next.coco" << endl;
			Info << "Time = " << runTime.timeName() << nl << endl; // Might be deleted when linked to CoCoNuT (which already outputs current time step)
		}
    	
    	if (exists("continue.coco"))
		{

    		// Define movement of the coupling interface
    		label patchWallID = mesh.boundaryMesh().findPatchID("mantle");
            const fvPatch& patchWallFaces = mesh.boundary()[patchWallID];

            // Info << "In Next" << nl << endl;

            // *** Set patch displacement for motion solver.*** //
            // Find the reference to the pointDisplacement field (this appears to work)
        	pointVectorField& PointDisplacement = const_cast<pointVectorField&>
            (
                mesh.objectRegistry::lookupObject<pointVectorField >
                (
                "pointDisplacement"
                )
            );

            //OFstream testfile(runTime.path()/"Example_pointDispFile");
            //testfile << PointDisplacement<< endl;
            //PointDisplacement.write();

            // Info << PointDisplacement.instance() << nl << endl; //Instance is a part of the path referring to the file that should be read and is updated (verified this by printing)

            //Get the vector field of the patch
            vectorField &pDisp=refCast<vectorField>(PointDisplacement.boundaryFieldRef()[patchWallID]);

            Info<< "Reading pointDisplacement\n" << endl;

            Info<< "prev_runTime" << prev_runTime << endl;

            pointVectorField pointDisplacement_temp_
            (
                IOobject
                (
                    "pointDisplacement_Next",
                    prev_runTime,
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                pointMesh::New(mesh)
            );

            pointVectorField& PointDisplacementTemp = const_cast<pointVectorField&>
            (
                pointDisplacement_temp_
            );

            vectorField &pDispTemp=refCast<vectorField>(PointDisplacementTemp.boundaryFieldRef()[patchWallID]);

            //Info << pDispTemp <<nl << endl;

            //Info << pointDisplacement_temp_ << nl <<endl;
            //Info << pointDisplacement_temp_.boundaryField()[patchWallID] << nl<< endl;
            //Info << "-------------------------------------------------------------------" << nl<< endl;
            //Info << PointDisplacement.boundaryField()[patchWallID] << nl<< endl;

            // Assign the new boundary displacements
            PointDisplacement.boundaryFieldRef()[patchWallID] ==  pDispTemp;

            //Info << PointDisplacement.boundaryField()[patchWallID] << nl<< endl;
    		
    		
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
            remove("continue.coco");
            // Return the coupling interface output

            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

            runTime.run();
            Info << "I get past the save" << nl << endl;
            OFstream outfile ("continue_ready.coco");
            outfile << "Joris says good job on continue.coco" << endl;
    		
		}
		
    	if (exists("save.coco"))
		{

    		runTime.write(); // OF-command: loops over all objects and requests writing - writing is done based on the specific settings of each variable (AUTO_WRITE, NO_WRITE)
    		remove("save.coco");
    		OFstream outfile ("save_ready.coco");
			outfile << "Joris says: good job on save.coco" << endl;
		}
    	
    	if (exists("stop.coco"))
		{
    		remove("stop.coco");
    		OFstream outfile ("stop_ready.coco"); 
    		outfile << "Joris says: good job on stop.coco" << endl;
    		break;
		}  
    }
    return 0;
}

// ************************************************************************* //
