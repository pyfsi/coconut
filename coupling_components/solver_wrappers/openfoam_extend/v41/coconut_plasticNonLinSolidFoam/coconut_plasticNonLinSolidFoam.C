/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    coconut_plasticNonLinSolidFoam

Description
    Finite volume structural solver for elastoplasticity employing an
    incremental strain updated Lagrangian approach (UL).

    Employing the Kirchhoff stress tensor (Cauchy stress tensor scaled to the
    reference volume). Note: when elastic deformations are small and plastic
    deformation is isochoric, the Kirchhoff and Cauchy tensors coincide.

    The stress is calculated by the mechanicalModel and can be any number of
    constitutive relations.

    Derivation of momentum using Nanson's formula to relate deformed area Sf and
    reference area Sf_0
    force =
    = sigmaCauchy & Sf
    = sigmaCauchy & (J*Finv.T() & Sf_0)
    = ((sigmaCauchy*J) & Finv.T()) & Sf_0
    = (tau & Finv.T()) & Sf_0
    = (Finv.T() & Sf_0) & tau

    Note: when using an updated Lagrangian formulation then we must scale by
    relJ instead of J so (a) becomes:
    force = (relFinvf.T() & Sf_u) & tau * (relJ/J)
    where relJ == J for total Lagrangian.

    This solver is modified to use in the numerical FSI coupling tool: CoCoNuT

Author
    Philip Cardiff UCD
    Peter De Jaeger Bekaert
    Michael Clancy UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mechanicalModel.H"
#include "newLeastSquaresVolPointInterpolation.H"
#include "twoDPointCorrector.H"
#include "symmetryPolyPatch.H"
//#include "solidContactFvPatchVectorField.H"
#include "thermalModel.H"
#include "thermalContactFvPatchScalarField.H"
//#include "thermalGeneralContactFvPatchScalarField.H"
#include "nonLinearGeometry.H"
#include "transformGeometricField.H"
#include "processorFvPatchFields.H"
#include "logVolFields.H"
#include "fvcInterpolateP.H"
#include "fvcReconstructNew.H"
#include "dynamicFvMesh.H"
#include "materialGgiFvPatchVectorField.H"
#include "materialCouplingFvPatchVectorField.H"
#include "regionCouplePolyPatch.H"
#include "regionCoupledFvPatchFields.H"
#include "pointGaussLeastSquaresGrad.H"
#include "interpolationTable.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include <stdlib.h>
#include <assert.h>
#include <set>
#include <string>
#include <map>
#include <sstream>
#include <unistd.h>
#include <fstream>
#include <iostream>

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "checkDynamicMeshDict.H"
    #include "createDynamicFvMesh.H"
    #include "createNonLinearGeometry.H"
    #include "createFields.H"
    #include "axisymmetricSolid.H"
    #include "readSolidMechanicsControls.H"
    #include "checkForGlobalFaceZones.H"
    #include "checkForGlobalFaceZones.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    runTime.run();
    word prevRunTime;
    pointField oldPoints_;

    int maxUIterReached = 0;
    int maxTIterReached = 0;
    unsigned int iteration = 0;
    int iCorr = 0;
    lduSolverPerformance solverPerf;
    bool converged = false;
    blockLduMatrix::debug = 0;

    IOdictionary controlDict(IOobject("controlDict",runTime.system(),mesh,IOobject::MUST_READ,IOobject::NO_WRITE));
    wordList boundary_names (controlDict.lookup("boundary_names"));

    while (true)
    {
          scalar relativeResidual = 1.0;
          scalar materialResidual = 0.0;

        usleep(10000);

        if (exists("next.coco"))
        {
            mesh.update();
            mechanical.updateYieldStress();

            // Update surface fields after a topoChange
            mu = mechanical.mu();
            lambda = mechanical.lambda();
            twoMuLambda = 2.0 * mu + lambda;
            twoMuLambdaf = fvc::interpolate(twoMuLambda, "twoMuLambda");
            RhieChowScaleFactor = mechanical.RhieChowScaleFactor();

            #include "setDeltaT.H"
            #include "createNonLinearGeometry.H"

            prevRunTime = runTime.timeName();
            oldPoints_ = mesh.points();
            runTime++;
            //runtime.write();//Navaneeth
            remove("next.coco");
            OFstream outfile ("next_ready.coco");
            outfile << "next.coco" << endl;
            Info<< "Time: " << runTime.timeName() << nl << endl;
            iteration = 0;
//          lduSolverPerformance solverPerf;
            converged = false;
            blockLduMatrix::debug = 0;

            // Optional predictor
            if (predictor)
                {
                    solve(fvm::laplacian(twoMuLambdaf, DU, "laplacian(DDU,DU)"));
                }

        }

        if (exists("continue.coco"))
        {

            pointField &oldPoints = const_cast<pointField&>(oldPoints_);
            mesh.movePoints(oldPoints);

//            Info<< "mesh_points_before" << mesh.points() << nl << endl;

//            curS == curS.oldTime();
//            DU == DU.oldTime();
//            cauchyTraction == cauchyTraction.oldTime();
            mechanical.resetYieldStress();

            gradDU = fvc::grad(DU);

            // Relative deformation gradient
            relF = I + gradDU.T();

            // Inverse relative deformation gradient
            relFinv = hinv(relF);

             // Total deformation gradient
            F = relF & F.oldTime();

            // Relative Jacobian (Jacobian of relative deformation gradient)
            relJ = det(relF);

            // Relative deformation gradient with volumetric strain removed
            relFbar = pow(relJ, -1.0/3.0)*relF;

             // Jacobian of deformation gradient
            J = det(F);

            // Update tau using material model
            mechanical.correct(tau);

            rho == rho.oldTime();
            U == U.oldTime();

            iteration++;
            Info<< "Coupling iteration = " << iteration << nl << endl;

            iCorr = 0;

            // Momentum loop
            do
                {
                // Store previous iteration of DU to allow under-relaxation and
                // calculation of a relative residual
                DU.storePrevIter();

                if (regionCoupled)
                {
                    // Attach region coupled patches
                    #include "attachPatches.H"
                }

                // Update Cauchy traction vectors
                #include "calcCauchyTraction.H"

                // Discretise the linear momentum equation using an updated
                // Lagrangian approach

                fvVectorMatrix DUEqn
                (
                    fvm::d2dt2(rho, DU)
                  + fvc::d2dt2(rho.oldTime(), U.oldTime())
                  ==fvm::laplacian(twoMuLambdaf, DU, "laplacian(DDU,DU)")
                  + fvc::div(cauchyTraction*mag(curS))
                  - fvc::laplacian(twoMuLambdaf, DU, "laplacian(DDU,DU)")
                  + RhieChowScaleFactor
                   *(
                        fvc::laplacian(twoMuLambdaf, DU, "laplacian(DDU,DU)")
                      - fvc::div(twoMuLambdaf*(mesh.Sf() & gradDUf))
                     )

                );

                // to stop the RVE 'drifting' off we fix the first cell
                if (representativeVolumeElementCase)
                {
                    DUEqn.setReference(0,vector::zero);
                }

                // Under-relax the linear system
                DUEqn.relax();

                // Solve the linear system
                solverPerf = DUEqn.solve();

                // Under-relax the DU field
                DU.relax();

                // Detach region coupled patches
                if (regionCoupled)
                {
                    #include "detachPatches.H"
                }

                // Interpolate DU from cell-centres to points
                //if (updatePointFields)
                if (materialGgi)
                {
                    volToPointInterp.interpolate(DU, pointDU);
                }

                // Calculate gradient of DU at cell-centres using pointDU field
                //gradDU = fvc::grad(DU, pointDU);
                gradDU = fvc::grad(DU);

                // Correct gradDU on materialGgi patches
                #include "correctDispGrad.H"

                // Relative deformation gradient
                relF = I + gradDU.T();

                // Add average total displacement contribution for RVE-style cases
                if (representativeVolumeElementCase)
                {
                    // Lookup new and old specified total average deformation gradients
                    const tensor avgF = avgFSeries(runTime.timeOutputValue());
                    const tensor avgFold =
                        avgFSeries(runTime.timeOutputValue() - runTime.deltaTValue());



                    // Add relative deformation contribution
                    relF += (avgF & inv(avgFold)) - I;
                }

                // Inverse relative deformation gradient
                relFinv = hinv(relF);

                // Total deformation gradient
                F = relF & F.oldTime();

                // Relative Jacobian (Jacobian of relative deformation gradient)
                if (stabilisePressure)
                {
                    // Reconstruct face gradients to relJ: this will avoid
                    // pressure oscillations
                    // This is just det(I + gradDU.T()) where the gradDU has
                    // been averaged from the face gradDU field
                    relJ =
                        det
                        (
                            I
                            + fvc::reconstructNew
                            (
                                mesh.magSf()*fvc::snGrad(DU)
                            )().T()
                        );
                }
                else
                {
                    relJ = det(relF);
                }

                if
                (
                    gMin(relJ.internalField()) < 0.01
                 || gMax(relJ.internalField()) > 100.0
                )
                {
                    forAll(relJ, cellI)
                    {
                        if (relJ[cellI] < SMALL)
                        {
                            Pout<< "Cell " << cellI
                                << " with centre " << mesh.C()[cellI]
                                << " has a become inverted!"
                                << " relJ: " << relJ[cellI] << endl;
                        }
                        else if (relJ[cellI] > 100.0)
                        {
                            Pout<< "Cell " << cellI
                                << " with centre " << mesh.C()[cellI]
                                << " has increased in volume by 100 times"
                                << endl;
                        }
                    }

                    FatalError
                        << "Cells have become inverted! see details above."
                        << abort(FatalError);
                }

                // Relative deformation gradient with volumetric strain removed
                relFbar = pow(relJ, -1.0/3.0)*relF;

                // Jacobian of deformation gradient
                J = det(F);

                // Store coefficient field as it may be used by the mechanical law
                // to calculate the pressure
                const volScalarField DUEqnA("DUEqnA", DUEqn.A());

                // Update tau using material model
                mechanical.correct(tau);

                // Calculate residuals and check convergence
                    #include "checkConvergence.H"
            }
            while (!converged && ++iCorr < nCorr);

            // Update density
            rho = rho.oldTime()/relJ;

            // Cauchy stress
            sigmaCauchy = tau/J;

            // Total displacement at cell-centres
            gradU = fvc::grad(U.oldTime() + DU);
            U = U.oldTime() + DU;
//            U.component(2)*= 0.0;
//            Info << 'Z_direction' << U << endl;

            V = DU/runTime.deltaT();
            if (thermalStress && runTime.value() > thermalStressStartTime)
            {
                #include "TEqn.H"
            }

            // Mesh update
            // Note: after remeshing, geometric surface field may be incorrect on
            // the new faces so we will re-created all surface fields here
            //mesh.update();

            // Let the mechanicalLaw know that we have reached the end of the
            // time-step so it may increment its accumulated values
            // Note: we call this after mesh.update() in case it writes a field to
            // disk
//            mechanical.updateYieldStress();
//
//            // Update surface fields after a topoChange
//            mu = mechanical.mu();
//            lambda = mechanical.lambda();
//            twoMuLambda = 2.0 * mu + lambda;
//            twoMuLambdaf = fvc::interpolate(twoMuLambda, "twoMuLambda");
//            RhieChowScaleFactor = mechanical.RhieChowScaleFactor();

            remove("continue.coco");
            //Return the coupling interface output

            Info<< "ExecutionTime = "<< runTime.elapsedCpuTime() << " s"
                << " clockTime = "<< runTime.elapsedClockTime() << " s"
                << nl << endl;

            runTime.run();
            runTime.write(); //debugging
            Info << "Coupling iteration " << iteration << " end" << nl << endl;
            OFstream outfile ("continue_ready.coco");

        }

        if (maxTIterReached > 0)
        {
            Warning
                << "Max iterations for the energy equation were reached in "
                << maxTIterReached << " time-steps" << endl;
        }

        if (exists("save.coco"))
        {
            #include "writeFields.H"
            //runTime.write(); // OF-command: loops over all objects and requests writing - writing is done based on the specific settings of each variable (AUTO_WRITE, NO_WRITE)
            remove("save.coco");
            OFstream outfile ("save_ready.coco");
        }

        if (exists("stop.coco"))
        {
            remove("stop.coco"); // should not be uncommented
            runTime.stopAt(Time::stopAtControls::saNoWriteNow);
            OFstream outfile ("stop_ready.coco");
            break;
        }
    }

    Info<< "\nEnd\n" << endl;

    return(0);
}


// ************************************************************************* //
