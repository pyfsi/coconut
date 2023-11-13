/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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

InClass
    standardPenalty

\*---------------------------------------------------------------------------*/

#include "standardPenalty.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(standardPenalty, 0);
  addToRunTimeSelectionTable(normalContactModel, standardPenalty, dictionary);


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //


void standardPenalty::calcPenaltyFactor() const
{
    if (penaltyFactor_ > -SMALL)
    {
        FatalErrorIn("void standardPenalty::calcPenaltyFactor() const")
            << "value already set!" << abort(FatalError);
    }

    // Set penalty factor

    // Approximate penaltyFactor from the mechanical properties
    // This can then be scaled using the penaltyScale

    const label masterPatchIndex =  masterPatchID();
    const label slavePatchIndex =  slavePatchID();

    // Lookup stiffness fields
    const volScalarField& mu = mesh_.lookupObject<volScalarField>("mu");
    const volScalarField& lambda = mesh_.lookupObject<volScalarField>("lambda");

    // Lookup slave side values
    const scalar slaveK =
        gAverage
        (
            lambda.boundaryField()[slavePatchIndex]
          + (2.0/3.0)*mu.boundaryField()[slavePatchIndex]
        );

    const scalar slaveMagSf =
        gAverage(mesh_.magSf().boundaryField()[slavePatchIndex]);

    const volScalarField::DimensionedInternalField& V = mesh_.V();
    scalarField slaveV(mesh_.boundary()[slavePatchIndex].size(), 0.0);
    {
        const unallocLabelList& faceCells =
            mesh_.boundary()[slavePatchIndex].faceCells();
        forAll(mesh_.boundary()[slavePatchIndex], facei)
        {
            slaveV[facei] = V[faceCells[facei]];
        }
    }

    // Lookup master side values
    if (masterPatchIndex != -1)
    {
        const scalar masterK =
            gAverage
            (
                lambda.boundaryField()[masterPatchIndex]
              + (2.0/3.0)*mu.boundaryField()[masterPatchIndex]
            );

        const scalar masterMagSf =
            gAverage(mesh_.magSf().boundaryField()[masterPatchIndex]);

        scalarField masterV(mesh_.boundary()[masterPatchIndex].size(), 0.0);
        {
            const unallocLabelList& faceCells =
                mesh_.boundary()[masterPatchIndex].faceCells();
            forAll(mesh_.boundary()[masterPatchIndex], facei)
            {
                masterV[facei] = V[faceCells[facei]];
            }
        }

        // Avarage contact patch bulk modulus
        const scalar bulkModulus = 0.5*(masterK + slaveK);

        // Average contact patch face area
        const scalar faceArea = 0.5*(masterMagSf + slaveMagSf);

        // Average contact patch cell volume
        const scalar cellVolume = 0.5*(gAverage(masterV) + gAverage(slaveV));

        // Approximate penalty factor based on Hallquist et al.
        // we approximate penalty factor for traction instead of force
        penaltyFactor_ = penaltyScale_*bulkModulus*faceArea/cellVolume;
    }
    else
    {
        // Avarage contact patch bulk modulus
        const scalar bulkModulus = slaveK;

        // Average contact patch face area
        const scalar faceArea = slaveMagSf;

        // Average contact patch cell volume
        const scalar cellVolume = gAverage(slaveV);

        // Approximate penalty factor based on Hallquist et al.
        // we approximate penalty factor for traction instead of force
        penaltyFactor_ = penaltyScale_*bulkModulus*faceArea/cellVolume;
    }

    Info<< "    normal penalty factor: " << penaltyFactor_ << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardPenalty::standardPenalty
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const standAlonePatch& masterFaceZonePatch,
    const standAlonePatch& slaveFaceZonePatch
)
:
    normalContactModel
    (
        name,
        patch,
        dict,
        masterPatchID,
        slavePatchID,
        masterFaceZonePatch,
        slaveFaceZonePatch
    ),
    normalContactModelDict_(dict.subDict(name + "NormalModelDict")),
    mesh_(patch.boundaryMesh().mesh()),
    slavePressureVolField_
    (
        IOobject
        (
            "slavePressure_" + mesh_.boundaryMesh()[slavePatchID].name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimPressure, vector::zero)
    ),
    areaInContactVolField_
    (
        IOobject
        (
            "areaInContact_" + mesh_.boundaryMesh()[slavePatchID].name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimArea, 0.0)
    ),
    penaltyFactor_(-1),
    penaltyScale_
    (
        normalContactModelDict_.lookupOrDefault<scalar>("penaltyScale", 1.0)
    ),
    relaxFac_
    (
        normalContactModelDict_.lookupOrDefault<scalar>
        (
            "relaxationFactor", 0.02
        )
    ),
    adjustPenaltyFactor_
    (
        normalContactModelDict_.lookupOrDefault<Switch>
        (
            "adjustPenaltyScale", false
        )
    ),
    adjustPenaltyFactorDict_
    (
        adjustPenaltyFactor_
      ? normalContactModelDict_.subDict("adjustPenaltyScaleDict")
      : dictionary()
    ),
    adjustPenaltyFactorMethod_
    (
        adjustPenaltyFactorDict_.lookupOrDefault<word>("method", "average")
    ),
    maxPenaltyScale_
    (
        adjustPenaltyFactor_
      ? readScalar(adjustPenaltyFactorDict_.lookup("maxPenaltyScale"))
      : 1.0
    ),
    minPenaltyScale_
    (
        adjustPenaltyFactor_
      ? readScalar(adjustPenaltyFactorDict_.lookup("minPenaltyScale"))
      : 1.0
    ),
    targetAvPen_
    (
        adjustPenaltyFactor_
      ? readScalar(adjustPenaltyFactorDict_.lookup("targetAveragePenetration"))
      : 1.0
    ),
    averagePenetration_(0),
    minPenetration_(0),
    checkSolverConvergence_
    (
        normalContactModelDict_.lookupOrDefault<Switch>
        (
            "checkSolverConvergence", Switch(false)
        )
    ),
    pressureMethod_
    (
        normalContactModelDict_.lookupOrDefault<word>
        (
            "pressureMethod", "gaussian"
        )
    ),
    epsilon0_
    (
        normalContactModelDict_.lookupOrDefault<scalar>("epsilon0", 1e-06)
    ),
    p0_
    (
        normalContactModelDict_.lookupOrDefault<scalar>("p0", 1e9)
    ),
    beta_(Foam::log(100.0)/(epsilon0_ + SMALL)),
    contactIterNum_(0)
{
    Info<< "        pressureMethod: " << pressureMethod_ << nl
        << "            epsilon0: " << epsilon0_ << nl
        << "            p0: " << p0_ << nl
        << "        adjustPenaltyFactor: " << adjustPenaltyFactor_ << endl;

    // If gaussian pressure is enabled then adjustPenaltyFactor should also be
    // enabled; this is because the target average penetration is used to scale
    // the pressure
    if (pressureMethod_ == "gaussian" && !adjustPenaltyFactor_)
    {
        FatalErrorIn
        (
            "standardPenalty::standardPenalty\n"
            "(\n"
            "    const word& name,\n"
            "    const fvPatch& patch,\n"
            "    const dictionary& dict,\n"
            "    const label masterPatchID,\n"
            "    const label slavePatchID,\n"
            "    const standAlonePatch& masterFaceZonePatch,\n"
            "    const standAlonePatch& slaveFaceZonePatch\n"
            ")"
        )   << "When the pressureMethod is gaussian, "
            << "then adjustPenaltyFactor should also be enabled!"
            << abort(FatalError);
    }
}


standardPenalty::standardPenalty(const standardPenalty& nm)
:
    normalContactModel(nm),
    normalContactModelDict_(nm.normalContactModelDict_),
    mesh_(nm.mesh_),
    slavePressureVolField_(nm.slavePressureVolField_),
    areaInContactVolField_(nm.areaInContactVolField_),
    penaltyFactor_(nm.penaltyFactor_),
    penaltyScale_(nm.penaltyScale_),
    relaxFac_(nm.relaxFac_),
    adjustPenaltyFactor_(nm.adjustPenaltyFactor_),
    adjustPenaltyFactorDict_(nm.adjustPenaltyFactorDict_),
    adjustPenaltyFactorMethod_(nm.adjustPenaltyFactorMethod_),
    maxPenaltyScale_(nm.maxPenaltyScale_),
    minPenaltyScale_(nm.minPenaltyScale_),
    targetAvPen_(nm.targetAvPen_),
    averagePenetration_(nm.averagePenetration_),
    minPenetration_(nm.minPenetration_),
    checkSolverConvergence_(nm.checkSolverConvergence_),
    pressureMethod_(nm.pressureMethod_),
    epsilon0_(nm.epsilon0_),
    p0_(nm.p0_),
    beta_(nm.beta_),
    contactIterNum_(nm.contactIterNum_)
{}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //


void standardPenalty::correct
(
    const vectorField& slavePatchFaceNormals,
    const scalarField& slavePointPenetration,
    const vectorField& slaveDU,
    const vectorField& masterDUInterpToSlave
)
{
    // Preliminaries
    const fvMesh& mesh = mesh_;
    const label slavePatchIndex = slavePatchID();

    // Calculate area in contact for slave patch

    const faceList& slavePatchLocalFaces =
        mesh.boundaryMesh()[slavePatchIndex].localFaces();

    const pointField& slavePatchLocalPoints =
        mesh.boundaryMesh()[slavePatchIndex].localPoints();

    scalarField slavePatchLocalFaceAreas(slavePatchLocalFaces.size(), 0.0);

    scalarField& areaInContact = this->areaInContact();
    forAll(slavePatchLocalFaces, faceI)
    {
        areaInContact[faceI] =
            slavePatchLocalFaces[faceI].areaInContact
            (
                slavePatchLocalPoints,
                slavePointPenetration
            );

        if (areaInContact[faceI] < -SMALL)
        {
            const labelList& labels = slavePatchLocalFaces[faceI];
            scalarField vertexValue(labels.size());
            forAll(labels, i)
            {
                vertexValue[i] = slavePointPenetration[labels[i]];
            }

            FatalErrorIn(type())
                << "areaInContact is less than zero!" << nl
                << "areaInContact[" << faceI << "] = " << areaInContact[faceI]
                << nl
                << "vertexValue = " << vertexValue << nl
                << endl;
        }

        slavePatchLocalFaceAreas[faceI] =
            mag(slavePatchLocalFaces[faceI].normal(slavePatchLocalPoints));
    }

    // Calculate the point pressures
    // We will also record the average and minium penetrations

    const scalar penaltyFac = penaltyFactor();
    scalarField totalSlavePointPressure(slavePointPenetration.size(), 0.0);
    averagePenetration_ = 0.0;
    minPenetration_ = 0.0;
    int nPointsInContact = 0;

    forAll(totalSlavePointPressure, pointI)
    {
        // Take copy of penetration
        const scalar d = slavePointPenetration[pointI];

        // Note: penetration is negative for points in contact
        if (pressureMethod_ == "standard")
        {
            // The force is linearly proportional the penetration, like a spring
            if (d < epsilon0_)
            {
                totalSlavePointPressure[pointI] =
                    max(penaltyFac*(epsilon0_ - d), 0.0);

                averagePenetration_ += d;
                minPenetration_ = min(minPenetration_, d);
                nPointsInContact++;
            }
            else
            {
                totalSlavePointPressure[pointI] = 0.0;
            }
        }
        else if (pressureMethod_ == "gaussian")
        {
            // Multiply slave pressure with a Gaussian to smooth the contact
            // for points which are in the contact
            if (d < 0.0)
            {
                totalSlavePointPressure[pointI] =
                   -penaltyFac*d*(1 - exp(-0.5*pow(d/(targetAvPen_/6.0), 2.0)));

                averagePenetration_ += d;
                minPenetration_ = min(minPenetration_, d);
                nPointsInContact++;
            }
            else
            {
                totalSlavePointPressure[pointI] = 0.0;
            }
        }
        else if (pressureMethod_ == "exponential")
        {
            if (d < epsilon0_)
            {
                totalSlavePointPressure[pointI] = p0_*Foam::exp(-beta_*d);

                averagePenetration_ += d;
                minPenetration_ = min(minPenetration_, d);
                nPointsInContact++;
            }
            else
            {
                totalSlavePointPressure[pointI] = 0.0;
            }
        }
        else if (pressureMethod_ == "smooth")
        {
            if (d <= 0.0)
            {
                // Linear region
                totalSlavePointPressure[pointI] =
                   - penaltyFac*d + 0.5*epsilon0_*penaltyFac;

                averagePenetration_ += d;
                minPenetration_ = min(minPenetration_, d);
                nPointsInContact++;
            }
            else if (d <= epsilon0_)
            {
                // Quadratic region
                totalSlavePointPressure[pointI] =
                    (0.5*penaltyFac/epsilon0_)*d*d
                  - penaltyFac*d + 0.5*penaltyFac*epsilon0_;

                averagePenetration_ += d;
                minPenetration_ = min(minPenetration_, d);
                nPointsInContact++;
            }
            else
            {
                totalSlavePointPressure[pointI] = 0.0;
            }
        }
        else
        {
            FatalErrorIn
            (
                "void standardPenalty::correct\n"
                "(\n"
                "    const vectorField& slavePatchFaceNormals,\n"
                "    const scalarField& slavePointPenetration,\n"
                "    const vectorField& slaveDU,\n"
                "    const vectorField& masterDUInterpToSlave\n"
                ")"
            )   << "Unknown slave point pressure method: " << pressureMethod_
                << abort(FatalError);
        }
    }

    // Find global minimum penetration
    // IB 11/2018
    reduce(minPenetration_, minOp<scalar>());

    // Update the average penetration
    reduce(averagePenetration_, sumOp<scalar>());
    reduce(nPointsInContact, sumOp<label>());
    if (nPointsInContact > 0)
    {
        averagePenetration_ /= nPointsInContact;
    }
    else
    {
        averagePenetration_ = 0.0;
    }


    // Interpolate point pressures to the faces

    // Create local patch interpolation: No need to interpolate using the entire
    // face zone patch
    primitivePatchInterpolation localSlaveInterpolator
    (
        mesh.boundaryMesh()[slavePatchIndex]
    );

    // Interpolate point pressures to the face centres and apply in the negative
    // normal direction
    vectorField newSlaveTraction =
        localSlaveInterpolator.pointToFaceInterpolate<scalar>
        (
            totalSlavePointPressure
        )*(-slavePatchFaceNormals);

    // Under-relax pressure/traction
    // Note: slavePressure_ is really a traction vector
    slavePressure() =
        relaxFac_*newSlaveTraction + (1.0 - relaxFac_)*slavePressure();
}


void standardPenalty::newTimeStep() const
{
    if (adjustPenaltyFactor_)
    {
        const scalar growthRate = 1.25;
        const scalar reductionRate = 0.4;

        scalar pen = 0.0;
        if (adjustPenaltyFactorMethod_ == "average")
        {
            pen = -averagePenetration_;
        }
        else if (adjustPenaltyFactorMethod_ == "minimum")
        {
            pen = -minPenetration_;
        }
        else if (adjustPenaltyFactorMethod_ == "mixed")
        {
            pen = 0.5*(-averagePenetration_ + -minPenetration_);
        }
        else
        {
            FatalErrorIn("void standardPenalty::newTimeStep() const")
                << "adjustPenaltyFactorMethod " << adjustPenaltyFactorMethod_
                << " is unknown!" << abort(FatalError);
        }

        // Check if the momentum loop converged for the last time-step in the
        // solver; if it did not then we will reduce the penalty factor: this is
        // because the penalty factor being too large is one of the main reasons
        // for lack of convergence in the solver
        bool maxUIterReached = false;

        if (checkSolverConvergence_)
        {
            maxUIterReached =
                mesh_.solutionDict().subDict
                (
                    "solidMechanics"
                ).lookupOrDefault<bool>
                (
                    "maxUIterReached", false
                );
        }

        if (pen > 0.0 && !maxUIterReached)
        {
            // If the averagePenetration_ is more negative than negative
            // edgeLength then we will increase the penalty scale
            const scalar scaleFac =
                max(min(pen/targetAvPen_, growthRate), reductionRate);

            penaltyScale_ =
                max
                (
                    min
                    (
                        scaleFac*penaltyScale_,
                        maxPenaltyScale_
                    ),
                    minPenaltyScale_
                );

            reduce(penaltyScale_, maxOp<scalar>());

            // Force recalculation of penaltyFactor
            penaltyFactor_ = -1;

            Info<< "    method: " << adjustPenaltyFactorMethod_ << nl
                << "    penetration: " << pen << nl
                << "    average penetration: " << mag(averagePenetration_) << nl
                << "    greatest penetration: " << mag(minPenetration_) << nl
                << "    target penetration: " << targetAvPen_ << endl;
        }
        else if (gMax(mag(slavePressure())) > SMALL)
        {
            // Reduce penalty factor is faces have a pressure but are not in
            // contact
            penaltyScale_ =
                max
                (
                    reductionRate*penaltyScale_,
                    minPenaltyScale_
                );

            // Force re-calculation of penaltyFactor
            penaltyFactor_ = -1;

            if (maxUIterReached)
            {
                Info<< "    Reducing the penalty factor because the momentum "
                    << "loop did not converge in the previous time-step"
                    << endl;
            }
        }

        Info<< "    penaltyScale: " << penaltyScale_ << endl;
    }
}


scalar standardPenalty::penaltyFactor() const
{
    if (penaltyFactor_ < -SMALL)
    {
        calcPenaltyFactor();
    }

    return penaltyFactor_;
}


scalar standardPenalty::updatePenaltyScale(const scalar previousPenaltyScale)
{
    if (previousPenaltyScale > 0.0)
    {
        // Lookup initial value for penaltyScale
        const scalar initialPenaltyScale =
            normalContactModelDict_.lookupOrDefault<scalar>
            (
                "penaltyScale", 1.0
            );

        if (mag(initialPenaltyScale - penaltyScale_) < SMALL)
        {
            // After a topo change, use the previous value of penalty scale
            penaltyScale_ = previousPenaltyScale;
        }
    }

    return penaltyScale_;
}


void standardPenalty::autoMap(const fvPatchFieldMapper& m)
{
    if (debug)
    {
        InfoIn
        (
            "void standardPenalty::autoMap(const fvPatchFieldMapper& m)"
        )   << "autoMap" << endl;
    }

    normalContactModel::autoMap(m);

    // The internal fields for the volFields should always be zero
    // We will reset them as they may not be zero after field advection
    slavePressureVolField_.internalField() = vector::zero;
    areaInContactVolField_.internalField() = 0.0;
}


void standardPenalty::writeDict(Ostream& os) const
{
    // Update the penalty scale in the dictionary
    normalContactModelDict_.set("penaltyScale", penaltyScale_);

    // Write the dictionary
    word keyword(name() + "NormalModelDict");
    os.writeKeyword(keyword)
        << normalContactModelDict_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
