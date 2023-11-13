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

Class
    solidGeneralContactFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solidGeneralContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"
#include "Switch.H"
#include "pointFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::solidGeneralContactFvPatchVectorField::movingMesh() const
{
    // If the deformation gradient "F" and the displacement increment DU" are
    // found then we can assume it is a moving mesh (updated Lagrangian) case
    if
    (
        db().foundObject<volVectorField>("DU")
     && db().foundObject<volTensorField>("F")
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::solidGeneralContactFvPatchVectorField::
moveZonesToDeformedConfiguration()
{
    // Only the master moves the zones
    if (!globalMaster())
    {
        return;
    }

    // Reference to mesh for tidiness
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Method
    // We will interpolate the patch face displacements to the patch vertices
    // and then add these vertex/point displacements to the initial patch
    // points
    // We need to take care in parallel, and also realise that the solidModel
    // might have a moving or stationary mesh

    // Shadow patch and zone indices
    const labelList& shadPatchIndices = shadowPatchIndices();
    const labelList& shadZoneIndices = shadowZoneIndices();

    forAll(shadPatchIndices, shadowI)
    {
        // Assemble the zone face displacement field to move the zones
        vectorField zoneD(zone().size(), vector::zero);
        vectorField shadowZoneD(shadowZone(shadowI).size(), vector::zero);

        // For a non-moving mesh, we will move the zones by the total
        // displacement, whereas for a moving mesh (updated Lagrangian), we will
        // move the zones by the displacement increment

        if (movingMesh())
        {
            // Updated Lagrangian, so we will move the zones by the displacement
            // increment

            // Lookup the current total displacement field
            const volVectorField& DD = db().lookupObject<volVectorField>("DU");

            // Take a reference to the patch face displacement increment field
            const vectorField& patchDD =
                DD.boundaryField()[patch().index()];
            const vectorField& shadowPatchDD =
                DD.boundaryField()[shadPatchIndices[shadowI]];

            zoneD =
                zoneField(zoneIndex(), patch().index(), patchDD);
            shadowZoneD =
                    zoneField
                    (
                        shadZoneIndices[shadowI],
                        shadPatchIndices[shadowI],
                        shadowPatchDD
                    );
        }
        else
        {
            // Non-moving mesh: we will move the zones by the total displacement

            // Lookup the current total displacement field
            const volVectorField& D = db().lookupObject<volVectorField>("U");

            // Take a reference to the patch face total displacement field
            const vectorField& patchD =
                D.boundaryField()[patch().index()];

            const vectorField& shadowPatchDD =
                D.boundaryField()[shadPatchIndices[shadowI]];

            zoneD =
                zoneField(zoneIndex(), patch().index(), patchD);
            shadowZoneD =
                zoneField
                (
                    shadZoneIndices[shadowI],
                    shadPatchIndices[shadowI],
                    shadowPatchDD
                );
        }

        // Interpolate the zone face field to the zone points
        const pointField zonePointD =
            zoneFaceToPointInterpolate(zoneIndex(), zoneD, -1);
        const pointField shadowZonePointD =
            zoneFaceToPointInterpolate
            (
                shadZoneIndices[shadowI], shadowZoneD, shadowI
            );

        // The zone deformed points are the initial position plus the
        // displacement
        const pointField zoneNewPoints =
            mesh.faceZones()[zoneIndex()]().localPoints()
          + zonePointD;
        const pointField shadowZoneNewPoints =
            mesh.faceZones()[shadZoneIndices[shadowI]]().localPoints()
          + shadowZonePointD;

        // Move the zones

        // Remove zones weights
        if (shadowI == 0)
        {
            zone().movePoints(zoneNewPoints);
        }
        shadowZone(shadowI).movePoints(shadowZoneNewPoints);

        // We need to use const_cast to move the standAlonePatch points as the
        // movePoints function only clears weights
        // Also, be careful to move the points are opposed to the localPoints
        if (shadowI == 0)
        {
            const_cast<pointField&>(zone().points()) = zoneNewPoints;
        }
        const_cast<pointField&>(shadowZone(shadowI).points()) =
            shadowZoneNewPoints;
    }
}


void Foam::solidGeneralContactFvPatchVectorField::calcGlobalMaster() const
{
    if (globalMasterPtr_)
    {
        FatalErrorIn
            (
                "void Foam::solidGeneralContactFvPatchVectorField::"
                "calcGlobalMaster() const"
            )   << "globalMasterPtr_ already set" << abort(FatalError);
    }

    // The global master is the first solidGeneralContact patch i.e. the one
    // with the lowest patch index

    if (globalMasterIndex() == patch().index())
    {
        globalMasterPtr_ = new bool(true);
    }
    else
    {
        globalMasterPtr_ = new bool(false);
    }
}


void Foam::solidGeneralContactFvPatchVectorField::calcGlobalMasterIndex() const
{
    if (globalMasterIndexPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solidGeneralContactFvPatchVectorField::"
            "calcGlobalMasterIndex() const"
        )   << "globalMasterIndexPtr_ already set" << abort(FatalError);
    }

    // The global master is the first solidGeneralContact patch i.e. the one
    // with the lowest patch index

    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    globalMasterIndexPtr_ = new label(-1);
    label& gMasterID = *globalMasterIndexPtr_;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == solidGeneralContactFvPatchVectorField::typeName
        )
        {
            gMasterID = patchI;

            break;
        }
    }

    // Check there is only one global master

    label GMasterID = returnReduce(gMasterID, maxOp<label>());

    if (gMasterID != GMasterID)
    {
        FatalErrorIn
        (
            "solidGeneralContactFvPatchVectorField::"
            "calcGlobalMasterIndex() const"
        )   << "There are multiple global masters" << abort(FatalError);
    }

    if (debug)
    {
        Pout<< nl << "The global master contact patch is "
            << patch().boundaryMesh()[gMasterID].name() << endl;
    }
}


void Foam::solidGeneralContactFvPatchVectorField::calcLocalSlave() const
{
    if (localSlavePtr_)
    {
        FatalErrorIn
        (
            "label Foam::solidGeneralContactFvPatchVectorField::"
            "calcLocalSlave() const"
        )   << "localSlavePtr_ already set" << abort(FatalError);
    }

    localSlavePtr_ = new boolList(shadowPatchNames().size(), false);

    boolList& localSlave = *localSlavePtr_;

    forAll(localSlave, shadowI)
    {
        if (patch().index() < shadowPatchIndices()[shadowI])
        {
            localSlave[shadowI] = true;

            Info<< "solidGeneralContact: "
                << shadowPatchNames()[shadowI] << " (master)" << " to "
                << patch().name() << " (slave)" << endl;
        }
    }
}


const Foam::boolList&
Foam::solidGeneralContactFvPatchVectorField::localSlave() const
{
    if (!localSlavePtr_)
    {
        calcLocalSlave();
    }

    return *localSlavePtr_;
}


void Foam::solidGeneralContactFvPatchVectorField::calcShadowPatchNames() const
{
    if (shadowPatchNamesPtr_ || shadowPatchIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solidGeneralContactFvPatchVectorField::"
            "calcShadowPatchNames() const"
        )   << "shadowPatchNames_ or shadowPatchIndices_ already set"
            << abort(FatalError);
    }

    // Add each solidGeneralContact patch in the order of increasing patch index

    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    // Count shadow patches

    label nShadPatches = 0;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == solidGeneralContactFvPatchVectorField::typeName
            && patchI != patch().index()
        )
        {
            nShadPatches++;
        }
    }

    shadowPatchNamesPtr_ = new wordList(nShadPatches);
    wordList& shadowPatchNames = *shadowPatchNamesPtr_;

    shadowPatchIndicesPtr_ = new labelList(nShadPatches);
    labelList& shadowPatchIndices = *shadowPatchIndicesPtr_;

    // Record shadow patch names

    label shadowI = 0;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == solidGeneralContactFvPatchVectorField::typeName
            && patchI != patch().index()
        )
        {
            shadowPatchNames[shadowI] = patch().boundaryMesh()[patchI].name();

            shadowPatchIndices[shadowI++] = patchI;
        }
    }
}


void Foam::solidGeneralContactFvPatchVectorField::calcShadowZoneNames() const
{
    if (shadowZoneNamesPtr_ || shadowZoneIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::solidGeneralContactFvPatchVectorField::"
            "calcShadowZoneNames() const"
        )   << "shadowZoneNames_ or shadowZoneIndices_ already set"
            << abort(FatalError);
    }

    const wordList& shadNames = shadowPatchNames();

    shadowZoneNamesPtr_ = new wordList(shadNames.size());
    wordList& shadowZoneNames = *shadowZoneNamesPtr_;

    shadowZoneIndicesPtr_ = new labelList(shadNames.size());
    labelList& shadowZoneIndices = *shadowZoneIndicesPtr_;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    forAll(shadNames, shadowI)
    {
        word zoneName = shadNames[shadowI] + "FaceZone";

        faceZoneID zone(zoneName, mesh.faceZones());

        if (!zone.active())
        {
            FatalErrorIn("solidGeneralContactFvPatchVectorField")
                << "Face zone name " << zoneName
                << " not found.  Please check your zone definition."
                << abort(FatalError);
        }

        shadowZoneNames[shadowI] = zoneName;

        shadowZoneIndices[shadowI] = zone.index();
    }
}


void Foam::solidGeneralContactFvPatchVectorField::calcNormalModels() const
{
    if (!normalModels_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::"
            "calcNormalModel() const"
        )   << "normalModels already set" << abort(FatalError);
    }

    normalModels_.setSize(shadowPatchNames().size());

    const boolList& locSlave = localSlave();

    forAll(normalModels_, shadowI)
    {
        // Only the local slave creates the contact model
        if (locSlave[shadowI])
        {
            // Calculate normal contact forces
            normalModels_.set
            (
                shadowI,
                normalContactModel::New
                (
                    word(dict().lookup("normalContactModel")),
                    patch().boundaryMesh()[shadowPatchIndices()[shadowI]],
                    dict(),
                    shadowPatchIndices()[shadowI], // master
                    patch().index(), // slave
                    shadowZone(shadowI), // master
                    zone() // slave
                )
            );
        }
    }
}


Foam::normalContactModel&
Foam::solidGeneralContactFvPatchVectorField::normalModel(const label shadowI)
{
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("normalModel(const label shadowI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (normalModels_.empty())
    {
        calcNormalModels();
    }

    return normalModels_[shadowI];
}


const Foam::normalContactModel&
Foam::solidGeneralContactFvPatchVectorField::normalModel
(
    const label shadowI
) const
{
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("normalModel(const label shadowI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (normalModels_.empty())
    {
        calcNormalModels();
    }

    return normalModels_[shadowI];
}


void Foam::solidGeneralContactFvPatchVectorField::calcFrictionModels() const
{
    if (!frictionModels_.empty())
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::"
            "calcFrictionModel(shadowI) const"
        )   << "frictionModelPtr_[shadowI] already set" << abort(FatalError);
    }

    frictionModels_.setSize(shadowPatchNames().size());

    const boolList& locSlave = localSlave();

    forAll(frictionModels_, shadowI)
    {
        if (locSlave[shadowI])
        {
            frictionModels_.set
                (
                    shadowI,
                    frictionContactModel::New
                    (
                        word(dict().lookup("frictionContactModel")),
                        patch().boundaryMesh()[shadowPatchIndices()[shadowI]],
                        dict(),
                        shadowPatchIndices()[shadowI], // master
                        patch().index() // slave
                    )
                );
        }
    }
}


Foam::frictionContactModel&
Foam::solidGeneralContactFvPatchVectorField::frictionModel(const label shadowI)
{
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("frictionModel(const label shadowI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (frictionModels_.empty())
    {
        calcFrictionModels();
    }

    return frictionModels_[shadowI];
}


const Foam::frictionContactModel&
Foam::solidGeneralContactFvPatchVectorField::frictionModel
(
    const label shadowI
) const
{
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("frictionModel(const label shadowI)")
            << "Only the local slave can call the contact model"
            << abort(FatalError);
    }

    if (frictionModels_.empty())
    {
        calcFrictionModels();
    }

    return frictionModels_[shadowI];
}


void Foam::solidGeneralContactFvPatchVectorField::calcZoneIndex() const
{
    if (zoneIndex_ != -1)
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::calcZoneIndex()"
            "const"
        )   << "zoneIndex_ already set" << abort(FatalError);
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    word zoneName = patch().name() + "FaceZone";

    faceZoneID zone(zoneName, mesh.faceZones());

    if (!zone.active())
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::calcZoneIndex()"
            "const"
        )   << "Face zone name " << zoneName
            << " not found.  Please check your zone definition." << nl
            << "Current faceZones are:" << mesh.faceZones().names()
            << abort(FatalError);
    }

    zoneIndex_ = zone.index();
}


void Foam::solidGeneralContactFvPatchVectorField::calcZone() const
{
    if (zonePtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::calcZone() const"
        )   << "zonePtr_ already set" << abort(FatalError);
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Note: the main mesh will either be in the initial configuration or the
    // updated configuration
    zonePtr_ =
        new standAlonePatch
        (
            mesh.faceZones()[zoneIndex()]().localFaces(),
            mesh.faceZones()[zoneIndex()]().localPoints()
        );
}


Foam::label Foam::solidGeneralContactFvPatchVectorField::zoneIndex() const
{
    if (zoneIndex_ == -1)
    {
        calcZoneIndex();
    }

    return zoneIndex_;
}


const Foam::standAlonePatch&
Foam::solidGeneralContactFvPatchVectorField::zone() const
{
    if (!zonePtr_)
    {
        calcZone();
    }

    return *zonePtr_;
}


Foam::standAlonePatch& Foam::solidGeneralContactFvPatchVectorField::zone()
{
    if (!zonePtr_)
    {
        calcZone();
    }

    return *zonePtr_;
}


const Foam::standAlonePatch&
Foam::solidGeneralContactFvPatchVectorField::shadowZone
(
    const label shadowI
) const
{
    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    const solidGeneralContactFvPatchVectorField& shadowPatchField =
        refCast<const solidGeneralContactFvPatchVectorField>
        (
            field.boundaryField()[shadowPatchIndices()[shadowI]]
        );

    return shadowPatchField.zone();
}


Foam::standAlonePatch&
Foam::solidGeneralContactFvPatchVectorField::shadowZone
(
    const label shadowI
)
{
    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    // Const cast away the const-ness
    solidGeneralContactFvPatchVectorField& shadowPatchField =
        const_cast<solidGeneralContactFvPatchVectorField&>
        (
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[shadowI]]
            )
        );

    return shadowPatchField.zone();
}


void Foam::solidGeneralContactFvPatchVectorField::calcZoneToZones() const
{
    // Create zone-to-zone interpolation
    if (!zoneToZones_.empty())
    {
        FatalErrorIn
        (
            "void solidGeneralContactFvPatchVectorField::calcZoneToZones()"
            "const"
        )   << "Zone to zone interpolation already calculated"
            << abort(FatalError);
    }

    zoneToZones_.setSize(shadowPatchNames().size());

    const boolList& locSlave = localSlave();

    forAll(zoneToZones_, shadowI)
    {
        // Only the local slave creates the interpolator
        if (locSlave[shadowI])
        {
            zoneToZones_.set
                (
                    shadowI,
                    new newGgiStandAlonePatchInterpolation
                    (
                        shadowZone(shadowI), // master
                        zone(), // slave
                        tensorField(0),
                        tensorField(0),
                        vectorField(0), // Slave-to-master separation.
                        true,           // global data
                        0,              // Non-overlapping face tolerances
                        0,              //
                        true,           // Rescale weighting factors.
                        newGgiInterpolation::AABB
                        //newGgiInterpolation::BB_OCTREE
                        //newGgiInterpolation::THREE_D_DISTANCE
                        //newGgiInterpolation::N_SQUARED
                    )
                );
        }
    }
}


const Foam::newGgiStandAlonePatchInterpolation&
Foam::solidGeneralContactFvPatchVectorField::zoneToZone
(
    const label shadowI
) const
{
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("zoneToZone(const label shadowI)")
            << "Only the local slave can call the zoneToZone interpolator"
            << abort(FatalError);
    }

    if (zoneToZones_.empty())
    {
        word zoneName =
            patch().boundaryMesh().mesh().faceZones()[zoneIndex()].name();

        if (debug)
        {
            Info<< "Initializing the GGI interpolators for " << zoneName
                << endl;
        }

        calcZoneToZones();
    }

    return zoneToZones_[shadowI];
}


Foam::newGgiStandAlonePatchInterpolation&
Foam::solidGeneralContactFvPatchVectorField::zoneToZone(const label shadowI)
{
    if (!localSlave()[shadowI])
    {
        FatalErrorIn("zoneToZone(const label shadowI)")
            << "Only the local slave can call the zoneToZone interpolator"
            << abort(FatalError);
    }

    if (zoneToZones_.empty())
    {
        word zoneName =
            patch().boundaryMesh().mesh().faceZones()[zoneIndex()].name();

        if (debug)
        {
            Info<< "Initializing the GGI interpolators for " << zoneName
                << endl;
        }

        calcZoneToZones();
    }

    return zoneToZones_[shadowI];
}


Foam::label Foam::solidGeneralContactFvPatchVectorField::findShadowID
(
    const label patchID
) const
{
    label shadowI = -1;

    const labelList shadowIDs = shadowPatchIndices();

    forAll(shadowIDs, I)
    {
        if (patchID == shadowIDs[I])
        {
            shadowI = I;
            break;
        }
    }

    if (shadowI == -1)
    {
        FatalErrorIn("findShadowID(const label patchID)")
            << "shadow patch not found!" << abort(FatalError);
    }

    return shadowI;
}


void Foam::solidGeneralContactFvPatchVectorField::makeCurPatchTractions() const
{
    if (curPatchTractionPtr_)
    {
        FatalErrorIn
        (
            "void Foam::solidGeneralContactFvPatchVectorField::"
            "makeCurPatchTractions() const"
        )   << "curPatchTractionPtr_ already set" << abort(FatalError);
    }

    curPatchTractionPtr_ =
        new List<vectorField>
        (
            shadowPatchNames().size(),
            vectorField(patch().size(), vector::zero)
        );
}


void Foam::solidGeneralContactFvPatchVectorField::calcQc() const
{
    if (QcPtr_)
    {
        FatalErrorIn("solidGeneralContactFvPatchVectorField::calcQc")
            << "QcPtr_ already set!" << abort(FatalError);
    }

    QcPtr_ = new scalarField(patch().size(), 0.0);

    scalarField& Qc = *QcPtr_;

    // For now, we assume traction is constant over time-step
    // Todo: use trapezoidal rule
    vectorField curTraction(Qc.size(), vector::zero);

    // sigma/sigmaCauchy is up-to-date as Qc is called after momentum loop
    // has converged and sigma has been updated and mesh moved
    if
    (
        db().objectRegistry::foundObject<volSymmTensorField>
        (
            "sigmaCauchy"
        )
    )
    {
        const symmTensorField& sigma =
            db().objectRegistry::lookupObject<volSymmTensorField>
            (
                "sigmaCauchy"
            ).boundaryField()[patch().index()];

        curTraction = patch().nf() & sigma;
    }
    else
    {
        const symmTensorField& sigma =
            db().objectRegistry::lookupObject<volSymmTensorField>
            (
                "sigma"
            ).boundaryField()[patch().index()];

        curTraction = patch().nf() & sigma;
    }


    // Accumulate Qc for from all shadows

    const scalar deltaT = patch().boundaryMesh().mesh().time().deltaT().value();

    const volVectorField& field =
        db().objectRegistry::lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    const boolList& locSlave = localSlave();

    forAll(locSlave, shadowI)
    {
        vectorField curPatchSlip(patch().size(), vector::zero);

        // Calculate slip
        if (locSlave[shadowI])
        {
            curPatchSlip = frictionModel(shadowI).slip();
        }
        else
        {
            const solidGeneralContactFvPatchVectorField& shadowPatchField =
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[shadowI]]
            );

            const label locShadowID =
                shadowPatchField.findShadowID(patch().index());

            vectorField shadowPatchSlip =
                shadowPatchField.frictionModel(locShadowID).slip();

            vectorField shadowZoneSlip =
                zoneField
                (
                    shadowZoneIndices()[shadowI],
                    shadowPatchIndices()[shadowI],
                    shadowPatchSlip
                );

            // Interpolate from shadow to the current patch
            // Face-to-face

            vectorField curZoneSlip =
                shadowPatchField.zoneToZone(locShadowID).slaveToMaster
                (
                    shadowZoneSlip
                );

            curPatchSlip =
                patchField
                (
                    patch().index(),
                    zoneIndex(),
                    curZoneSlip
                );
        }

        // Heat flux rate: rate of dissipated frictional energy
        // The dot product of the traction vectors and the slip vectors
        // gives the dissipated frictional energy rate per unit area, which
        // is always positive
        Qc += mag(curTraction & (curPatchSlip/deltaT));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::solidGeneralContactFvPatchVectorField::
solidGeneralContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(-1),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    rigidMaster_(false),
    dict_(NULL),
    normalModels_(0),
    frictionModels_(0),
    zonePtr_(NULL),
    zoneToZones_(0),
    alg_(Foam::intersection::VISIBLE),
    dir_(Foam::intersection::CONTACT_SPHERE),
    curTimeIndex_(-1),
    curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(0.0)
{
    if (debug)
    {
        InfoIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const fvPatch& p,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        ) << endl;
    }
}


Foam::solidGeneralContactFvPatchVectorField::
solidGeneralContactFvPatchVectorField
(
    const solidGeneralContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(ptf, p, iF, mapper),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(ptf.zoneIndex_),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    zonePtr_(NULL),
    zoneToZones_(0),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(ptf.bbOffset_)
{
    if (debug)
    {
        InfoIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf,"
            "    const fvPatch& p,"
            "    const DimensionedField<vector, volMesh>& iF,"
            "    const fvPatchFieldMapper& mapper"
            ")"
        ) << endl;
    }

    // Copy pointer objects

    if (ptf.globalMasterPtr_)
    {
        globalMasterPtr_ = new bool(*ptf.globalMasterPtr_);
    }

    if (ptf.globalMasterIndexPtr_)
    {
        globalMasterIndexPtr_ = new label(*ptf.globalMasterIndexPtr_);
    }

    if (ptf.localSlavePtr_)
    {
        localSlavePtr_ = new boolList(*ptf.localSlavePtr_);
    }

    if (ptf.shadowPatchNamesPtr_)
    {
        shadowPatchNamesPtr_ = new wordList(*ptf.shadowPatchNamesPtr_);
    }

    if (ptf.shadowPatchIndicesPtr_)
    {
        shadowPatchIndicesPtr_ = new labelList(*ptf.shadowPatchIndicesPtr_);
    }

    if (ptf.shadowZoneNamesPtr_)
    {
        shadowZoneNamesPtr_ = new wordList(*ptf.shadowZoneNamesPtr_);
    }

    if (ptf.shadowZoneIndicesPtr_)
    {
        shadowZoneIndicesPtr_ = new labelList(*ptf.shadowZoneIndicesPtr_);
    }

    if (ptf.zonePtr_)
    {
        zonePtr_ = new standAlonePatch(*ptf.zonePtr_);
    }

    if (!ptf.zoneToZones_.empty())
    {
        // I will not copy the GGI interpolators
        // They can be re-created when required
        WarningIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        )   << "solidGeneralContact: zoneToZone GGI interpolators not mapped"
            << endl;
    }

    if (ptf.curPatchTractionPtr_)
    {
        curPatchTractionPtr_ =
            new List<vectorField>(*ptf.curPatchTractionPtr_);
    }

    if (ptf.QcPtr_)
    {
        QcPtr_ = new scalarField(*ptf.QcPtr_);
    }

    if (ptf.QcsPtr_)
    {
        QcsPtr_ = new List<scalarField>(*ptf.QcsPtr_);
    }
}


Foam::solidGeneralContactFvPatchVectorField::
solidGeneralContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:   solidTractionFvPatchVectorField(p, iF),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(-1),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    rigidMaster_(dict.lookupOrDefault<Switch>("rigidMaster", false)),
    dict_(dict),
    normalModels_(0),
    frictionModels_(0),
    zonePtr_(0),
    zoneToZones_(0),
    alg_(Foam::intersection::VISIBLE),
    dir_(Foam::intersection::CONTACT_SPHERE),
    curTimeIndex_(-1),
    curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(0.0)
{
    Info<< "Creating " << solidGeneralContactFvPatchVectorField::typeName
        << " patch" << endl;

    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        Field<vector>::operator=
        (
            patchInternalField() + gradient()/patch().deltaCoeffs()
        );
    }
}


Foam::solidGeneralContactFvPatchVectorField::
solidGeneralContactFvPatchVectorField
(
    const solidGeneralContactFvPatchVectorField& ptf
)
:
    solidTractionFvPatchVectorField(ptf),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(ptf.zoneIndex_),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    zonePtr_(NULL),
    zoneToZones_(0),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    curPatchTractionPtr_(NULL),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(ptf.bbOffset_)
{
    if (debug)
    {
        InfoIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf"
            ")"
        ) << endl;
    }

    // Copy pointer objects

    if (ptf.globalMasterPtr_)
    {
        globalMasterPtr_ = new bool(*ptf.globalMasterPtr_);
    }

    if (ptf.globalMasterIndexPtr_)
    {
        globalMasterIndexPtr_ = new label(*ptf.globalMasterIndexPtr_);
    }

    if (ptf.localSlavePtr_)
    {
        localSlavePtr_ = new boolList(*ptf.localSlavePtr_);
    }

    if (ptf.shadowPatchNamesPtr_)
    {
        shadowPatchNamesPtr_ = new wordList(*ptf.shadowPatchNamesPtr_);
    }

    if (ptf.shadowPatchIndicesPtr_)
    {
        shadowPatchIndicesPtr_ = new labelList(*ptf.shadowPatchIndicesPtr_);
    }

    if (ptf.shadowZoneNamesPtr_)
    {
        shadowZoneNamesPtr_ = new wordList(*ptf.shadowZoneNamesPtr_);
    }

    if (ptf.shadowZoneIndicesPtr_)
    {
        shadowZoneIndicesPtr_ = new labelList(*ptf.shadowZoneIndicesPtr_);
    }

    if (ptf.zonePtr_)
    {
        zonePtr_ = new standAlonePatch(*ptf.zonePtr_);
    }

    if (!ptf.zoneToZones_.empty())
    {
        // I will not copy the GGI interpolators
        // They can be re-created when required
        WarningIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        )   << "solidGeneralContact: zoneToZone GGI interpolators not mapped"
            << endl;
    }

    if (ptf.curPatchTractionPtr_)
    {
        curPatchTractionPtr_ =
            new List<vectorField>(*ptf.curPatchTractionPtr_);
    }

    if (ptf.QcPtr_)
    {
        QcPtr_ = new scalarField(*ptf.QcPtr_);
    }

    if (ptf.QcsPtr_)
    {
        QcsPtr_ = new List<scalarField>(*ptf.QcsPtr_);
    }
}


Foam::solidGeneralContactFvPatchVectorField::
solidGeneralContactFvPatchVectorField
(
    const solidGeneralContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(ptf, iF),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(ptf.zoneIndex_),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    rigidMaster_(ptf.rigidMaster_),
    dict_(ptf.dict_),
    normalModels_(ptf.normalModels_),
    frictionModels_(ptf.frictionModels_),
    zonePtr_(NULL),
    zoneToZones_(0),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    QcPtr_(NULL),
    QcsPtr_(NULL),
    bbOffset_(ptf.bbOffset_)
{
    if (debug)
    {
        InfoIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        ) << endl;
    }

    // Copy pointer objects

    if (ptf.globalMasterPtr_)
    {
        globalMasterPtr_ = new bool(*ptf.globalMasterPtr_);
    }

    if (ptf.globalMasterIndexPtr_)
    {
        globalMasterIndexPtr_ = new label(*ptf.globalMasterIndexPtr_);
    }

    if (ptf.localSlavePtr_)
    {
        localSlavePtr_ = new boolList(*ptf.localSlavePtr_);
    }

    if (ptf.shadowPatchNamesPtr_)
    {
        shadowPatchNamesPtr_ = new wordList(*ptf.shadowPatchNamesPtr_);
    }

    if (ptf.shadowPatchIndicesPtr_)
    {
        shadowPatchIndicesPtr_ = new labelList(*ptf.shadowPatchIndicesPtr_);
    }

    if (ptf.shadowZoneNamesPtr_)
    {
        shadowZoneNamesPtr_ = new wordList(*ptf.shadowZoneNamesPtr_);
    }

    if (ptf.shadowZoneIndicesPtr_)
    {
        shadowZoneIndicesPtr_ = new labelList(*ptf.shadowZoneIndicesPtr_);
    }

    if (ptf.zonePtr_)
    {
        zonePtr_ = new standAlonePatch(*ptf.zonePtr_);
    }

    if (!ptf.zoneToZones_.empty())
    {
        // I will not copy the GGI interpolators
        // They can be re-created when required
        WarningIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "solidGeneralContactFvPatchVectorField"
            "("
            "    const solidGeneralContactFvPatchVectorField& ptf,"
            "    const DimensionedField<vector, volMesh>& iF"
            ")"
        )   << "solidGeneralContact: zoneToZone GGI interpolators not mapped"
            << endl;
    }

    if (ptf.curPatchTractionPtr_)
    {
        curPatchTractionPtr_ =
            new List<vectorField>(*ptf.curPatchTractionPtr_);
    }

    if (ptf.QcPtr_)
    {
        QcPtr_ = new scalarField(*ptf.QcPtr_);
    }

    if (ptf.QcsPtr_)
    {
        QcsPtr_ = new List<scalarField>(*ptf.QcsPtr_);
    }
}


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::solidGeneralContactFvPatchVectorField::
~solidGeneralContactFvPatchVectorField()
{
    if (debug)
    {
        InfoIn
        (
            "Foam::solidGeneralContactFvPatchVectorField::"
            "~solidGeneralContactFvPatchVectorField()"
        ) << endl;
    }

    deleteDemandDrivenData(globalMasterPtr_);
    deleteDemandDrivenData(globalMasterIndexPtr_);
    deleteDemandDrivenData(localSlavePtr_);
    deleteDemandDrivenData(shadowPatchNamesPtr_);
    deleteDemandDrivenData(shadowPatchIndicesPtr_);
    deleteDemandDrivenData(shadowZoneNamesPtr_);
    deleteDemandDrivenData(shadowZoneIndicesPtr_);

    normalModels_.clear();
    frictionModels_.clear();

    deleteDemandDrivenData(zonePtr_);

    zoneToZones_.clear();

    deleteDemandDrivenData(curPatchTractionPtr_);
    deleteDemandDrivenData(QcPtr_);
    deleteDemandDrivenData(QcsPtr_);
}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //


bool Foam::solidGeneralContactFvPatchVectorField::globalMaster() const
{
    if (!globalMasterPtr_)
    {
        calcGlobalMaster();
    }

    return *globalMasterPtr_;
}


Foam::label Foam::solidGeneralContactFvPatchVectorField::globalMasterIndex()
const
{
    if (!globalMasterIndexPtr_)
    {
        calcGlobalMasterIndex();
    }

    return *globalMasterIndexPtr_;
}


const Foam::List<Foam::word>&
Foam::solidGeneralContactFvPatchVectorField::shadowPatchNames() const
{
    if (!shadowPatchNamesPtr_)
    {
        calcShadowPatchNames();
    }

    return *shadowPatchNamesPtr_;
}


const Foam::List<Foam::label>&
Foam::solidGeneralContactFvPatchVectorField::shadowPatchIndices() const
{
    if (!shadowPatchIndicesPtr_)
    {
        calcShadowPatchNames();
    }

    return *shadowPatchIndicesPtr_;
}


const Foam::List<Foam::word>&
Foam::solidGeneralContactFvPatchVectorField::shadowZoneNames() const
{
    if (!shadowZoneNamesPtr_)
    {
        calcShadowZoneNames();
    }

    return *shadowZoneNamesPtr_;
}


const Foam::List<Foam::label>&
Foam::solidGeneralContactFvPatchVectorField::shadowZoneIndices() const
{
    if (!shadowZoneIndicesPtr_)
    {
        calcShadowZoneNames();
    }

    return *shadowZoneIndicesPtr_;
}


// Map from self
void Foam::solidGeneralContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    FatalErrorIn
    (
        "void Foam::solidGeneralContactFvPatchVectorField::autoMap"
        "("
        "    const fvPatchFieldMapper& m"
        ")"
    )   << "member mapping not implemented" << endl;

    solidTractionFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::solidGeneralContactFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    FatalErrorIn
    (
        "void Foam::solidGeneralContactFvPatchVectorField::rmap"
        "("
        "    const fvPatchField<vector>& ptf,"
        "    const labelList& addr"
        ")"
    )   << "member mapping not implemented" << endl;

    solidTractionFvPatchVectorField::rmap(ptf, addr);
}


void Foam::solidGeneralContactFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    boolList activeContactPairs(shadowPatchNames().size(), false);

    // if it is a new time step then reset iCorr
    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        // Delete friction heat rate to force its recalculation when thermal
        // boundaries ask for it
        deleteDemandDrivenData(QcPtr_);
        deleteDemandDrivenData(QcsPtr_);

        if (globalMaster())
        {
            forAll(activeContactPairs, shadowI)
            {
                // Let the contact models know that it is a new time-step, in
                // case they need to update anything
                normalModel(shadowI).newTimeStep();
                frictionModel(shadowI).newTimeStep();
            }
        }
    }

    // Method
    // Move all global face zones to the deformed configuration
    // Clear interpolator weights
    // Perform quick check to find potential contacting pairs
    // Call normal and frction contact models for active contacting pairs
    // Accumulate contact force contributions for all active contact pairs


    if (rigidMaster_)
    {
        // Set to master to traction free to mimic a rigid patch
        traction() = vector::zero;
    }
    else
    {
        // Move all global face zones to the deformed configuration
        if (globalMaster())
        {
           // Move the master and slave zone to the deformed configuration
            moveZonesToDeformedConfiguration();
        }


        // Clear interpolator weights

        forAll(activeContactPairs, slaveI)
        {
            if (localSlave()[slaveI])
            {
                zoneToZone(slaveI).movePoints
                (
                    tensorField(0), tensorField(0), vectorField(0)
                );
            }
        }


        // Accumulated traction for the current patch
        vectorField curPatchTraction(patch().size(), vector::zero);

        // Only the local masters calculates the contact force and the local
        // master interpolates this force
        const boolList& locSlave = localSlave();

        // Create master bounding box used for quick check
        boundBox masterBb(zone().localPoints(), false);

        // The BB may have zero thickness in one of the directions e.g. for a
        // flat patch, so we will check for this and, if found, create an offset
        const scalar bbOff = bbOffset();
        if (masterBb.minDim() < bbOff)
        {
            const vector bbDiag = masterBb.max() - masterBb.min();

            if (bbDiag.x() < bbOff)
            {
                vector offset(bbOff, 0, 0);
                masterBb.min() -= offset;
                masterBb.max() += offset;
            }
            else if (bbDiag.y() < bbOff)
            {
                vector offset(0, bbOff, 0);
                masterBb.min() -= offset;
                masterBb.max() += offset;
            }
            else if (bbDiag.z() < bbOff)
            {
                vector offset(0, 0, bbOff);
                masterBb.min() -= offset;
                masterBb.max() += offset;
            }
        }

        forAll(activeContactPairs, shadowI)
        {
            // Perform quick check to find potential contacting pairs
            // The quick check is based on the bounding box (BB) of the contact
            // pairs: if the BBs of pair intersect then we will designate the
            // pair as active.

            // Create shadow bounding box
            boundBox shadowBb(shadowZone(shadowI).localPoints(), false);

            // Check for a zero dimension in the shadowBb
            if (shadowBb.minDim() < bbOff)
            {
                const vector bbDiag = shadowBb.max() - shadowBb.min();

                if (bbDiag.x() < bbOff)
                {
                    vector offset(bbOff, 0, 0);
                    shadowBb.min() -= offset;
                    shadowBb.max() += offset;
                }
                else if (bbDiag.y() < bbOff)
                {
                    vector offset(0, bbOff, 0);
                    shadowBb.min() -= offset;
                    shadowBb.max() += offset;
                }
                else if (bbDiag.z() < bbOff)
                {
                    vector offset(0, 0, bbOff);
                    shadowBb.min() -= offset;
                    shadowBb.max() += offset;
                }
            }

            if (masterBb.overlaps(shadowBb))
            {
                activeContactPairs[shadowI] = true;
            }

            // Call normal and frction contact models for active contacting
            // pairs
            // Accumulate contact force contributions for all active contact
            // pairs

            if (activeContactPairs[shadowI])
            {
                if (locSlave[shadowI])
                {
                    // Correct normal and friction contact models for the
                    // current contact pair

                    // Calculate the slave patch face unit normals as they are
                    // units by both the normal and friction models
                    const vectorField shadowPatchFaceNormals =
                        patchField
                        (
                            shadowPatchIndices()[shadowI],
                            shadowZoneIndices()[shadowI],
                            shadowZone(shadowI).faceNormals()
                        );

                    // Interpolate the master displacement increment to the
                    // slave patch as it is required by specific normal and
                    // friction contact models

                    vectorField patchDD(patch().size(), vector::zero);
                    vectorField shadowPatchDD
                    (
                        patch().boundaryMesh()
                        [
                            shadowPatchIndices()[shadowI]
                        ].size(),
                        vector::zero
                    );

                    if (movingMesh())
                    {
                        // Updated Lagrangian, we will directly lookup the
                        // displacement increment

                        const volVectorField& DD =
                            db().lookupObject<volVectorField>("DU");

                        patchDD = DD.boundaryField()[patch().index()];
                        shadowPatchDD =
                            DD.boundaryField()[shadowPatchIndices()[shadowI]];
                    }
                    else
                    {
                        // We will lookup the total displacement and old total
                        // displacement

                        const volVectorField& D =
                            db().lookupObject<volVectorField>("U");

                        patchDD =
                            D.boundaryField()[patch().index()]
                          - D.oldTime().boundaryField()[patch().index()];
                        shadowPatchDD =
                            D.boundaryField()[shadowPatchIndices()[shadowI]]
                          - D.oldTime().boundaryField()
                            [
                                shadowPatchIndices()[shadowI]
                            ];
                    }

                    // Master zone DD
                    const vectorField zoneDD =
                        zoneField
                        (
                            zoneIndex(),
                            patch().index(),
                            patchDD
                        );

                    // Master patch DD interpolated to the slave patch
                    const vectorField patchDDInterpToShadowPatch =
                        patchField
                        (
                            shadowPatchIndices()[shadowI],
                            shadowZoneIndices()[shadowI],
                            zoneToZone(shadowI).masterToSlave(zoneDD)()
                        );

                    FatalError
                        << "Disabled: use jasakSolidContact" << abort(FatalError);
                    // normalModel(shadowI).correct
                    // (
                    //     shadowPatchFaceNormals,
                    //     zoneToZone(shadowI),
                    //     shadowPatchDD,
                    //     patchDDInterpToShadowPatch
                    // );

                    frictionModel(shadowI).correct
                    (
                        normalModel(shadowI).slavePressure(),
                        shadowPatchFaceNormals,
                        normalModel(shadowI).areaInContact(),
                        shadowPatchDD,
                        patchDDInterpToShadowPatch
                    );

                    // Accumulate traction

                    curPatchTractions(shadowI) =
                        frictionModel(shadowI).slaveTraction()
                        + normalModel(shadowI).slavePressure();

                    curPatchTraction += curPatchTractions(shadowI);
                }
                else // local master
                {
                    // Get traction from local slave

                    const volVectorField& field =
                        db().lookupObject<volVectorField>
                        (
                            dimensionedInternalField().name()
                        );

                    const solidGeneralContactFvPatchVectorField&
                        localMasterField =
                        refCast<const solidGeneralContactFvPatchVectorField>
                        (
                            field.boundaryField()
                            [
                                shadowPatchIndices()[shadowI]
                            ]
                        );

                    const label masterShadowI =
                        localMasterField.findShadowID(patch().index());

                    vectorField shadowPatchTraction =
                        -localMasterField.frictionModel
                        (
                            masterShadowI
                        ).slaveTractionForMaster()
                        -localMasterField.normalModel
                        (
                            masterShadowI
                        ).slavePressure();

                    vectorField shadowZoneTraction =
                        zoneField
                        (
                            shadowZoneIndices()[shadowI],
                            shadowPatchIndices()[shadowI],
                            shadowPatchTraction
                        );

                    // Face-to-face
                    vectorField masterZoneTraction =
                        localMasterField.zoneToZone
                        (
                            masterShadowI
                        ).slaveToMaster(shadowZoneTraction);

                    // We store master patch traction as thermalGeneralContact
                    // uses it
                    curPatchTractions(shadowI) =
                        patchField
                        (
                            patch().index(),
                            zoneIndex(),
                            masterZoneTraction
                        );

                    curPatchTraction += curPatchTractions(shadowI);
                }
            } // if contact pair is active
        } // forAll contact pairs

        // Set master gradient based on accumulated traction
        traction() = curPatchTraction;
    }

    solidTractionFvPatchVectorField::updateCoeffs();
}


const Foam::scalarField& Foam::solidGeneralContactFvPatchVectorField::Qc() const
{
    if (!QcPtr_)
    {
        calcQc();
    }

    return *QcPtr_;
}


//- Increment of dissipated energy due to friction for each pair
const Foam::scalarField& Foam::solidGeneralContactFvPatchVectorField::Qc
(
    const label shadowI
) const
{
    if (!QcsPtr_)
    {
        calcQcs();
    }

    return (*QcsPtr_)[shadowI];
}


void Foam::solidGeneralContactFvPatchVectorField::calcQcs() const
{
    const boolList& locSlave = localSlave();

    QcsPtr_ = new List<scalarField>(locSlave.size());

    bool sigmaCauchyFound =
        db().foundObject<volSymmTensorField>
        (
            "sigmaCauchy"
        );

    const scalar deltaT = patch().boundaryMesh().mesh().time().deltaT().value();

    const volVectorField& field =
        db().lookupObject<volVectorField>
        (
            dimensionedInternalField().name()
        );

    forAll(locSlave, shadowI)
    {
        scalarField& Qc = (*QcsPtr_)[shadowI];
        Qc.setSize(patch().size(), 0.0);

        // For now, we assume traction is constant over time-step
        // Todo: use trapezoidal rule
        vectorField curTraction(Qc.size(), vector::zero);

        // sigma/sigmaCauchy is up-to-date as Qc is called after momentum loop
        // has converged and sigma has been updated and mesh moved
        if (sigmaCauchyFound)
        {
            const symmTensorField& sigma =
                db().lookupObject<volSymmTensorField>
                (
                    "sigmaCauchy"
                ).boundaryField()[patch().index()];

            curTraction = patch().nf() & sigma;
        }
        else
        {
            const symmTensorField& sigma =
            db().lookupObject<volSymmTensorField>
            (
                "sigma"
            ).boundaryField()[patch().index()];

            curTraction = patch().nf() & sigma;
        }

        // Calculate Qc for shadowI

        vectorField curPatchSlip(Qc.size(), vector::zero);

        // Calculate slip
        if (locSlave[shadowI])
        {
            curPatchSlip = frictionModel(shadowI).slip();
        }
        else
        {
            const solidGeneralContactFvPatchVectorField& shadowPatchField =
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()[shadowPatchIndices()[shadowI]]
            );

            const label locShadowID =
            shadowPatchField.findShadowID(patch().index());

            vectorField shadowPatchSlip =
            shadowPatchField.frictionModel(locShadowID).slip();

            vectorField shadowZoneSlip =
            zoneField
            (
                shadowZoneIndices()[shadowI],
                shadowPatchIndices()[shadowI],
                shadowPatchSlip
            );

            // Interpolate from shadow to the current patch
            // Face-to-face

            vectorField curZoneSlip =
            shadowPatchField.zoneToZone(locShadowID).slaveToMaster
            (
                shadowZoneSlip
            );

            curPatchSlip =
            patchField
            (
                patch().index(),
                zoneIndex(),
                curZoneSlip
            );
        }

        // Heat flux rate: rate of dissipated frictional energy
        // The dot product of the traction vectors and the slip vectors
        // gives the dissipated frictional energy rate per unit area, which
        // is always positive
        Qc = mag(curTraction & (curPatchSlip/deltaT));
    }
}


void Foam::solidGeneralContactFvPatchVectorField::calcBbOffset() const
{
    if (bbOffset_ != 0)
    {
        FatalErrorIn
        (
            "Foam::scalar Foam::solidGeneralContactFvPatchVectorField::"
            "calcBbOffset() const"
        )   << "already set" << abort(FatalError);
    }

    // We will set the BB offset to five times the average dimension of the
    // smallest face on the zone

    scalar minDim = GREAT;

    if (patch().size() > 0)
    {
        minDim = min(sqrt(patch().magSf()));
    }

    bbOffset_ = 5.0*returnReduce(minDim, minOp<scalar>());

    if (debug)
    {
        Info<< nl << "The bbOffset is " << bbOffset_ << endl;
    }
}



Foam::scalar Foam::solidGeneralContactFvPatchVectorField::bbOffset() const
{
    if (bbOffset_ == 0)
    {
        calcBbOffset();
    }

    return bbOffset_;
}


const Foam::vectorField&
Foam::solidGeneralContactFvPatchVectorField::curPatchTractions
(
    const label shadowI
) const
{
    if (!curPatchTractionPtr_)
    {
        makeCurPatchTractions();
    }

    return (*curPatchTractionPtr_)[shadowI];
}


Foam::vectorField&
Foam::solidGeneralContactFvPatchVectorField::curPatchTractions
(
    const label shadowI
)
{
    if (!curPatchTractionPtr_)
    {
        makeCurPatchTractions();
    }

    return (*curPatchTractionPtr_)[shadowI];
}


// Write
void Foam::solidGeneralContactFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);

    os.writeKeyword("rigidMaster")
        << rigidMaster_ << token::END_STATEMENT << nl;

    // Write the dict from the first contact model

    const label shadowI = 0;

    if (localSlave()[shadowI])
    {
        os.writeKeyword("normalContactModel")
            << normalModel(shadowI).type() << token::END_STATEMENT << nl;
        normalModel(shadowI).writeDict(os);

        os.writeKeyword("frictionContactModel")
            << frictionModel(shadowI).type() << token::END_STATEMENT << nl;
        frictionModel(shadowI).writeDict(os);
    }
    else
    {
        const volVectorField& field =
            db().lookupObject<volVectorField>
            (
                dimensionedInternalField().name()
            );

        const solidGeneralContactFvPatchVectorField& localSlaveField =
            refCast<const solidGeneralContactFvPatchVectorField>
            (
                field.boundaryField()
                [
                    shadowPatchIndices()[shadowI]
                ]
            );

        const label localSlaveID =
            localSlaveField.findShadowID(patch().index());

        os.writeKeyword("normalContactModel")
            << localSlaveField.normalModel(localSlaveID).type()
            << token::END_STATEMENT << nl;
        localSlaveField.normalModel(localSlaveID).writeDict(os);

        os.writeKeyword("frictionContactModel")
            << localSlaveField.frictionModel(localSlaveID).type()
            << token::END_STATEMENT << nl;
        localSlaveField.frictionModel(localSlaveID).writeDict(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField, solidGeneralContactFvPatchVectorField
    )
}


// ************************************************************************* //
