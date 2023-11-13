/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "newThermalContactFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "polyPatchID.H"
#include "ZoneIDs.H"

#include "solidContactFvPatchVectorField.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::newThermalContactFvPatchScalarField::shadowPatchIndex() const
{
    if (shadowPatchIndex_ == -1)
    {
        setShadowPatchIndex();
    }

    return shadowPatchIndex_;
}


void Foam::newThermalContactFvPatchScalarField::setShadowPatchIndex() const
{
    if (shadowPatchIndex_ != -1)
    {
        FatalErrorIn
        (
            "void Foam::newThermalContactFvPatchScalarField::"
            "setShadowPatchIndex() const"
        ) << "index already set!" << abort(FatalError);
    }

    // Shadow patch index
    polyPatchID shadow(shadowPatchName_, patch().patch().boundaryMesh());

    if (!shadow.active())
    {
        FatalErrorIn("newThermalContactFvPatchScalarField")
            << "Shadow patch name " << shadowPatchName_ << " not found."
            << abort(FatalError);
    }

    shadowPatchIndex_ = shadow.index();
}


Foam::label Foam::newThermalContactFvPatchScalarField::zoneIndex() const
{
    if (zoneIndex_ == -1)
    {
        setZoneIndex();
    }

    return zoneIndex_;
}


void Foam::newThermalContactFvPatchScalarField::setZoneIndex() const
{
    if (zoneIndex_ != -1)
    {
        FatalErrorIn
        (
            "void Foam::newThermalContactFvPatchScalarField::setZoneIndex() const"
        ) << "index already set!" << abort(FatalError);
    }

    // Zone index
    word zoneName = patch().name() + "FaceZone";

    faceZoneID zone(zoneName, patch().boundaryMesh().mesh().faceZones());

    if (!zone.active())
    {
        FatalErrorIn("newThermalContactFvPatchScalarField")
            << "Face zone name " << zoneName
            << " not found.  Please check your zone definition."
            << abort(FatalError);
    }

    zoneIndex_ = zone.index();
}


Foam::label Foam::newThermalContactFvPatchScalarField::shadowZoneIndex() const
{
    if (shadowZoneIndex_ == -1)
    {
        setShadowZoneIndex();
    }

    return shadowZoneIndex_;
}


void Foam::newThermalContactFvPatchScalarField::setShadowZoneIndex() const
{
    if (shadowZoneIndex_ != -1)
    {
        FatalErrorIn
        (
            "void Foam::thermalContactFvZoneScalarField::"
            "setShadowZoneIndex() const"
        ) << "index already set!" << abort(FatalError);
    }

    // Shadow zone index
    word shadowZoneName = shadowPatchName_ + "FaceZone";

    faceZoneID shadowZone
    (
        shadowZoneName,
        patch().boundaryMesh().mesh().faceZones()
    );

    if (!shadowZone.active())
    {
        FatalErrorIn("newThermalContactFvPatchScalarField")
            << "Face zone name " << shadowZoneName
            << " not found.  Please check your zone definition."
            << abort(FatalError);
    }

    shadowZoneIndex_ = shadowZone.index();
}


void Foam::newThermalContactFvPatchScalarField::calcZoneToZone() const
{
    // Create zone-to-zone interpolation
    if (zoneToZonePtr_)
    {
        FatalErrorIn
        (
            "void newThermalContactFvPatchScalarField::calcZoneToZone() const"
        )
            << "Zone to zone interpolation already calculated"
            << abort(FatalError);
    }


    // Check master and slave patch
    const volScalarField& field =
        this->db().objectRegistry::lookupObject<volScalarField>
        (
            this->dimensionedInternalField().name()
        );

    const newThermalContactFvPatchScalarField& shadowPatchField =
        refCast<const newThermalContactFvPatchScalarField>
        (
            field.boundaryField()[shadowPatchIndex()]
        );

    if (master_)
    {
        if (shadowPatchField.master() == true)
        {
            FatalErrorIn("newThermalContactFvPatchScalarField")
                << "There are two master patches"
                    << abort(FatalError);
        }
    }
    else
    {
        if (shadowPatchField.master() == false)
        {
            FatalErrorIn("newThermalContactFvPatchScalarField")
                << "There is no master patch"
                    << abort(FatalError);
        }
    }

    if (master())
    {
        // Create interpolation for patches
        zoneToZonePtr_ =
            // new ggiZoneInterpolation
            new newGgiZoneInterpolation
            (
                patch().boundaryMesh().mesh().faceZones()[zoneIndex()](),
                patch().boundaryMesh().mesh().faceZones()[shadowZoneIndex()](),
                tensorField(0),
                tensorField(0),
                vectorField(0), // Slave-to-master separation. Bug fix
                true,           // global data
                0,              // Non-overlapping face tolerances
                0,              // HJ, 24/Oct/2008
                true,           // Rescale weighting factors.  Bug fix, MB.
                newGgiInterpolation::AABB
                //newGgiInterpolation::BB_OCTREE
            );
    }
    else
    {
        FatalErrorIn
        (
            "void newThermalContactFvPatchScalarField::calcZoneToZone() const"
        )
            << "Attempting to create GGIInterpolation on a slave"
            << abort(FatalError);
    }
}


void Foam::newThermalContactFvPatchScalarField::calcContact() const
{
    // Calculating contact based on the contact area from normalContact
    // model, i.e. there is contact if contact area is greater then 0.
    // Storing contact area into contact variable makes sense, since face can be
    // in a partial contact with other face, i.e. part of the face can be
    // uncovered.

    // Create contact indicator field
    if (contactPtr_)
    {
        FatalErrorIn
        (
            "void newThermalContactFvPatchScalarField::calcContact() const"
        )
            << "Contact indicator interpolation already calculated"
            << abort(FatalError);
    }

    contactPtr_ = new scalarField(patch().size(), 0);
    scalarField& contact = *contactPtr_;

    // Calculate current deformed normals
    const fvPatchField<symmTensor>& sigmaCauchy =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaCauchy");

    // Relative deformation gradient inverse
    const fvPatchField<tensor>& Finv =
        patch().lookupPatchField<volTensorField, tensor>("relFinv");

    // Calculate current deformed face normals
    // nCurrent should be unit vectors but we will force normalisation
    // to remove any errors
    const vectorField n = patch().nf();
    vectorField nCurrent = Finv.T() & n;
    nCurrent /= mag(nCurrent);

    // Calcuate contact pressure
    scalarField contactPressure = -nCurrent & (nCurrent & sigmaCauchy);

    // Careful: normalModels()[shadowZoneIndex()] is not allowed:
    // shadowZoneIndex is not related to the normalModels list
    // Lookup areaInContact from normalContactModel
    // const volVectorField& dispField =
    //     this->db().objectRegistry::lookupObject<volVectorField>
    //     (
    //         DUName_
    //     );

    // if (!master())
    // {
    //     const solidContactFvPatchVectorField& DUpatch =
    //         refCast<const solidContactFvPatchVectorField>
    //         (
    //             dispField.boundaryField()[patch().index()]
    //         );

    //     const scalarField& areaInContact =
    //         DUpatch.normalModels()[zoneIndex()].areaInContact();

    //     contact = areaInContact;
    // }
    // else
    // {
    //     const solidContactFvPatchVectorField& DUshadowPatch =
    //         refCast<const solidContactFvPatchVectorField>
    //         (
    //             dispField.boundaryField()[shadowPatchIndex()]
    //         );
    //     const scalarField& areaInContactShadowPatch =
    //         DUshadowPatch.normalModels()[shadowZoneIndex()].areaInContact();

    //     scalarField curZoneAreaInContact =
    //         zoneField
    //         (
    //             shadowZoneIndex(),
    //             shadowPatchIndex(),
    //             areaInContactShadowPatch
    //         );

    //     contact =
    //         patchField
    //         (
    //             patch().index(),
    //             zoneIndex(),
    //             zoneToZone().slaveToMaster(curZoneAreaInContact)()
    //         );
    // }
}


void Foam::newThermalContactFvPatchScalarField::calcHc() const
{
    // Contact conductance dependent on contact pressure
    HcPtr_ = new scalarField(patch().size(), 0.0);
    scalarField& h = *HcPtr_;

    //scalarField contactPressure(patch().size(), 0.0);

    // Pressure dependent contact conductance is expressed as:
    // h = hRef*(p/H)^beta
    // where:
    // p = contact pressure [Pa]
    // H = Vicker's hardness of softer material [Pa]
    // hRef = reference contact conductance
    // Rc = reference contact resistance (Rc = 1/hRef)
    // beta = experimentally determined coefficient describing sensitivity of
    //        contact conductance on pressure
    //
    //        beta > 1 very sensitive
    //        beta < 0.01 very insensitive

    // Lookup slavePressure from normalContactModel
    // const volVectorField& dispField =
    //     this->db().objectRegistry::lookupObject<volVectorField>
    //     (
    //         DUName_
    //     );

    // Calculate current deformed normals
    const fvPatchField<symmTensor>& sigmaCauchy =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaCauchy");

    // Relative deformation gradient inverse
    const fvPatchField<tensor>& Finv =
        patch().lookupPatchField<volTensorField, tensor>("relFinv");

    // Calculate current deformed face normals
    // nCurrent should be unit vectors but we will force normalisation
    // to remove any errors
    const vectorField n = patch().nf();
    vectorField nCurrent = Finv.T() & n;
    nCurrent /= mag(nCurrent);

    // Calcuate contact pressure
    scalarField contactPressure = -nCurrent & (nCurrent & sigmaCauchy);

    // Limit contact pressure to avoid FPE in pow
    contactPressure = max(contactPressure, SMALL);

    // Be careful: zoneIndex and shadowZoneIndex are not related to the
    // addressing of the normalModels i.e. normalModels()[shadowZoneIndex()]
    // does not make sense.
    // The idea of directly looking up the contactPressure field makes sense;
    // however to do it correctly we would need to accumulated the slavePressure
    // field for each contact; remember the master stores the normalModels and
    // there can be multiple; for now, let's revert to the previous method.
    // if (!master())
    // {
    //     const solidContactFvPatchVectorField& DUpatch =
    //         refCast<const solidContactFvPatchVectorField>
    //         (
    //             dispField.boundaryField()[patch().index()]
    //         );

    //     const vectorField& slavePressure =
    //         DUpatch.normalModels()[zoneIndex()].slavePressure();

    //     contactPressure = mag(slavePressure);
    // }
    // else
    // {
    //     const solidContactFvPatchVectorField& DUshadowPatch =
    //         refCast<const solidContactFvPatchVectorField>
    //         (
    //             dispField.boundaryField()[shadowPatchIndex()]
    //         );
    //     const scalarField magSlavePressureShadowPatch =
    //         mag
    //         (
    //          DUshadowPatch.normalModels()[shadowZoneIndex()].slavePressure()
    //         );

    //     scalarField curZoneMagSlavePressure =
    //         zoneField
    //         (
    //             shadowZoneIndex(),
    //             shadowPatchIndex(),
    //             magSlavePressureShadowPatch
    //         );

    //     contactPressure =
    //         patchField
    //         (
    //             patch().index(),
    //             zoneIndex(),
    //             zoneToZone().slaveToMaster(curZoneMagSlavePressure)()
    //         );
    // }


    // Vickers hardness is defined by:
    // H ~ 0.3*sigmaU, where:
    // H [kgf/mm2] - Vickers hardness
    // sigmaU [MPa] - ultimate tensile strength
    //
    // In order to calculate hardness in [MPa], the Vickers hardness in
    // [kgf/mm2] should be multiplied by gravitational constant (9.80665).
    // The final formula which gives Vickers hardness in [Pa] is:
    // H[Pa] = 0.3*sigmaU[Pa]*9.80665;
    //h = pow(contactPressure/(0.3*UTS()*9.80665), beta())/Rc();
    // Philip Note: in the Bathe paper, they are scaling "h" by
    // "contactPressure/H" where they say H is the hardness
    // We know that this scaling parameter should be dimensionless, which means
    // that H should be in Pa (if contactPressure is in Pa)
    // So then then question is how do we get a measure of hardness in Pa?
    // A natural estimate of the hardness of material in Pa is the yield
    // strength (or ultimate tensile strength). So it seems to me that the
    // intention by Bathe for this relation was that "H" is a measure of the
    // strength of the material e.g sigmaY or UTS or in my case I used 0.3*UTS
    // as a measure of Vicker's hardness in Pa
    // Vanja, I understand the unit conversions you used but in effect this
    // means that H is approximately equal to 10*UTS which I do not beleive is
    // in keeping with the intention of the formulation. Either way, this is
    // just an emprical scaling factor so it may not be critical
    h = pow(contactPressure/UTS(), beta())/Rc();

    if (debug)
    {
        Info << "contact conductance:" << nl
             << "    max: " << max(h) << nl
             << "    min: " << min(h) << nl
             << "contact resistance: " << nl
             << "    max: " << 1.0/(max(h) + SMALL) << nl
             << "    min: " << 1.0/(min(h) + SMALL) << endl;
    }
}


const Foam::scalarField&
Foam::newThermalContactFvPatchScalarField::curPatchK() const
{
    // Return thermal conductivity of a current patch
    const volScalarField& k =
        this->db().objectRegistry::lookupObject<volScalarField>
        (
            "k"
        );

    return k.boundaryField()[patch().index()];
}


const Foam::tmp<Foam::scalarField>
Foam::newThermalContactFvPatchScalarField::
shadowPatchKOnCurPatch() const
{
    // Return interpolated thermal conductivity of a shadow patch onto current
    // patch
    tmp<scalarField> tShadowPatchKOnCurPatch
    (
        new scalarField(patch().size(), 0.0)
    );
    scalarField& shadowPatchKOnCurPatch = tShadowPatchKOnCurPatch();

    const volScalarField& k =
        this->db().objectRegistry::lookupObject<volScalarField>
        (
            "k"
        );
    const scalarField& shadowPatchK = k.boundaryField()[shadowPatchIndex()];

    scalarField shadowZoneK =
        zoneField(shadowZoneIndex(), shadowPatchIndex(), shadowPatchK);

    scalarField shadowZoneKOnCurPatch(patch().size(), 0.0);

    if (master())
    {
        shadowZoneKOnCurPatch = zoneToZone().slaveToMaster(shadowZoneK);
    }
    else
    {
        shadowZoneKOnCurPatch = zoneToZone().masterToSlave(shadowZoneK);
    }

    shadowPatchKOnCurPatch =
        patchField(patch().index(), zoneIndex(), shadowZoneKOnCurPatch);

    return tShadowPatchKOnCurPatch;
}


const Foam::scalarField&
Foam::newThermalContactFvPatchScalarField::curPatchRhoC() const
{
    // Return rho*C of a current patch
    const volScalarField& rhoC =
        this->db().objectRegistry::lookupObject<volScalarField>("(rho*C)");

    return rhoC.boundaryField()[patch().index()];
}


const Foam::tmp<Foam::scalarField>
Foam::newThermalContactFvPatchScalarField::
shadowPatchRhoCOnCurPatch() const
{
    // Return interpolated rho*C of a shadow patch onto current patch
    tmp<scalarField> tShadowPatchRhoCOnCurPatch
    (
        new scalarField(patch().size(), 0.0)
    );
    scalarField& shadowPatchRhoCOnCurPatch = tShadowPatchRhoCOnCurPatch();

    const volScalarField& rhoC =
        this->db().objectRegistry::lookupObject<volScalarField>("(rho*C)");

    const scalarField& shadowPatchRhoC =
        rhoC.boundaryField()[shadowPatchIndex()];

    scalarField shadowZoneRhoC =
        zoneField(shadowZoneIndex(), shadowPatchIndex(), shadowPatchRhoC);

    scalarField shadowZoneRhoCOnCurPatch(patch().size(), 0.0);

    if (master())
    {
        shadowZoneRhoCOnCurPatch = zoneToZone().slaveToMaster(shadowZoneRhoC);
    }
    else
    {
        shadowZoneRhoCOnCurPatch = zoneToZone().masterToSlave(shadowZoneRhoC);
    }

    shadowPatchRhoCOnCurPatch =
        patchField(patch().index(), zoneIndex(), shadowZoneRhoCOnCurPatch);

    return tShadowPatchRhoCOnCurPatch;
}


const Foam::tmp<Foam::scalarField>
Foam::newThermalContactFvPatchScalarField::
shadowPatchTempOnCurPatch() const
{
    // Return interpolated temperature field of a shadow patch onto current
    // patch
    tmp<scalarField> tShadowPatchTempOnCurPatch
    (
        new scalarField(patch().size(), 0.0)
    );
    scalarField& shadowPatchTempOnCurPatch = tShadowPatchTempOnCurPatch();

    const volScalarField& field =
        this->db().objectRegistry::lookupObject<volScalarField>
        (
            this->dimensionedInternalField().name()
        );

    const newThermalContactFvPatchScalarField& shadowPatchField =
        refCast<const newThermalContactFvPatchScalarField>
        (
            field.boundaryField()[shadowPatchIndex()]
        );

    scalarField shadowZoneTemp =
        zoneField(shadowZoneIndex(), shadowPatchIndex(), shadowPatchField);

    scalarField shadowZoneTempOnCurPatch(patch().size(), 0.0);

    if (master())
    {
        shadowZoneTempOnCurPatch = zoneToZone().slaveToMaster(shadowZoneTemp);
    }
    else
    {
        shadowZoneTempOnCurPatch = zoneToZone().masterToSlave(shadowZoneTemp);
    }

    shadowPatchTempOnCurPatch =
        patchField(patch().index(),zoneIndex(), shadowZoneTempOnCurPatch);

    return tShadowPatchTempOnCurPatch;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::newThermalContactFvPatchScalarField::
newThermalContactFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    master_(false),
    shadowPatchName_("undefined"),
    shadowPatchIndex_(-1),
    zoneIndex_(-1),
    shadowZoneIndex_(-1),
    underRelaxation_(1),
    alpha_(0),
    Tinf_(0),
    beta_(0),
    Rc_(0),
    UTS_(0),
    Qc_(p.size(), 0),
    useFrictionQc_(false),
    DUName_("undefined"),
    curTimeIndex_(-1),
    zoneToZonePtr_(NULL),
    HcPtr_(NULL),
    contactPtr_(NULL)
{}


Foam::newThermalContactFvPatchScalarField::
newThermalContactFvPatchScalarField
(
    const newThermalContactFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    master_(ptf.master_),
    shadowPatchName_(ptf.shadowPatchName_),
    shadowPatchIndex_(ptf.shadowPatchIndex_),
    zoneIndex_(ptf.zoneIndex_),
    shadowZoneIndex_(ptf.shadowZoneIndex_),
    underRelaxation_(ptf.underRelaxation_),
    alpha_(ptf.alpha_),
    Tinf_(ptf.Tinf_),
    beta_(ptf.beta_),
    Rc_(ptf.Rc_),
    UTS_(ptf.UTS_),
    Qc_(ptf.Qc_, mapper),
    useFrictionQc_(ptf.useFrictionQc_),
    DUName_(ptf.DUName_),
    curTimeIndex_(ptf.curTimeIndex_),
    zoneToZonePtr_(NULL),
    HcPtr_(NULL),
    contactPtr_(NULL)
{
    // copy pointer objects

    // if (ptf.zoneToZonePtr_)
    // {
    //     // zoneToZone can be recalculated when needed
    // }

    if (ptf.HcPtr_)
    {
        HcPtr_ = new scalarField(*ptf.HcPtr_);
    }

    if (ptf.contactPtr_)
    {
        contactPtr_ = new scalarField(*ptf.contactPtr_);
    }
}


Foam::newThermalContactFvPatchScalarField::
newThermalContactFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    master_(dict.lookupOrDefault<Switch>("master", false)),
    shadowPatchName_(dict.lookup("shadowPatch")),
    shadowPatchIndex_(-1),
    zoneIndex_(-1),
    shadowZoneIndex_(-1),
    underRelaxation_(1),
    alpha_(readScalar(dict.lookup("alpha"))),
    Tinf_(readScalar(dict.lookup("Tinf"))),
    beta_(0.0),
    Rc_(0.0),
    UTS_(0.0),
    Qc_(),
    useFrictionQc_(true),
    DUName_(dict.lookupOrDefault<word>("DUName", "DU")),
    curTimeIndex_(-1),
    zoneToZonePtr_(NULL),
    HcPtr_(NULL),
    contactPtr_(NULL)
{
    // Read only on master
    if (master())
    {
        underRelaxation_ = readScalar(dict.lookup("underRelaxation"));

        Rc_ = readScalar(dict.lookup("Rc"));

        beta_ = dict.lookupOrDefault<scalar>("beta", 0.1);

        UTS_ = dict.lookupOrDefault<scalar>("UTS", 2e9);

        if (Rc_ < SMALL)
        {
            Warning
                << "Contact conductivity resistance cannot be exactly zero"
                << nl << "Rc = " << SMALL << " is assumed" << endl;

            Rc_ = SMALL;
        }

        if (dict.found("Qc"))
        {
            Qc_ = scalarField("Qc", dict, p.size());

            // if Qc is found and is larger than zero then we will use this
            // instead of frictional energy
            if (max(Qc_) > SMALL)
            {
                Warning
                    << "As Qc is found in the thermal boundary and is greater "
                    << "then zero, it will be used instead of frictional energy"
                    << endl;
                useFrictionQc_ = false;
            }
        }
    }

    if (dict.found("gradient"))
    {
        this->gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        this->gradient() = 0.0;
    }

    if (dict.found("value"))
    {
        Field<scalar>::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        Field<scalar>::operator=
        (
            this->patchInternalField()
          + this->gradient()/this->patch().deltaCoeffs()
        );
    }
}


Foam::newThermalContactFvPatchScalarField::
newThermalContactFvPatchScalarField
(
    const newThermalContactFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    master_(tppsf.master_),
    shadowPatchName_(tppsf.shadowPatchName_),
    shadowPatchIndex_(tppsf.shadowPatchIndex_),
    zoneIndex_(tppsf.zoneIndex_),
    shadowZoneIndex_(tppsf.shadowZoneIndex_),
    underRelaxation_(tppsf.underRelaxation_),
    alpha_(tppsf.alpha_),
    Tinf_(tppsf.Tinf_),
    beta_(tppsf.beta_),
    Rc_(tppsf.Rc_),
    UTS_(tppsf.UTS_),
    Qc_(tppsf.Qc_),
    useFrictionQc_(tppsf.useFrictionQc_),
    DUName_(tppsf.DUName_),
    curTimeIndex_(tppsf.curTimeIndex_),
    zoneToZonePtr_(NULL),
    HcPtr_(NULL),
    contactPtr_(NULL)
{
    // copy pointer objects

    // if (tppsf.zoneToZonePtr_)
    // {
    //     // zoneToZone can be recalculated when needed
    // }

    if (tppsf.HcPtr_)
    {
        HcPtr_ = new scalarField(*tppsf.HcPtr_);
    }

    if (tppsf.contactPtr_)
    {
        contactPtr_ = new scalarField(*tppsf.contactPtr_);
    }
}


Foam::newThermalContactFvPatchScalarField::
newThermalContactFvPatchScalarField
(
    const newThermalContactFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    master_(tppsf.master_),
    shadowPatchName_(tppsf.shadowPatchName_),
    shadowPatchIndex_(tppsf.shadowPatchIndex_),
    zoneIndex_(tppsf.zoneIndex_),
    shadowZoneIndex_(tppsf.shadowZoneIndex_),
    underRelaxation_(tppsf.underRelaxation_),
    alpha_(tppsf.alpha_),
    Tinf_(tppsf.Tinf_),
    beta_(tppsf.beta_),
    Rc_(tppsf.Rc_),
    UTS_(tppsf.UTS_),
    Qc_(tppsf.Qc_),
    useFrictionQc_(tppsf.useFrictionQc_),
    DUName_(tppsf.DUName_),
    curTimeIndex_(tppsf.curTimeIndex_),
    zoneToZonePtr_(NULL),
    HcPtr_(NULL),
    contactPtr_(NULL)
{
    // copy pointer objects

    // if (tppsf.zoneToZonePtr_)
    // {
    //     // zoneToZone can be recalculated when needed
    // }

    if (tppsf.HcPtr_)
    {
        HcPtr_ = new scalarField(*tppsf.HcPtr_);
    }

    if (tppsf.contactPtr_)
    {
        contactPtr_ = new scalarField(*tppsf.contactPtr_);
    }
}



// * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * * //

Foam::newThermalContactFvPatchScalarField::
~newThermalContactFvPatchScalarField()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::newThermalContactFvPatchScalarField::underRelaxation() const
{
    if (master())
    {
        return underRelaxation_;
    }
    else
    {
        const volScalarField& field =
            this->db().objectRegistry::lookupObject<volScalarField>
            (
                this->dimensionedInternalField().name()
            );

        const newThermalContactFvPatchScalarField& shadowPatchField =
            refCast<const newThermalContactFvPatchScalarField>
            (
                field.boundaryField()[shadowPatchIndex()]
            );

        return shadowPatchField.underRelaxation();
    }
}


Foam::scalar Foam::newThermalContactFvPatchScalarField::Rc() const
{
    if (master())
    {
        return Rc_;
    }
    else
    {
        const volScalarField& field =
            this->db().objectRegistry::lookupObject<volScalarField>
            (
                this->dimensionedInternalField().name()
            );

        const newThermalContactFvPatchScalarField& shadowPatchField =
            refCast<const newThermalContactFvPatchScalarField>
            (
                field.boundaryField()[shadowPatchIndex()]
            );

        return shadowPatchField.Rc();
    }
}


Foam::scalar Foam::newThermalContactFvPatchScalarField::beta() const
{
    if (master())
    {
        return beta_;
    }
    else
    {
        const volScalarField& field =
            this->db().objectRegistry::lookupObject<volScalarField>
            (
                this->dimensionedInternalField().name()
            );

        const newThermalContactFvPatchScalarField& shadowPatchField =
            refCast<const newThermalContactFvPatchScalarField>
            (
                field.boundaryField()[shadowPatchIndex()]
            );

        return shadowPatchField.beta();
    }
}


Foam::scalar Foam::newThermalContactFvPatchScalarField::UTS() const
{
    if (master())
    {
        return UTS_;
    }
    else
    {
        const volScalarField& field =
            this->db().objectRegistry::lookupObject<volScalarField>
            (
                this->dimensionedInternalField().name()
            );

        const newThermalContactFvPatchScalarField& shadowPatchField =
            refCast<const newThermalContactFvPatchScalarField>
            (
                field.boundaryField()[shadowPatchIndex()]
            );

        return shadowPatchField.UTS();
    }
}


Foam::tmp<Foam::scalarField>
Foam::newThermalContactFvPatchScalarField::Hc() const
{
    if (!HcPtr_)
    {
        calcHc();
    }

    return *HcPtr_;
}


Foam::tmp<Foam::scalarField>
Foam::newThermalContactFvPatchScalarField::Qc() const
{
    tmp<scalarField> tQc(new scalarField());

    // Initialise frixFlux field
    scalarField fricFlux(patch().size(), 0.0);


    // If displacement field exists and if frictional energy should be
    // used for frictional flux, then read it from the normal contact model
    if
    (
        this->db().objectRegistry::foundObject<volVectorField>(DUName_)
        && useFrictionQc()
    )
    {
        // Lookup displacement field
        const volVectorField& DU =
            this->db().objectRegistry::lookupObject<volVectorField>
            (
                DUName_
            );

        // Now, check if the current patch is master or slave
        if (master())
        {
            // There may be multiple contacts so we will directly lookup master
            // DU
            const label masterID = patch().index();
            if
            (
                DU.boundaryField()[masterID].type()
                != solidContactFvPatchVectorField::typeName
            )
            {
                FatalErrorIn("newThermalContactFvPatchScalarField::Qc()")
                    << "DU patch " << patch().boundaryMesh()[masterID].name()
                    << " should be of type solidContact "
                    << abort(FatalError);
            }
            else
            {
                const solidContactFvPatchVectorField& contactDU =
                    refCast<const solidContactFvPatchVectorField>
                    (
                        DU.boundaryField()[masterID]
                    );

                // Lookup friction flux from normal contact model
                fricFlux = contactDU.Qc();
            }
        }
        else
        {
            // Lookup shadow patch (master DU)
            if
            (
                DU.boundaryField()[shadowPatchIndex()].type()
                != solidContactFvPatchVectorField::typeName
            )
            {
                FatalErrorIn("newThermalContactFvPatchScalarField::Qc()")
                    << "DU patch "
                    << patch().boundaryMesh()[shadowPatchIndex()].name()
                    << " should be of type solidContact "
                    << abort(FatalError);
            }
            else
            {
                const solidContactFvPatchVectorField& shadowPatchContactDU =
                    refCast<const solidContactFvPatchVectorField>
                    (
                        DU.boundaryField()[shadowPatchIndex()]
                    );

                // Lookup friction flux from normal contact model
                scalarField masterZoneQc =
                    zoneField
                    (
                        shadowZoneIndex(),
                        shadowPatchIndex(),
                        shadowPatchContactDU.Qc()()
                    );

                fricFlux =
                    patchField
                    (
                        patch().index(),
                        zoneIndex(),
                        zoneToZone().masterToSlave(masterZoneQc)()
                    );
            }
        }
    }
    else
    {
        if (master())
        {
            fricFlux = Qc_;
        }
        else
        {
            // Interpolate friction flux (Qc) from master to slave
            const volScalarField& field =
                this->db().objectRegistry::lookupObject<volScalarField>
                (
                    this->dimensionedInternalField().name()
                );

            const newThermalContactFvPatchScalarField& shadowPatchField =
                refCast<const newThermalContactFvPatchScalarField>
                (
                    field.boundaryField()[shadowPatchIndex()]
                );

            scalarField masterZoneQc =
                zoneField
                (
                    shadowZoneIndex(),
                    shadowPatchIndex(),
                    shadowPatchField.Qc()()
                );

            fricFlux =
                patchField
                (
                    patch().index(),
                    zoneIndex(),
                    zoneToZone().masterToSlave(masterZoneQc)()
                );
        }
    }

    // In the contact region:
    // masterFlux + slaveFlux - fricFlux = 0
    // And
    // masterFlux = hBar*(slaveTemp - masterTemp) + masterFricFlux
    // Where
    // masterFricFlux = fricFlux*(masterK/masterKsi)
    // masterKsi = sqrt(slaveRhoC/masterRhoC)*sqrt(masterK*slaveK) + masterK
    // K = conductivity
    // rho = density
    // C = specific heat capacity
    // fricFlux = lossCoeff*shearTraction*relativeVelocity

    // Calculate current patch flux
    scalarField curPatchKsi =
        sqrt
        (
            (shadowPatchRhoCOnCurPatch()()/curPatchRhoC())
           *curPatchK()*shadowPatchKOnCurPatch()()
        )
        + curPatchK();

    tQc() = fricFlux*(curPatchK()/curPatchKsi);

    return tQc;
}


const Foam::scalarField&
Foam::newThermalContactFvPatchScalarField::contact() const
{
    if (!contactPtr_)
    {
        calcContact();
    }

    return *contactPtr_;
}


bool Foam::newThermalContactFvPatchScalarField::useFrictionQc() const
{
    if (master())
    {
        return useFrictionQc_;
    }

    // If slave, then read it from master
    const volScalarField& field =
        this->db().objectRegistry::lookupObject<volScalarField>
        (
            this->dimensionedInternalField().name()
        );

    const newThermalContactFvPatchScalarField& shadowPatchField =
        refCast<const newThermalContactFvPatchScalarField>
        (
            field.boundaryField()[shadowPatchIndex()]
        );

    return shadowPatchField.useFrictionQc();
}


const Foam::newGgiZoneInterpolation&
Foam::newThermalContactFvPatchScalarField::zoneToZone() const
{
    if (master())
    {
        if (!zoneToZonePtr_)
        {
            word zoneName =
                patch().boundaryMesh().mesh()
               .faceZones()[zoneIndex()].name();

            word shadowZoneName =
                patch().boundaryMesh().mesh()
               .faceZones()[shadowZoneIndex()].name();

            if (debug)
            {
                Info<< "Initializing the GGI interpolator between "
                    << "master/shadow zones: "
                    << zoneName << "/" << shadowZoneName
                    << endl;
            }

            calcZoneToZone();
        }

        return *zoneToZonePtr_;
    }
    else
    {
        const volScalarField& field =
            this->db().objectRegistry::lookupObject<volScalarField>
            (
                this->dimensionedInternalField().name()
            );

        const newThermalContactFvPatchScalarField& shadowPatchField =
            refCast<const newThermalContactFvPatchScalarField>
            (
                field.boundaryField()[shadowPatchIndex()]
            );

        return shadowPatchField.zoneToZone();
    }
}


void Foam::newThermalContactFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
}


void Foam::newThermalContactFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const newThermalContactFvPatchScalarField& tiptf =
        refCast<const newThermalContactFvPatchScalarField>(ptf);

    Qc_.rmap(tiptf.Qc_, addr);
}


void Foam::newThermalContactFvPatchScalarField::clearOut()
{
    shadowPatchIndex_ = -1;
    zoneIndex_ = -1;
    shadowZoneIndex_ = -1;
    deleteDemandDrivenData(zoneToZonePtr_);
    deleteDemandDrivenData(HcPtr_);
    deleteDemandDrivenData(contactPtr_);
}


void Foam::newThermalContactFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Delete ggi interpolation at the begining of each time step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        curTimeIndex_ = this->db().time().timeIndex();

        deleteDemandDrivenData(contactPtr_);

        if (master())
        {
            if (debug)
            {
                Info<< "Deleting the GGI interpolation" << endl;
            }
            deleteDemandDrivenData(zoneToZonePtr_);
        }
    }


    // In the contact region:
    // masterFlux + slaveFlux - fricFlux = 0
    // And
    // masterFlux = hBar*(slaveTemp - masterTemp) + masterFricFlux
    // Where
    // masterFricFlux = fricFlux*(masterK/masterKsi)
    // masterKsi = sqrt(slaveRhoC/masterRhoC)*sqrt(masterK*slaveK) + masterK
    // K = conductivity
    // rho = density
    // C = specific heat capacity
    // fricFlux = lossCoeff*shearTraction*relativeVelocity

    // For now, we will calculate the master and slave fluxes using the
    // method
    // We could consider calculating the slaveFlux and then interpolating this
    // to the master to calculate the slave flux.

    // Calculate patch flux
    scalarField curPatchPatchFlux =
        - Hc()()*(shadowPatchTempOnCurPatch()() - *this)
        - Qc()();

    scalarField curPatchPatchSnGrad = -curPatchPatchFlux/curPatchK();

    // Calculate flux for faces not in contact i.e. thermal convection boundary
    // condition
    // k*snGradT = -alpha*(T - Tinf)
    // i.e. heat flux within solid == heat flux due to convection at the surface
    //
    // Therefore:
    // snGradT = -(alpha/k)*(T - Tinf)
    curPatchPatchSnGrad =
        contact()*curPatchPatchSnGrad
      - (1.0 - contact())*(alpha_/curPatchK())*(*this - Tinf_);

    // Set gradient using under-relaxation
    this->gradient() =
        underRelaxation()*curPatchPatchSnGrad
        + (1.0 - underRelaxation())*this->gradient();

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::newThermalContactFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const fvPatchField<vector>& gradT =
        patch().lookupPatchField<volVectorField, vector>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        );

    vectorField n = patch().nf();
    vectorField delta = patch().delta();

    //- non-orthogonal correction vectors
    vectorField k = delta - n*(n&delta);

    Field<scalar>::operator=
    (
        this->patchInternalField()
      + (k & gradT.patchInternalField())
      + gradient()/this->patch().deltaCoeffs()
    );

    fvPatchField<scalar>::evaluate();
}


void Foam::newThermalContactFvPatchScalarField::write(Ostream& os) const
{
    os.writeKeyword("master")
        << master_ << token::END_STATEMENT << nl;
    os.writeKeyword("shadowPatch") << shadowPatchName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("alpha")
        << alpha_ << token::END_STATEMENT << nl;
    os.writeKeyword("Tinf")
        << Tinf_ << token::END_STATEMENT << nl;

    if (master())
    {
        os.writeKeyword("underRelaxation") << underRelaxation_
            << token::END_STATEMENT << nl;
        os.writeKeyword("Rc")
            << Rc_ << token::END_STATEMENT << nl;
        os.writeKeyword("beta")
            << beta_ << token::END_STATEMENT << nl;
        os.writeKeyword("UTS")
            << UTS_ << token::END_STATEMENT << nl;

        if (!useFrictionQc_)
        {
            Qc_.writeEntry("Qc", os);
        }
    }

    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        newThermalContactFvPatchScalarField
    );
}

// ************************************************************************* //
