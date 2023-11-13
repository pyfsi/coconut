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
    thermalGeneralContactFvPatchScalarField

\*---------------------------------------------------------------------------*/

#include "thermalGeneralContactFvPatchScalarField.H"
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


void Foam::thermalGeneralContactFvPatchScalarField::calcGlobalMaster() const
{
    if (globalMasterPtr_)
    {
        FatalErrorIn
        (
            "void Foam::thermalGeneralContactFvPatchScalarField::"
            "calcGlobalMaster() const"
        )   << "globalMasterPtr_ already set" << abort(FatalError);
    }

    // The global master is the first thermalGeneralContact patch i.e. the one
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


void Foam::thermalGeneralContactFvPatchScalarField
::calcGlobalMasterIndex() const
{
    if (globalMasterIndexPtr_)
    {
        FatalErrorIn
        (
            "label Foam::thermalGeneralContactFvPatchScalarField::"
            "calcGlobalMasterIndex() const"
        )   << "globalMasterIndexPtr_ already set" << abort(FatalError);
    }

    // The global master is the first thermalGeneralContact patch i.e. the one
    // with the lowest patch index

    const volScalarField& field =
        this->db().objectRegistry::lookupObject<volScalarField>
        (
            this->dimensionedInternalField().name()
        );

    globalMasterIndexPtr_ = new label(-1);
    label& gMasterID = *globalMasterIndexPtr_;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == thermalGeneralContactFvPatchScalarField::typeName
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
            "thermalGeneralContactFvPatchScalarField::"
            "calcGlobalMasterIndex() const"
        )   << "There are multiple global masters" << abort(FatalError);
    }

    if (debug)
    {
        Pout<< nl << "The global master contact patch is "
            << patch().boundaryMesh()[gMasterID].name() << endl;
    }
}


void Foam::thermalGeneralContactFvPatchScalarField::calcLocalSlave() const
{
    if (localSlavePtr_)
    {
        FatalErrorIn
        (
            "label Foam::thermalGeneralContactFvPatchScalarField::"
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

            Info<< "thermalGeneralContact: "
                << shadowPatchNames()[shadowI] << " (master)" << " to "
                << patch().name() << " (slave)" << endl;
        }
    }
}


const Foam::boolList&
Foam::thermalGeneralContactFvPatchScalarField::localSlave() const
{
    if (!localSlavePtr_)
    {
        calcLocalSlave();
    }

    return *localSlavePtr_;
}


void Foam::thermalGeneralContactFvPatchScalarField::calcShadowPatchNames() const
{
    if (shadowPatchNamesPtr_ || shadowPatchIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::thermalGeneralContactFvPatchScalarField::"
            "calcShadowPatchNames() const"
        )   << "shadowPatchNames_ or shadowPatchIndices_ already set"
            << abort(FatalError);
    }

    // Add each thermalGeneralContact patch in the order of increasing patch
    // index

    const volScalarField& field =
        this->db().objectRegistry::lookupObject<volScalarField>
        (
            this->dimensionedInternalField().name()
        );

    // Count shadow patches

    label nShadPatches = 0;

    forAll(field.boundaryField(), patchI)
    {
        if
        (
            field.boundaryField()[patchI].type()
            == thermalGeneralContactFvPatchScalarField::typeName
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
            == thermalGeneralContactFvPatchScalarField::typeName
            && patchI != patch().index()
        )
        {
            shadowPatchNames[shadowI] = patch().boundaryMesh()[patchI].name();

            shadowPatchIndices[shadowI++] = patchI;
        }
    }
}


void Foam::thermalGeneralContactFvPatchScalarField::calcShadowZoneNames() const
{
    if (shadowZoneNamesPtr_ || shadowZoneIndicesPtr_)
    {
        FatalErrorIn
        (
            "label Foam::thermalGeneralContactFvPatchScalarField::"
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
            FatalErrorIn("thermalGeneralContactFvPatchScalarField")
                << "Face zone name " << zoneName
                << " not found.  Please check your zone definition."
                << abort(FatalError);
        }

        shadowZoneNames[shadowI] = zoneName;

        shadowZoneIndices[shadowI] = zone.index();
    }
}


// const Foam::extendedGgiZoneInterpolation&
// Foam::thermalGeneralContactFvPatchScalarField::zoneToZone
// (
//     const label shadowI
// ) const
// {
//     if (!localSlave()[shadowI])
//     {
//         FatalErrorIn("zoneToZone(const label shadowI)")
//             << "Only the local slave can call the zoneToZone interpolator"
//             << abort(FatalError);
//     }

//     if (zoneToZones_.empty())
//     {
//         if (debug)
//         {
//             word zoneName =
//                patch().boundaryMesh().mesh().faceZones()[zoneIndex()].name();

//             Info<< "Initializing the GGI interpolators for " << zoneName
//                 << endl;
//         }

//         calcZoneToZones();
//     }

//     return zoneToZones_[shadowI];
// }


// Foam::extendedGgiZoneInterpolation&
// Foam::thermalGeneralContactFvPatchScalarField::zoneToZone
//(const label shadowI)
// {
//     if (!localSlave()[shadowI])
//     {
//         FatalErrorIn("zoneToZone(const label shadowI)")
//             << "Only the local slave can call the zoneToZone interpolator"
//             << abort(FatalError);
//     }

//     if (zoneToZones_.empty())
//     {
//         if (debug)
//         {
//             word zoneName =
//               patch().boundaryMesh().mesh().faceZones()[zoneIndex()].name();

//             Info<< "Initializing the GGI interpolators for " << zoneName
//                 << endl;
//         }

//         calcZoneToZones();
//     }

//     return zoneToZones_[shadowI];
// }


void Foam::thermalGeneralContactFvPatchScalarField::calcContact() const
{
    // Create contact indicator field
    if (contactPtr_)
    {
        FatalErrorIn
        (
            "void thermalGeneralContactFvPatchScalarField::calcContact() const"
        )   << "Contact indicator interpolation already calculated"
            << abort(FatalError);
    }

    contactPtr_ = new scalarField(patch().size(), 0);
    scalarField& c = *contactPtr_;


    // Accumulate contact for each pair

    forAll(localSlave(), shadowI)
    {
        c += contact(shadowI);
    }

    // Limit contact to 1.0
    c = min(c, 1.0);

    /*
    word sigmaName("sigma");
    bool foundSigma =
        this->db().objectRegistry::foundObject<volSymmTensorField>(sigmaName);

    if (!foundSigma)
    {
        sigmaName = "sigmaCauchy";
        foundSigma =
            this->db().objectRegistry::foundObject<volSymmTensorField>
            (
                sigmaName
            );
    }

    if (foundSigma)
    {
        const volSymmTensorField& sigma =
            this->db().objectRegistry::lookupObject<volSymmTensorField>
            (
                sigmaName
            );

        // Independently calculate the master and slave contacts

        vectorField n = patch().nf();

        // Faces with a positive contact pressure (i.e. negative normal
        // traction) are deemed to be in contact.
        // Problem: sometimes uncovered faces think they are in contact and
        // this causes numerical cooling.
        scalarField normalTraction =
            n & (n & sigma.boundaryField()[patch().index()]);

        // Use mac pressure as a normalisation factor
        // maxSigma will be used to normalise contact tractions
        scalar normFactor = gMax(-normalTraction) + SMALL;
        normalTraction /= normFactor;

        // Only consider faces with a pressure greater than 10% of the
        // max to be in contact, for thermal contact calculations
        const scalar contactRatioLimit = 0.1;

        forAll(normalTraction, faceI)
        {
            if (normalTraction[faceI] < -contactRatioLimit)
            {
                contact[faceI] = 1;
            }
        }
    }
    else
    {
        FatalErrorIn("calcContact() const")
            << "sigma/sigmaCauchy not field found!" << abort(FatalError);
        }*/
}


const Foam::scalarField& Foam::thermalGeneralContactFvPatchScalarField::
contact() const
{
    if (!contactPtr_)
    {
        calcContact();
    }

    return *contactPtr_;
}


const Foam::scalarField& Foam::thermalGeneralContactFvPatchScalarField::contact
(
    const label shadowI
) const
{
    if (!contactsPtr_)
    {
        calcContacts();
    }

    return (*contactsPtr_)[shadowI];
}


void Foam::thermalGeneralContactFvPatchScalarField::calcContacts() const
{
    const boolList& locSlave = localSlave();

    contactsPtr_ = new List<scalarField>(locSlave.size());

    const volVectorField& DU =
        this->db().objectRegistry::lookupObject<volVectorField>
        (
            DUName_
        );

    const solidGeneralContactFvPatchVectorField& DUPatch =
        refCast<const solidGeneralContactFvPatchVectorField>
        (
            DU.boundaryField()[patch().index()]
        );

    // We will designate faces with a traction greater than a certain tolerance
    // to be in contact, for the thermal contact boundaries.
    // But what is the best choice for this tolerance...
    // We will normalise it relative to the average stiffness of the
    // materials: we will use the shear modulus mu.

    const volScalarField& muField =
        this->db().objectRegistry::lookupObject<volScalarField>("mu");

    forAll(locSlave, shadowI)
    {
        const scalarField masterMagTrac =
            mag(DUPatch.curPatchTractions(shadowI));

        scalarField& contact = (*contactsPtr_)[shadowI];
        contact.setSize(masterMagTrac.size(), 0.0);

        const label shadowID = shadowPatchIndices()[shadowI];

        const scalar avMu =
            gAverage
            (
                muField.boundaryField()[patch().index()]
            )
            + gAverage
            (
                muField.boundaryField()[shadowID]
            );

        const scalar tol = 1e-5*avMu;

        forAll(masterMagTrac, faceI)
        {
            if (masterMagTrac[faceI] > tol)
            {
                contact[faceI] = 1.0;
            }
        }
    }
}


void Foam::thermalGeneralContactFvPatchScalarField::calcHc() const
{
    HcPtr_ = new scalarField(patch().size(), 0.0);

    scalarField& h = *HcPtr_;

    // Contact resistance dependent on contact pressure

    // Pressure dependent contact conductance eg:
    // h = hRef*((p/H)**beta)
    // where
    // p = contact pressure
    // H = Vicker's hardness of softer material
    // hRef, beta = experimentally fit coefficients
    // beta determines pressure sensitivity
    // beta > 1 very sensitive
    // beta < 0.01 very insensitive

    // Contact resistance is the reciprocal of contact conductance

    word sigmaName("sigma");
    bool foundSigma =
        this->db().objectRegistry::foundObject<volSymmTensorField>
        (
            sigmaName
        );

    if (!foundSigma)
    {
        sigmaName = "sigmaCauchy";
        foundSigma =
            this->db().objectRegistry::foundObject<volSymmTensorField>
            (
                sigmaName
            );
    }

    if (!foundSigma)
    {
        FatalErrorIn("thermalGeneralContactFvPatchScalarField::Hc()")
            << "sigma or sigmaCauchy field not found"
            << abort(FatalError);
    }

    const volSymmTensorField& sigma =
        this->db().objectRegistry::lookupObject<volSymmTensorField>
        (
            sigmaName
        );

    vectorField n = patch().nf();

    scalarField contactPressure =
        -n & (n & sigma.boundaryField()[patch().index()]);

    // Limit minimum contact pressure to avoid FPE in pow
    contactPressure = max(contactPressure, SMALL);

    // Hmnn formula says use Vicker's hardness, but surely it should be
    // Vicker's hardness by 1e6
    // as Vicker's hardness = 0.3*UTS in MPa

    h = Foam::pow(contactPressure/(0.3*UTS()), beta())/Rc();

    if (debug)
    {
        Info<< "contact conductance:" << nl
            << "    max: " << max(h) << nl
            << "    min: " << min(h) << nl
            << "contact resistance: " << nl
            << "    max: " << 1.0/(max(h) + SMALL) << nl
            << "    min: " << 1.0/(min(h) + SMALL) << endl;
    }
}


Foam::label Foam::thermalGeneralContactFvPatchScalarField::findShadowID
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::thermalGeneralContactFvPatchScalarField::
thermalGeneralContactFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    fieldName_("undefined"),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(-1),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    underRelaxation_(1.0),
    thermalConductivityName_("undefined"),
    alpha_(p.size(), 0),
    Tinf_(0),
    Rc_(0),
    beta_(0),
    UTS_(0),
    Qc_(p.size(), 0),
    useFrictionQc_(false),
    DUName_("undefined"),
    curTimeIndex_(-1),
    //zoneToZones_(0),
    HcPtr_(NULL),
    contactPtr_(NULL),
    contactsPtr_(NULL)
{
    if (debug)
    {
        InfoIn
        (
            "Foam::thermalGeneralContactFvPatchScalarField::"
            "thermalGeneralContactFvPatchScalarField"
            "("
            "    const fvPatch& p,"
            "    const DimensionedField<scalar, volMesh>& iF"
            ")"
        ) << endl;
    }
}


Foam::thermalGeneralContactFvPatchScalarField::
thermalGeneralContactFvPatchScalarField
(
    const thermalGeneralContactFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(ptf.zoneIndex_),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    underRelaxation_(ptf.underRelaxation_),
    thermalConductivityName_(ptf.thermalConductivityName_),
    alpha_(ptf.alpha_),
    Tinf_(ptf.Tinf_),
    Rc_(ptf.Rc_),
    beta_(ptf.beta_),
    UTS_(ptf.UTS_),
    Qc_(ptf.Qc_),
    useFrictionQc_(ptf.useFrictionQc_),
    DUName_(ptf.DUName_),
    curTimeIndex_(ptf.curTimeIndex_),
    //zoneToZones_(0),
    HcPtr_(NULL),
    contactPtr_(NULL),
    contactsPtr_(NULL)
{
    if (debug)
    {
        InfoIn
        (
            "Foam::thermalGeneralContactFvPatchScalarField::"
            "thermalGeneralContactFvPatchScalarField"
            "("
            "    const thermalGeneralContactFvPatchScalarField& ptf,"
            "    const fvPatch& p,"
            "    const DimensionedField<scalar, volMesh>& iF,"
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

    // if (!ptf.zoneToZones_.empty())
    // {
    //     // I will not copy the GGI interpolators
    //     // They can be re-created when required
    // }

    if (ptf.HcPtr_)
    {
        HcPtr_ = new scalarField(*ptf.HcPtr_);
    }

    if (ptf.contactPtr_)
    {
        contactPtr_ = new scalarField(*ptf.contactPtr_);
    }

    if (ptf.contactsPtr_)
    {
        contactsPtr_ = new List<scalarField>(*ptf.contactsPtr_);
    }
}


Foam::thermalGeneralContactFvPatchScalarField::
thermalGeneralContactFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:   fixedGradientFvPatchScalarField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(-1),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    underRelaxation_(1),
    thermalConductivityName_("undefined"),
    alpha_("alpha", dict, p.size()),
    Tinf_(readScalar(dict.lookup("Tinf"))),
    Rc_(0.0),
    beta_(0.0),
    UTS_(0.0),
    Qc_(),
    useFrictionQc_(true),
    DUName_(dict.lookupOrDefault<word>("DUName", "DU")),
    curTimeIndex_(-1),
    //zoneToZones_(0),
    HcPtr_(NULL),
    contactPtr_(NULL),
    contactsPtr_(NULL)
{
    Info<< "Creating " << thermalGeneralContactFvPatchScalarField::typeName
        << " patch" << endl;

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    // Zone index
    word zoneName = patch().name() + "FaceZone";

    faceZoneID zone(zoneName, mesh.faceZones());

    if (!zone.active())
    {
        FatalErrorIn("thermalGeneralContactFvPatchScalarField")
            << "Face zone name " << zoneName
            << " not found.  Please check your zone definition."
            << abort(FatalError);
    }

    zoneIndex_ = zone.index();

    underRelaxation_ = readScalar(dict.lookup("underRelaxation"));

    thermalConductivityName_ =
        dict.lookupOrDefault<word>("thermalConductivityName", "k");

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
            WarningIn
            (
                "Foam::thermalGeneralContactFvPatchScalarField::"
                "thermalGeneralContactFvPatchScalarField"
            )   << "As Qc is found in the thermal boundary and is greater "
                << "then zero, it will be used instead of frictional energy"
                << endl;

            useFrictionQc_ = false;
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


Foam::thermalGeneralContactFvPatchScalarField::
thermalGeneralContactFvPatchScalarField
(
    const thermalGeneralContactFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    fieldName_(ptf.fieldName_),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(ptf.zoneIndex_),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    underRelaxation_(ptf.underRelaxation_),
    thermalConductivityName_(ptf.thermalConductivityName_),
    alpha_(ptf.alpha_),
    Tinf_(ptf.Tinf_),
    Rc_(ptf.Rc_),
    beta_(ptf.beta_),
    UTS_(ptf.UTS_),
    Qc_(ptf.Qc_),
    useFrictionQc_(ptf.useFrictionQc_),
    DUName_(ptf.DUName_),
    curTimeIndex_(ptf.curTimeIndex_),
    //zoneToZones_(0),
    HcPtr_(NULL),
    contactPtr_(NULL),
    contactsPtr_(NULL)
{
    if (debug)
    {
        InfoIn
        (
            "Foam::thermalGeneralContactFvPatchScalarField::"
            "thermalGeneralContactFvPatchScalarField"
            "("
            "    const thermalGeneralContactFvPatchScalarField& ptf"
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

    if (ptf.HcPtr_)
    {
        HcPtr_ = new scalarField(*ptf.HcPtr_);
    }

    if (ptf.contactPtr_)
    {
        contactPtr_ = new scalarField(*ptf.contactPtr_);
    }

    if (ptf.contactsPtr_)
    {
        contactsPtr_ = new List<scalarField>(*ptf.contactsPtr_);
    }
}


Foam::thermalGeneralContactFvPatchScalarField::
thermalGeneralContactFvPatchScalarField
(
    const thermalGeneralContactFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    fieldName_(ptf.fieldName_),
    globalMasterPtr_(NULL),
    globalMasterIndexPtr_(NULL),
    localSlavePtr_(NULL),
    shadowPatchNamesPtr_(NULL),
    shadowPatchIndicesPtr_(NULL),
    zoneIndex_(ptf.zoneIndex_),
    shadowZoneNamesPtr_(NULL),
    shadowZoneIndicesPtr_(NULL),
    underRelaxation_(ptf.underRelaxation_),
    thermalConductivityName_(ptf.thermalConductivityName_),
    alpha_(ptf.alpha_),
    Tinf_(ptf.Tinf_),
    Rc_(ptf.Rc_),
    beta_(ptf.beta_),
    UTS_(ptf.UTS_),
    Qc_(ptf.Qc_),
    useFrictionQc_(ptf.useFrictionQc_),
    DUName_(ptf.DUName_),
    curTimeIndex_(ptf.curTimeIndex_),
    //zoneToZones_(0),
    HcPtr_(NULL),
    contactPtr_(NULL),
    contactsPtr_(NULL)
{
    if (debug)
    {
        InfoIn
        (
            "Foam::thermalGeneralContactFvPatchScalarField::"
            "thermalGeneralContactFvPatchScalarField"
            "("
            "    const thermalGeneralContactFvPatchScalarField& ptf,"
            "    const DimensionedField<scalar, volMesh>& iF"
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

    if (ptf.HcPtr_)
    {
        HcPtr_ = new scalarField(*ptf.HcPtr_);
    }

    if (ptf.contactPtr_)
    {
        contactPtr_ = new scalarField(*ptf.contactPtr_);
    }

    if (ptf.contactsPtr_)
    {
        contactsPtr_ = new List<scalarField>(*ptf.contactsPtr_);
    }
}


// * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * * //


Foam::thermalGeneralContactFvPatchScalarField::
~thermalGeneralContactFvPatchScalarField()
{
    if (debug)
    {
        InfoIn
        (
            "Foam::thermalGeneralContactFvPatchScalarField::"
            "~thermalGeneralContactFvPatchScalarField()"
        ) << endl;
    }

    deleteDemandDrivenData(globalMasterPtr_);
    deleteDemandDrivenData(globalMasterIndexPtr_);
    deleteDemandDrivenData(localSlavePtr_);
    deleteDemandDrivenData(shadowPatchNamesPtr_);
    deleteDemandDrivenData(shadowPatchIndicesPtr_);
    deleteDemandDrivenData(shadowZoneNamesPtr_);
    deleteDemandDrivenData(shadowZoneIndicesPtr_);
    deleteDemandDrivenData(HcPtr_);
    deleteDemandDrivenData(contactPtr_);
    deleteDemandDrivenData(contactsPtr_);
}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //


bool Foam::thermalGeneralContactFvPatchScalarField::globalMaster() const
{
    if (!globalMasterPtr_)
    {
        calcGlobalMaster();
    }

    return *globalMasterPtr_;
}


Foam::label Foam::thermalGeneralContactFvPatchScalarField::globalMasterIndex()
const
{
    if (!globalMasterIndexPtr_)
    {
        calcGlobalMasterIndex();
    }

    return *globalMasterIndexPtr_;
}


const Foam::List<Foam::word>&
Foam::thermalGeneralContactFvPatchScalarField::shadowPatchNames() const
{
    if (!shadowPatchNamesPtr_)
    {
        calcShadowPatchNames();
    }

    return *shadowPatchNamesPtr_;
}


const Foam::List<Foam::label>&
Foam::thermalGeneralContactFvPatchScalarField::shadowPatchIndices() const
{
    if (!shadowPatchIndicesPtr_)
    {
        calcShadowPatchNames();
    }

    return *shadowPatchIndicesPtr_;
}


const Foam::List<Foam::word>&
Foam::thermalGeneralContactFvPatchScalarField::shadowZoneNames() const
{
    if (!shadowZoneNamesPtr_)
    {
        calcShadowZoneNames();
    }

    return *shadowZoneNamesPtr_;
}


const Foam::List<Foam::label>&
Foam::thermalGeneralContactFvPatchScalarField::shadowZoneIndices() const
{
    if (!shadowZoneIndicesPtr_)
    {
        calcShadowZoneNames();
    }

    return *shadowZoneIndicesPtr_;
}


const Foam::scalarField&
Foam::thermalGeneralContactFvPatchScalarField::Hc() const
{
    if (!HcPtr_)
    {
        calcHc();
    }

    return *HcPtr_;
}


Foam::tmp<Foam::scalarField>
Foam::thermalGeneralContactFvPatchScalarField::Qc() const
{
    tmp<scalarField> tQc(new scalarField(patch().size(), 0.0));

    scalarField& Qc = tQc();

    if
    (
        this->db().objectRegistry::foundObject<volVectorField>(DUName_)
        && useFrictionQc_
    )
    {
        const volVectorField& DU =
            this->db().objectRegistry::lookupObject<volVectorField>
            (
                DUName_
            );

        if
        (
            DU.boundaryField()[patch().index()].type()
            != solidGeneralContactFvPatchVectorField::typeName
        )
        {
            FatalErrorIn("thermalGeneralContactFvPatchScalarField::Qc()")
                << "DU patch " << patch().boundaryMesh()[patch().index()].name()
                << " should be of type solidContact " << abort(FatalError);
        }
        else
        {
            const solidGeneralContactFvPatchVectorField& contactDU =
                refCast<const solidGeneralContactFvPatchVectorField>
                (
                    DU.boundaryField()[patch().index()]
                );

            // In the contact region:
            // masterFlux + slaveFlux - fricFlux = 0
            // And
            // masterFlux = hBar*(slaveTemp - masterTemp) + masterFricFlux
            // Where
            // masterFricFlux = fricFlux*(masterK/masterKsi)
            // masterKsi =
            //     sqrt(slaveRhoC/masterRhoC)*sqrt(masterK*slaveK) + masterK
            // K = conductivity
            // rho = density
            // C = specific heat capacity
            // fricFlux = lossCoeff*shearTraction*relativeVelocity


            const volScalarField& k =
                this->db().objectRegistry::lookupObject<volScalarField>
                (
                    thermalConductivityName()
                );

            const scalarField& curPatchK = k.boundaryField()[patch().index()];

            const volScalarField& rhoC =
                this->db().objectRegistry::lookupObject<volScalarField>
                (
                    "(rho*C)"
                );

            const scalarField& curPatchRhoC =
                rhoC.boundaryField()[patch().index()];

            const boolList& locSlave = localSlave();

            forAll(locSlave, shadowI)
            {
                const scalarField fricFlux = contactDU.Qc(shadowI);

                const scalarField& shadowPatchK =
                    k.boundaryField()[shadowPatchIndices()[shadowI]];

                const scalarField& shadowPatchRhoC =
                    rhoC.boundaryField()[shadowPatchIndices()[shadowI]];


                // Interpolate shadow fields to the curPatch

                scalarField shadowZoneK =
                    zoneField
                    (
                        shadowZoneIndices()[shadowI],
                        shadowPatchIndices()[shadowI],
                        shadowPatchK
                    );

                scalarField shadowZoneRhoC =
                    zoneField
                    (
                        shadowZoneIndices()[shadowI],
                        shadowPatchIndices()[shadowI],
                        shadowPatchRhoC
                    );

                scalarField shadowZoneKOnCurPatch
                (
                    patch().boundaryMesh().mesh().faceZones()
                    [
                        shadowZoneIndices()[shadowI]
                    ].size(),
                    0.0
                );

                scalarField shadowZoneRhoCOnCurPatch
                (
                    patch().boundaryMesh().mesh().faceZones()
                    [
                        shadowZoneIndices()[shadowI]
                    ].size(),
                    0.0
                );

                if (locSlave[shadowI])
                {
                    shadowZoneKOnCurPatch =
                        contactDU.zoneToZone
                        (
                            shadowI
                        ).masterToSlave(shadowZoneK);

                    shadowZoneRhoCOnCurPatch =
                        contactDU.zoneToZone
                        (
                            shadowI
                        ).masterToSlave(shadowZoneRhoC);
                }
                else
                {
                    const solidGeneralContactFvPatchVectorField&
                        shadowContactDU =
                        refCast<const solidGeneralContactFvPatchVectorField>
                        (
                            DU.boundaryField()[shadowPatchIndices()[shadowI]]
                        );

                    const label locShadowID =
                        shadowContactDU.findShadowID(patch().index());

                    shadowZoneKOnCurPatch =
                        shadowContactDU.zoneToZone
                        (
                            locShadowID
                        ).slaveToMaster(shadowZoneK);

                    shadowZoneRhoCOnCurPatch =
                        shadowContactDU.zoneToZone
                        (
                            locShadowID
                        ).slaveToMaster(shadowZoneRhoC);
                }

                scalarField shadowPatchKOnCurPatch =
                    patchField
                    (
                        patch().index(),
                        zoneIndex(),
                        shadowZoneKOnCurPatch
                    );

                scalarField shadowPatchRhoCOnCurPatch =
                    patchField
                    (
                        patch().index(),
                        zoneIndex(),
                        shadowZoneRhoCOnCurPatch
                    );
                scalarField curPatchKsi =
                    Foam::sqrt((shadowPatchRhoCOnCurPatch/curPatchRhoC)
                    *curPatchK*shadowPatchKOnCurPatch)
                    + curPatchK;

                Qc += fricFlux*(curPatchK/curPatchKsi);
            }
        }
    }

    return tQc;
}


// Map from self
void Foam::thermalGeneralContactFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    FatalErrorIn
    (
        "void Foam::thermalGeneralContactFvPatchScalarField::autoMap"
        "("
        "    const fvPatchFieldMapper& m"
        ")"
    )   << "member mapping not implemented" << endl;

    fixedGradientFvPatchScalarField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void Foam::thermalGeneralContactFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    FatalErrorIn
    (
        "void Foam::thermalGeneralContactFvPatchScalarField::rmap"
        "("
        "    const fvPatchField<scalar>& ptf,"
        "    const labelList& addr"
        ")"
    )   << "member mapping not implemented" << endl;

    fixedGradientFvPatchScalarField::rmap(ptf, addr);
}


void Foam::thermalGeneralContactFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // if it is a new time step then reset iCorr
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        curTimeIndex_ = this->db().time().timeIndex();

        deleteDemandDrivenData(HcPtr_);
        deleteDemandDrivenData(contactPtr_);
        deleteDemandDrivenData(contactsPtr_);
    }


    // Loop through all shadow patches

    const boolList& locSlave = localSlave();

    const volScalarField& field =
        this->db().objectRegistry::lookupObject<volScalarField>
        (
            this->dimensionedInternalField().name()
        );

    const volScalarField& k =
        this->db().objectRegistry::lookupObject<volScalarField>
        (
            thermalConductivityName()
        );

    const scalarField& curPatchK = k.boundaryField()[patch().index()];

    // Accumulated snGrad (i.e. -flux/k) for the current patch
    scalarField curPatchPatchSnGrad(patch().size(), 0.0);

    boolList activeContactPairs(shadowPatchNames().size(), true);

    forAll(activeContactPairs, shadowI)
    {
        const thermalGeneralContactFvPatchScalarField& shadowPatchField =
            refCast<const thermalGeneralContactFvPatchScalarField>
            (
                field.boundaryField()[shadowPatchIndices()[shadowI]]
            );

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

        // For now, we will calculate the master and slave fluxes using this
        // method
        // We could consider calculating the slaveFlux and then interpolating
        // this to the master to calculate the slave flux.


        // Interpolate shadow temperature field to the curPatch

        scalarField shadowZoneTemp =
            zoneField
            (
                shadowZoneIndices()[shadowI],
                shadowPatchIndices()[shadowI],
                shadowPatchField
            );

        scalarField shadowZoneTempOnCurPatch(patch().size(), 0.0);

        const volVectorField& DU =
            this->db().objectRegistry::lookupObject<volVectorField>
            (
                DUName_
            );

        if (locSlave[shadowI])
        {
            const solidGeneralContactFvPatchVectorField& slaveContactDU =
                refCast<const solidGeneralContactFvPatchVectorField>
                (
                    DU.boundaryField()[patch().index()]
                );

            shadowZoneTempOnCurPatch =
                slaveContactDU.zoneToZone
                //zoneToZone
                (
                    shadowI
                ).masterToSlave(shadowZoneTemp);
        }
        else
        {
            const solidGeneralContactFvPatchVectorField& slaveContactDU =
                refCast<const solidGeneralContactFvPatchVectorField>
                (
                    DU.boundaryField()[shadowPatchIndices()[shadowI]]
                );

            const label locShadowID =
                //shadowPatchField.findShadowID(patch().index());
                slaveContactDU.findShadowID(patch().index());

            shadowZoneTempOnCurPatch =
                slaveContactDU.zoneToZone
                //shadowPatchField.zoneToZone
                (
                    locShadowID
                ).slaveToMaster(shadowZoneTemp);
        }

        scalarField shadowPatchTempOnCurPatch =
            patchField
            (
                patch().index(),
                zoneIndex(),
                shadowZoneTempOnCurPatch
            );

        // Calculate curPatchH
        // We will include pressure dependent contact conductance:
        // hBar = hRef*((p/H)**beta)
        // where
        // p = contact pressure
        // H = Vicker's hardness of softer material
        // hRef, beta = experimentally fit coefficients
        scalarField curPatchH = Hc();

        scalarField curPatchPatchFlux =
            -curPatchH*(shadowPatchTempOnCurPatch - *this) - Qc();

        curPatchPatchSnGrad -= contact(shadowI)*curPatchPatchFlux/curPatchK;
    }

    // Calculate flux for faces not in contact i.e. thermal convection
    // boundary condition
    // k*snGradT = alpha*(T - Tinf)
    // i.e. heat flux within solid == heat flux due to convection at the
    // surface
    // Therefore:
    // snGradT = (alpha/k)*(T - Tinf)

    curPatchPatchSnGrad =
        contact()*curPatchPatchSnGrad
        + (1.0 - contact())*(alpha_/curPatchK)*(*this - Tinf_);

    // Set gradient using under-relaxation
    this->gradient() =
        underRelaxation()*curPatchPatchSnGrad
        + (1.0 - underRelaxation())*this->gradient();

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::thermalGeneralContactFvPatchScalarField::evaluate
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

// Foam::tmp<Foam::Field<Foam::vector> >
// Foam::thermalGeneralContactFvPatchScalarField::snGrad() const
// {
//     vectorField pif = this->patchInternalField();

//     vectorField normalValue = transform(valueFraction(), refValue());

//     const fvPatchField<tensor>& gradField =
//         patch().lookupPatchField<volTensorField, tensor>
//         (
//             "grad(" + fieldName_ + ")"
//         );
//     vectorField n = patch().nf();
//     vectorField delta = patch().delta();
//     //- correction vector
//     vectorField k = delta - n*(n&delta);

//     vectorField gradValue =
//       this->patchInternalField()
//       + (k&gradField.patchInternalField())
//       + refGrad()/this->patch().deltaCoeffs();

//     vectorField transformGradValue =
//       transform(I - valueFraction(), gradValue);

//     return
//       (
//        (normalValue + transformGradValue)
//        - (pif + (k&gradField.patchInternalField()))
//        )*this->patch().deltaCoeffs();
// }


// Write
void Foam::thermalGeneralContactFvPatchScalarField::write(Ostream& os) const
{
    os.writeKeyword("underRelaxation") << underRelaxation_
        << token::END_STATEMENT << nl;

    os.writeKeyword("thermalConductivityName")
        << thermalConductivityName_ << token::END_STATEMENT << nl;

    os.writeKeyword("Rc") << Rc_ << token::END_STATEMENT << nl;

    os.writeKeyword("beta") << beta_ << token::END_STATEMENT << nl;

    os.writeKeyword("UTS") << UTS_ << token::END_STATEMENT << nl;

    if (!useFrictionQc_)
    {
        Qc_.writeEntry("Qc", os);
    }

    alpha_.writeEntry("alpha", os);

    os.writeKeyword("Tinf") << Tinf_ << token::END_STATEMENT << nl;

    fixedGradientFvPatchScalarField::write(os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField(fvPatchScalarField, thermalGeneralContactFvPatchScalarField)
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//} // End namespace Foam

// ************************************************************************* //
