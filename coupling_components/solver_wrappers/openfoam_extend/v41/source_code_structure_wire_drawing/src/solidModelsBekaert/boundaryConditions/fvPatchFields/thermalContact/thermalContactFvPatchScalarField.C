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

#include "thermalContactFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //


void Foam::thermalContactFvPatchScalarField::checkConsistentMaster() const
{
    if (master_ != solidContactPatch().master())
    {
        FatalErrorIn("void checkConsistentMaster() const")
            << "The solidContact master patch should be the same as the "
            << "thermalContact master patch!" << nl
            << abort(FatalError);
    }
}


const Foam::thermalContactFvPatchScalarField&
Foam::thermalContactFvPatchScalarField::shadowPatchField() const
{
    const labelList& shadowPatchIndices =
        solidContactPatch().shadowPatchIndices();

    if (shadowPatchIndices.size() != 1)
    {
        FatalErrorIn
        (
            "const Foam::thermalContactFvPatchScalarField&\n"
            "Foam::thermalContactFvPatchScalarField::shadowPatchField() const"
        )   << "This function can only be called for a patch with 1 shadow "
            << "patch; this patch has " << shadowPatchIndices.size()
            << " shadow patches!" << abort(FatalError);
    }

    return shadowPatchField(0);
}


const Foam::thermalContactFvPatchScalarField&
Foam::thermalContactFvPatchScalarField::shadowPatchField
(
    const label shadI
) const
{
    if (shadI < 0)
    {
        FatalErrorIn
        (
            "const Foam::thermalContactFvPatchScalarField&\n"
            "Foam::thermalContactFvPatchScalarField::"
            "shadowPatchField(const label shadI) const"
        )   << "shadI must be a non-negative number!"
            << abort(FatalError);
    }

    const labelList& shadowPatchIndices =
        solidContactPatch().shadowPatchIndices();

    if (shadI >= shadowPatchIndices.size())
    {
        FatalErrorIn
        (
            "const Foam::thermalContactFvPatchScalarField&\n"
            "Foam::thermalContactFvPatchScalarField::"
            "shadowPatchField(const label shadI) const"
        )   << "shadI is " << shadI << " but this patch has only "
            << shadowPatchIndices.size() << " shadow patches!"
            << abort(FatalError);
    }

    const volScalarField& field =
        db().lookupObject<volScalarField>(dimensionedInternalField().name());

    return
        refCast<const thermalContactFvPatchScalarField>
        (
            field.boundaryField()[shadowPatchIndices[shadI]]
        );
}


const Foam::solidContactFvPatchVectorField&
Foam::thermalContactFvPatchScalarField::solidContactPatch() const
{
    return
        refCast<const solidContactFvPatchVectorField>
        (
            patch().lookupPatchField<volVectorField, vector>(DUName_)
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalContactFvPatchScalarField::thermalContactFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    master_(false),
    dict_(),
    underRelaxation_(1),
    alpha_(p.size(), 0),
    Tinf_(0),
    Rc_(0),
    beta_(0),
    UTS_(0),
    DUName_("undefined"),
    curTimeIndex_(-1)
{
    fvPatchScalarField::operator=(patchInternalField());
    gradient() = 0.0;
}


Foam::thermalContactFvPatchScalarField::thermalContactFvPatchScalarField
(
    const thermalContactFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    master_(ptf.master_),
    dict_(ptf.dict_),
    underRelaxation_(ptf.underRelaxation_),
    alpha_(ptf.alpha_, mapper),
    Tinf_(ptf.Tinf_),
    Rc_(ptf.Rc_),
    beta_(ptf.beta_),
    UTS_(ptf.UTS_),
    DUName_(ptf.DUName_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


Foam::thermalContactFvPatchScalarField::thermalContactFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    master_(dict.lookupOrDefault<Switch>("master", false)),
    dict_(dict),
    underRelaxation_(1),
    alpha_("alpha", dict, p.size()),
    Tinf_(readScalar(dict.lookup("Tinf"))),
    Rc_(0.0),
    beta_(0.0),
    UTS_(0.0),
    DUName_(dict.lookupOrDefault<word>("DUName", "DU")),
    curTimeIndex_(-1)
{
    if (debug)
    {
        Info<< patch().name() << ": " << type() << endl;
    }

    // Read only on master
    if (master())
    {
        underRelaxation_ = readScalar(dict.lookup("underRelaxation"));

        Rc_ = readScalar(dict.lookup("Rc"));

        beta_ = dict.lookupOrDefault<scalar>("beta", 0.1);

        UTS_ = dict.lookupOrDefault<scalar>("UTS", 2e9);

        if (Rc_ < SMALL)
        {
            WarningIn(type() + "::" +  type() + "(...)")
                << "Contact conductivity resistance cannot be exactly zero"
                << nl << "Rc = " << SMALL << " is assumed" << endl;

            Rc_ = SMALL;
        }
    }

    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        gradient() = 0.0;
    }

    if (dict.found("value"))
    {
        Field<scalar>::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        Field<scalar>::operator=
        (
            patchInternalField() + gradient()/patch().deltaCoeffs()
        );
    }
}


Foam::thermalContactFvPatchScalarField::thermalContactFvPatchScalarField
(
    const thermalContactFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    master_(ptf.master_),
    dict_(ptf.dict_),
    underRelaxation_(ptf.underRelaxation_),
    alpha_(ptf.alpha_),
    Tinf_(ptf.Tinf_),
    Rc_(ptf.Rc_),
    beta_(ptf.beta_),
    UTS_(ptf.UTS_),
    DUName_(ptf.DUName_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


Foam::thermalContactFvPatchScalarField::thermalContactFvPatchScalarField
(
    const thermalContactFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    master_(ptf.master_),
    dict_(ptf.dict_),
    underRelaxation_(ptf.underRelaxation_),
    alpha_(ptf.alpha_),
    Tinf_(ptf.Tinf_),
    Rc_(ptf.Rc_),
    beta_(ptf.beta_),
    UTS_(ptf.UTS_),
    DUName_(ptf.DUName_),
    curTimeIndex_(ptf.curTimeIndex_)
{}



// * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * * //

Foam::thermalContactFvPatchScalarField::~thermalContactFvPatchScalarField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::wordList&
Foam::thermalContactFvPatchScalarField::shadowPatchNames() const
{
    return solidContactPatch().shadowPatchNames();
}


const Foam::labelList&
Foam::thermalContactFvPatchScalarField::shadowPatchIndices() const
{
    return solidContactPatch().shadowPatchIndices();
}


const Foam::scalarField& Foam::thermalContactFvPatchScalarField::contact() const
{
    return solidContactPatch().contact();
}


Foam::scalar Foam::thermalContactFvPatchScalarField::underRelaxation() const
{
    if (master())
    {
        return underRelaxation_;
    }
    else
    {
        return shadowPatchField().underRelaxation();
    }
}


Foam::scalar Foam::thermalContactFvPatchScalarField::Rc() const
{
    if (master())
    {
        return Rc_;
    }
    else
    {
        return shadowPatchField().Rc();
    }
}


Foam::scalar Foam::thermalContactFvPatchScalarField::beta() const
{
    if (master())
    {
        return beta_;
    }
    else
    {
        return shadowPatchField().beta();
    }
}


Foam::scalar Foam::thermalContactFvPatchScalarField::UTS() const
{
    if (master())
    {
        return UTS_;
    }
    else
    {
        return shadowPatchField().UTS();
    }
}


Foam::tmp<Foam::scalarField> Foam::thermalContactFvPatchScalarField::Hc() const
{
    checkConsistentMaster();

    tmp<scalarField> tHc
    (
        new scalarField(patch().size(), 0.0)
    );

    scalarField& Hc = tHc();

    // Contact resistance dependent on contact pressure

    // Pressure dependent contact conductance, where
    // h = hRef*((p/H)**beta)
    // where
    // p = contact pressure
    // H = a measure of the hardness of softer material in Pascal
    // hRef, beta = experimentally fit coefficients
    // beta determines pressure sensitivity
    // beta > 1 very sensitive
    // beta < 0.01 very insensitive
    // See reference: P. Wriggers, Computational Contact Mechanics, Wiley, 2002.

    // Contact resistance is the reciprocal of contact conductance

    // Lookup contact field from DU solidContact patch

    // Unit normals
    // Note: the mesh should be in the deformed position so these should be the
    // deformed configuration normals
    const vectorField n = patch().nf();

    // Calculate contact pressure and limit to avoid division by zero in pow
    const scalarField contactPressure =
        max(-n & solidContactPatch().traction(), SMALL);

    // Hmnn formula says use Vicker's hardness, but surely it should be
    // Vicker's hardness by 1e6
    // as Vicker's hardness = 0.3*UTS in MPa

    Hc = Foam::pow(contactPressure/(0.3*UTS()), beta())/Rc();

    return tHc;
}


Foam::tmp<Foam::scalarField>
Foam::thermalContactFvPatchScalarField::frictionFluxRateForThisPatch() const
{
    checkConsistentMaster();

    tmp<scalarField> tfrictionFluxRateForThisPatch
    (
        new scalarField(patch().size(), 0.0)
    );

    // Check the DU patch is of type solidContact
    if
    (
        solidContactPatch().type()
     != solidContactFvPatchVectorField::typeName
    )
    {
        FatalErrorIn
        (
            "thermalContactFvPatchScalarField::frictionFluxRateForThisPatch()"
        )   << "DU patch " << patch().boundaryMesh()[patch().index()].name()
            << " should be of type solidContact "
            << abort(FatalError);
    }
    else
    {
        // Heat flux generated due to friction, stored on the master surface
        scalarField fricFlux(patch().size(), 0.0);

        if (master_)
        {
            fricFlux = solidContactPatch().frictionHeatRate();
        }
        else
        {
            // Heat flux on the master global patch
            const scalarField masterZoneFricFlux =
                solidContactPatch().shadowPatchField().zone().patchFaceToGlobal
                (
                    solidContactPatch().shadowPatchField().frictionHeatRate()()
                );

            // Interpolate the heat flux from the master global patch to the
            // shadow global patch
            const scalarField shadowZoneFricFlux =
                solidContactPatch().zoneToZoneForThisSlave().masterToSlave
                (
                    masterZoneFricFlux
                );

            // Convert the global shadow patch to the shadow patch
            // Note: solidContactPatch().zone() always returns the master zone
            // and solidContactPatch().shadowZones() always returns the slave
            // zones
            fricFlux =
                solidContactPatch().zoneForThisSlave().globalFaceToPatch
                (
                    shadowZoneFricFlux
                );
        }

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
        // Here we assume the lossCoeff is zero

        // Lookup the thermal conductivity field
        const volScalarField& k =
            db().lookupObject<volScalarField>("k");

        // Lookup density times specific heat field
        const volScalarField& rhoC =
            db().lookupObject<volScalarField>("(rho*C)");

        // Current patch index
        const label currentPatchID = patch().index();

        // Get k on the current patch
        const scalarField& curPatchK = k.boundaryField()[currentPatchID];

        // Get rhoC on the current patch
        const scalarField& curPatchRhoC = rhoC.boundaryField()[currentPatchID];

        // For the slave, shadowPatchNames will be length 1, whereas it could be
        // longer for the master
        const wordList shadowPatchNames = this->shadowPatchNames();

        // Reset to zero as we will accumulate frictionFluxRateForThisPatch for
        // all shadow patch
        tfrictionFluxRateForThisPatch() = 0.0;

        forAll(shadowPatchNames, shadI)
        {
            // Get k on the shadow patch
            const scalarField& shadowPatchK =
                k.boundaryField()[shadowPatchIndices()[shadI]];

            // Get rhoC on the shadow patch
            const scalarField& shadowPatchRhoC =
                rhoC.boundaryField()[shadowPatchIndices()[shadI]];

            // Interpolate shadow fields to the curPatch

            // Note: solidContactPatch().zone() always returns the master zone
            // and solidContactPatch().shadowZones() always returns the slave
            // zones

            scalarField shadowPatchKOnCurPatch;
            scalarField shadowPatchRhoCOnCurPatch;

            if (master_)
            {
                const scalarField shadowZoneK =
                    solidContactPatch().shadowZones()[shadI].patchFaceToGlobal
                    (
                        shadowPatchK
                    );

                const scalarField shadowZoneRhoC =
                    solidContactPatch().shadowZones()[shadI].patchFaceToGlobal
                    (
                        shadowPatchRhoC
                    );

                const scalarField shadowZoneKOnCurPatch =
                    solidContactPatch().zoneToZones()[shadI].slaveToMaster
                    (
                        shadowZoneK
                    );

                const scalarField shadowZoneRhoCOnCurPatch =
                    solidContactPatch().zoneToZones()[shadI].slaveToMaster
                    (
                        shadowZoneRhoC
                    );
                shadowPatchKOnCurPatch =
                    solidContactPatch().zone().globalFaceToPatch
                    (
                        shadowZoneKOnCurPatch
                    );

                shadowPatchRhoCOnCurPatch =
                    solidContactPatch().zone().globalFaceToPatch
                    (
                        shadowZoneRhoCOnCurPatch
                    );
            }
            else
            {
                const scalarField shadowZoneK =
                    solidContactPatch().zone().patchFaceToGlobal
                    (
                        shadowPatchK
                    );

                const scalarField shadowZoneRhoC =
                    solidContactPatch().zone().patchFaceToGlobal
                    (
                        shadowPatchRhoC
                    );

                const scalarField shadowZoneKOnCurPatch =
                    solidContactPatch().zoneToZoneForThisSlave().masterToSlave
                    (
                        shadowZoneK
                    );

                const scalarField shadowZoneRhoCOnCurPatch =
                    solidContactPatch().zoneToZoneForThisSlave().masterToSlave
                    (
                        shadowZoneRhoC
                    );

                shadowPatchKOnCurPatch =
                    solidContactPatch().shadowZones()[shadI].globalFaceToPatch
                    (
                        shadowZoneKOnCurPatch
                    );

                shadowPatchRhoCOnCurPatch =
                    solidContactPatch().shadowZones()[shadI].globalFaceToPatch
                    (
                        shadowZoneRhoCOnCurPatch
                    );
            }

            const scalarField curPatchKsi =
                Foam::sqrt((shadowPatchRhoCOnCurPatch/curPatchRhoC)
                *curPatchK*shadowPatchKOnCurPatch)
              + curPatchK;

            tfrictionFluxRateForThisPatch() +=
                solidContactPatch().contactPerShadow()[shadI]
                *fricFlux*(curPatchK/curPatchKsi);
        }
    }

    return tfrictionFluxRateForThisPatch;
}


void Foam::thermalContactFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);

    alpha_.autoMap(m);
}


void Foam::thermalContactFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const thermalContactFvPatchScalarField& dmptf =
        refCast<const thermalContactFvPatchScalarField>(ptf);

    alpha_.rmap(dmptf.alpha_, addr);
}


void Foam::thermalContactFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Delete ggi interpolation at the begining of each time step
    if (curTimeIndex_ != db().time().timeIndex())
    {
        checkConsistentMaster();
        curTimeIndex_ = db().time().timeIndex();
    }

    // Lookup the T field
    // const volScalarField& field =
    //     db().lookupObject<volScalarField>
    //     (
    //         this->dimensionedInternalField().name()
    //     );

    // Lookup the thermal conductivity field
    const volScalarField& k = db().lookupObject<volScalarField>("k");

    // K on the current patch
    const scalarField& curPatchK = k.boundaryField()[patch().index()];

    // Loop through all shadow patches
    // Note: a slave patch will have only one shadow patch (i.e. the master),
    // whereas the master can have multiple shadow patches (i.e. multiple
    // slaves)
    const wordList& shadowPatchNames = this->shadowPatchNames();

    // Accumulate the normal gradient field in the contact areas
    scalarField curPatchSnGradInContactArea(patch().size(), 0.0);

    forAll(shadowPatchNames, shadI)
    {
        // Create the shadow zone temperature field
        scalarField shadowPatchTOnCurPatch;
        if (master())
        {
            // Note: solidContactPatch().zone() always returns the master zone
            // and solidContactPatch().shadowZones() always returns the slave
            // zones
            const scalarField shadowZoneT =
                solidContactPatch().shadowZones()[shadI].patchFaceToGlobal
                (
                    shadowPatchField(shadI)
                );

            // Interpolate shadow temperature field to the current zone
            const scalarField shadowZoneTOnCurPatch =
                solidContactPatch().zoneToZones()[shadI].slaveToMaster
                (
                    shadowZoneT
                );

            // Create the temperature field for this patch
            shadowPatchTOnCurPatch =
                solidContactPatch().zone().globalFaceToPatch
                (
                    shadowZoneTOnCurPatch
                );
        }
        else
        {
            // Note: zone() is the master global patch
            // Note 2: the slave always has only one master i.e.
            // one shadowPatchField
            const scalarField shadowZoneT =
                solidContactPatch().zone().patchFaceToGlobal
                (
                    shadowPatchField()
                );

            // Interpolate shadow temperature field to the current zone
            const scalarField shadowZoneTOnCurPatch =
                solidContactPatch().zoneToZoneForThisSlave().masterToSlave
                (
                    shadowZoneT
                );

            // Create the temperature field for this patch
            shadowPatchTOnCurPatch =
                solidContactPatch().shadowZones()[shadI].globalFaceToPatch
                (
                    shadowZoneTOnCurPatch
                );
        }

        // Calculate current contact conductance
        const scalarField curPatchH = Hc();

        // Calculate the heat flux through the current patch (in the contact
        // area)
        const scalarField curPatchFluxInContactArea =
          - curPatchH*(shadowPatchTOnCurPatch - *this)
          - frictionFluxRateForThisPatch();

        // Get the contact indicator field for this contact pair
        const scalarField contactPerShadow =
            solidContactPatch().contactPerShadow()[shadI];

        // Convert flux to normal gradient
        // fluxInContactArea == k*snGradTInContactArea
        // therefore
        // snGradTInContactArea == fluxInContactArea/k

        // Accumulate for all contact pairs
        curPatchSnGradInContactArea +=
            -contactPerShadow*curPatchFluxInContactArea/curPatchK;
    }

    // Next, we will add the thermal convections contributions for regions not
    // in contact:
    // k*snGradTNotInContact = -alpha*(T - Tinf)
    // i.e. heat flux within solid == heat flux due to convection at the
    // surface, therefore:
    // snGradTsnGradTNotInContact = -(alpha/k)*(T - Tinf)
    // Info<< nl
    //     << solidContactPatch().contact().size() << nl
    //     << alpha_.size() << nl
    //     << this->size() << nl
    //     << curPatchK.size() << endl;
    const scalarField curPatchPatchSnGradNotInContact =
        (1.0 - solidContactPatch().contact())
       *(-alpha_/curPatchK)*(*this - Tinf_);

    // Merge contributions for "in contact" and "not in contact" regions
    const scalarField curPatchPatchSnGrad =
        curPatchSnGradInContactArea + curPatchPatchSnGradNotInContact;

    // Set gradient using under-relaxation
    gradient() =
        underRelaxation()*curPatchPatchSnGrad
      + (1.0 - underRelaxation())*gradient();

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::thermalContactFvPatchScalarField::evaluate
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

    // Unit normals
    const vectorField n = patch().nf();

    // Delta vectors
    const vectorField delta = patch().delta();

    // Correction vectors
    const vectorField k = (I - sqr(n)) & delta;

    Field<scalar>::operator=
    (
        patchInternalField()
      + (k & gradT.patchInternalField())
      + gradient()/patch().deltaCoeffs()
    );

    fvPatchField<scalar>::evaluate();
}


void Foam::thermalContactFvPatchScalarField::write(Ostream& os) const
{
    os.writeKeyword("master")
        << master_ << token::END_STATEMENT << nl;
    alpha_.writeEntry("alpha", os);
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
        thermalContactFvPatchScalarField
    );
}

// ************************************************************************* //
