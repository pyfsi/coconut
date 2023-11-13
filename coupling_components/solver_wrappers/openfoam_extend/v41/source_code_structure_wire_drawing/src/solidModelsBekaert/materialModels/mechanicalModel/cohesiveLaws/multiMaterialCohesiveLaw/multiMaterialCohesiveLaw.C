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

\*---------------------------------------------------------------------------*/

#include "multiMaterialCohesiveLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvc.H"
#include "crackerFvMesh.H"
#include "multiMaterial.H"
#include "mechanicalModel.H"
#include "cohesiveFvPatch.H"
#include "cohesiveLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiMaterialCohesiveLaw, 0);
    addToRunTimeSelectionTable
    (
        cohesiveLaw, multiMaterialCohesiveLaw, dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::multiMaterialCohesiveLaw::indicator
(
    const label i
) const
{
    const scalarField& mat = materials_.internalField();

    tmp<scalarField> tresult(new scalarField(mat.size(), 0.0));
    scalarField& result = tresult();

    forAll (mat, matI)
    {
        if (mat[matI] > i - SMALL && mat[matI] < i + 1 - SMALL)
        {
            result[matI] = 1.0;
        }
    }

    return tresult;
}


Foam::scalar Foam::multiMaterialCohesiveLaw::indicator
(
    const label index,
    const label cellID
) const
{
    const scalar mat = materials_.internalField()[cellID];
    scalar result = 0.0;

    if (mat > index - SMALL && mat < index + 1 - SMALL)
      {
    result = 1.0;
      }

    return result;
}

Foam::scalar Foam::multiMaterialCohesiveLaw::interfaceID
(
    const label mat1,
    const label mat2
) const
{
    word mat1name = (*this)[mat1].name();
    word mat2name = (*this)[mat2].name();

    word interfaceName("interface_"+mat1name+"_"+mat2name);
    label interfaceLawID = -1;
    forAll(interfaceCohesiveLaws_, lawI)
    {
        if (interfaceCohesiveLaws_[lawI].name() == interfaceName)
        {
            interfaceLawID = lawI;
            break;
        }
    }
    if (interfaceLawID == -1)
    {
        // flip name
        interfaceName = word("interface_"+mat2name+"_"+mat1name);
        forAll(interfaceCohesiveLaws_, lawI)
        {
            if (interfaceCohesiveLaws_[lawI].name() == interfaceName)
            {
                interfaceLawID = lawI;
                break;
            }
        }
        if (interfaceLawID == -1)
        {
            FatalError
                << "Cannot find cohesive interfaceLaw "
                    << word("interface_"+mat1name+"_"+mat2name) << " or "
                    << interfaceName << nl
                    << "One of these should be defined!"
                    << exit(FatalError);
        }
    }

    return interfaceLawID;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::multiMaterialCohesiveLaw::multiMaterialCohesiveLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    cohesiveLaw(name, mesh, dict),
    PtrList<cohesiveLaw>(),
    materials_
    (
        mesh.objectRegistry::lookupObject<volScalarField>("materials")
    ),
    interfaceCohesiveLaws_()
{
    PtrList<cohesiveLaw>& laws = *this;

    PtrList<entry> lawEntries(dict.lookup("laws"));
    laws.setSize(lawEntries.size());

    forAll (laws, lawI)
    {
        laws.set
        (
            lawI,
            cohesiveLaw::New
            (
                lawEntries[lawI].keyword(),
                mesh,
                lawEntries[lawI].dict()
            )
        );
    }

    if
    (
        min(materials_).value() < 0
     || max(materials_).value() > laws.size() + SMALL
    )
    {
        FatalErrorIn
        (
            "multiMaterialCohesiveLaw::multiMaterialCohesiveLaw\n"
            "(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Invalid definition of material indicator field.  "
            << "Number of materials: " << laws.size()
            << " max index: " << max(materials_)
            << ".  Should be " << laws.size() - 1
            << abort(FatalError);
    }

    // cohesive laws must be the same size as mechanical laws
    const mechanicalModel& conModel =
        mesh.objectRegistry::lookupObject<mechanicalModel>
        (
            "mechanicalProperties"
        );
    if (conModel.law().type() == "multiMaterial")
    {
        const multiMaterial& mulMatLaw =
            refCast<const multiMaterial>(conModel.law());

        if (laws.size() != mulMatLaw.size())
        {
            FatalError
                << "There should be the same number of cohesive laws "
                << "as mechanical laws" << nl
                << "Currently there are " << mulMatLaw.size()
                << " mechanical laws "
                << "and " << laws.size() << " cohesive laws!"
                << abort(FatalError);
        }
    }

    // set interfaceCohesiveLaws_
    PtrList<entry> interfaceLawEntries(dict.lookup("interfaceLaws"));
    // if (interfaceLawEntries.size() != int(laws.size()*(laws.size()-1)/2))
    if
    (
        mag(interfaceLawEntries.size() - (laws.size()*(laws.size()-1)/2))
        > SMALL
    )
    {
        // number of interfaces is a triangular number of number of materials
        // ((n)*(n-1)/2)
        FatalError
            << "There are " << interfaceLawEntries.size()
            << " interface cohesive"
            << " laws defined, but there should be "
            << (laws.size()*(laws.size()-1)/2)
            << "," << nl
            << "as there are " << laws.size()
            << " materials in cohesive laws" << nl
            << abort(FatalError);
    }
    interfaceCohesiveLaws_.setSize(interfaceLawEntries.size());
    forAll (interfaceCohesiveLaws_, lawI)
    {
        interfaceCohesiveLaws_.set
        (
            lawI,
           cohesiveLaw::New
       (
        interfaceLawEntries[lawI].keyword(),
        mesh,
        interfaceLawEntries[lawI].dict()
            )
       );
      }


    // Set materialsSurf
    // materialsSurf_ = fvc::interpolate(materials_);
    // forAll(mesh().boundary(), patchi)
    //   {
    //      materialsSurf_.boundaryField()[patchi] =
    //        materials_.boundaryField()[patchi].patchInternalField();
    //   }
    // // Fix interface values
    // const labelList owner = mesh().owner();
    // const labelList neighbour = mesh().neighbour();
    // forAll (materialsSurf_.internalField(), faceI)
    //   {
    //      // round value to integer and check difference
    //      // if it is small then the face is not on a multi-material
    //      // interface
    //      scalar matIDscalar = materialsSurf_.internalField()[faceI];
    //      label matID = int(matIDscalar);
    //      if (mag(matIDscalar - matID) > SMALL)
    //        {
    //          // find which interface it is on
    //          const label own = owner[faceI];
    //          const label nei = neighbour[faceI];

    //          materialsSurf_.internalField()[faceI] =
    //            interfaceID(materials_[own], materials_[nei]) + laws.size();
    //        }
    //   }

    // philipc
    // processor boundaries are meant to hold the patchNeighbourField
    // but the proc boundaries are interpolated by decomposePar
    // so we must correct them. Now the proc boundaries hold the
    // patchNiehgbourField
    //materials_.correctBoundaryConditions(); // done by mechanicalLaw
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiMaterialCohesiveLaw::~multiMaterialCohesiveLaw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::multiMaterialCohesiveLaw::materials() const
{
    return materials_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiMaterialCohesiveLaw::sigmaMax() const
{
    tmp<surfaceScalarField> tresult
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "sigmaMax",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimForce/dimArea, 0)
            )
        );
    surfaceScalarField& result = tresult();

    // Accumulate data for all fields
    const PtrList<cohesiveLaw>& laws = *this;
    volScalarField indic
        (
            IOobject
            (
                "indic",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );

    // We interpolate indicate field and then fix interface fields after
    forAll (laws, lawI)
    {
        indic.internalField() = indicator(lawI)();
        surfaceScalarField indicatorSurf = fvc::interpolate(indic);

        // Fix boundary fields
        forAll(mesh().boundary(), patchi)
        {
            indicatorSurf.boundaryField()[patchi] =
                indic.boundaryField()[patchi].patchInternalField();
        }

        result += indicatorSurf*laws[lawI].sigmaMax()();
    }

    // Fix interfaces
    // forAll surfaces check if surface is a material interface
    // i.e. material indicator should read non-integer
    // Get the two materials it is an interface of and
    // look up value of sigmaMax in dictionary
    // Overwrite existing value on surface with new value

    const labelList& owner = mesh().owner();
    const labelList& neighbour = mesh().neighbour();

    forAll(result.internalField(), faceI)
    {
        label ownerMat = label(materials_[owner[faceI]]);
        label neighbourMat = label(materials_[neighbour[faceI]]);

        if (ownerMat != neighbourMat)
        {
            result.internalField()[faceI] =
                interfaceSigmaMax(ownerMat, neighbourMat);
        }
    }

    forAll(mesh().boundary(), patchI)
    {
        // Faces in the cohesive patch may be on an interface so we must check
        if (mesh().boundary()[patchI].type() == cohesiveFvPatch::typeName)
        {
            const crackerFvMesh& crackerMesh =
                refCast<const crackerFvMesh>(mesh());

            // We must use owner values
            scalarField localPatchMaterials =
                materials_.boundaryField()[patchI].patchInternalField();
            scalarField globalPatchMaterials =
                crackerMesh.globalCrackField(localPatchMaterials);
            // crackerMesh.globalCrackField(materials_.boundaryField()[patchI]);
            const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
            label globalIndex = crackerMesh.localCrackStart();

            forAll(mesh().boundaryMesh()[patchI], facei)
            {
                label ownerMat = label(globalPatchMaterials[globalIndex]);
                label neighbourMat =
                    label(globalPatchMaterials[gcfa[globalIndex]]);

                if (ownerMat != neighbourMat)
                {
                    result.boundaryField()[patchI][facei] =
                        interfaceSigmaMax(ownerMat, neighbourMat);
                }
                globalIndex++;
            }
        }
        else if (mesh().boundary()[patchI].coupled())
        {
            const scalarField ownerMatField =
                materials_.boundaryField()[patchI].patchInternalField();
            const scalarField neighbourMatField =
                materials_.boundaryField()[patchI].patchNeighbourField();

            forAll(mesh().boundaryMesh()[patchI], facei)
            {
                label ownerMat = label(ownerMatField[facei]);
                label neighbourMat = label(neighbourMatField[facei]);

                if (ownerMat != neighbourMat)
                {
                    result.boundaryField()[patchI][facei] =
                        interfaceSigmaMax(ownerMat, neighbourMat);
                }
            }
        }
    }

    return tresult;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiMaterialCohesiveLaw::tauMax() const
{
    tmp<surfaceScalarField> tresult
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "tauMax",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimForce/dimArea, 0)
            )
        );
    surfaceScalarField& result = tresult();

    const PtrList<cohesiveLaw>& laws = *this;
    volScalarField indic
        (
            IOobject
            (
                "indic",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );
    forAll (laws, lawI)
    {
        // interface fields will not be correct but we fix them after
        indic.internalField() = indicator(lawI)();
        surfaceScalarField indicatorSurf = fvc::interpolate(indic);
        // fix boundary fields
        forAll(mesh().boundary(), patchi)
        {
            indicatorSurf.boundaryField()[patchi] =
                indic.boundaryField()[patchi].patchInternalField();
        }
        result += indicatorSurf*laws[lawI].tauMax()();
    }

    // forAll surfaces check if surface is a material interface
    // material indicator should read non integer
    // Get the two materials it is an interface of
    // Look up value of tauMaxf in dictionary
    // Overwrite existing value on surface with new value

    const labelList& owner = mesh().owner();
    const labelList& neighbour = mesh().neighbour();

    forAll(result.internalField(), faceI)
    {
        label ownerMat = label(materials_[owner[faceI]]);
        label neighbourMat = label(materials_[neighbour[faceI]]);

        if (ownerMat != neighbourMat)
        {
            result.internalField()[faceI] =
                interfaceTauMax(ownerMat, neighbourMat);
        }
    }

    forAll(mesh().boundary(), patchI)
    {

        if (mesh().boundary()[patchI].type() == cohesiveFvPatch::typeName)
        {
            //label size = (mesh().boundary()[patchI].size())/2;
            const crackerFvMesh& crackerMesh =
                refCast<const crackerFvMesh>(mesh());
            // we must use owner values
            scalarField localPatchMaterials =
                materials_.boundaryField()[patchI].patchInternalField();
            scalarField globalPatchMaterials =
                crackerMesh.globalCrackField(localPatchMaterials);
            // crackerMesh.globalCrackField(materials_.boundaryField()[patchI]);
            const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
            label globalIndex = crackerMesh.localCrackStart();

            forAll(mesh().boundaryMesh()[patchI], facei)
            {
                label ownerMat = label(globalPatchMaterials[globalIndex]);
                label neighbourMat =
                    label(globalPatchMaterials[gcfa[globalIndex]]);

                if (ownerMat != neighbourMat)
                {
                    result.boundaryField()[patchI][facei] =
                        interfaceTauMax(ownerMat, neighbourMat);
                }
                globalIndex++;
            }
        }
        else if (mesh().boundary()[patchI].coupled())
        {
            const scalarField ownerMatField =
                materials_.boundaryField()[patchI].patchInternalField();
            const scalarField neighbourMatField =
                materials_.boundaryField()[patchI].patchNeighbourField();

            forAll(mesh().boundaryMesh()[patchI], facei)
            {
                label ownerMat = label(ownerMatField[facei]);
                label neighbourMat = label(neighbourMatField[facei]);

                if (ownerMat != neighbourMat)
                {
                    result.boundaryField()[patchI][facei] =
                        interfaceTauMax(ownerMat, neighbourMat);
                }
            }
        }
    }

    return tresult;
}

Foam::tmp<Foam::surfaceScalarField> Foam::multiMaterialCohesiveLaw::GIc() const
{
    tmp<surfaceScalarField> tresult
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "GIc",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimForce*dimLength/dimArea, 0)
            )
        );
    surfaceScalarField& result = tresult();

    const PtrList<cohesiveLaw>& laws = *this;
    volScalarField indic
        (
            IOobject
            (
                "indic",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );
    forAll (laws, lawI)
    {
        // interface fields will not be correct but we fix them after
        indic.internalField() = indicator(lawI)();
        surfaceScalarField indicatorSurf = fvc::interpolate(indic);
        // fix boundary fields
        forAll(mesh().boundary(), patchi)
        {
            indicatorSurf.boundaryField()[patchi] =
                indic.boundaryField()[patchi].patchInternalField();
        }
        result += indicatorSurf*laws[lawI].GIc()();
    }

    const labelList& owner = mesh().owner();
    const labelList& neighbour = mesh().neighbour();

    forAll(result.internalField(), faceI)
    {
        label ownerMat = label(materials_[owner[faceI]]);
        label neighbourMat = label(materials_[neighbour[faceI]]);

        if (ownerMat != neighbourMat)
        {
            result.internalField()[faceI] =
                interfaceGIc(ownerMat, neighbourMat);
        }
    }

    forAll(mesh().boundary(), patchI)
    {

        if (mesh().boundary()[patchI].type() == cohesiveFvPatch::typeName)
        {
            //label size = (mesh().boundary()[patchI].size())/2;
            // const labelList& fCells = mesh().boundary()[patchI].faceCells();
            scalarField localPatchMaterials =
                materials_.boundaryField()[patchI].patchInternalField();
            const crackerFvMesh& crackerMesh =
                refCast<const crackerFvMesh>(mesh());
            scalarField globalPatchMaterials =
                crackerMesh.globalCrackField(localPatchMaterials);
            const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
            label globalIndex = crackerMesh.localCrackStart();

            forAll(mesh().boundaryMesh()[patchI], facei)
            {
                label ownerMat = label(globalPatchMaterials[globalIndex]);
                label neighbourMat =
                    label(globalPatchMaterials[gcfa[globalIndex]]);

                if (ownerMat != neighbourMat)
                {
                    result.boundaryField()[patchI][facei] =
                        interfaceGIc(ownerMat, neighbourMat);
                }
                globalIndex++;
            }
        }
        else if (mesh().boundary()[patchI].coupled())
        {
            const scalarField ownerMatField =
                materials_.boundaryField()[patchI].internalField();
            const scalarField neighbourMatField =
                materials_.boundaryField()[patchI].patchNeighbourField();

            forAll(mesh().boundaryMesh()[patchI], facei)
            {
                label ownerMat = label(ownerMatField[facei]);
                label neighbourMat = label(neighbourMatField[facei]);

                if (ownerMat != neighbourMat)
                {
                    //result.boundaryField()[patchI][facei] = iterGIc();
                    result.boundaryField()[patchI][facei] =
                        interfaceGIc(ownerMat, neighbourMat);
                }
            }
        }
    }

    return tresult;
}

Foam::tmp<Foam::surfaceScalarField> Foam::multiMaterialCohesiveLaw::GIIc() const
{
    tmp<surfaceScalarField> tresult
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "GIIc",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimForce*dimLength/dimArea, 0)
            )
        );
    surfaceScalarField& result = tresult();

    const PtrList<cohesiveLaw>& laws = *this;
    volScalarField indic
        (
            IOobject
            (
                "indic",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );
    forAll (laws, lawI)
    {
        // interface fields will not be correct but we fix them after
        indic.internalField() = indicator(lawI)();
        surfaceScalarField indicatorSurf = fvc::interpolate(indic);
        // fix boundary fields
        forAll(mesh().boundary(), patchi)
        {
            indicatorSurf.boundaryField()[patchi] =
                indic.boundaryField()[patchi].patchInternalField();
        }
        result += indicatorSurf*laws[lawI].GIIc()();
    }

    const labelList& owner = mesh().owner();
    const labelList& neighbour = mesh().neighbour();

    forAll(result.internalField(), faceI)
    {
        label ownerMat = label(materials_[owner[faceI]]);
        label neighbourMat = label(materials_[neighbour[faceI]]);

        if (ownerMat != neighbourMat)
        {
            result.internalField()[faceI] =
                interfaceGIIc(ownerMat, neighbourMat);
        }
    }

    forAll(mesh().boundary(), patchI)
    {

        if (mesh().boundary()[patchI].type() == cohesiveFvPatch::typeName)
        {
            //label size = (mesh().boundary()[patchI].size())/2;
            // const labelList& fCells = mesh().boundary()[patchI].faceCells();
            scalarField localPatchMaterials =
                materials_.boundaryField()[patchI].patchInternalField();
            const crackerFvMesh& crackerMesh =
                refCast<const crackerFvMesh>(mesh());
            scalarField globalPatchMaterials =
                crackerMesh.globalCrackField(localPatchMaterials);
            const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
            label globalIndex = crackerMesh.localCrackStart();

            forAll(mesh().boundaryMesh()[patchI], facei)
            {
                label ownerMat = label(globalPatchMaterials[globalIndex]);
                label neighbourMat =
                    label(globalPatchMaterials[gcfa[globalIndex]]);

                if (ownerMat != neighbourMat)
                {
                    result.boundaryField()[patchI][facei] =
                        interfaceGIIc(ownerMat, neighbourMat);
                }
                globalIndex++;
            }
        }
        else if (mesh().boundary()[patchI].coupled())
        {
            const scalarField ownerMatField =
                materials_.boundaryField()[patchI].internalField();
            const scalarField neighbourMatField =
                materials_.boundaryField()[patchI].patchNeighbourField();

            forAll(mesh().boundaryMesh()[patchI], facei)
            {
                label ownerMat = label(ownerMatField[facei]);
                label neighbourMat = label(neighbourMatField[facei]);

                if (ownerMat != neighbourMat)
                {
                    result.boundaryField()[patchI][facei] =
                        interfaceGIIc(ownerMat, neighbourMat);
                }
            }
        }

    }

    return tresult;
}

Foam::scalar Foam::multiMaterialCohesiveLaw::interfaceSigmaMax
(
    double mat1,
    double mat2
) const
{
  label interfaceLawID = interfaceID(mat1, mat2);

  return interfaceCohesiveLaws_[interfaceLawID].sigmaMaxValue();
}

Foam::scalar Foam::multiMaterialCohesiveLaw::interfaceTauMax
(
    double mat1,
    double mat2
) const
{
  label interfaceLawID = interfaceID(mat1, mat2);

  return interfaceCohesiveLaws_[interfaceLawID].tauMaxValue();
}

Foam::scalar Foam::multiMaterialCohesiveLaw::interfaceGIc
(
    double mat1,
    double mat2
) const
{
  label interfaceLawID = interfaceID(mat1, mat2);

  return interfaceCohesiveLaws_[interfaceLawID].GIcValue();
}

Foam::scalar Foam::multiMaterialCohesiveLaw::interfaceGIIc
(
    double mat1,
    double mat2
) const
{
  label interfaceLawID = interfaceID(mat1, mat2);

  return interfaceCohesiveLaws_[interfaceLawID].GIIcValue();
}


void Foam::multiMaterialCohesiveLaw::damageTractions
(
    scalar& tN,
    scalar& tS,
    const scalar deltaN,
    const scalar deltaS,
    const scalar GI,
    const scalar GII,
    const label faceID,
    const scalarField& globalPatchMaterials
 ) const
{
    // Find out which cohesive law does the face belong to
    const crackerFvMesh& crackerMesh = refCast<const crackerFvMesh>(mesh());
    const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
    label ownerMat = label(globalPatchMaterials[faceID]);
    label neighbourMat = label(globalPatchMaterials[gcfa[faceID]]);

    if (ownerMat != neighbourMat)
    {
        label matID = interfaceID(ownerMat, neighbourMat);

        // face is on multi-material interface
        interfaceCohesiveLaws_[matID].damageTractions
            (tN, tS, deltaN, deltaS, GI, GII, faceID, globalPatchMaterials);
    }
    else
    {
        // face is within one material
        // call material law function
        label matID = ownerMat;
        (*this)[matID].damageTractions
            (tN, tS, deltaN, deltaS, GI, GII, faceID, globalPatchMaterials);
    }
}


//- Cohesive tractions from implicit damage procedure
void Foam::multiMaterialCohesiveLaw::damageTractions
(
    scalar& tN,
    scalar& tS,
    const scalar deltaN,
    const scalar deltaS,
    const scalar GI,
    const scalar GII,
    const label internalFaceID
) const
{
    // Find out which cohesive law does the face belong to
    label ownerMat = materials_[mesh().owner()[internalFaceID]];
    label neighbourMat = materials_[mesh().neighbour()[internalFaceID]];

    if (ownerMat != neighbourMat)
    {
        label matID = interfaceID(ownerMat, neighbourMat);

        // face is on multi-material interface
        interfaceCohesiveLaws_[matID].damageTractions
            (tN, tS, deltaN, deltaS, GI, GII, internalFaceID);
    }
    else
    {
        // face is within one material
        // call material law function
        label matID = ownerMat;
        (*this)[matID].damageTractions
            (tN, tS, deltaN, deltaS, GI, GII, internalFaceID);
    }
}


void Foam::multiMaterialCohesiveLaw::correct()
{
    PtrList<cohesiveLaw>& laws = *this;

    forAll (laws, lawI)
    {
        laws[lawI].correct();
    }
}


// ************************************************************************* //
