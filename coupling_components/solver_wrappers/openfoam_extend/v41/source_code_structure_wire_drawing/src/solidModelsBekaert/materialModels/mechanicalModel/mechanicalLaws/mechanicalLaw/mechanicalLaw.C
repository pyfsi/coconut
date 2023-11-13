/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    mechanicalLaw

Description
    Material mechanical for solids.

\*---------------------------------------------------------------------------*/

#include "mechanicalLaw.H"
#include "volFields.H"
#include "fvc.H"
#include "mechanicalModel.H"
#include "ggiPolyPatch.H"
#include "regionCouplePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(mechanicalLaw, 0);
defineRunTimeSelectionTable(mechanicalLaw, dictionary);


// * * * * * * * * * * * Private Member functions * * * * * * * * * * * * * * //

void mechanicalLaw::calcCurMaterial() const
{
    if (curMaterialPtr_)
    {
        FatalErrorIn
        (
            "const volScalarField& plasticityStressReturn::curMaterial() const"
        )   << "pointer already set" << abort(FatalError);
    }

    const fvMesh& mesh = this->mesh();

    if (mesh.foundObject<volScalarField>("materials"))
    {
        curMaterialPtr_ =
            new volScalarField
            (
                IOobject
                (
                    // PC, 21Jul18. remove '_' to avoid confusion with oldTime
                    // fields
                    "curMaterial" + Foam::name(curMatIndex_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0", dimless, 0.0)
            );

        volScalarField& curMaterial = *curMaterialPtr_;

        scalarField& curMaterialI = curMaterial.internalField();

        // Check the current material index has been set
        if (curMatIndex_ < 0)
        {
            FatalErrorIn("void mechanicalLaw::calcCurMaterial() const")
                << "The current material index has not been set"
                << abort(FatalError);
        }

        // Lookup materials index field
        const volScalarField& materials =
            mesh.objectRegistry::lookupObject<volScalarField>("materials");
        const scalarField& materialsI = materials.internalField();

        forAll(materialsI, cellI)
        {
            if (label(materialsI[cellI]) == curMatIndex_)
            {
                curMaterialI[cellI] = 1.0;
            }
        }

        forAll(materials.boundaryField(), patchI)
        {
            if
            (
                !materials.boundaryField()[patchI].coupled()
             && mesh.boundaryMesh()[patchI].type() != emptyPolyPatch::typeName
            )
            {
                scalarField& pCurMaterial =
                    curMaterial.boundaryField()[patchI];
                const labelList& faceCells =
                    mesh.boundaryMesh()[patchI].faceCells();

                forAll(pCurMaterial, faceI)
                {
                    if (label(materialsI[faceCells[faceI]]) == curMatIndex_)
                    {
                        pCurMaterial[faceI] = 1.0;
                    }
                }
            }
        }

        curMaterial.correctBoundaryConditions();

        // ZT, 16/09/16
        // Correct material indicator for ggi and region couple patches
        // This is required for materialGgi and materialCoupling fvPatchField-s
        forAll(materials.boundaryField(), patchI)
        {
            if
            (
                isA<ggiPolyPatch>(mesh.boundaryMesh()[patchI])
             || isA<regionCouplePolyPatch>(mesh.boundaryMesh()[patchI])
            )
            {
                scalarField& pCurMaterial =
                    curMaterial.boundaryField()[patchI];
                const labelList& faceCells =
                    mesh.boundaryMesh()[patchI].faceCells();

                forAll(pCurMaterial, faceI)
                {
                    if (label(materialsI[faceCells[faceI]]) == curMatIndex_)
                    {
                        pCurMaterial[faceI] = 1.0;
                    }
                    else
                    {
                        pCurMaterial[faceI] = 0.0;
                    }
                }
            }
        }
    }
    else
    {
        curMaterialPtr_ =
            new volScalarField
            (
                IOobject
                (
                    "curMaterial",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("1.0", dimless, 1.0),
                zeroGradientFvPatchScalarField::typeName
            );
    }

    // Info << curMatIndex_ << ", "
    //      << curMaterial().boundaryField() << endl;
}


void mechanicalLaw::calcCurMaterialf() const
{
    if (curMaterialfPtr_)
    {
        FatalErrorIn
        (
            "void mechanicalLaw::calcCurMaterialf() const"
        )   << "pointer already set" << abort(FatalError);
    }

    // Interpolate the cell-field
    curMaterialfPtr_ =
        new surfaceScalarField
        (
            "curMaterialf_" + Foam::name(curMatIndex_),
            fvc::interpolate(curMaterial())
        );

    WarningIn("void mechanicalLaw::calcCurMaterialf() const")
        << "interpolating curMaterial field to faces" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mechanicalLaw::mechanicalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const label lawIndex
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    curMaterialPtr_(NULL),
    curMaterialfPtr_(NULL),
    curMatIndex_(lawIndex)
{}


// * * * * * * * * * * * * * * * * * Destructors  * * * * * * * * * * * * * * //

mechanicalLaw::~mechanicalLaw()
{
    deleteDemandDrivenData(curMaterialPtr_);
    deleteDemandDrivenData(curMaterialfPtr_);
}

// * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * //


tmp<volScalarField> Foam::mechanicalLaw::bulkModulus() const
{
    notImplemented(type() + "::bulkModulus()");

    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
              "bulkModulus",
              mesh().time().timeName(),
              mesh(),
              IOobject::NO_READ,
              IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    return tresult;
}


tmp<volScalarField> Foam::mechanicalLaw::RhieChowScaleFactor() const
{
    // Each mechanicalLaw may override this function to define their
    // own scale factor
    // Defaults to 1.0 (full stabilisation)
    tmp<volScalarField> tresult
    (
        new volScalarField
        (
            IOobject
            (
              "RhieChowScaleFactor",
              mesh().time().timeName(),
              mesh(),
              IOobject::NO_READ,
              IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    return tresult;
}


const Foam::volScalarField& Foam::mechanicalLaw::curMaterial() const
{
    if (!curMaterialPtr_)
    {
        calcCurMaterial();
    }

    return *curMaterialPtr_;
}


Foam::volScalarField& Foam::mechanicalLaw::curMaterial()
{
    if (!curMaterialPtr_)
    {
        calcCurMaterial();
    }

    return *curMaterialPtr_;
}


const Foam::surfaceScalarField& Foam::mechanicalLaw::curMaterialf() const
{
    if (!curMaterialfPtr_)
    {
        calcCurMaterialf();
    }

    return *curMaterialfPtr_;
}


Foam::surfaceScalarField& Foam::mechanicalLaw::curMaterialf()
{
    if (!curMaterialfPtr_)
    {
        calcCurMaterialf();
    }

    return *curMaterialfPtr_;
}


Foam::scalar Foam::mechanicalLaw::residual()
{
    // Calculate residual based on change in stress tensor

    if (mesh().foundObject<volSymmTensorField>("tau"))
    {
        const volSymmTensorField& sigma =
            mesh().lookupObject<volSymmTensorField>("tau");

        return gMax
        (
            mag
            (
                sigma.internalField()
              - sigma.prevIter().internalField()
            )
        )/gMax(SMALL + mag(sigma.prevIter().internalField()));
    }
    else if (mesh().foundObject<volSymmTensorField>("sigma"))
    {
        const volSymmTensorField& sigma =
            mesh().lookupObject<volSymmTensorField>("sigma");

        return gMax
        (
            mag
            (
                sigma.internalField()
              - sigma.prevIter().internalField()
            )
        )/gMax(SMALL + mag(sigma.prevIter().internalField()));
    }

    return 0;
}


void Foam::mechanicalLaw::updateYieldStress()
{}

void Foam::mechanicalLaw::resetYieldStress()//MATHIEU
{}

Foam::scalar Foam::mechanicalLaw::newDeltaT() const
{
    // Default to a large number
    return mesh_.time().endTime().value();
}


} // End namespace Foam

// ************************************************************************* //
