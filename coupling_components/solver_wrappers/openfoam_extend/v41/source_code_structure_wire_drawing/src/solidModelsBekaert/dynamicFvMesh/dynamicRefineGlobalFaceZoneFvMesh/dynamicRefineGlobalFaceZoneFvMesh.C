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

\*---------------------------------------------------------------------------*/

#include "dynamicRefineGlobalFaceZoneFvMesh.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "resErrorLaplacian.H"
#include "solidContactFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicRefineGlobalFaceZoneFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh, dynamicRefineGlobalFaceZoneFvMesh, IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::dynamicRefineGlobalFaceZoneFvMesh::dynamicRefineGlobalFaceZoneFvMesh
(
    const IOobject& io
)
:
    dynamicRefineFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict("dynamicRefineFvMeshCoeffs")
    ),
    residualError_
    (
        IOobject
        (
            "residualError",
            time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar("zero", dimLength, 0.0)
        //dimensionedScalar("zero", dimLength/dimTime, 0.0)
    ),
    bb_
    (
        dict_.lookupOrDefault<boundBox>
        (
            "boundBox",
            boundBox(vector::min, vector::max)
        )
    ),
    gfz_(*this, bb_),
    indicator_
    (
        IOobject
        (
            "globalFaceZones",
            time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    updateIndicator_
    (
        dict_.lookupOrDefault<Switch>("updateIndicatorField", true)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicRefineGlobalFaceZoneFvMesh::~dynamicRefineGlobalFaceZoneFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicRefineGlobalFaceZoneFvMesh::update()
{
    // Lookup fields from solver
    const volVectorField& DU = lookupObject<volVectorField>("DU");
    const volScalarField& twoMuLambda =
        lookupObject<volScalarField>("twoMuLambda");
    const surfaceScalarField& twoMuLambdaf =
        lookupObject<surfaceScalarField>("twoMuLambdaf");
    const volScalarField& relJ = lookupObject<volScalarField>("relJ");
    const volScalarField& J = lookupObject<volScalarField>("J");
    const volTensorField& relFinv = lookupObject<volTensorField>("relFinv");
    const volTensorField& gradDU = lookupObject<volTensorField>("grad(DU)");
    const volSymmTensorField& tau =
        lookupObject<volSymmTensorField>("tauKirchhoff");

    // Estimate spatial discretisation error
    errorEstimate<vector> DUEqnError
    (
      - resError::laplacian(twoMuLambdaf, DU)
      - fvc::div((relJ/J)*tau & relFinv.T(), "div(sigma)")
      + fvc::div(twoMuLambda*gradDU, "div(sigma)")
    );

    // Absolute error
    residualError_ = mag(DUEqnError.error());
    // Normalised error
    //residualError_ = mag(DUEqnError.residual());

    const bool topoChange = dynamicRefineFvMesh::update();

    if (topoChange)
    {
        // Remove all current face zones
        gfz_.removeFaceZones();

        // Re-create global face zones from contact patches
        forAll(DU.boundaryField(), patchI)
        {
            const word pType = DU.boundaryField()[patchI].type();

            if (pType == solidContactFvPatchVectorField::typeName)
            {
                gfz_.createFaceZone(patchI);

                solidContactFvPatchVectorField& contactPatch =
                    const_cast<solidContactFvPatchVectorField&>
                    (
                        refCast<const solidContactFvPatchVectorField>
                        (
                            DU.boundaryField()[patchI]
                        )
                    );

                // Force contact patch addressing to be re-calculated
                contactPatch.clearOut();
            }
        }

        // Reorder each face zone so that the faces are in the same order on all
        // procs; this allows parallel syncing of zone data
        forAll(this->faceZones(), zoneI)
        {
            gfz_.reorderFaceZone(zoneI);
        }

        // Set indicator field for visualisation of current zones
        if (updateIndicator_)
        {
            // Reset field to zero
            indicator_ = dimensionedScalar("zero", dimless, 0.0);

            const label nActiveFaces = this->nFaces();
            forAll(this->faceZones(), zoneI)
            {
                const faceZone& curZone = this->faceZones()[zoneI];

                forAll(curZone, fI)
                {
                    const label faceID = curZone[fI];

                    if (faceID < nActiveFaces)
                    {
                        const label patchID =
                            this->boundaryMesh().whichPatch(faceID);

                        if (patchID != -1)
                        {
                            const label start =
                                this->boundaryMesh()[patchID].start();
                            indicator_.boundaryField()[patchID]
                            [
                                faceID - start
                            ] = 1.0;
                        }
                    }
                }
            }

            // Set internal cell values
            scalarField& indicatorI = indicator_.internalField();
            forAll(indicator_.boundaryField(), patchI)
            {
                const word pType = indicator_.boundaryField()[patchI].type();

                if (!this->boundaryMesh()[patchI].coupled() && pType != "empty")
                {
                    scalarField& pIndicator =
                        indicator_.boundaryField()[patchI];
                    const unallocLabelList& faceCells =
                        this->boundaryMesh()[patchI].faceCells();

                    forAll(faceCells, faceI)
                    {
                        indicatorI[faceCells[faceI]] =
                            max
                            (
                                indicatorI[faceCells[faceI]],
                                pIndicator[faceI]
                            );
                    }
                }
            }
        }
    }

    return topoChange;
}


// ************************************************************************* //
