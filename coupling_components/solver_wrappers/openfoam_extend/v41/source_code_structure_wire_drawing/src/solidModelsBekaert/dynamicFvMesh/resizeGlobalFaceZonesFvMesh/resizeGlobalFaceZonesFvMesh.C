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

#include "resizeGlobalFaceZonesFvMesh.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "solidContactFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(resizeGlobalFaceZonesFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh, resizeGlobalFaceZonesFvMesh, IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::resizeGlobalFaceZonesFvMesh::resizeGlobalFaceZonesFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io),
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
        ).subDict(typeName + "Coeffs")
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
{
    Info<< "Creating resizeGlobalFaceZonesFvMesh" << nl
        << "    boundBox: " << bb_ << nl
        << "    updateIndicatorField: " << updateIndicator_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::resizeGlobalFaceZonesFvMesh::~resizeGlobalFaceZonesFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::resizeGlobalFaceZonesFvMesh::update()
{
    // We will remove all current face zones
    // This is to avoid problems with faces being present in multiple face zones
    gfz_.removeFaceZones();

    // Lookup DU field
    const volVectorField& DU = this->lookupObject<volVectorField>("DU");

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
                        indicator_.boundaryField()[patchID][faceID - start] =
                            1.0;
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
                scalarField& pIndicator = indicator_.boundaryField()[patchI];
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

    return true;
}


// ************************************************************************* //
