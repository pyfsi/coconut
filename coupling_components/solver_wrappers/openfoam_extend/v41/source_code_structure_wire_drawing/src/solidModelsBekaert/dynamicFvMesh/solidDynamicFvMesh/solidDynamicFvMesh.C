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

#include "solidDynamicFvMesh.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "symmetryPolyPatch.H"
#include "twoDPointCorrector.H"
#include "materialGgiFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidDynamicFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, solidDynamicFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solidDynamicFvMesh::makePointDU() const
{
    if (pointDUPtr_)
    {
        FatalErrorIn("void Foam::solidDynamicFvMesh::makePointDU() const")
            << "pointer not set" << endl;
    }

    wordList types
    (
        pMesh_.boundary().size(),
        calculatedFvPatchVectorField::typeName
    );

    // forAll(types, patchI)
    // {
    //     // pMesh boundary contains global patches
    //     if (patchI < this->boundary().size())
    //     {
    //         if
    //         (
    //             DU.boundaryField()[patchI].type()
    //          == "fixedDisplacementZeroShear"
    //         )
    //         {
    //             types[patchI] = fixedValueFvPatchVectorField::typeName;
    //         }
    //     }
    // }

    pointDUPtr_ =
        new pointVectorField
        (
            IOobject
            (
                "pointDU",
                time().timeName(),
                *this,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            pMesh_,
            dimensionedVector("zero", dimLength, vector::zero),
            types
        );
}


Foam::pointVectorField& Foam::solidDynamicFvMesh::pointDU()
{
    if (!pointDUPtr_)
    {
        makePointDU();
    }

    return *pointDUPtr_;
}


void Foam::solidDynamicFvMesh::calcGlobalFaceZones() const
{
    // Find global face zones
    if (globalFaceZonesPtr_)
    {
        FatalErrorIn
        (
            "void solidModel::calcGlobalFaceZones() const"
        )
            << "Global face zones already found"
            << abort(FatalError);
    }

    if (Pstream::parRun())
    {
        SLList<label> globalFaceZonesSet;

        // Previous method
        // const faceZoneMesh& faceZones = mesh().faceZones();
        // forAll(faceZones, zoneI)
        // {
        //     const faceZone& curFaceZone = faceZones[zoneI];
        //     forAll(curFaceZone, faceI)
        //     {
        //         // If unused face exist
        //         if (curFaceZone[faceI] >= mesh().nFaces())
        //         {
        //             globalFaceZonesSet.insert(zoneI);
        //             break;
        //         }
        //     }
        // }

        // New method: directly lookup globalFaceZones from decomposeParDict


        // For FSI cases, we need to look in a different location for the dict

        word decompDictName = "system/decomposeParDict";

        if
        (
            isDir
            (
                this->time().rootPath()/this->time().caseName()
                /"../system/solid"
            )
        )
        {
            decompDictName = "../system/solid/decomposeParDict";
        }

        Info<< "Reading decomposeParDict " << decompDictName << endl;

        IOdictionary decompDict
        (
            IOobject
            (
                decompDictName,
                this->time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        if (decompDict.found("globalFaceZones"))
        {
            wordList globalFaceZoneNames(decompDict.lookup("globalFaceZones"));

            const faceZoneMesh& faceZones = this->faceZones();

            forAll(globalFaceZoneNames, nameI)
            {
                const label zoneID =
                    faceZones.findZoneID(globalFaceZoneNames[nameI]);

                if (zoneID == -1)
                {
                    FatalErrorIn(type() + "::findGlobalFaceZones")
                        << "Cannot find globalFaceZone:"
                        << " " << globalFaceZoneNames[nameI]
                        << abort(FatalError);
                }

                globalFaceZonesSet.insert(zoneID);
            }

            globalFaceZonesPtr_ = new labelList(globalFaceZonesSet);
        }
        else
        {
            globalFaceZonesPtr_ = new labelList(0);
        }
    }
    else
    {
        globalFaceZonesPtr_ = new labelList(0);
    }
}


void Foam::solidDynamicFvMesh::calcGlobalToLocalFaceZonePointMap() const
{
    // Find global face zones
    if (globalToLocalFaceZonePointMapPtr_)
    {
        FatalErrorIn
        (
            "void solidModel::calcGlobalToLocalFaceZonePointMap() const"
        )   << "Global to local face zones point map already exists"
            << abort(FatalError);
    }

    globalToLocalFaceZonePointMapPtr_ =
        new labelListList(globalFaceZones().size());

    labelListList& globalToLocalFaceZonePointMap =
        *globalToLocalFaceZonePointMapPtr_;

    const labelList& globalFaceZones = this->globalFaceZones();

    forAll(globalFaceZones, zoneI)
    {
        const label curZoneID = globalFaceZones[zoneI];

        Info<< "Creating faceMap for globalFaceZones "
            << this->faceZones()[curZoneID].name()<< endl;

        labelList curMap(this->faceZones()[curZoneID]().nPoints(), -1);

        vectorField fzGlobalPoints =
            this->faceZones()[curZoneID]().localPoints();

        // Set all slave points to zero because only the master order is used
        if(!Pstream::master())
        {
            fzGlobalPoints = vector::zero;
        }

        // Pass points to all procs
        reduce(fzGlobalPoints, sumOp<vectorField>());

        // Now every proc has the master's list of FZ points
        // every proc must now find the mapping from their local FZ points to
        // the global FZ points

        const vectorField& fzLocalPoints =
            this->faceZones()[curZoneID]().localPoints();

        const edgeList& fzLocalEdges =
            this->faceZones()[curZoneID]().edges();

        const labelListList& fzPointEdges =
            this->faceZones()[curZoneID]().pointEdges();

        scalarField minEdgeLength(fzLocalPoints.size(), GREAT);

        forAll(minEdgeLength, pI)
        {
            const labelList& curPointEdges = fzPointEdges[pI];

            forAll(curPointEdges, eI)
            {
                const scalar Le =
                    fzLocalEdges[curPointEdges[eI]].mag(fzLocalPoints);

                if (Le < minEdgeLength[pI])
                {
                    minEdgeLength[pI] = Le;
                }
            }
        }

        forAll(fzGlobalPoints, globalPointI)
        {
            boolList visited(fzLocalPoints.size(), false);

            forAll(fzLocalPoints, procPointI)
            {
                if (!visited[procPointI])
                {
                    visited[procPointI] = true;

                    label nextPoint = procPointI;

                    scalar curDist =
                        mag
                        (
                            fzLocalPoints[nextPoint]
                          - fzGlobalPoints[globalPointI]
                        );

                    if (curDist < 1e-4*minEdgeLength[nextPoint])
                    {
                        curMap[globalPointI] = nextPoint;
                        break;
                    }

                    label found = false;

                    while (nextPoint != -1)
                    {
                        const labelList& nextPointEdges =
                            fzPointEdges[nextPoint];

                        scalar minDist = GREAT;
                        label index = -1;
                        forAll(nextPointEdges, edgeI)
                        {
                            label curNgbPoint =
                                fzLocalEdges[nextPointEdges[edgeI]]
                               .otherVertex(nextPoint);

                            if (!visited[curNgbPoint])
                            {
                                visited[curNgbPoint] = true;

                                scalar curDist =
                                    mag
                                    (
                                        fzLocalPoints[curNgbPoint]
                                      - fzGlobalPoints[globalPointI]
                                    );

                                if (curDist < 1e-4*minEdgeLength[curNgbPoint])
                                {
                                    curMap[globalPointI] = curNgbPoint;
                                    found = true;
                                    break;
                                }
                                else if (curDist < minDist)
                                {
                                    minDist = curDist;
                                    index = curNgbPoint;
                                }
                            }
                        }

                        nextPoint = index;
                    }

                    if (found)
                    {
                        break;
                    }
                }
            }
        }

        forAll(curMap, globalPointI)
        {
            if (curMap[globalPointI] == -1)
            {
                FatalErrorIn
                (
                    "solidModel::calcGlobalToLocalFaceZonePointMap()"
                )   << "local to global face zone point map is not correct"
                    << abort(FatalError);
            }
        }

        globalToLocalFaceZonePointMap[zoneI] = curMap;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::solidDynamicFvMesh::solidDynamicFvMesh
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
    pMesh_(*this),
    pointDUPtr_(NULL),
    volToPointInterp_(*this),
    globalFaceZonesPtr_(NULL),
    globalToLocalFaceZonePointMapPtr_(NULL),
    meshSmootherPtr_(),
    meshSmoothFreq_(1)
{
    globalToLocalFaceZonePointMap();

    // Create mesh smoother, if specified
    if (dict_.lookupOrDefault<Switch>("meshSmoothing", false))
    {
        Info<< "Creating mesh smoother" << endl;

        dictionary& smoothDict = dict_.subDict("meshSmoothingCoeffs");

        meshSmootherPtr_ = meshSmoother::New(*this, smoothDict);
        meshSmoothFreq_ =
            smoothDict.lookupOrDefault<int>("smoothingFrequency", 1);

        Info<< "    smoothing frequency: " << meshSmoothFreq_ << endl;
    }
    else
    {
        Info<< "No mesh smoother specified" << endl;
    }

    // If materialGgi is used, pointDU field is required
    // befor first update() is called
    pointDU();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidDynamicFvMesh::~solidDynamicFvMesh()
{
    deleteDemandDrivenData(pointDUPtr_);
    deleteDemandDrivenData(globalFaceZonesPtr_);
    deleteDemandDrivenData(globalToLocalFaceZonePointMapPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidDynamicFvMesh::update()
{
    // Lookup displacement increment field from the solver
    const volVectorField& DU = lookupObject<volVectorField>("DU");

    // Interpolate volfield to points
    volToPointInterp_.interpolate(DU, pointDU());

    pointDU().correctBoundaryConditions();
    vectorField& pointDUI = pointDU().internalField();

#   include "correctGgiPointDU.H"

    // Correct symmetryPlane points

    forAll(boundaryMesh(), patchI)
    {
        if (isA<symmetryPolyPatch>(boundaryMesh()[patchI]))
        {
            const labelList& meshPoints = boundaryMesh()[patchI].meshPoints();

            vector avgN =
                gAverage(boundaryMesh()[patchI].pointNormals());

            vector i(1, 0, 0);
            vector j(0, 1, 0);
            vector k(0, 0, 1);

            if (mag(avgN & i) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDUI[meshPoints[pI]].x() = 0;
                }
            }
            else if (mag(avgN & j) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDUI[meshPoints[pI]].y() = 0;
                }
            }
            else if (mag(avgN & k) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDUI[meshPoints[pI]].z() = 0;
                }
            }
        }
        else if (isA<emptyPolyPatch>(boundaryMesh()[patchI]))
        {
            const labelList& meshPoints =
                boundaryMesh()[patchI].meshPoints();

            vector avgN = gAverage(boundaryMesh()[patchI].pointNormals());

            vector k(0, 0, 1);

            if (mag(avgN & k) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointDUI[meshPoints[pI]].z() = 0;
                }
            }
        }
    }

    pointField newPoints = allPoints();

    forAll (pointDUI, pointI)
    {
        newPoints[pointI] += pointDUI[pointI];
    }

#   include "updateGlobalFaceZoneNewPoints.H"

    twoDPointCorrector twoDCorrector(*this);
    twoDCorrector.correctPoints(newPoints);
    twoDCorrector.correctPoints(pointDU().internalField());

    movePoints(newPoints);

    V00();
    moving(false);
    changing(false);
    setPhi().writeOpt() = IOobject::NO_WRITE;

    // Optional mesh smoothing at specified frequency
    // Note: we could perform smoothing based on some mesh quality criteria, but
    // using a fixed frequency should be fine
    if (meshSmootherPtr_.valid())
    {
        if
        (
            (time().timeIndex() != 0)
         && (time().timeIndex() % meshSmoothFreq_ == 0)
        )
        {
            // Clear out previous addressing
            meshSmootherPtr_().clearOut();

            // Store the old points before smoothing
            meshSmootherPtr_().storeOldPoints(points());

            // Smooth the mesh
            meshSmootherPtr_().smooth();

            // Map the fields, if they have not been advected
            if (!meshSmootherPtr_().fieldsAdvected())
            {
                meshSmootherPtr_().mapFields();
            }

            // Update the field values on the wire downstream patch
            // We will extrapolate the values from the internal field using zero
            // gradient extrapolation
            // updateVolFieldsOnPatch<scalar>(wireOutletPatchID_);
            // updateVolFieldsOnPatch<vector>(wireOutletPatchID_);
            // updateVolFieldsOnPatch<tensor>(wireOutletPatchID_);
            // updateVolFieldsOnPatch<symmTensor>(wireOutletPatchID_);
            // updateVolFieldsOnPatch<diagTensor>(wireOutletPatchID_);
            // updateVolFieldsOnPatch<sphericalTensor>(wireOutletPatchID_);

            // Note: no need to update surface and point fields like curS
            // and gradDUf because they get recalculated every iteration anyway

            // Comment on the conservativeMeshToMesh class used by the
            // meshSmoother: by default it does not interpolate the boundary
            // patches, it just does a direct map: this is not OK! So I have
            // used patchToPatch to fix this, which perform a first order
            // inverse distance interpolation (2nd order would be better but
            // this is OK for now).
            // Also, by default conservativeMeshToMesh only interpolated scalars
            // and vectors but I have modified this to now interpolate tensors

            // Flag that the mesh has changed
            //globalMeshChanged = true;
        }
    }

    return true;
}


const Foam::labelList& Foam::solidDynamicFvMesh::globalFaceZones() const
{
    if (!globalFaceZonesPtr_)
    {
        calcGlobalFaceZones();
    }

    return *globalFaceZonesPtr_;
}


const Foam::labelListList&
Foam::solidDynamicFvMesh::globalToLocalFaceZonePointMap() const
{
    if (!globalToLocalFaceZonePointMapPtr_)
    {
        calcGlobalToLocalFaceZonePointMap();
    }

    return *globalToLocalFaceZonePointMapPtr_;
}


// ************************************************************************* //
