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

#include "effectiveTractionLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "crackerFvMesh.H"
#include "mechanicalModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(effectiveTractionLaw, 0);
    addToRunTimeSelectionTable
    (
        faceBreakerLaw, effectiveTractionLaw, dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::effectiveTractionLaw::calcCohesivePatchID() const
{
    if (cohesivePatchIDPtr_)
    {
        FatalErrorIn
        (
            "void Foam::effectiveTractionLaw::calcCohesivePatchID() const"
        ) << "pointer already set" << abort(FatalError);
    }

    const fvMesh& mesh = this->mesh();

    cohesivePatchIDPtr_ = new label(-1);
    label& cohesivePatchID = *cohesivePatchIDPtr_;

    forAll (mesh.boundaryMesh(), patchI)
    {
        if (mesh.boundaryMesh()[patchI].type() == "cohesive")
        {
            cohesivePatchID = patchI;
            break;
        }
    }

    if (cohesivePatchID == -1)
    {
        FatalErrorIn
        (
            "void Foam::effectiveTractionLaw::calcCohesivePatchID() const"
        )   << "boundary patch of type cohesive not found" << abort(FatalError);
    }
}


Foam::label Foam::effectiveTractionLaw::cohesivePatchID() const
{
    if (!cohesivePatchIDPtr_)
    {
        calcCohesivePatchID();
    }

    return *cohesivePatchIDPtr_;
}


void Foam::effectiveTractionLaw::calcAllFacesToBreak() const
{
    if (facesToBreakPtr_ || coupledFacesToBreakPtr_ || facesToBreakFlipPtr_)
    {
        FatalErrorIn
        (
            "void Foam::effectiveTractionLaw::calcAllFacesToBreak() const"
        ) << "pointer already set" << abort(FatalError);
    }

    // First, we check if any internal faces need to be broken, then we will
    // check if any coupled (processor boundary) faces need to be broken.
    // If we find more than one face with a traction fraction greater than 1.0,
    // we will select the face with the highest value

    int nFacesToBreak = 0;
    int nCoupledFacesToBreak = 0;
    //bool topoChange = false;

    // Cast the mesh to a crackerFvMesh

    if (!isA<crackerFvMesh>(this->mesh()))
    {
        FatalErrorIn("Foam::label Foam::effectiveTractionLaw::updateMesh()")
            << "Mesh should be of type: " << crackerFvMesh::typeName
            << abort(FatalError);
    }

    const crackerFvMesh& mesh = refCast<const crackerFvMesh>(this->mesh());


    // Clear out demand driven data
    // clearOut();


    // Lookup cohesive tensile and shear strengths from constitutive model
    // Note: we re-create these fields every time because it is not straight-
    // forward to map fields after topological changes in multi-material cases

    const mechanicalModel& mechanical =
        mesh.lookupObject<mechanicalModel>("mechanicalProperties");

    const surfaceScalarField sigmaMax = mechanical.cohLaw().sigmaMax();
    const surfaceScalarField tauMax = mechanical.cohLaw().tauMax();

    // Lookup the traction field from the solver: this should be updated before
    // calling updateMesh
    const surfaceVectorField& traction =
        mesh.lookupObject<surfaceVectorField>("traction");

    // Calculate normal and shear tractions, where only tensile normal traction
    // is considered i.e. a compressive traction does cannot solely cause
    // failure

    // Face unit normals
    const surfaceVectorField n = mesh.Sf()/mesh.magSf();

    const surfaceScalarField normalTrac =
        max(dimensionedScalar("zero", dimPressure, 0.0), n & traction);

    const surfaceScalarField shearTrac = mag((I - sqr(n)) & traction);

    // Calculate effective traction fraction
    surfaceScalarField effTracFrac =
        Foam::sqrt
        (
            pow(normalTrac/sigmaMax, 2.0) + pow(shearTrac/tauMax, 2.0)
        );

    // Apply crack path limiting if specified
    if (pathLimiterPtr_.valid())
    {
        pathLimiterPtr_().clearOut();
        effTracFrac *= pathLimiterPtr_().facesAllowedToBreak();
    }

    const scalar maxEffTracFrac = gMax(effTracFrac.internalField());

    Info<< nl << "Max traction fraction: " << maxEffTracFrac << endl;


    label faceToBreakIndex = -1;
    scalar faceToBreakEffTracFrac = 0.0;

    if (maxEffTracFrac > 1.0)
    {
        // Find face with maximum traction fraction

        forAll(effTracFrac, faceI)
        {
            if (mag(effTracFrac[faceI] - maxEffTracFrac) < SMALL)
            {
                faceToBreakIndex = faceI;
                faceToBreakEffTracFrac = maxEffTracFrac;
                break;
            }
        }

        nFacesToBreak = 1;

        if (Pstream::parRun())
        {
            // Find processor with greatest traction fraction

            bool procHasFaceToBreak = false;

            if (mag(faceToBreakEffTracFrac - maxEffTracFrac) < SMALL)
            {
                procHasFaceToBreak = true;
            }

            // Check if maximum is present on more then one processors

            label procID = Pstream::nProcs();

            if (procHasFaceToBreak)
            {
                procID = Pstream::myProcNo();
            }

            label minProcID = returnReduce<label>(procID, minOp<label>());

            if (procID != minProcID)
            {
                nFacesToBreak = 0;
            }
        }
    }

    // Check coupled (processor) patches

    label coupledFaceToBreakIndex = -1;
    scalar coupledFaceToBreakEffTracFrac = 0.0;
    scalar maxCoupledEffTracFrac = 0.0;

    if (allowCoupledFaces_ && Pstream::parRun())
    {
        forAll(mesh.boundary(), patchI)
        {
            // Find coupled face with the largest traction fraction
            if (mesh.boundary()[patchI].coupled())
            {
                const scalarField& pEffTracFrac =
                    effTracFrac.boundaryField()[patchI];
                const label start = mesh.boundaryMesh()[patchI].start();

                forAll(pEffTracFrac, faceI)
                {
                    if (pEffTracFrac[faceI] > coupledFaceToBreakEffTracFrac)
                    {
                        coupledFaceToBreakIndex = faceI + start;
                        coupledFaceToBreakEffTracFrac = pEffTracFrac[faceI];
                    }
                }
            }
        }

        // Get global coupled traction fraction
        maxCoupledEffTracFrac =
            returnReduce(coupledFaceToBreakEffTracFrac, maxOp<scalar>());

        Info<< "Max coupled traction fraction: " << maxCoupledEffTracFrac
            << endl;

        if (maxCoupledEffTracFrac > 1.0)
        {
            nCoupledFacesToBreak = 1;

            // Find processor with the greatest traction fraction on a processor
            // face

            bool procHasCoupledFaceToBreak = false;

            if
            (
                mag(maxCoupledEffTracFrac - coupledFaceToBreakEffTracFrac)
              < SMALL
            )
            {
                procHasCoupledFaceToBreak = true;
            }

            // Check if maximum is present on more then one processors

            label procID = Pstream::nProcs();

            if (procHasCoupledFaceToBreak)
            {
                procID = Pstream::myProcNo();
            }

            label minProcID = returnReduce<label>(procID, minOp<label>());

            if (procID != minProcID)
            {
                nCoupledFacesToBreak = 0;
            }
        }
    }

    // If both an internal and coupled face want to break then we pick the face
    // with the greatest traction fraction

    if (maxCoupledEffTracFrac > maxEffTracFrac)
    {
        // Break coupled face
        nFacesToBreak = 0;
    }
    else
    {
        // Break internal face
        nCoupledFacesToBreak = 0;
    }


    // Make sure that coupled faces are broken in pairs

    labelList ngbProc(Pstream::nProcs(), -1);
    labelList index(Pstream::nProcs(), -1);

    if (nCoupledFacesToBreak)
    {
        label patchID =
            mesh.boundaryMesh().whichPatch(coupledFaceToBreakIndex);

        if (patchID == -1)
        {
            FatalErrorIn("calcAllFacesToBreak()")
                << "something is wrong: patchID is -1 for coupled face"
                << abort(FatalError);
        }

        label start = mesh.boundaryMesh()[patchID].start();
        label localIndex = coupledFaceToBreakIndex - start;

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(mesh.boundaryMesh()[patchID]);

        label ngbProcNo = procPatch.neighbProcNo();

        ngbProc[Pstream::myProcNo()] = ngbProcNo;
        index[Pstream::myProcNo()] = localIndex;
    }

    if (returnReduce(nCoupledFacesToBreak, maxOp<label>()))
    {
        reduce(ngbProc, maxOp<labelList>());
        reduce(index, maxOp<labelList>());

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                if (ngbProc[procI] == Pstream::myProcNo())
                {
                    forAll(mesh.boundaryMesh(), patchI)
                    {
                        if
                        (
                            mesh.boundaryMesh()[patchI].type()
                            == processorPolyPatch::typeName
                        )
                        {
                            const processorPolyPatch& procPatch =
                                refCast<const processorPolyPatch>
                                (
                                    mesh.boundaryMesh()[patchI]
                                );

                            label ngbProcNo = procPatch.neighbProcNo();

                            if (ngbProcNo == procI)
                            {
                                label start =
                                    mesh.boundaryMesh()[patchI].start();
                                coupledFaceToBreakIndex = start + index[procI];
                                nCoupledFacesToBreak = 1;
                            }
                        }
                    }
                }
            }
        }
    }


    // Note: It is not necessary to scale the tractions as the cohesive boundary
    // condition will do this

    if (returnReduce(nCoupledFacesToBreak, sumOp<label>()) > 2)
    {
        FatalErrorIn("effectiveTractionLaw::calcAllFacesToBreak()")
            << "More than two processors are trying to break a coupled face"
            << abort(FatalError);
    }


    // Now we break the face

    if (nFacesToBreak)
    {
        facesToBreakPtr_ = new labelList(nFacesToBreak, faceToBreakIndex);

        facesToBreakFlipPtr_ = new boolList(nFacesToBreak, false);

        coupledFacesToBreakPtr_ = new labelList(0);
    }
    else if (nCoupledFacesToBreak)
    {
        facesToBreakPtr_ = new labelList(0);

        facesToBreakFlipPtr_ = new boolList(0);

        coupledFacesToBreakPtr_ =
            new labelList(nCoupledFacesToBreak, coupledFaceToBreakIndex);
    }
    else
    {
        facesToBreakPtr_ = new labelList(0);

        facesToBreakFlipPtr_ = new boolList(0);

        coupledFacesToBreakPtr_ = new labelList(0);
    }

    if (debug)
    {
        const labelList& f = *facesToBreakPtr_;
        forAll(f, faceI)
        {
            Pout<< "coupled face to break: " << mesh.Cf()[f[faceI]] << endl;
        }

        const labelList& cf = *coupledFacesToBreakPtr_;

        forAll(cf, faceI)
        {
            const label patchID = mesh.boundaryMesh().whichPatch(cf[faceI]);
            const label start = mesh.boundaryMesh()[patchID].start();

            Pout<< "coupled face to break: "
                << mesh.boundaryMesh()[patchID].faceCentres()[cf[faceI] - start]
                << endl;
        }
    }
}


void Foam::effectiveTractionLaw::clearOut()
{
    deleteDemandDrivenData(cohesivePatchIDPtr_);
    deleteDemandDrivenData(facesToBreakPtr_);
    deleteDemandDrivenData(facesToBreakFlipPtr_);
    deleteDemandDrivenData(coupledFacesToBreakPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from dictionary
Foam::effectiveTractionLaw::effectiveTractionLaw
(
    const word& name,
    fvMesh& mesh,
    const dictionary& dict
)
:
    faceBreakerLaw(name, mesh, dict),
    cohesivePatchIDPtr_(NULL),
    allowCoupledFaces_(dict.lookupOrDefault<Switch>("allowCoupledFaces", true)),
    facesToBreakPtr_(NULL),
    facesToBreakFlipPtr_(NULL),
    coupledFacesToBreakPtr_(NULL),
    pathLimiterPtr_(NULL)
{
    if (!allowCoupledFaces_)
    {
        Warning
            << name << ": allowCoupledFaces is false" << endl;
    }

    // Create crackPathLimiter if specified
    if (dict.found("crackPathLimiter"))
    {
        pathLimiterPtr_ =
            crackPathLimiter::New
            (
                "law", mesh, dict.subDict("crackPathLimiter")
            );
    }
    else
    {
        Info<< "crackPathLimiter not specified" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::effectiveTractionLaw::~effectiveTractionLaw()
{
    clearOut();
}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //


const Foam::labelList& Foam::effectiveTractionLaw::facesToBreak() const
{
    if (!facesToBreakPtr_)
    {
        calcAllFacesToBreak();
    }

    return *facesToBreakPtr_;
}


const Foam::boolList& Foam::effectiveTractionLaw::facesToBreakFlip() const
{
    if (!facesToBreakFlipPtr_)
    {
        calcAllFacesToBreak();
    }

    return *facesToBreakFlipPtr_;
}


const Foam::labelList& Foam::effectiveTractionLaw::coupledFacesToBreak() const
{
    if (!facesToBreakPtr_)
    {
        calcAllFacesToBreak();
    }

    return *coupledFacesToBreakPtr_;
}


// ************************************************************************* //
