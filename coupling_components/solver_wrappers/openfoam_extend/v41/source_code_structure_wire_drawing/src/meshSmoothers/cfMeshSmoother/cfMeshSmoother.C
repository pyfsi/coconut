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

#include "cfMeshSmoother.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "twoDPointCorrector.H"
#include "newLeastSquaresVolPointInterpolation.H"
#include "ZoneID.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "fvcDiv.H"
#include "fixedGradientFvPatchFields.H"
#include "cellQuality.H"
#include "volPointInterpolation.H"
#include "faMesh.H"
#include "areaMesh.H"
#include "fac.H"
#include "pointNormalsConsistent.H"

// cfMesh headers
#include "polyMeshGenModifier.H"
#include "meshOptimizer.H"
#include "meshSurfaceEngine.H"
#include "meshSurfaceOptimizer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cfMeshSmoother, 0);
    addToRunTimeSelectionTable
    (
        meshSmoother, cfMeshSmoother, dictionary
    );


// * * * * * * * * * * * * Private Data Members  * * * * * * * * * * * * * * //

fvMesh& cfMeshSmoother::subMeshToSmooth()
{
    return subsetMeshToSmooth().subMesh();
}


const fvMesh& cfMeshSmoother::subMeshToSmooth() const
{
    return subsetMeshToSmooth().subMesh();
}


newFvMeshSubset& cfMeshSmoother::subsetMeshToSmooth()
{
    if (subMeshToSmooth_.empty())
    {
        makeSubMeshToSmooth();
    }

    return subMeshToSmooth_();
}


const newFvMeshSubset&
cfMeshSmoother::subsetMeshToSmooth() const
{
    if (subMeshToSmooth_.empty())
    {
        makeSubMeshToSmooth();
    }

    return subMeshToSmooth_();
}


void cfMeshSmoother::makeSubMeshToSmooth() const
{
    if (debug)
    {
        InfoIn(type() + "::makeSubMeshToSmooth() const")
            << endl;
    }

    if (!subMeshToSmooth_.empty())
    {
        FatalErrorIn(type() + "::makeSubMeshToSmooth() const")
            << "pointer already set!" << abort(FatalError);
    }

    // Take a reference to the mesh for efficiency
    const fvMesh& mesh = this->mesh();

    // Read the list of cell zones to smooth from the dict; if none are
    // specified then we will smooth the entire mesh

    // List indicating cells to be smoothed (1) and cells to not smooth (0)
    labelList cellsToSmooth(mesh.nCells(), 0);

    if (dict().found("cellZones"))
    {
        // Read the cell zones to smooth
        const wordList cellZonesNames = wordList(dict().lookup("cellZones"));

        forAll(cellZonesNames, czI)
        {
            // Check the cell zone exist
            const ZoneID<cellZone> cellZoneID =
                ZoneID<cellZone>(cellZonesNames[czI], mesh.cellZones());

            if (!cellZoneID.active())
            {
                FatalErrorIn(type() + "::smooth()")
                    << "cellZone " << cellZonesNames[czI] << " not found!"
                    << abort(FatalError);
            }

            // Set cells in zone to be smoothed
            Info<< "Cell zone " << cellZonesNames[czI] << " will be smoothed"
                << endl;
            const labelList& curCellZone = mesh.cellZones()[cellZoneID.index()];
            forAll(curCellZone, cI)
            {
                cellsToSmooth[curCellZone[cI]] = 1;
            }
        }
    }
    else
    {
        // Smooth all the cells
        Info<< "All cells in the mesh will be smoothed"
            << endl;
        cellsToSmooth = 1.0;
    }


    // Create a subMesh containing only the cells to be smoothed
    subMeshToSmooth_.set
    (
        new newFvMeshSubset
        (
            IOobject
            (
                "cellsToSmooth",
                mesh.time().constant(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );

    if (debug)
    {
        InfoIn(type() + "::makeSubMeshToSmooth() const")
            << "Creating the fvMeshSubset" << endl;
    }

    // Select cells with cellsToSmooth set to 1
    subMeshToSmooth_().setLargeCellSubset(cellsToSmooth, 1);
}


void cfMeshSmoother::calculatePointMotionCourantNo
(
    const pointVectorField& pointMotionD,
    const fvMesh& mesh
)
{
    // To calculate a characteristic length field, we will interpolate the
    // cell volume field to the points and take the cubed-root of it
    // This may be a bit liberal e.g. for high aspect ratio cells

    // Create a cell volume field
    wordList patchTypes(mesh.boundaryMesh().size(), "zeroGradient");
    volScalarField cellVolumes
    (
        IOobject
        (
            "cellVolumes",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimVolume, 0.0),
        patchTypes
    );
    cellVolumes.internalField() = mesh.V();
    cellVolumes.correctBoundaryConditions();

    // Interpolate cell volumes to the points
    volPointInterpolation volToPoint(mesh);
    const pointScalarField pointVolumes = volToPoint.interpolate(cellVolumes);

    // Calculate the local characteristic length
    scalarField pointDelta(pointVolumes.size(), 0.0);
    const scalarField& pointVolumesI = pointVolumes.internalField();
    forAll(pointDelta, pI)
    {
        pointDelta[pI] = Foam::cbrt(pointVolumesI[pI]);
    }

    Info<< "Min pointVolume = " << gMin(pointVolumes) << nl
        << "Max pointDelta = " << gMax(pointDelta) << nl
        << "Min pointDelta = " << gMin(pointDelta) << endl;

    // Calculate maximum Courant number
    const scalar maxCo =
        gMax(mag(pointMotionD.internalField()/pointDelta));

    Info<< "Maximum point motion Courant number = " << maxCo << endl;
}

void cfMeshSmoother::smoothSurface(polyMeshGen& pmg, fvMesh& mesh) const
{
    // Store mesh points before smoothing
    const pointField oldPoints = mesh.points();

    // Store old patch point normals
    // Take care that the point normals are calculated consistently in parallel
    PtrList<pointField> oldPatchPointNormals(mesh.boundaryMesh().size());
    forAll(mesh.boundaryMesh(), patchI)
    {
        oldPatchPointNormals.set
        (
            patchI,
            new pointField(mesh.boundaryMesh()[patchI].nPoints(), vector::zero)
        );
    }
    pointNormalsConsistent(oldPatchPointNormals, mesh);

    // Outer mesh smoothing loop
    int oCorr = 0;
    const int nOCorr = iMQCoeffs_.nSurfaceIterations_;
    do
    {
        Info<< nl << "Outer corr = " << oCorr << endl;

        // First, smooth surface only
        meshSurfaceEngine mse(pmg);
        meshSurfaceOptimizer msOpt(mse);
        //msOpt.optimizeSurface(nSurfaceIterations);
        msOpt.optimizeSurface(1);

        // Calculate point displacement
        const pointMesh& pMesh = pointMesh::New(mesh);
        pointVectorField pointDispField
        (
            IOobject
            (
                "pointDisp",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pMesh,
            dimensionedVector("zero", dimLength, vector::zero)
        );
        pointField& pointDisp = pointDispField.internalField();
        pointDisp = pmg.points() - oldPoints;

        // Set disp to zero on "fixedPatches" and enforce flat patches
        correctFixedPatchesMotion(pointDispField);
        correctFlatPatchesMotion(pointDispField);

        // At this point, pmg.points() are smoothed but geometry is not
        // preserved
        // We will perform a correction to the motion of the surface in an
        // attempt preserve the geometry

        int iCorr = 0;
        scalar res = GREAT;
        scalar prevMaxMagDisp = GREAT;
        const scalar tol = 1e-10;
        const int nCorr = 10;

        do
        {
            // Move the fvMesh to be the same as the pmg
            mesh.movePoints(pmg.points());

            // Optionally: we will not project the points on the final iteration
            // In this way, the geometry will be slightly smoothed, in an
            // attempt to remove any sucface oscillations
            if
            (
                oCorr != (nOCorr - 1)
             || !iMQCoeffs_.noSurfaceProjectionOnFinalIteration_
            )
            {
                // Calculate new point normals
                PtrList<pointField> newPatchPointNormals
                (
                    mesh.boundaryMesh().size()
                );
                forAll(mesh.boundaryMesh(), patchI)
                {
                    newPatchPointNormals.set
                    (
                        patchI,
                        new pointField
                        (
                            mesh.boundaryMesh()[patchI].nPoints(), vector::zero
                        )
                    );
                }
                pointNormalsConsistent(newPatchPointNormals, mesh);

                // For all patches (except for coupled patches), remove
                // surface-normal component of displacement
                forAll(mesh.boundaryMesh(), patchI)
                {
                    if (mesh.boundaryMesh()[patchI].coupled())
                    {
                        continue;
                    }

                    const polyPatch& ppatch = mesh.boundaryMesh()[patchI];
                    const labelList& meshPoints = ppatch.meshPoints();

                    // Define average point normals
                    const vectorField nAverage =
                        0.5
                       *(
                            newPatchPointNormals[patchI]
                          + oldPatchPointNormals[patchI]
                        );

                    // Calculate new point displacement for points on this patch
                    forAll(nAverage, pI)
                    {
                        const label pointID = meshPoints[pI];
                        pointDisp[pointID] =
                            ((I - sqr(nAverage[pI])) & pointDisp[pointID]);
                    }
                }
            }
            else
            {
                Info<< "    final iteration: not projecting surface" << endl;
            }

            // Update pmg points
            pointFieldPMG& points =
                const_cast<pointFieldPMG&>(pmg.points());
            pointDispField.correctBoundaryConditions();
            forAll(points, pointI)
            {
                points[pointI] = oldPoints[pointI] + pointDisp[pointI];
            }

            // Define residual
            if (iCorr == 0)
            {
                res = 1.0;
            }
            else
            {
                res = mag
                (
                    gMax(mag(pointDisp)) - prevMaxMagDisp
                )/(prevMaxMagDisp + SMALL);
            }
            prevMaxMagDisp = gMax(mag(pointDisp));

            Info<< "i = " << iCorr << ", res = " << res << endl;
        }
        while (res > tol && ++iCorr < nCorr);
    }
    while (++oCorr < nOCorr);
}


// void cfMeshSmoother::calculateSweptAreas
// (
//     PtrList<edgeScalarField>& sweptAreas,
//     const PtrList<faMesh>& aMeshes,
//     const fvMesh& mesh,
//     const pointField& oldPoints
// ) const
// {
//     forAll(sweptAreas, patchI)
//     {
//         if (mesh.boundaryMesh()[patchI].type() == "empty")
//         {
//             continue;
//         }

//         // Take a reference to the polyPatch
//         const polyPatch& ppatch = mesh.boundaryMesh()[patchI];

//         // List of local edges
//         const edgeList& edges = ppatch.edges();

//         // List of local deformed points
//         const pointField& newPatchPoints = ppatch.localPoints();

//         // Mesh points
//         const labelList& meshPoints = ppatch.meshPoints();

//         // Edge bi-normals
//         const edgeVectorField& Le = aMeshes[patchI].Le();
//         const vectorField& LeI = Le.internalField();

//         // Initialise swept areas field
//         sweptAreas.set
//         (
//             patchI,
//             new edgeScalarField
//             (
//                 IOobject
//                 (
//                     "sweptAreas_patch" + mesh.boundaryMesh()[patchI].name(),
//                     mesh.time().timeName(),
//                     mesh,
//                     IOobject::NO_READ,
//                     IOobject::NO_WRITE
//                 ),
//                 aMeshes[patchI],
//                 dimensionedScalar("zero", dimArea, 0.0)
//             )
//         );
//         scalarField& sweptArea = sweptAreas[patchI].internalField();

//         forAll(edges, edgeI)
//         {
//             if (ppatch.isInternalEdge(edgeI))
//             {
//                 const edge& curEdge = edges[edgeI];

//                 // Create a face using the two old points and two new points
//                 face curFace(4);
//                 curFace[0] = 0;
//                 curFace[1] = 1;
//                 curFace[2] = 2;
//                 curFace[3] = 3;

//                 pointField pointsInFace(4);
//                 pointsInFace[0] = oldPoints[meshPoints[curEdge.start()]];
//                 pointsInFace[1] = oldPoints[meshPoints[curEdge.end()]];
//                 pointsInFace[2] = newPatchPoints[curEdge.end()];
//                 pointsInFace[3] = newPatchPoints[curEdge.start()];

//                 // Vector in the direction of the edge motion
//                 const vector X =
//                 (
//                     pointsInFace[3] + pointsInFace[2]
//                   - pointsInFace[1] - pointsInFace[0]
//                 );

//                 // Swept area can be positive or negative
//                 sweptArea[edgeI] =
//                     sign(LeI[edgeI] & X)*curFace.mag(pointsInFace);
//             }
//         }
//     }
// }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
cfMeshSmoother::cfMeshSmoother
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    meshSmoother(mesh, dict),
    subMeshToSmooth_(),
    iMQCoeffs_
    (
        meshSmoother::dict().lookupOrDefault<int>("nLoops", 100),
        meshSmoother::dict().lookupOrDefault<int>("nIterations", 1000),
        meshSmoother::dict().lookupOrDefault<int>("nSurfaceIterations", 0),
        meshSmoother::dict().lookupOrDefault<scalar>("qualityThreshold", 0.1),
        meshSmoother::dict().lookupOrDefault<word>
        (
            "constrainedCellsSet", word()
        ),
        meshSmoother::dict().lookupOrDefault<int>
        (
            "nIterationsNearBoundaries", 100
        ),
        meshSmoother::dict().lookupOrDefault<int>("nLayers", 20),
        meshSmoother::dict().lookupOrDefault<Switch>("smoothSurface", false),
        meshSmoother::dict().lookupOrDefault<Switch>("optimizeMeshFV", false),
        meshSmoother::dict().lookupOrDefault<Switch>
        (
            "noSurfaceProjectionOnFinalIteration", false
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfMeshSmoother::~cfMeshSmoother()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void cfMeshSmoother::clearOut()
{
    subMeshToSmooth_.clear();
}


scalar cfMeshSmoother::smooth()
{
    if (debug)
    {
        InfoIn("scalar cfMeshSmoother::smooth()")
            << endl;
    }

    Info<< "Smoothing the mesh" << endl;

    // Get the sub-mesh containing the cells to be smoothed
    // This could be the whole mesh or specified cells zones
    fvMesh& subMesh = subMeshToSmooth();

    // Create the point mesh
    const pointMesh& pMesh = pointMesh::New(subMesh);

    // Create the point motion field
    pointVectorField pointMotionD
    (
        IOobject
        (
            "pointMotionD",
            subMesh.time().timeName(),
            subMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    // Perform smoothing
    blockLduMatrix::debug = 0;

    // Create lists needed to make a polyMeshGen object
    wordList patchNames(subMesh.boundaryMesh().size(), "");
    wordList patchTypes(subMesh.boundaryMesh().size(), "");
    labelList patchStarts(subMesh.boundaryMesh().size(), -1);
    labelList nFacesInPatch(subMesh.boundaryMesh().size(), -1);
    labelList procPatchMyProcNo(subMesh.boundaryMesh().size(), -1);
    labelList procPatchNeighbProcNo(subMesh.boundaryMesh().size(), -1);
    forAll(patchStarts, patchI)
    {
        patchNames[patchI] = subMesh.boundaryMesh()[patchI].name();
        patchTypes[patchI] = subMesh.boundaryMesh()[patchI].type();
        patchStarts[patchI] = subMesh.boundaryMesh()[patchI].start();
        nFacesInPatch[patchI] = subMesh.boundaryMesh()[patchI].size();

        if (patchTypes[patchI] == processorPolyPatch::typeName)
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>
                (
                    subMesh.boundaryMesh()[patchI]
                );

            procPatchMyProcNo[patchI] = procPatch.myProcNo();
            procPatchNeighbProcNo[patchI] = procPatch.neighbProcNo();
        }
    }

    // Create a cfMesh polyMeshGen object from the current polyMesh
    Info<< endl;
    if (debug)
    {
        Info<< "Creating polyMeshGen" << endl;
    }
    polyMeshGen pmg
    (
        subMesh.time(),
        subMesh.points(),
        subMesh.faces(),
        subMesh.cells(),
        patchNames,
        patchTypes,
        patchStarts,
        nFacesInPatch,
        procPatchMyProcNo,
        procPatchNeighbProcNo,
        subMesh.instance(),
        subMesh.meshDir()
    );

    // Construct the smoother
    if (debug)
    {
        Info<< "Creating meshOptimizer" << endl;
    }
    meshOptimizer mOpt(pmg);
    if (debug)
    {
        Info<< "Creating meshOptimizer: done" << endl;
    }

    // PC: we could use this constrained method instead of using a
    // subMesh...
    // For now, we will not but will leave the option
    // if (!iMQCoeffs_.constrainedCellsSet_.empty())
    // {
    //     //- lock cells in constrainedCellsSet
    //     mOpt.lockCellsInSubset
    //     (
    //         iMQCoeffs_.constrainedCellsSet_
    //     );

    //     //- find boundary faces which shall be locked
    //     labelLongList lockedBndFaces, selectedCells;

    //     const label sId =
    //         pmg.cellSubsetIndex(iMQCoeffs_.constrainedCellsSet_);
    //     pmg.cellsInSubset(sId, selectedCells);

    //     boolList activeCell(pmg.cells().size(), false);
    //     forAll(selectedCells, i)
    //     {
    //         activeCell[selectedCells[i]] = true;
    //     }
    // }

    // Clear geometry information before volume smoothing
    pmg.clearAddressingData();

    // Perform optimisation using the laplace smoother on entire mesh
    // This tends to produce poor cells near the boundary
    if (iMQCoeffs_.optimizeMeshFV_)
    {
        if (debug)
        {
            Info<< "Running optimizeMeshFV" << endl;
        }
        mOpt.optimizeMeshFV
        (
            iMQCoeffs_.nLoops_,
            iMQCoeffs_.nLoops_,
            iMQCoeffs_.nIterations_,
            iMQCoeffs_.nSurfaceIterations_
        );
    }

    // Optionally, smooth the surface
    if (iMQCoeffs_.smoothSurface_)
    {
        if (debug)
        {
            Info<< "Smoothing the surface" << endl;
        }

        // We will make sure to reset the subMesh afterwards as it also gets
        // moved
        const pointField oldPoints = subMesh.points();
        smoothSurface(pmg, subMesh);
        subMesh.movePoints(oldPoints);
    }

    // Optimize mesh near boundaries
    // This smoother does a good job of smoothing cells near the boundary but it
    // can be slow
    Info<< nl << "Smoothing cells near boundary" << endl;
    mOpt.optimizeMeshNearBoundaries
    (
        iMQCoeffs_.nIterationsNearBoundaries_,
        iMQCoeffs_.nLayers_
    );

    // Calculate point motion
    WarningIn(type() + "::smooth()")
        << "To-do: limit mesh motion to the mesh Courant number!" << endl;
    vectorField& pointMotionDI = pointMotionD.internalField();
    pointMotionDI = pmg.points() - subMesh.points();
    Info<< nl << "Maximum point motion is " << gMax(mag(pointMotionDI))
        << endl;
    calculatePointMotionCourantNo(pointMotionD, subMesh);

    // Sub-mesh new points
    pointField newPointsSubMesh = subMesh.points() + pointMotionDI;

    // Apply corrections for 2-D meshes
    twoDPointCorrector twoDCorrectorSubMesh(subMesh);
    twoDCorrectorSubMesh.correctPoints(newPointsSubMesh);

    const scalarField oldVSubMesh = subMesh.V();

    // It is important to set the old points otherwise the swept volumes
    // will be wrong
    subMesh.setOldPoints(subMesh.points());

    // Move the subMesh
    subMesh.movePoints(newPointsSubMesh);

    // Disable flags saying that the mesh is moved as we do not need to
    // store the mesh motion fluxes
    subMesh.moving(false);
    subMesh.changing(false);
    subMesh.setPhi().writeOpt() = IOobject::NO_WRITE;

    // Take a reference to the main mesh: be careful not to confuse the
    // sub-mesh and main mesh
    fvMesh& mesh = this->mesh();

    // Map the point mesh displacements from the sub-mesh to the main mesh
    const labelList& pointMap = subsetMeshToSmooth().pointMap();

    // Apply correction at symmetry planes
    correctFlatPatchesMotion(pointMotionD);
    setNormalDisplacementToZeroOnSymmetryPlanes(pointMotionD);

    // Set newPoints field by mapping the sub-mesh point displacement to the
    // main mesh
    pointField newPoints = mesh.points();
    forAll(pointMap, subMeshPI)
    {
        const label pointID = pointMap[subMeshPI];
        newPoints[pointID] += pointMotionDI[subMeshPI];
    }

    // Apply corrections for 2-D meshes
    twoDPointCorrector twoDCorrector(mesh);
    twoDCorrector.correctPoints(newPoints);

    // Store the old cell volumes
    const scalarField volOld = scalarField(mesh.V());

    // It is important to set the old points otherwise the swept volumes
    // will be wrong
    mesh.setOldPoints(mesh.points());

    // Take a copy of the old points
    const pointField oldPoints = mesh.points();

    // Patches field may be advected using the finiteArea method, so we need to
    // neccessary fields
    // PtrList<faMesh> aMeshes(mesh.boundaryMesh().size());
    // PtrList<areaScalarField> areasOld(mesh.boundaryMesh().size());
    // forAll(areasOld, patchI)
    // {
    //     if (mesh.boundaryMesh()[patchI].type() == "empty")
    //     {
    //         continue;
    //     }

    //     aMeshes.set(patchI, new faMesh(mesh, patchI));

    //     areasOld.set
    //     (
    //         patchI,
    //         new areaScalarField
    //         (
    //             IOobject
    //             (
    //                 "areaOldPatch" + mesh.boundaryMesh()[patchI].name(),
    //                 mesh.time().timeName(),
    //                 mesh,
    //                 IOobject::NO_READ,
    //                 IOobject::NO_WRITE
    //             ),
    //             aMeshes[patchI],
    //             dimensionedScalar("zero", dimArea, 0.0),
    //             "zeroGradient"
    //         )
    //     );

    //     areasOld[patchI].internalField() = mesh.boundary()[patchI].magSf();
    //     areasOld[patchI].correctBoundaryConditions();
    // }

    // Move the mesh and store the swept face volumes
    const scalarField sweptVol = mesh.movePoints(newPoints);

    // Store new areas for finiteArea mesh
    // PtrList<areaScalarField> areasNew(mesh.boundaryMesh().size());
    // forAll(areasNew, patchI)
    // {
    //     if (mesh.boundaryMesh()[patchI].type() == "empty")
    //     {
    //         continue;
    //     }

    //     areasNew.set
    //     (
    //         patchI,
    //         new areaScalarField
    //         (
    //             IOobject
    //             (
    //                 "areaNewPatch" + mesh.boundaryMesh()[patchI].name(),
    //                 mesh.time().timeName(),
    //                 mesh,
    //                 IOobject::NO_READ,
    //                 IOobject::NO_WRITE
    //             ),
    //             aMeshes[patchI],
    //             dimensionedScalar("zero", dimArea, 0.0),
    //             "zeroGradient"
    //         )
    //     );

    //     areasNew[patchI].internalField() = mesh.boundary()[patchI].magSf();
    //     areasNew[patchI].correctBoundaryConditions();
    // }

    if (Switch(dict().lookup("advectFields")))
    {
        // Store the new cell volumes
        const scalarField volNew = scalarField(mesh.V());

        // Calculate areas swept by the edges on the patch
        //PtrList<edgeScalarField> sweptAreas(mesh.boundaryMesh().size());
        //calculateSweptAreas(sweptAreas, aMeshes, mesh, oldPoints);

        // Advect the volFields
       //advectFields(sweptVol, volOld, volNew, areasOld, areasNew, sweptAreas);
        advectFields(sweptVol, volOld, volNew);
    }

    // Disable flags saying that the mesh is moved as we do not need to
    // store the mesh motion fluxes
    mesh.moving(false);
    mesh.changing(false);
    mesh.setPhi().writeOpt() = IOobject::NO_WRITE;

    // Clear mesh data
    mesh.clearOut();
    subMesh.clearOut();

    blockLduMatrix::debug = 0;

    return 0;
}

// ************************************************************************* //

} // end of namespace foam
