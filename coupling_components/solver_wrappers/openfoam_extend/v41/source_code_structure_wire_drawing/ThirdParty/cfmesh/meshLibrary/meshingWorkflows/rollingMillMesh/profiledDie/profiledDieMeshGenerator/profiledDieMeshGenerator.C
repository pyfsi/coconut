/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "profiledDieMeshGenerator.H"
#include "demandDrivenData.H"
#include "meshSurfaceEngine.H"
#include "rollingMillMesh.H"
#include "profiledDieGeometryInterpolator.H"
#include "meshSurfaceOptimizer.H"
#include "meshOptimizer.H"

# ifdef ExtendSpecific
#include "polyMesh.H"
#include "tetPolyMesh.H"
#include "tetFem.H"
#include "tetFemMatrix.H"
# endif

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

const meshSurfaceEngine& profiledDieMeshGenerator::surfaceEngine() const
{
    if( !msePtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
        {
            FatalErrorIn
            (
                "const meshSurfaceEngine& profiledDieMeshGenerator"
                "::surfaceEngine() const"
            ) << "Constructing meshSurfaceEngine in parallel region"
              << abort(FatalError);
        }
        # endif

        msePtr_ = std::make_shared<meshSurfaceEngine>(mesh_);
    }

    return *msePtr_;
}

const profiledDieGeometryInterpolator&
profiledDieMeshGenerator::interpolator() const
{
    if( !interpolatorPtr_ )
    {
        interpolatorPtr_ =
            std::make_shared<profiledDieGeometryInterpolator>
            (
                dieGeom_,
                patchNamesHandler_,
                mesh_
            );
    }

    return *interpolatorPtr_;
}

void profiledDieMeshGenerator::detectLockedBndPoints
(
    labelLongList& lockedBndPoints
) const
{
    lockedBndPoints.clear();

    const meshSurfaceEngine& mse = surfaceEngine();
    const pointFieldPMG& points = mse.points();
    const faceList::subList& bFaces = mse.boundaryFaces();
    const labelLongList& bp = mse.bp();
    const labelLongList& facePatch = mse.boundaryFacePatches();

    //- detect the range of x coordinates
    scalar xMin(VGREAT), xMax(-VGREAT);

    forAll(points, pointI)
    {
        const point& p = points[pointI];

        xMin = min(xMin, p.x());
        xMax = max(xMax, p.x());
    }

    const word upstreamName =
        patchNamesHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIEUPSTREAM
        );
    const word downstreamName =
        patchNamesHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIEDOWNSTREAM
        );

    const word frontName =
        patchNamesHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIEFRONT
        );
    const word backName =
        patchNamesHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIEBACK
        );

    const word symmYName =
        patchNamesHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIESYMMY
        );
    const word symmZName =
        patchNamesHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::DIESYMMZ
        );

    label upstreamPatchId(-1);
    label downstreamPatchId(-1);
    label frontPatchId(-1);
    label backPatchId(-1);
    label symmYPatchId(-1);
    label symmZPatchId(-1);

    forAll(mesh_.boundaries(), patchI)
    {
        const word pName = mesh_.boundaries()[patchI].patchName();

        if( pName.find(upstreamName) != word::npos )
        {
            upstreamPatchId = patchI;
        }
        else if( pName.find(downstreamName) != word::npos )
        {
            downstreamPatchId = patchI;
        }
        else if( pName.find(symmYName) != word::npos )
        {
            symmYPatchId = patchI;
        }
        else if( pName.find(symmZName) != word::npos )
        {
            symmZPatchId = patchI;
        }
        else if( pName.find(frontName) != word::npos )
        {
            frontPatchId = patchI;
        }
        else if( pName.find(backName) != word::npos )
        {
            backPatchId = patchI;
        }
    }

    //- lock vertices belonging to faces whose
    labelHashSet locked;
    forAll(bFaces, bfI)
    {
        const label patchI = facePatch[bfI];

        //- skip patches that may be moved
        if
        (
            patchI == symmYPatchId ||
            patchI == symmZPatchId ||
            patchI == frontPatchId ||
            patchI == backPatchId
        )
        {
            continue;
        }

        const face& bf = bFaces[bfI];

        if( patchI == upstreamPatchId || patchI == downstreamPatchId )
        {
            //- lock vertices that are not located in the front and back plane
            const point fc = help::faceCentre(points, bf);

            if( fc.x() > (xMin+SMALL) && fc.x() < (xMax-SMALL) )
            {
                //- lock points
                forAll(bf, pI)
                {
                    locked.insert(bp[bf[pI]]);
                }
            }
        }
        else
        {
            //- lock points at patches that shall not be moved
            forAll(bf, pI)
            {
                locked.insert(bp[bf[pI]]);
            }
        }
    }

    lockedBndPoints.clear();
    forAllConstIter(labelHashSet, locked, it)
        lockedBndPoints.append(it.key());

    const labelLongList& bPoints = mse.boundaryPoints();
    const label lockedPointsId = mesh_.addPointSubset("lockedPoints");
    forAll(lockedBndPoints, i)
        mesh_.addPointToSubset(lockedPointsId, bPoints[lockedBndPoints[i]]);
}

void profiledDieMeshGenerator::calculateDisplacements()
{
    const pointFieldPMG& points = mesh_.points();

    //- construct cross section at all x coordinates found in the mesh
    const profiledDieGeometryInterpolator& geomCrossSections = interpolator();

    Info << "Calculating new positions" << endl;
    vectorLongList pointDisplacements;
    geomCrossSections.calculateDisplacementsAll(pointDisplacements);

    polyMeshGenModifier meshModifier(mesh_);
    forAll(pointDisplacements, pI)
        meshModifier.movePoint(pI, points[pI] + pointDisplacements[pI]);

    meshModifier.clearAll();
}

void profiledDieMeshGenerator::updateMeshPointsFVM
(
    const std::map<label, vector>& bndDisplacements
)
{
    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();
    const labelLongList& owner = mesh_.owner();
    const labelLongList& neighbour = mesh_.neighbour();

    const polyMeshGenAddressing& addressingData = mesh_.addressingData();
    const edgeLongList& edges = addressingData.edges();
    const VRWGraph& pointEdges = addressingData.pointEdges();
    const VRWGraph& faceEdges = addressingData.faceEdges();
    const vectorLongList& cellCentres = addressingData.cellCentres();

    //- detect points at fixed diameter
    labelLongList pointsAtFixedDiameter;
    detectLockedBndPoints(pointsAtFixedDiameter);
    labelHashSet fixedDiameterPoints(pointsAtFixedDiameter.size());
    forAll(pointsAtFixedDiameter, i)
        fixedDiameterPoints.insert(pointsAtFixedDiameter[i]);
    pointsAtFixedDiameter.setSize(0);

    scalarLongList edgeCoefficients(edges.size());

    vectorLongList displacements(points.size());
    vectorLongList newDisplacements(points.size());

    polyMeshGenModifier meshModifier(mesh_);

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static) nowait
        # endif
        forAll(edgeCoefficients, i)
            edgeCoefficients[i] = 0.0;

        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(displacements, i)
        {
            if( bndDisplacements.find(i) != bndDisplacements.end() )
            {
                displacements[i] = bndDisplacements.find(i)->second;
            }
            else
            {
                displacements[i] = vector::zero;
            }
        }

        std::map<label, scalar> localCoefficients;

        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(faceEdges, faceI)
        {
            //- calculate centres of adjacent cells. Use face centres
            //- for boundary faces
            const point& cOwn = cellCentres[owner[faceI]];
            point cNei;
            if( neighbour[faceI] >= 0 )
            {
                cNei = cellCentres[neighbour[faceI]];
            }
            else
            {
                cNei = help::faceAreaVector(points, faces[faceI]);
            }

            //- create a local area of a triangle formed by the edge centre
            //- and the centres of adjacent cells
            forAllRow(faceEdges, faceI, feI)
            {
                const label edgeI = faceEdges(faceI, feI);

                if( localCoefficients.find(edgeI) == localCoefficients.end() )
                    localCoefficients[edgeI] = 0.0;

                const edge& e = edges[edgeI];

                triangle<point, point> tri
                (
                    e.centre(points),
                    cOwn,
                    cNei
                );

                localCoefficients[edgeI] += tri.mag();
            }
        }

        //- update global coefficients
        # ifdef USE_OMP
        # pragma omp critical
        # endif
        {
            for
            (
                auto it=localCoefficients.begin();
                it!=localCoefficients.end();
                ++it
            )
            {
                edgeCoefficients[it->first] += it->second;
            }
        }

        localCoefficients.clear();

        //- normalize coefficients by the magnitude of the edge
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(edgeCoefficients, edgeI)
        {
            edgeCoefficients[edgeI] /= edges[edgeI].mag(points);

            const edge& e = edges[edgeI];
            const point& p0 = points[e.start()];
            const point& p1 = points[e.end()];

            const scalar r1Sq = magSqr(p0 - point(p0.x(), 0.0, 0.0));
            const scalar r2Sq = magSqr(p1 - point(p1.x(), 0.0, 0.0));

            edgeCoefficients[edgeI] /= (max(r1Sq, r2Sq) + VSMALL);
        }

        //- start updating displacements
        for(label iter=0;iter<100;++iter)
        {
            # ifdef USE_OMP
            # pragma omp barrier
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(displacements, pI)
            {
                //- skip updating points with fixed displacement
                if( bndDisplacements.find(pI) != bndDisplacements.end() )
                {
                    newDisplacements[pI] = displacements[pI];
                    continue;
                }

                //- update displacement
                vector& disp = newDisplacements[pI];
                disp = vector::zero;

                scalar sum = 0.0;

                forAllRow(pointEdges, pI, peI)
                {
                    const label edgeI = pointEdges(pI, peI);

                    const edge& e = edges[edgeI];

                    const label opI = e.otherVertex(pI);

                    disp += edgeCoefficients[edgeI] * displacements[opI];
                    sum += edgeCoefficients[edgeI];
                }

                //- calculate displacement
                disp /= (sum + VSMALL);

                if( fixedDiameterPoints.found(pI) )
                {
                    const point& p = points[pI];

                    const scalar r = sqrt(sqr(p.y()) + sqr(p.z()));

                    const point pNew = p + disp;

                    const scalar rNew = sqrt(sqr(pNew.y()) + sqr(pNew.z()));

                    vector v = pNew - point(pNew.x(), 0.0, 0.0);

                    disp += (r - rNew) * (v / (mag(v)+VSMALL));
                }

                //- do not move the point in x direction
                disp.x() = 0.0;

                //- update new displacement
                newDisplacements[pI] = 0.8 * disp;
            }

            //- update the original displacement
            # ifdef USE_OMP
            # pragma omp for schedule(dynamic, 100)
            # endif
            forAll(newDisplacements, i)
            {
                displacements[i] = newDisplacements[i];
            }
        }

        //- update displacements to mesh points
        scalar maxDisp = 0.0;
        scalar maxBndDisp = 0.0;

        # ifdef USE_OMP
        # pragma omp barrier
        # pragma omp for schedule(dynamic, 100)
        # endif
        forAll(displacements, pI)
        {
            if( bndDisplacements.find(pI) == bndDisplacements.end() )
            {
                maxDisp = max(maxDisp, mag(displacements[pI]));
            }
            else
            {
                maxBndDisp = max(maxBndDisp, mag(displacements[pI]));
            }

            meshModifier.movePoint(pI, points[pI] + displacements[pI]);
        }

        Info << "Max inner displacement " << maxDisp << endl;
        Info << "Max bnd displacement " << maxBndDisp << endl;
    }
}

void profiledDieMeshGenerator::updateMeshPointsFEM
(
    const std::map<label, vector>& displacements
)
{
    # ifdef ExtendSpecific
    //- write the mesh to disk and read it as polyMesh
    mesh_.write();

    //- create polyMesh from polyMeshGen
    Info << "reading polyMesh" << endl;
    polyMesh mesh
    (
        IOobject
        (
            "",
            mesh_.returnTime().constant(),
            mesh_.returnTime(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info << "Constructing tetPolyMesh" << endl;
    tetPolyMesh tetMesh(mesh);

    Info << "Constructing motion displacements" << endl;
    tetPointVectorField motionDisp
    (
        IOobject
        (
            "motionDisp",
            mesh_.returnTime().timeName(),
            mesh_.returnTime(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tetMesh,
        dimensionedVector("disp", dimLength, vector::zero),
        "fixedValue"
    );

    //- apply displacement to given points
    for(auto it=displacements.begin();it!=displacements.end();++it)
        motionDisp[it->first] = it->second;

    Info << "Creating motion matrix" << endl;
    tetFemMatrix<vector> motionMatrix
    (
        tetFem::laplacian(motionDisp)
    );

    Info << "Solving" << endl;
    dictionary d;
    d.add("motionDisp", dictionary());
    d.add("solver", "GAMG");
    d.add("smoother", "GaussSeidel");
    d.add("agglomerator", "faceAreaPair");
    d.add("nCellsInCoarsestLevel", 20);
    d.add("mergeLevels", 1);
    d.add("tolerance", 1e-6);
    motionMatrix.solve(d);
    # else
    # warning Cannot use FEM motion with this version of OpenFOAM
    # endif
}

void profiledDieMeshGenerator::updateMeshPoints()
{
    //- perform optimisation of the mesh to ensure there are no poor quality
    //- cells in the mesh. the procedure starts by optimising the surface
    //- of the mesh, followed by the vertices inside the mesh.
    const meshSurfaceEngine& mse = surfaceEngine();

    {
        meshSurfaceEngineModifier(mse).updateGeometry();

        //- optimize surface vertices
        //- lock boundary vertices except upstream, downstream and symmetry
        //- patches
        labelLongList lockedBndPoints;
        detectLockedBndPoints(lockedBndPoints);

        meshSurfaceOptimizer surfOptimizer(mse);
        surfOptimizer.lockFeatureEdges();
        surfOptimizer.lockBoundaryPoints(lockedBndPoints);

        surfOptimizer.untangleSurface(0, 10, 20, 0);
        surfOptimizer.optimizeSurface(10);

        mesh_.addressingData().updateGeometry();
    }

    meshOptimizer mOpt(mesh_);
    mOpt.optimizeMeshNearBoundaries(2, 1000);
    mOpt.optimizeMeshFV(2, 10, 50, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

profiledDieMeshGenerator::profiledDieMeshGenerator
(
    polyMeshGen& mesh,
    const dieGeometryInfo& dieGeom,
    const rollingMillPatchNamesHandler& patchHandler
)
:
    mesh_(mesh),
    dieGeom_(dieGeom),
    patchNamesHandler_(patchHandler),
    interpolatorPtr_(),
    msePtr_(),
    bndDisplacements_()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

profiledDieMeshGenerator::~profiledDieMeshGenerator()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void profiledDieMeshGenerator::generateProfiledDie()
{
    //- calculate displacements of boundary vertices
    calculateDisplacements();

    //- calculate new positions of internal points and modify the mesh
//    updateMeshPoints();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
