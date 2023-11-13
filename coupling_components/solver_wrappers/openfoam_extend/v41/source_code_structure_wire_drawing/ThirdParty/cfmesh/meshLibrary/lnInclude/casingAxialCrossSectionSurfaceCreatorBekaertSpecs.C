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

#include "casingAxialCrossSectionSurfaceCreatorBekaertSpecs.H"
#include "rollingMillMesh.H"
#include "demandDrivenData.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "polyMeshGenModifier.H"
#include "helperFunctions.H"
#include "addToRunTimeSelectionTable.H"

//#define DEBUGDie

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(casingAxialCrossSectionSurfaceCreatorBekaertSpecs, 0);
addToRunTimeSelectionTable
(
    casingAxialCrossSectionSurfaceCreator,
    casingAxialCrossSectionSurfaceCreatorBekaertSpecs,
    dictionary
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void casingAxialCrossSectionSurfaceCreatorBekaertSpecs::createLine
(
    const point& firstPoint,
    const point& secondPoint,
    const label code,
    LongList<point>& points,
    edgeLongList& edges,
    labelLongList& edgePatchCode,
    labelLongList& featurePoints
) const
{
    const label nPoints = points.size();

    points.append(firstPoint);
    points.append(secondPoint);

    edges.append(edge(nPoints, nPoints+1));

    edgePatchCode.append(code);

    featurePoints.append(nPoints);
    featurePoints.append(nPoints+1);
}

void casingAxialCrossSectionSurfaceCreatorBekaertSpecs::createArc
(
    const point& firstPoint,
    const point& secondPoint,
    const point& centre,
    const label code,
    LongList<point>& points,
    edgeLongList& edges,
    labelLongList& edgePatchCode,
    labelLongList& featurePoints
) const
{
    vector v0 = firstPoint - centre;
    const scalar r = mag(v0);
    v0 /= (r + VSMALL);

    vector v1 = secondPoint - centre;
    v1 /= (r + VSMALL);

    //- calculate angle step
    const scalar angle = ::asin((v0 ^ v1).z());

    scalar angleTol = Foam::acos(1.0 - geometryTol_ / (r + VSMALL));

    const label nDivisions = max(ceil(angle / angleTol), 2);

    angleTol = angle / nDivisions;

    //- create edges across the arc
    vector v = (vector(0., 0., 1.0) ^ v0);
    v /= (mag(v) + VSMALL);

    label nPoints = points.size();
    featurePoints.append(nPoints);
    points.append(firstPoint);
    scalar currAngle(0.0);

    for(label i=0;i<nDivisions;++i)
    {
        if( i == (nDivisions-1) )
        {
            currAngle = angle;
        }
        else
        {
            currAngle += angleTol;
        }

        //- create next point
        points.append
        (
            centre + r * (v0 * cos(currAngle) + v * sin(currAngle))
        );

        //- create new edge
        edges.append(edge(nPoints, nPoints+1));
        ++nPoints;

        //- append new patch code
        edgePatchCode.append(code);
    }

    featurePoints.append(nPoints);
}

void casingAxialCrossSectionSurfaceCreatorBekaertSpecs::
createTriangulatedSurface
(
    const scalar diameter,
    const scalar height,
    const scalar BC,
    const scalar TC,
    const scalar EX_OD,
    const scalar N_P,
    const scalar EX_ID,
    const scalar N_X,
    const scalar EN_ID,
    const scalar EN_OD
)
{
    const scalar sqrt2 = sqrt(2.0);
    FixedList<point, 9> profilePoints(vector::zero);

    // Point 1
    profilePoints[0].x() = N_P;
    profilePoints[0].y() = 0.5 * EX_OD;

    // Point 2
    profilePoints[1].x() = N_P;
    profilePoints[1].y() = 0.5 * (diameter - BC * sqrt2);

    // Point 3
    profilePoints[2].x() = N_P - BC * sqrt2 / 2.0;
    profilePoints[2].y() = 0.5 * diameter;

    // Point 4
    profilePoints[3].x() = N_P - (height - TC * sqrt2 / 2.0);
    profilePoints[3].y() = 0.5 * diameter;

    // Point 5
    profilePoints[4].x() = N_P - height;
    profilePoints[4].y() = 0.5 * (diameter - TC * sqrt2);

    // Point 6
    profilePoints[5].x() = N_P - height;
    profilePoints[5].y() = 0.5 * EN_OD;

    // Point 7
    profilePoints[6].x() = -N_X;
    profilePoints[6].y() = 0.5 * EN_ID;

    // Point 8
    profilePoints[7].y() = 0.5 * EN_ID;

    // Point 9
    profilePoints[8].y() = 0.5 * EX_ID;

    # ifdef DEBUGDie
    forAll(profilePoints, i)
    {
        Info << "Profile point " << i << " has coordinates "
             << profilePoints[i] << endl;
    }
    # endif

    //- create surface segments
    LongList<point> points;
    edgeLongList edges;
    labelLongList edgePatchCodes;
    labelLongList featurePoints;

    //- line 0-1
    createLine
    (
        profilePoints[0],
        profilePoints[1],
        rollingMillPatchNamesHandler::CASINGDOWNSTREAM,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 1-2
    createLine
    (
        profilePoints[1],
        profilePoints[2],
        rollingMillPatchNamesHandler::CASINGTOOUTSIDE,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 2-3
    createLine
    (
        profilePoints[2],
        profilePoints[3],
        rollingMillPatchNamesHandler::CASINGTOOUTSIDE,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 3-4
    createLine
    (
        profilePoints[3],
        profilePoints[4],
        rollingMillPatchNamesHandler::CASINGTOOUTSIDE,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 4-5
    createLine
    (
        profilePoints[4],
        profilePoints[5],
        rollingMillPatchNamesHandler::CASINGUPSTREAM,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 5-6
    createLine
    (
        profilePoints[5],
        profilePoints[6],
        rollingMillPatchNamesHandler::CASINGENTRANCE,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 6-7
    createLine
    (
        profilePoints[6],
        profilePoints[7],
        rollingMillPatchNamesHandler::CASINGTODIERADIAL,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 7-8
    createLine
    (
        profilePoints[7],
        profilePoints[8],
        rollingMillPatchNamesHandler::CASINGTODIEAXIAL,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 8-0
    createLine
    (
        profilePoints[8],
        profilePoints[0],
        rollingMillPatchNamesHandler::CASINGEXIT,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- generate the surface of a cross section
    triSurf surf;
    triSurfModifier surfMod(surf);
    pointField& pts = surfMod.pointsAccess();

    //- copy surface points
    pts.setSize(points.size());
    forAll(pts, i)
        pts[i] = points[i];

    //- create edge subsets
    std::map<label, label> patchCodeToId;

    forAll(edgePatchCodes, i)
        patchCodeToId[edgePatchCodes[i]] = -1;

    for
    (
        std::map<label, label>::iterator it=patchCodeToId.begin();
        it!=patchCodeToId.end();
        ++it
    )
    {
        const word pName = patchHandler_.patchNameForCasing
        (
            rollingMillPatchNamesHandler::casingPatchKeys(it->first)
        );

        it->second = surf.addEdgeSubset(pName);
    }

    //- add edges to the surface
    forAll(edges, eI)
    {
        surf.appendFeatureEdge(edges[eI]);
        surf.addEdgeToSubset(patchCodeToId[edgePatchCodes[eI]], eI);
    }

    //- add feature points
    const label cornerId = surf.addPointSubset("_corners_");
    forAll(featurePoints, pI)
        surf.addPointToSubset(cornerId, featurePoints[pI]);

    //- extrude edges into a 2D surface
    triSurfaceExtrude2DEdges extruder(surf);
    const triSurf* extrudedPtr = extruder.extrudeSurface();

    surfPtr_ = new triSurf(*extrudedPtr);

    deleteDemandDrivenData(extrudedPtr);

    # ifdef DEBUGDie
    surfPtr_->writeSurface("casingBeakertSpecs.stl");
    Info << "Casing geometry written" << endl;
    ::exit(0);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

casingAxialCrossSectionSurfaceCreatorBekaertSpecs::
casingAxialCrossSectionSurfaceCreatorBekaertSpecs
(
    const rollingMillPatchNamesHandler& patchHandler,
    const dictionary& dict,
    const scalar tol
)
:
    casingAxialCrossSectionSurfaceCreator(patchHandler, dict, tol)
{
    // casing geometry parameters accoding Bekaert specification
    //- casing diameter
    scalar diameter;
    if( !dict.readIfPresent("diameter", diameter) )
    {
        FatalError << "diameter is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- height of the casing
    scalar height;
    if( !dict.readIfPresent("height", height) )
    {
        FatalError << "height is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- BC chamfer
    scalar BC;
    if( !dict.readIfPresent("B_C", BC) )
    {
        FatalError << "B_C is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- TC chamfer
    scalar TC;
    if( !dict.readIfPresent("T_C", TC) )
    {
        FatalError << "T_C is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- EX outer diameter
    scalar EX_OD;
    if( !dict.readIfPresent("EX_OD", EX_OD) )
    {
        FatalError << "EX_OD is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- N_P point
    scalar N_P;
    if( !dict.readIfPresent("N_P", N_P) )
    {
        FatalError << "N_P is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- EX inner diameter
    scalar EX_ID;
    if( !dict.readIfPresent("EX_ID", EX_ID) )
    {
        FatalError << "EX_ID is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- N_X
    scalar N_X;
    if( !dict.readIfPresent("N_X", N_X) )
    {
        FatalError << "N_X is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- EN_ID inner diameter
    scalar EN_ID;
    if( !dict.readIfPresent("EN_ID", EN_ID) )
    {
        FatalError << "EN_ID is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- EN_OD inner diameter
    scalar EN_OD;
    if( !dict.readIfPresent("EN_OD", EN_OD) )
    {
        FatalError << "EN_OD is not present in dictionary " << dict
            << exit(FatalError);
    }

    createTriangulatedSurface
    (
        diameter,
        height,
        BC,
        TC,
        EX_OD,
        N_P,
        EX_ID,
        N_X,
        EN_ID,
        EN_OD
    );
}

casingAxialCrossSectionSurfaceCreatorBekaertSpecs::
casingAxialCrossSectionSurfaceCreatorBekaertSpecs
(
    const casingAxialCrossSectionSurfaceCreatorBekaertSpecs& die
)
:
    casingAxialCrossSectionSurfaceCreator
    (
        die.patchHandler_,
        die.dict_,
        die.geometryTol_
    )
{
    surfPtr_ = new triSurf(die.surface());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
