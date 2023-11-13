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

#include "dieSurfaceCreatorBekaertSpecs.H"
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

defineTypeNameAndDebug(dieSurfaceCreatorBekaertSpecs, 0);
addToRunTimeSelectionTable
(
    dieSurfaceCreator,
    dieSurfaceCreatorBekaertSpecs,
    dictionary
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dieSurfaceCreatorBekaertSpecs::createLine
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

void dieSurfaceCreatorBekaertSpecs::createArc
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

void dieSurfaceCreatorBekaertSpecs::createTriangulatedSurface
(
    const scalar H,
    const scalar W,
    const scalar Con,
    const scalar TC,
    const scalar gamma2,
    const scalar BC,
    const scalar gamma1,
    const scalar CR,
    const scalar CRA,
    const scalar EBR,
    const scalar EBL,
    const scalar ECA,
    const scalar ECL,
    const scalar RCA,
    const scalar RCL,
    const scalar RCR,
    const scalar D,
    const scalar BL,
    const scalar BRA,
    const scalar BRL,
    const scalar BRR,
    const scalar XCA,
    const scalar XCL,
    const scalar XCR
)
{
    const scalar ConR = Con * M_PI / 180.0;
    const scalar gamma2R = gamma2 * M_PI / 180.0;
    const scalar gamma1R = gamma1 * M_PI / 180.0;
    const scalar CRAR = CRA * M_PI / 180.0;
    const scalar ECAR = ECA * M_PI / 180.0;
    const scalar RCAR = RCA * M_PI / 180.0;
    const scalar BRAR = BRA * M_PI / 180.0;
    const scalar XCAR = XCA * M_PI / 180.0;

    // Calculate distances between 2 neighbouring points A and B in
    // x-direction. This is represented as dx_A_B
//    const scalar dx_11_12 = 0.;
    const scalar dx_10_11 = XCR - XCR*sin(XCAR/2.);
    const scalar dx_9_10 = XCL - dx_10_11;
    const scalar dx_7_8 = BRR*sin(BRAR/2.);
    const scalar dx_8_9 = BRL - dx_7_8;
    const scalar dx_6_7 = BL;
    const scalar dx_5_6 = RCR*sin(RCAR/2.);
    const scalar dx_4_5 = RCL - dx_5_6;
    const scalar dx_3_4 = ECL;
    const scalar dx_1_2 = CR*sin(CRAR);
    const scalar dx_2_3 = EBL - dx_1_2;
    const scalar dx_12_13 = BC*sin(gamma1R);
    const scalar dx_14_15 = TC*sin(gamma2R);
    const scalar dx_13_14 = H - dx_12_13 - dx_14_15;
//    const scalar dx_15_1 = 0.;

    // Calculate distances between 2 neighbouring points A and B in
    // y-direction. This is represented as dy_A_B
    const scalar dy_10_11 = XCR*cos(XCAR/2.);
    const scalar dy_9_10 = dx_9_10*tan(XCAR/2.);
    const scalar dy_8_9 = dx_8_9*tan(BRAR/2.);
    const scalar dy_7_8 = BRR - BRR*cos(BRAR/2.);
//    const scalar dy_6_7 = 0;
    const scalar dy_5_6 = RCR - RCR*cos(RCAR/2.);
    const scalar dy_4_5 = dx_4_5*tan(RCAR/2.);
    const scalar dy_3_4 = ECL*tan(ECAR/2.);
    const scalar dy_2_3_arg =
        sqr(2.*EBR*sin(M_PI/2.-CRAR-ECAR/2.)) - sqr(dx_2_3);
    if( dy_2_3_arg < 0.0 )
    {
        FatalError
          << "pow(2.*EBR*sin(M_PI/2. - CRAR - ECAR/2.), 2) - pow(dx_2_3, 2)"
          << " is negative. " << dy_2_3_arg
          << " Please update your settings " << exit(FatalError);
    }
    const scalar dy_2_3 = sqrt(dy_2_3_arg);
    const scalar dy_1_2 = CR*cos(CRAR);
    const scalar dy_12_13 = BC*cos(gamma1R);
    const scalar dy_13_14 = dx_13_14*tan(ConR);
    const scalar dy_14_15 = TC*cos(gamma2R);
    const scalar dy_15_1 = W/2. - D/2. - dy_14_15 - dy_1_2 - dy_2_3 - dy_3_4 -
                           dy_4_5 - dy_5_6 - dx_14_15*tan(ConR);
    const scalar dy_11_12 = W/2. - D/2. - dy_13_14 - dy_12_13 - dy_7_8 -
                            dy_8_9 - dy_9_10 - dy_10_11 - dx_14_15*tan(ConR);

    FixedList<point, 15> profilePoints(vector::zero);

    // Point 1
    profilePoints[0].x() = -(dx_1_2 + dx_2_3 + dx_3_4 + dx_4_5 + dx_5_6 +
                              dx_6_7 + dx_7_8 + dx_8_9 + dx_9_10 + dx_10_11);
    profilePoints[0].y() = dy_5_6 + dy_4_5 + dy_3_4 + dy_2_3 + dy_1_2 + D/2.;

    // Point 2
    profilePoints[1].x() = -(dx_2_3 + dx_3_4 + dx_4_5 + dx_5_6 + dx_6_7 +
                              dx_7_8 + dx_8_9 + dx_9_10 + dx_10_11);
    profilePoints[1].y() = dy_5_6 + dy_4_5 + dy_3_4 + dy_2_3 + D/2.;

    // Point 3
    profilePoints[2].x() = -(dx_3_4 + dx_4_5 + dx_5_6 + dx_6_7 + dx_7_8 +
                              dx_8_9 + dx_9_10 + dx_10_11);
    profilePoints[2].y() = dy_5_6 + dy_4_5 + dy_3_4 + D/2.;

    // Point 4
    profilePoints[3].x() = -(dx_4_5 + dx_5_6 + dx_6_7 + dx_7_8 + dx_8_9 +
                              dx_9_10 + dx_10_11);
    profilePoints[3].y() = dy_5_6 + dy_4_5 + D/2.;

    // Point 5
    profilePoints[4].x() = -(dx_5_6 + dx_6_7 + dx_7_8 + dx_8_9 + dx_9_10 +
                              dx_10_11);
    profilePoints[4].y() = dy_5_6 + D/2.;

    // Point 6
    profilePoints[5].x() = -(dx_6_7 + dx_7_8 + dx_8_9 + dx_9_10 + dx_10_11);
    profilePoints[5].y() =  D/2.;

    // Point 7
    profilePoints[6].x() = -(dx_7_8 + dx_8_9 + dx_9_10 + dx_10_11);
    profilePoints[6].y() =  D/2.;

    // Point 8
    profilePoints[7].x() = -(dx_8_9 + dx_9_10 + dx_10_11);
    profilePoints[7].y() = dy_7_8 + D/2.;

    // Point 9
    profilePoints[8].x() = -(dx_9_10 + dx_10_11);
    profilePoints[8].y() = dy_7_8 + dy_8_9 + D/2.;

    // Point 10
    profilePoints[9].x() = -(dx_10_11);
    profilePoints[9].y() = dy_7_8 + dy_8_9 + dy_9_10 + D/2.;

    // Point 11
    profilePoints[10].x() = 0;
    profilePoints[10].y() = dy_7_8 + dy_8_9 + dy_9_10 + dy_10_11 + D/2.;

    // Point 12
    profilePoints[11].x() = 0;
    profilePoints[11].y() = dy_7_8 + dy_8_9 + dy_9_10 + dy_10_11 + dy_11_12 +
                               D/2.;

    // Point 13
    profilePoints[12].x() = -(dx_12_13);
    profilePoints[12].y() = dy_7_8 + dy_8_9 + dy_9_10 + dy_10_11 + dy_11_12 +
                               dy_12_13 + D/2.;

    // Point 14
    profilePoints[13].x() = -(dx_12_13 + dx_13_14);
    profilePoints[13].y() = dy_7_8 + dy_8_9 + dy_9_10 + dy_10_11 + dy_11_12 +
                               dy_12_13 + dy_13_14 + D/2.;

    // Point 15
    profilePoints[14].x() = -(dx_12_13 + dx_13_14 + dx_14_15);
    profilePoints[14].y() = dy_5_6 + dy_4_5 + dy_3_4 + dy_2_3 + dy_1_2 +
                               dy_15_1 + D/2.;

    //- set inlet diameter
    inletDiameter_ = 2. * profilePoints[0].y();

    //- set outlet diameter
    outletDiameter_ = 2. * profilePoints[10].y();

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
        rollingMillPatchNamesHandler::DIEENTRANCECONE,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- arc 1-2
    const scalar dist_1_2 = mag(profilePoints[2] - profilePoints[1]);

    const point eCentre
    (
        0.5 * (profilePoints[1] + profilePoints[2]) +
        sqrt
        (
            EBR * EBR - (0.5 * dist_1_2) * (0.5 * dist_1_2)
        ) *
        (
            vector(0., 0., 1.0) ^ (profilePoints[2] - profilePoints[1])
        ) /
        mag(profilePoints[2] - profilePoints[1])
    );

    createArc
    (
        profilePoints[1],
        profilePoints[2],
        eCentre,
        rollingMillPatchNamesHandler::DIEENTRANCECONE,
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
        rollingMillPatchNamesHandler::DIEENTRANCECONE,
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
        rollingMillPatchNamesHandler::DIEWIRE,
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
        rollingMillPatchNamesHandler::DIEWIRE,
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
        rollingMillPatchNamesHandler::DIEWIRE,
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
        rollingMillPatchNamesHandler::DIEWIRE,
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
        rollingMillPatchNamesHandler::DIEWIRE,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 8-9
    createLine
    (
        profilePoints[8],
        profilePoints[9],
        rollingMillPatchNamesHandler::DIEEXITCONE,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 9-10
    createLine
    (
        profilePoints[9],
        profilePoints[10],
        rollingMillPatchNamesHandler::DIEEXITCONE,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 10-11
    createLine
    (
        profilePoints[10],
        profilePoints[11],
        rollingMillPatchNamesHandler::DIEDOWNSTREAM,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 11-12
    createLine
    (
        profilePoints[11],
        profilePoints[12],
        rollingMillPatchNamesHandler::DIEDOWNSTREAM,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 12-13
    createLine
    (
        profilePoints[12],
        profilePoints[13],
        rollingMillPatchNamesHandler::DIEHOUSING,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 13-14
    createLine
    (
        profilePoints[13],
        profilePoints[14],
        rollingMillPatchNamesHandler::DIEUPSTREAM,
        points,
        edges,
        edgePatchCodes,
        featurePoints
    );

    //- line 14-0
    createLine
    (
        profilePoints[14],
        profilePoints[0],
        rollingMillPatchNamesHandler::DIEUPSTREAM,
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
        const word pName = patchHandler_.patchNameForDie
        (
            rollingMillPatchNamesHandler::diePatchKeys(it->first)
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

/*    pts.setSize(15);

    pts[0] = profilePoints[0];
    pts[1] = profilePoints[1];
    pts[2] = profilePoints[2];
    pts[3] = profilePoints[3];
    pts[4] = profilePoints[4];
    pts[5] = profilePoints[5];
    pts[6] = profilePoints[6];
    pts[7] = profilePoints[7];
    pts[8] = profilePoints[8];
    pts[9] = profilePoints[9];
    pts[10] = profilePoints[10];
    pts[11] = profilePoints[11];
    pts[12] = profilePoints[12];
    pts[13] = profilePoints[13];
    pts[14] = profilePoints[14];

    //- create edges of the cross-section. Last edge between point 15 and
    // point 1 on upstream side is not created here.
    //- entrance cone
    surfPtr_->appendFeatureEdge(edge(0, 1));
    surfPtr_->appendFeatureEdge(edge(1, 2));
    surfPtr_->appendFeatureEdge(edge(2, 3));

    //- contact region
    surfPtr_->appendFeatureEdge(edge(3, 4));
    surfPtr_->appendFeatureEdge(edge(4, 5));
    surfPtr_->appendFeatureEdge(edge(5, 6));
    surfPtr_->appendFeatureEdge(edge(6, 7));
    surfPtr_->appendFeatureEdge(edge(7, 8));

    //- exit cone
    surfPtr_->appendFeatureEdge(edge(8, 9));
    surfPtr_->appendFeatureEdge(edge(9, 10));

    //- die downstream
    surfPtr_->appendFeatureEdge(edge(10, 11));
    surfPtr_->appendFeatureEdge(edge(11, 12));

    //- die to housing
    surfPtr_->appendFeatureEdge(edge(12, 13));

    //- die upstream
    surfPtr_->appendFeatureEdge(edge(13, 14));
    surfPtr_->appendFeatureEdge(edge(14, 0));

    //- extrude edges into a 2D surface
    triSurfaceExtrude2DEdges extruder(*surfPtr_);
    const triSurf* extrudedPtr = extruder.extrudeSurface();

    deleteDemandDrivenData(surfPtr_);
    surfPtr_ = new triSurf(*extrudedPtr);
    deleteDemandDrivenData(extrudedPtr);

    //- add feature edges
    for(label eI=0;eI<15;++eI)
        surfPtr_->appendFeatureEdge(edge(eI, eI+15));
*/

    # ifdef DEBUGDie
    surfPtr_->writeSurface("conicalDie.stl");
    Info << "Conical die written" << endl;
    ::exit(0);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dieSurfaceCreatorBekaertSpecs::
dieSurfaceCreatorBekaertSpecs
(
    const rollingMillPatchNamesHandler& patchHandler,
    const dictionary& dict,
    const scalar tol
)
:
    dieSurfaceCreator(patchHandler, dict, tol),
    inletDiameter_(0.0),
    outletDiameter_(0.0)
{
    // die geometry parameters accoding Bekaert specification
    //- die heigth, i.e. axial length
    scalar H;
    if( !dict.readIfPresent("H", H) )
    {
        FatalError << "H is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- die width, i.e., the outer diameter
    scalar W;
    if( !dict.readIfPresent("W", W) )
    {
        FatalError << "W is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- die outer surface conicity
    scalar Con;
    if( !dict.readIfPresent("Con", Con) )
    {
         FatalError << "Con is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- top corner chamger length
    scalar TC;
    if( !dict.readIfPresent("TC", TC) )
    {
         FatalError << "TC is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- top corner chamfer angle
    scalar gamma2;
    if( !dict.readIfPresent("gamma2", gamma2) )
    {
         FatalError << "gamma2 is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- bottom corner chamfer length
    scalar BC;
    if( !dict.readIfPresent("BC", BC) )
    {
         FatalError << "BC is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- bottom corner chamfer angle
    scalar gamma1;
    if( !dict.readIfPresent("gamma1", gamma1) )
    {
         FatalError << "gamma1 is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- entrance bell chamfer length
    scalar CR;
    if( !dict.readIfPresent("C_R", CR) )
    {
         FatalError << "C_R is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- entrance bell chamfer length
    scalar CRA;
    if( !dict.readIfPresent("CR_A", CRA) )
    {
         FatalError << "CR_A is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- entrance bell radius
    scalar EBR;
    if( !dict.readIfPresent("EB_R", EBR) )
    {
         FatalError << "EB_R is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- entrance cone length in x-dir, consisting of
    //- entrance bell and chamfer x-dir length
    scalar EBL;
    if( !dict.readIfPresent("EB_L", EBL) )
    {
         FatalError << "EB_L is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- entrance cone angle
    scalar ECA;
    if( !dict.readIfPresent("EC_A", ECA) )
    {
         FatalError << "EC_A is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- entrance cone length
    scalar ECL;
    if( !dict.readIfPresent("EC_L", ECL) )
    {
         FatalError << "EC_L is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- reduction cone angle
    scalar RCA;
    if( !dict.readIfPresent("RC_A", RCA) )
    {
         FatalError << "RC_A is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- reduction come length in x-dir
    scalar RCL;
    if( !dict.readIfPresent("RC_L", RCL) )
    {
         FatalError << "EC_L is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- rounding radius between reduction come and bearing
    scalar RCR;
    if( !dict.readIfPresent("RC_R", RCR) )
    {
         FatalError << "RC_R is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- die inner diameter
    scalar D;
    if( !dict.readIfPresent("D", D) )
    {
         FatalError << "D is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- bearing length
    scalar BL;
    if( !dict.readIfPresent("BL", BL) )
    {
         FatalError << "BL is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- back relief angle
    scalar BRA;
    if( !dict.readIfPresent("BR_A", BRA) )
    {
         FatalError << "BR_A is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- back relief length
    scalar BRL;
    if( !dict.readIfPresent("BR_L", BRL) )
    {
         FatalError << "BR_L is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- rounding radius between bearing and back relief
    scalar BRR;
    if( !dict.readIfPresent("BR_R", BRR) )
    {
         FatalError << "BR_R is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- exit cone angle
    scalar XCA;
    if( !dict.readIfPresent("XC_A", XCA) )
    {
         FatalError << "XC_A is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- exit cone length in x-dir
    scalar XCL;
    if( !dict.readIfPresent("XC_L", XCL) )
    {
         FatalError << "XC_L is not present in dictionary " << dict
            << exit(FatalError);
    }

    //- rounding radius between exit cone and downstream surface
    scalar XCR;
    if( !dict.readIfPresent("XC_R", XCR) )
    {
         FatalError << "XCR is not present in dictionary " << dict
            << exit(FatalError);
    }

    createTriangulatedSurface
    (
        H,
        W,
        Con,
        TC,
        gamma2,
        BC,
        gamma1,
        CR,
        CRA,
        EBR,
        EBL,
        ECA,
        ECL,
        RCA,
        RCL,
        RCR,
        D,
        BL,
        BRA,
        BRL,
        BRR,
        XCA,
        XCL,
        XCR
    );
}

dieSurfaceCreatorBekaertSpecs::
dieSurfaceCreatorBekaertSpecs
(
    const dieSurfaceCreatorBekaertSpecs& die
)
:
    dieSurfaceCreator
    (
        die.patchHandler_,
        die.dict_,
        die.geometryTol_
    ),
    inletDiameter_(die.inletDiameter_),
    outletDiameter_(die.outletDiameter_)
{
    surfPtr_ = new triSurf(die.surface());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar dieSurfaceCreatorBekaertSpecs::inletDiameter() const
{
    return inletDiameter_;
}

scalar dieSurfaceCreatorBekaertSpecs::outletDiameter() const
{
    return outletDiameter_;
}

scalar dieSurfaceCreatorBekaertSpecs::outerDiameter() const
{
    const pointField& points = surfPtr_->points();

    scalar d(0.0);

    dict_.readIfPresent("outerDiameter", d);

    forAll(points, i)
        d = max(d, 2. * points[i].y());

    return d;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
