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

#include "rollerSurfaceCreatorFlat.H"
#include "demandDrivenData.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "IFstream.H"
#include "boundBox.H"

#include "addToRunTimeSelectionTable.H"

//#define DEBUGSurfaceCreator

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rollerSurfaceCreatorFlat, 0);
addToRunTimeSelectionTable
(
    rollerSurfaceCreator,
    rollerSurfaceCreatorFlat,
    dictionary
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rollerSurfaceCreatorFlat::createTriangulatedSurface()
{
    //- create surface mesh
    triSurf surf;
    triSurfModifier sMod(surf);

    //- create surface points
    pointField& pts = sMod.pointsAccess();

    const label contactId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERCONTACT
            )
        );

    const label frontId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERFRONT
            )
        );

    const label backId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERBACK
            )
        );

    const label axisId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERAXIS
            )
        );

    if( contactWidth_ > SMALL )
    {
        const label rollerToAirId =
            surf.addEdgeSubset
            (
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERTOAIR
                )
            );

        //- use localised contact patch
        pts.setSize(6);

        contactWidth_ = min(contactWidth_, rollWidth_);

        pts[0] = vector(0.0, 0.5 * innerDiameter_, -0.5 * rollWidth_);
        pts[1] = vector(0.0, 0.5 * innerDiameter_, 0.5 * rollWidth_);
        pts[2] = vector(0.0, 0.5 * outerDiameter_, 0.5 * rollWidth_);
        pts[3] = vector(0.0, 0.5 * outerDiameter_, 0.5 * contactWidth_);
        pts[4] = vector(0.0, 0.5 * outerDiameter_, -0.5 * contactWidth_);
        pts[5] = vector(0.0, 0.5 * outerDiameter_, -0.5 * rollWidth_);

        //- add contact region
        surf.appendFeatureEdge(edge(3, 4));
        surf.addEdgeToSubset(contactId, 0);

        //- roller-to-air region
        surf.appendFeatureEdge(edge(2, 3));
        surf.addEdgeToSubset(rollerToAirId, 1);

        surf.appendFeatureEdge(edge(4, 5));
        surf.addEdgeToSubset(rollerToAirId, 2);

        //- add axis region
        surf.appendFeatureEdge(edge(0, 1));
        surf.addEdgeToSubset(axisId, 3);

        //- add front region
        surf.appendFeatureEdge(edge(1, 2));
        surf.addEdgeToSubset(frontId, 4);

        //- add back region
        surf.appendFeatureEdge(edge(0, 5));
        surf.addEdgeToSubset(backId, 5);
    }
    else
    {
        //- large contact patch
        pts.setSize(4);

        pts[0] = vector(0.0, 0.5 * innerDiameter_, -0.5 * rollWidth_);
        pts[1] = vector(0.0, 0.5 * innerDiameter_, 0.5 * rollWidth_);
        pts[2] = vector(0.0, 0.5 * outerDiameter_, 0.5 * rollWidth_);
        pts[3] = vector(0.0, 0.5 * outerDiameter_, -0.5 * rollWidth_);

        //- add contact region
        surf.appendFeatureEdge(edge(2, 3));
        surf.addEdgeToSubset(contactId, 0);

        //- add axis region
        surf.appendFeatureEdge(edge(0, 1));
        surf.addEdgeToSubset(axisId, 1);

        //- add front region
        surf.appendFeatureEdge(edge(1, 2));
        surf.addEdgeToSubset(frontId, 2);

        //- add back region
        surf.appendFeatureEdge(edge(0, 3));
        surf.addEdgeToSubset(backId, 3);
    }

    const label cornersId = surf.addPointSubset("_corners_");
    forAll(pts, i)
        surf.addPointToSubset(cornersId, i);

    //- transform points into the x-y plane
    forAll(pts, pI)
        pts[pI] = (transformationMatrix_ & pts[pI]);

    //- extrude feature edges to a 2D surface
    triSurfaceExtrude2DEdges extruder(surf);

    surfPtr_ = new triSurf();
    extruder.extrudeSurface(*surfPtr_);

    //- create patches
    triSurfModifier surfMod(*surfPtr_);

    const tensor inverseTransformation = inv(transformationMatrix_);
    pointField& surfPoints = surfMod.pointsAccess();
    forAll(surfPoints, pI)
    {
        point& p = surfPoints[pI];

        p = (inverseTransformation & p);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollerSurfaceCreatorFlat::rollerSurfaceCreatorFlat
(
    const word rollPosition,
    const direction symm,
    const rollingMillPatchNamesHandler& patchHandler,
    const dictionary& dict,
    const scalar tol
)
:
    rollerSurfaceCreator
    (
        rollPosition,
        symm,
        patchHandler,
        dict,
        tol
    )
{}

rollerSurfaceCreatorFlat::rollerSurfaceCreatorFlat
(
    const rollerSurfaceCreatorFlat& rollerCreator
)
:
    rollerSurfaceCreator(rollerCreator)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rollerSurfaceCreatorFlat::generateGeometry()
{
    createTriangulatedSurface();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
