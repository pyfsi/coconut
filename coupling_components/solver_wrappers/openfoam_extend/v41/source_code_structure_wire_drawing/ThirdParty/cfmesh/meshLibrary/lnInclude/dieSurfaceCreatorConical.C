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

#include "dieSurfaceCreatorConical.H"
#include "rollingMillMesh.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "addToRunTimeSelectionTable.H"

//#define DEBUGDie

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(dieSurfaceCreatorConical, 0);
addToRunTimeSelectionTable
(
    dieSurfaceCreator,
    dieSurfaceCreatorConical,
    dictionary
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dieSurfaceCreatorConical::createTriangulatedSurface()
{
    triSurf surf;
    triSurfModifier sMod(surf);

    pointField& pts = sMod.pointsAccess();

    pts.setSize(4);

    pts[0] = point(0.0, outletRadius_, 0.0);
    pts[1] = point(0.0, outerRadius_, 0.0);
    pts[2] = point(-axialLength_, outerRadius_, 0.0);
    pts[3] = point(-axialLength_, inletRadius_, 0.0);

    //- generate feature edge that will extruded into a ribbon
    surf.appendFeatureEdge(edge(0, 1));
    surf.appendFeatureEdge(edge(1, 2));
    surf.appendFeatureEdge(edge(2, 3));
    surf.appendFeatureEdge(edge(3, 0));

    const label downstreamId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEDOWNSTREAM
            )
        );

    const label outerId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEHOUSING
            )
        );

    const label inletId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEUPSTREAM
            )
        );

    const label contactId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEWIRE
            )
        );

    //- add edges to subset
    surf.addEdgeToSubset(downstreamId, 0);
    surf.addEdgeToSubset(outerId, 1);
    surf.addEdgeToSubset(inletId, 2);
    surf.addEdgeToSubset(contactId, 3);

    surfPtr_ = new triSurf();

    triSurfaceExtrude2DEdges extruder(surf);
    extruder.extrudeSurface(*surfPtr_);

    geometricSurfacePatchList& sPatches =
        triSurfModifier(*surfPtr_).patchesAccess();

    forAll(sPatches, patchI)
        sPatches[patchI].geometricType() = "patch";

    Info << "Number of surface triangles " << surfPtr_->size() << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dieSurfaceCreatorConical::
dieSurfaceCreatorConical
(
    const rollingMillPatchNamesHandler& patchHandler,
    const dictionary& dict,
    const scalar tol
)
:
    dieSurfaceCreator(patchHandler, dict, tol),
    inletRadius_(0.0),
    outletRadius_(0.0),
    axialLength_(0.0),
    outerRadius_(0.0)
{
    if( !dict.readIfPresent("inletDiameter", inletRadius_) )
    {
        FatalError << "inletDiameter no found in dictionary "
                << dict << exit(FatalError);
    }

    inletRadius_ *= 0.5;

    if( !dict.readIfPresent("outletDiameter", outletRadius_) )
    {
        FatalError << "outletDiameter no found in dictionary "
                << dict << exit(FatalError);
    }

    outletRadius_ *= 0.5;

    if( !dict.readIfPresent("axialLength", axialLength_) )
    {
        FatalError << "axialLength no found in dictionary "
                << dict << exit(FatalError);
    }

    if( !dict.readIfPresent("outerDiameter", outerRadius_) )
    {
        FatalError << "outerDiameter no found in dictionary "
                << dict << exit(FatalError);
    }

    outerRadius_ *= 0.5;

    createTriangulatedSurface();
}

dieSurfaceCreatorConical::
dieSurfaceCreatorConical
(
    const dieSurfaceCreatorConical& die
)
:
    dieSurfaceCreator
    (
        die.patchHandler_,
        die.dict_,
        die.geometryTol_
    ),
    inletRadius_(die.inletRadius_),
    outletRadius_(die.outletRadius_),
    axialLength_(die.axialLength_),
    outerRadius_(die.outerRadius_)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar dieSurfaceCreatorConical::inletDiameter() const
{
    return 2. * inletRadius_;
}

scalar dieSurfaceCreatorConical::outletDiameter() const
{
    return 2. * outletRadius_;
}

scalar dieSurfaceCreatorConical::outerDiameter() const
{
    return 2.0 * outerRadius_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
