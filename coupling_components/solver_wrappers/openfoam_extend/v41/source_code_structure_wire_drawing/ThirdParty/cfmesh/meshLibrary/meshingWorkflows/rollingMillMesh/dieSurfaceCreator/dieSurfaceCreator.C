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

#include "dieSurfaceCreator.H"
#include "rollingMillMesh.H"
#include "demandDrivenData.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "polyMeshGenModifier.H"
#include "splineBase.H"
#include "cubicBSpline.H"
#include "IFstream.H"
#include "boundBox.H"

//#define DEBUGDie

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

std::shared_ptr<triSurf> dieSurfaceCreator::createCircularProfile
(
    const scalar d
) const
{
    //- calculate the number of subdivisions
    const scalar r = 0.5 * d;
    scalar angle = 2.0 * Foam::acos(1.0 - geometryTol_ / r);

    const label nDivisions = 4 * ceil(2.0 * M_PI / angle);

    angle = 2.0 * M_PI / nDivisions;

    //- allocate the surface mesh
    std::shared_ptr<triSurf> surfPtr = std::make_shared<triSurf>();

    triSurfModifier sMod(*surfPtr);

    pointField& pts = sMod.pointsAccess();
    pts.setSize(nDivisions);

    const label sId =
        surfPtr->addEdgeSubset
        (
            patchHandler_.patchNameForDie(rollingMillPatchNamesHandler::DIEWIRE)
        );

    forAll(pts, pI)
    {
        pts[pI] = point(r * cos(angle * pI), r * sin(angle * pI), 0.0);

        surfPtr->appendFeatureEdge(edge(pI, pts.fcIndex(pI)));
        surfPtr->addEdgeToSubset(sId, pI);
    }

    //- extrude edges into the surface
    triSurfaceExtrude2DEdges extruder(*surfPtr);
    const triSurf* extrudedPtr = extruder.extrudeSurface();

    surfPtr.reset(const_cast<triSurf*>(extrudedPtr));

    triSurfModifier ptsMod(*surfPtr);
    pointField& pPts = ptsMod.pointsAccess();
    forAll(pPts, pI)
    {
        point& p = pPts[pI];

        const scalar x = p.x();
        p.x() = p.z();
        p.z() = x;
    }

    return surfPtr;
}

void dieSurfaceCreator::positionCrossSectionAtOrigin(triSurf& surf) const
{
    //- check if the user wishes to centre a die at the origin
    bool dieProfileCentering(false);
    dict_.readIfPresent("dieProfileCentering", dieProfileCentering);

    if( dieProfileCentering )
    {
        boundBox bb(surf.points());

        triSurfModifier sMod(surf);

        vector disp = bb.midpoint();
        disp.x() = 0.0;

        Info << "Bounding box centre " << disp << endl;

        sMod.pointsAccess() -= disp;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(dieSurfaceCreator, 0);
defineRunTimeSelectionTable(dieSurfaceCreator, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dieSurfaceCreator::dieSurfaceCreator
(
    const rollingMillPatchNamesHandler& patchNames,
    const dictionary& dict,
    const scalar tol
)
:
    dict_(dict),
    patchHandler_(patchNames),
    geometryTol_(tol),
    surfPtr_(NULL),
    crossSections_()
{}

dieSurfaceCreator::~dieSurfaceCreator()
{
    deleteDemandDrivenData(surfPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const triSurf& dieSurfaceCreator::surface() const
{
    return *surfPtr_;
}

const DynList<std::pair<scalar, std::shared_ptr<triSurf> > >&
dieSurfaceCreator::crossSections() const
{
    return crossSections_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
