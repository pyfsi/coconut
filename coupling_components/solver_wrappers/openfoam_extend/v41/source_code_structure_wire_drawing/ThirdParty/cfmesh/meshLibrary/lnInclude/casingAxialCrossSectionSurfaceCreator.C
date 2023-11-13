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

#include "casingAxialCrossSectionSurfaceCreator.H"
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

defineTypeNameAndDebug(casingAxialCrossSectionSurfaceCreator, 0);
defineRunTimeSelectionTable(casingAxialCrossSectionSurfaceCreator, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

casingAxialCrossSectionSurfaceCreator::casingAxialCrossSectionSurfaceCreator
(
    const rollingMillPatchNamesHandler& patchNames,
    const dictionary& dict,
    const scalar tol
)
:
    dict_(dict),
    patchHandler_(patchNames),
    geometryTol_(tol),
    surfPtr_(NULL)
{}

casingAxialCrossSectionSurfaceCreator::~casingAxialCrossSectionSurfaceCreator()
{
    deleteDemandDrivenData(surfPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const triSurf& casingAxialCrossSectionSurfaceCreator::surface() const
{
    return *surfPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
