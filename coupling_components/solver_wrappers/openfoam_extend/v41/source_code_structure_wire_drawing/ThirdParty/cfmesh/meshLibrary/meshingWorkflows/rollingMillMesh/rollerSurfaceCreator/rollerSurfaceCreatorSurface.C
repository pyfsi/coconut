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

#include "rollerSurfaceCreatorSurface.H"
#include "demandDrivenData.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "polyMeshGenModifier.H"
#include "splineBase.H"
#include "cubicBSpline.H"
#include "IFstream.H"
#include "boundBox.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "triSurfaceCleanupDuplicates.H"

#include "addToRunTimeSelectionTable.H"

#include <map>

//#define DEBUGSurfaceCreator

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rollerSurfaceCreatorSurface, 0);
addToRunTimeSelectionTable
(
    rollerSurfaceCreator,
    rollerSurfaceCreatorSurface,
    dictionary
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rollerSurfaceCreatorSurface::readGeometry()
{
    surfPtr_ = new triSurf(fName_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollerSurfaceCreatorSurface::rollerSurfaceCreatorSurface
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

rollerSurfaceCreatorSurface::rollerSurfaceCreatorSurface
(
    const rollerSurfaceCreatorSurface& creator
)
:
    rollerSurfaceCreator(creator)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rollerSurfaceCreatorSurface::generateGeometry()
{
    notImplemented("void rollerSurfaceCreatorSurface::generateGeometry()");

    if( !dict_.readIfPresent("surfaceFile", fName_) )
    {
        FatalErrorIn
        (
            "void rollerSurfaceCreatorSurface::generateGeometry()"
        ) << "dxfFile is not present in dictionary "
          << dict_ << exit(FatalError);
    }

    //- geometry is given as a dxf file
    if( rollPosition_.empty() )
    {
        FatalErrorIn
        (
            "void rollerSurfaceCreatorSurface::generateGeometry()"
        ) << "Invalid roll position " << rollPosition_ << exit(FatalError);
    }

    readGeometry();

    detectIndependentRegions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
