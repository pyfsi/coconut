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

#include "triSurfaceDetectFeatureEdges.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"
#include "triSurfModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurfaceDetectFeatureEdges::triSurfaceDetectFeatureEdges(triSurf& surface)
:
    surf_(surface),
    featureEdges_(surf_.edges().size(), direction(0)),
    angleTolerance_(-1.0),
    useDanglingEdges_(false)
{
    if( Pstream::parRun() )
        FatalError << "Feature edges detection does not run in parallel"
            << exit(FatalError);
}

triSurfaceDetectFeatureEdges::~triSurfaceDetectFeatureEdges()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceDetectFeatureEdges::setAngleTolerance(const scalar angleTol)
{
    angleTolerance_ = angleTol;
}

void triSurfaceDetectFeatureEdges::keepDanglingEdges()
{
    useDanglingEdges_ = true;
}

void triSurfaceDetectFeatureEdges::detectFeatureEdges()
{
    const edgeLongList& edges = surf_.edges();
    triSurfModifier surfMod(surf_);
    edgeLongList& featureEdges = surfMod.featureEdgesAccess();

    //- preserve existing feature edges
    const VRWGraph& pointEdges = surf_.pointEdges();
    forAll(featureEdges, eI)
    {
        const edge& e = featureEdges[eI];

        forAllRow(pointEdges, e.start(), peI)
        {
            const label edgeI = pointEdges(e.start(), peI);

            if( edges[edgeI] == e )
            {
                featureEdges_[edgeI] |= 1;
                break;
            }
        }
    }

    featureEdges.clear();

    //- detect feature edges
    detectFeatureEdgesAngleCriterion();

    //- remove edges with open-end corners
    //- they are not captured by the mesher
    removeDanglingEdges();

    //- create feature edges
    forAll(featureEdges_, eI)
    {
        if( featureEdges_[eI] )
            featureEdges.append(edges[eI]);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
