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

Class
    collapseBoundaryLayers

Description
    Refine existing boundary layers

SourceFiles
    collapseBoundaryLayers.C
    collapseBoundaryLayersFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef collapseBoundaryLayers_H
#define collapseBoundaryLayers_H

#include "polyMeshGenModifier.H"
#include "DynList.H"
#include "labelLongList.H"
#include "labelPair.H"

#include <map>
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class refineBoundaryLayers Declaration
\*---------------------------------------------------------------------------*/

class collapseBoundaryLayers
{
    // Private data
        //- Reference to the mesh
        polyMeshGen& mesh_;

        //- pointer to mesh surface engine
        mutable meshSurfaceEngine* msePtr_;

    // Private member functions
        //- Return reference to meshSurfaceEngine
        const meshSurfaceEngine& surfaceEngine() const;

        //- detect low quality hairs and modify faces
        void detectLowQualityHairsAndCollapseFaces();

        //- modify cells in the layer
        void modifyLayerCells();

        //- Disallow bitwise copy construct
        collapseBoundaryLayers(const collapseBoundaryLayers&);

        //- Disallow bitwise assignment
        void operator=(const collapseBoundaryLayers&);

public:

    // Constructors

        //- Construct from mesh reference
        collapseBoundaryLayers(polyMeshGen& mesh);

    // Destructor
        ~collapseBoundaryLayers();

    // Public member functions


    // Worker
        //- collapse layers based on the given settings
        void collapseLayers();

    // Static member functions
        //- read the settings from dictionary
        static void readSettings(const dictionary&, collapseBoundaryLayers&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
