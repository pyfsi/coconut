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
    meshSurfaceTopologyImprove

Description
    Checks

SourceFiles
    meshSurfaceTopologyImprove.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceTopologyImprove_H
#define meshSurfaceTopologyImprove_H

#include "polyMeshGenModifier.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshSurfaceEngine;
class meshOctree;

/*---------------------------------------------------------------------------*\
                Class meshSurfaceTopologyImprove Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceTopologyImprove
{
    // Private data
        //- Reference to polyMeshGen
        polyMeshGen& mesh_;

        //- const reference to meshOctree
        const meshOctree& octree_;

        //- pointer to meshSurfaceEngine
        mutable meshSurfaceEngine* meshSurfacePtr_;

        //- types of problems detected at this boundary face
        List<direction> problemAtFace_;

        //- new patch for a face
        labelLongList newFacePatch_;

    // Private member functions
        //- creates meshSurfaceEngine
        void createMeshSurface() const;

        //- returns mesh surface
        inline const meshSurfaceEngine& meshSurface() const
        {
            if( !meshSurfacePtr_ )
                createMeshSurface();

            return *meshSurfacePtr_;
        }

        //- deletes meshSurfaceEngine
        void clearMeshSurface();

        //- optimize the volume of the mesh
        void optimizeCells();

        //- detect faces at the surface of the mesh that cause
        //- a poor quality mesh
        label detectProblematicFaces();

        //- collapse or decompose deformed faces
        bool modifyDeformedFaces();

        //- move highly-rotated faces to new patches
        bool updateFacePatches();

        //- analyse new patches and ensure consistency
        bool analyseUpdatedPatches();

        //- Disallow default bitwise copy construct
        meshSurfaceTopologyImprove
        (
            const meshSurfaceTopologyImprove&
        );

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceTopologyImprove&);

    // Private enumerators
        enum faceProblemType
        {
            FACEOK = 0,
            BADVISIBILITY = 1,
            CONCAVEEDGE = 2,
            TWISTEDATEDGE = 4,
            TWISTEDATCORNER = 8
        };

public:

    // Constructors

        //- Construct from polyMeshGen and meshOctree
        meshSurfaceTopologyImprove(polyMeshGen& mesh, const meshOctree&);

    // Destructor

        ~meshSurfaceTopologyImprove();

    // Member Functions

        //- modify the mesh to better represent the geometry
        bool improveTopology();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
