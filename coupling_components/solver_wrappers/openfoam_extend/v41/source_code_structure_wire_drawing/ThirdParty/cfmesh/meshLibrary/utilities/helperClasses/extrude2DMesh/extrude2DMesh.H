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
    extrude2DMesh

Description
    Revolves a 2D mesh around pre-defined axis.

SourceFiles
    extrude2DMesh.C

\*---------------------------------------------------------------------------*/

#ifndef extrude2DMesh_H
#define extrude2DMesh_H

#include "polyMeshGenModifier.H"
#include "FixedList.H"

#include <map>

namespace Foam
{
// Forward declarations
class edgeMesh;

/*---------------------------------------------------------------------------*\
                        Class extrude2DMesh Declaration
\*---------------------------------------------------------------------------*/

class extrude2DMesh
{
    // Private data
        //- reference to the mesh
        polyMeshGen& mesh_;

        //- name of a patch that shall be extruded
        word extrusionPatch_;

        //- name of a subset containing created cells
        word cellSubsetName_;

        //- extrusion axis
        vector extrusionVector_;

        //- extrusion distance
        scalar extrusionDistance_;

        //- grading factor
        scalar gradingFactor_;

        //- number of subdivisions
        label nSubdivisions_;

        //- extrude in the normal direction of a patch
        bool extrudeInPatchNormalDirection_;

        //- remove remaining original cells from the mesh
        bool removeOriginalCells_;

    // Private member functions
        //- check whether the mesh can be extruded or not
        bool checkSetup() const;

        //- generate a extruded 3D mesh
        void extrudeMesh();

        //- disallow bitwise assignment
        void operator=(const extrude2DMesh&);

        //- disallow copy construct
        extrude2DMesh(const extrude2DMesh&);

public:

    // Constructors
        //- Construct from polyMeshGen
        extrude2DMesh(polyMeshGen&);

    // Destructor
        ~extrude2DMesh();

    // Public member functions
        //- provide the desired settings
        //- set translation vector
        void setExtrusionVector(const vector);

        //- length of extrusion path
        void extrusionLength(const scalar distance);

        //- set the patch that shall be extruded
        void setExtrusionPatch(const word&);

        //- set the number of subdivions
        void setNumberOfSubdivisions(const label);

        //- set the grading factor
        void setGradingFactor(const scalar);

        //- subset name for created cells
        void setCellSubsetName(const word&);

        //- shall it extrude in the normal direction of a patch
        //- extrusion vector does not need to be specified in this case
        void extrudeInPatchNormalDirection(const bool);

        //- shall the cells that were in the mesh before be deleted from the mesh
        //- the default is yes
        void removeOriginalCells(const bool);

        //- perform generation of a 3D mesh
        void generateExtrudedMesh();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
