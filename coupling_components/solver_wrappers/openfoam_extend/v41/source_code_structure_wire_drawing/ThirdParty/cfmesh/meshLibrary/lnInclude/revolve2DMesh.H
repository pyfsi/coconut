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
    revolve2DMesh

Description
    Revolves a 2D mesh around pre-defined axis.

SourceFiles
    revolve2DMesh.C

\*---------------------------------------------------------------------------*/

#ifndef revolve2DMesh_H
#define revolve2DMesh_H

#include "polyMeshGenModifier.H"
#include "FixedList.H"

#include <map>

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class revolve2DMesh Declaration
\*---------------------------------------------------------------------------*/

class revolve2DMesh
{
    // Private data
        //- reference to the mesh
        polyMeshGen& mesh_;

        //- name of a patch that shall be revolved into a volume mesh
        word revolvingPatch_;

        //- rotation axis
        vector rotationAxis_;

        //- origin of the rotation axis
        point origin_;

        //- maximum expansion ratio in the circumferential direction
        scalar maxExpansionRatio_;

        //- angle intervals and the applicable resolution over the interval
        typedef std::pair<scalar, scalar> ssPair;
        typedef std::map<ssPair, label> angleMap;

        angleMap angleIntervals_;

        //- angle interval and patches
        typedef std::map<label, label> llMap;
        typedef std::map<ssPair, llMap> patchesMap;
        patchesMap rplPatches_;

        //- start angle of a segment. angleIntervals are relative to it
        scalar startAngle_;

        //- perform a complete 360 degrees revolution or not
        bool completeRevolution_;

        //- create subsets of mesh points in each revolution plane
        bool generatePointSubsets_;

        //- does it generate a wedge
        bool isWedge_;

    // Private member functions
        //- check whether the mesh can be extruded or not
        bool checkSetup() const;

        //- analyse the angle setup
        void analyseRevolvingAngles();

        //- generate a revolved 3D mesh
        void revolveMesh();

        //- disallow bitwise assignment
        void operator=(const revolve2DMesh&);

        //- disallow copy construct
        revolve2DMesh(const revolve2DMesh&);

public:

    // Constructors
        //- Construct from polyMeshGen
        revolve2DMesh(polyMeshGen&);

    // Destructor
        ~revolve2DMesh();

    // Public member functions
        //- provide the desired settings
        //- set rotation axis
        void setRotationAxis(const vector&);

        //- set the name of the patch that shall be revolved into a volume mesh
        void setRevolvingPatch(const word&);

        //- set origin point at the rotation axis
        void setOrigin(const point&);

        //- set maximum expansion ratio in the circumferential direction
        void setMaxRatio(const scalar);

        //- clear angle intervals
        void clearAngleIntervals();

        //- set the desired resolution in the circumeferential direction
        //- all angles are provided in radians
        void setCircumResolution(const label nDivisions);
        void setCircumResolution(const scalar resolutionAngle);

        //- set the desired resolution over a given interval
        //- all angles are provided in radians
        void setIntervalResolution
        (
            const scalar startAngle,
            const scalar endAngle,
            const label nDivisions
        );
        void setIntervalResolution
        (
            const scalar startAngle,
            const scalar endAngle,
            const scalar resolutionAngle
        );

        void setIntervalPatches
        (
            const scalar startAngle,
            const scalar endAngle,
            const word origPatchName,
            const word newPatchName,
            const word newPatchType = "patch"
        );

        //- shall the revolution create a wedge boundary conditions
        void setWedgeBoundaryConditions();

        //- create points subsets for each circumferential plane
        void createPointSubsets();

        //- perform generation of a 3D mesh
        void generateRevolvedMesh();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
