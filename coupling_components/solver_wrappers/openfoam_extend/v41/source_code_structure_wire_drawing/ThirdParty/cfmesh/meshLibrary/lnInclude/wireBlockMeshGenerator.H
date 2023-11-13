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
    wireBlockMeshGenerator

Description

    Generates a block mesh in the domain of circular shape

SourceFiles
    wireBlockMeshGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef wireBlockMeshGenerator_H
#define wireBlockMeshGenerator_H

#include "polyMeshGenModifier.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMeshGen;
class wireGeometryInfo;
class rollingMillPatchNamesHandler;

/*---------------------------------------------------------------------------*\
                        Class blockWriter Declaration
\*---------------------------------------------------------------------------*/

class blockWriter
{
    // Private data
        //- type of shape
        word type_;

        //- shape points
        DynList<label> nodes_;

        //- resolution in all directions
        Vector<label> resolution_;

        //- grading
        vector grading_;

public:

    // Constructors
        //- Null construct
        blockWriter();

        //- Construct from components
        blockWriter
        (
            const word& shapeName,
            const DynList<label>& nodes,
            const Vector<label>& resolution,
            const vector& grading
        );

    // Destructor

        ~blockWriter();

    // Public member functions

    // Public operators

        void operator=(const blockWriter&);
        bool operator!=(const blockWriter&) const;

    // Friend functions
        friend Ostream& operator<<(Ostream&, const blockWriter&);
};

/*---------------------------------------------------------------------------*\
                        Class arcEdgeWriter Declaration
\*---------------------------------------------------------------------------*/

class arcEdgeWriter
{
    // Private data
        //- type of edge
        word type_;

        //- start point
        label start_;

        //- end point
        label end_;

        //- grading
        point midPoint_;

public:

    // Constructors
        //- Null construct
        arcEdgeWriter();

        //- Construct from components
        arcEdgeWriter
        (
            const word& shapeName,
            const label start,
            const label end,
            const point& midPoint
        );

    // Destructor

        ~arcEdgeWriter();

    // Public member functions

    // Public operators

        void operator=(const arcEdgeWriter&);
        bool operator!=(const arcEdgeWriter&) const;

    // Friend functions
        friend Ostream& operator<<(Ostream&, const arcEdgeWriter&);
};

/*---------------------------------------------------------------------------*\
                Class namedDictionary Declaration
\*---------------------------------------------------------------------------*/

class namedDictionary
{
    // Private data
        //- name of the dictionary
        word name_;

        //- dictionary entries
        dictionary dict_;

public:

    // Constructors
        //- Null construct
        namedDictionary();

        //- Construct from components
        namedDictionary
        (
            const word& name,
            const dictionary& dict
        );

    // Destructor

        ~namedDictionary();

    // Public member functions

    // Public operators

        void operator=(const namedDictionary&);
        bool operator!=(const namedDictionary&) const;

    // Friend functions
        friend Ostream& operator<<(Ostream&, const namedDictionary&);
};

/*---------------------------------------------------------------------------*\
                Class wireBlockMeshGenerator Declaration
\*---------------------------------------------------------------------------*/

class wireBlockMeshGenerator
{
    // Private data
        //- reference to runTime
        const Time& runTime_;

        //- const reference to patch names handler
        const rollingMillPatchNamesHandler& patchHandler_;

        //- radius of the wire
        const scalar radius_;

        //- region name
        const word regionName_;

        //- symmetry type
        const direction symmetryType_;

        //- geometric tolerance
        const scalar geometricTolerance_;

        //- const reference to the wire geometry
        const wireGeometryInfo& wireGeom_;

        //- block mesh dictionary
        IOdictionary blockMeshDict_;

    // Private member functions
        //- create a mesh dictionary for the axi-symmetric case
        void createAxiSymmetricMeshDictionary();

        //- assemble the dictionary for the block mesher
        void createMeshDictionary();

        //- Disallow default bitwise copy construct
        wireBlockMeshGenerator(const wireBlockMeshGenerator&);

        //- Disallow default bitwise assignment
        void operator=(const wireBlockMeshGenerator&);

public:

    // Constructors

        //- Construct from Time, diameter, symmetry type and tolerance
        wireBlockMeshGenerator
        (
            const Time& runTime,
            const rollingMillPatchNamesHandler& patchHandler,
            const wireGeometryInfo& wireGeom,
            const word regionName,
            const scalar tol=1e-6
        );


    // Destructor

        ~wireBlockMeshGenerator();


    // Member Functions

        //- generates a 2D mesh of a circular cross-section
        polyMeshGen* generateCrossSectionMesh();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
