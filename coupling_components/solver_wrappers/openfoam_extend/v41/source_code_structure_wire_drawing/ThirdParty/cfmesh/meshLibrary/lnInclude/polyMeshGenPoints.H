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
    polyMeshGenPoints

Description
    Mesh with geometry information

SourceFiles
    polyMeshGenPoints.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenPoints_H
#define polyMeshGenPoints_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenPointZones.H"
#include "meshSubset.H"
#include "pointFieldPMG.H"
#include "DynList.H"

# ifdef ExtendSpecific
#include "foamTime.H"
# else
#include "Time.H"
# endif

#include <map>

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenPoints
: public polyMeshGenPointZones
{
protected:

    // Protected data
        //- reference to the Time registry
        const Time& runTime_;

        //- mesh instance
        fileName instance_;

        //- directory with mesh files
        fileName meshDir_;

        //- list of vertices
        pointFieldPMG points_;

        //- backup of vertices
        std::map<label, point> origPoints_;

        //- store information about locked points
        std::set<label> lockedPoints_;

        //- map of point subsets
        std::map<label, meshSubset> pointSubsets_;

    // Disallow bitwise assignment
        void operator=(const polyMeshGenPoints&);

        polyMeshGenPoints(const polyMeshGenPoints&);

public:

    friend class polyMeshGenModifier;

    // Constructors
        //- Construct for Time
        polyMeshGenPoints
        (
            const Time&,
            const fileName instance = "constant",
            const fileName meshDir = "polyMesh"
        );

        //- Construct from components without the boundary
        polyMeshGenPoints
        (
            const Time&,
            const pointField& points,
            const fileName instance = "constant",
            const fileName meshDir = "polyMesh"
        );

    // Destructor
        ~polyMeshGenPoints();

    // Member functions
        //- access to Time
        inline const Time& returnTime() const;

        //- access to instance
        inline const fileName& instance() const;

        //- non-const access to instance
        inline fileName& instance();

        //- access to meshDir
        inline const fileName& meshDir() const;

        //- non-const access to meshDir
        inline fileName& meshDir();

        //- access to points
        inline const pointFieldPMG& points() const;

        //- append a vertex to the end of the list
        inline void appendVertex(const point& p);

        //- have the points been backed up or not
        inline bool hasPointsBackup() const;

        //- check if the given points has its backup
        inline bool hasBackup(const label) const;

        //- get coordinates of a point at the time when backup was performed
        //- in case the backup does not exist it returns the current coordinates
        inline bool getOrigPoint(const label pointI, point& p) const;

        //- return the set of locked points
        inline const std::set<label>& lockedPoints() const;

        //- has the mesh any locked points or not
        inline bool hasLockedPoints() const;

        //- check if the point is locked
        inline bool isLockedPoint(const label pointI) const;

        //- point subsets
        label addPointSubset(const word&);
        void removePointSubset(const label);
        word pointSubsetName(const label) const;
        label pointSubsetIndex(const word&) const;
        inline void addPointToSubset(const label, const label);
        inline void removePointFromSubset(const label, const label);
        inline void pointInSubsets(const label, DynList<label>&) const;
        inline void pointSubsetIndices(DynList<label>&) const;
        template<class ListType>
        inline void pointsInSubset(const label, ListType&) const;
        template<class ListType>
        inline void updatePointSubsets(const ListType&);

    // Read mesh
        void read();

        void readFromLatestTime();

    // Write mesh
        void write() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenPointsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
