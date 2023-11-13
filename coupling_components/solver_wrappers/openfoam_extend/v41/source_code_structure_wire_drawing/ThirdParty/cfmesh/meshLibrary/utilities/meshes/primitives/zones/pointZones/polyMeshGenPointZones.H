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
    polyMeshGenPointZones

Description
    This class implements zonal information for points. A point cannot be
    in more than one zone.

SourceFiles
    polyMeshGenPointZones.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenPointZones_H
#define polyMeshGenPointZones_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointFieldPMG.H"
#include "labelLongList.H"
#include "DynList.H"

#include <map>

namespace Foam
{

class VRWGraph;
class Time;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenPointZones
{
    // Private data
        //- const reference to cells
        const pointFieldPMG& points_;

        //- zone index for every point
        //- -1 for points not belonging to any zone
        labelLongList pointInZone_;

        //- name to index map
        std::map<word, label> nameToIndex_;

        //- index to name map
        std::map<label, word> indexToName_;

        //- counter of zone indices
        label index_;

    // Private member functions
        //- update the size of the cellInZone list
        void updateSize();

    // Disallow bitwise assignment
        void operator=(const polyMeshGenPointZones&);

        polyMeshGenPointZones(const polyMeshGenPointZones&);

public:

    // Constructors
        //- Construct from pointFieldPMG
        polyMeshGenPointZones(const pointFieldPMG&);

    // Destructor
        ~polyMeshGenPointZones();

    // Member functions
        //- add point zone
        inline label addPointZone(const word& name);

        //- remove point zone
        inline label removePointZone(const word& name);

        //- available point zones
        inline void pointZoneIndices(DynList<label>&) const;
        inline void pointZoneNames(DynList<word>&) const;

        //- add point to zone
        inline void addPointToZone(const label zoneId, const label pointI);

        //- remove point from zone
        inline void removePointFromZone(const label pointI);

        //- return the index of a point zone with a given name
        //- returns -1 if the zones does not exist
        inline label pointZoneIndex(const word&) const;

        //- return the name of a point zone with a given index
        inline word pointZoneName(const label) const;

        //- check if the point is in a zone
        inline bool isPointInZone(const word& name, const label pointI) const;
        inline bool isPointInZone(const label zoneId, const label pointI) const;

        //- return a zone index of a point
        inline label pointZone(const label pointI) const;

        //- update point zones from a graph-type information
        void updatePointZones(const VRWGraph&);

        //- update point zones from a list-type information
        template<class ListType>
        void updatePointZones(const ListType&);

    // Read zones
        void readPointZones
        (
            const Time&,
            const fileName& instance,
            const fileName& meshDir
        );

    // Write zones
        void writePointZones
        (
            const Time&,
            const fileName& instance,
            const fileName& meshDir
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenPointZonesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
