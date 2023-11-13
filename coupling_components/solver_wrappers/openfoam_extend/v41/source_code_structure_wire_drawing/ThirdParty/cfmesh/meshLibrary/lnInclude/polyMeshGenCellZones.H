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
    polyMeshGenCellZones

Description
    This class implements zonal information for cells. A cell cannot be
    in more than one zone.

SourceFiles
    polyMeshGenCellZones.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenCellZones_H
#define polyMeshGenCellZones_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellListPMG.H"
#include "labelLongList.H"
#include "DynList.H"

#include "IOList.H"
#include "Istream.H"

#include <map>

namespace Foam
{

class VRWGraph;
class Time;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenCellZones
{
    // Private data
        //- const reference to cells
        const cellListPMG& cells_;

        //- zone index for every cell
        //- -1 for cells not belonging to any zone
        labelLongList cellInZone_;

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
        void operator=(const polyMeshGenCellZones&);

        polyMeshGenCellZones(const polyMeshGenCellZones&);

public:

    // Constructors
        //- Construct from cellListPMG
        polyMeshGenCellZones(const cellListPMG&);

    // Destructor
        ~polyMeshGenCellZones();

    // Member functions
        //- add cell zone
        inline label addCellZone(const word& name);

        //- remove cell zone
        inline label removeCellZone(const word& name);

        //- available cell zones
        inline void cellZoneIndices(DynList<label>&) const;
        inline void cellZoneNames(DynList<word>&) const;

        //- add cell to zone
        inline void addCellToZone(const label zoneId, const label cellI);

        //- remove cell from zone
        inline void removeCellFromZone(const label cellI);

        //- return the index of a cell zone with a given name
        //- returns -1 if the zones does not exist
        inline label cellZoneIndex(const word&) const;

        //- return the name of a cell zone with a given index
        inline word cellZoneName(const label) const;

        //- check if the cell is in a zone
        inline bool isCellInZone(const word& name, const label cellI) const;
        inline bool isCellInZone(const label zoneId, const label cellI) const;

        //- return a zone index of a cell
        inline label cellZone(const label cellI) const;

        //- update cell zones from a graph-type information
        void updateCellZones(const VRWGraph&);

        //- update cell zones from a list-type information
        template<class ListType>
        void updateCellZones(const ListType&);

    // Read zones
        void readCellZones
        (
            const Time&,
            const fileName& instance,
            const fileName& meshDir
        );

    // Write zones
        void writeCellZones
        (
            const Time&,
            const fileName& instance,
            const fileName& meshDir
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenCellZonesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
