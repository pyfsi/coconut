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
    polyMeshGenFaceZones

Description
    This class implements zonal information for faces. A face cannot be
    in more than one zone.

SourceFiles
    polyMeshGenFaceZones.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenFaceZones_H
#define polyMeshGenFaceZones_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceListPMG.H"
#include "labelLongList.H"
#include "DynList.H"
#include "boolList.H"

#include <map>

namespace Foam
{

class VRWGraph;
class Time;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenFaceZones
{
    // Private data
        //- const reference to faces
        const faceListPMG& faces_;

        //- zone index for every face
        //- -1 for faces not belonging to any zone
        labelLongList faceInZone_;

        //- is a face is zone flipped or not
        boolList flipFace_;

        //- name to index map
        std::map<word, label> nameToIndex_;

        //- index to name map
        std::map<label, word> indexToName_;

        //- counter of zone indices
        label index_;

    // Private member functions
        //- update the size of faceInZone list
        void updateSize();

    // Disallow bitwise assignment
        void operator=(const polyMeshGenFaceZones&);

        polyMeshGenFaceZones(const polyMeshGenFaceZones&);

public:

    // Constructors
        //- Construct from faceListPMG
        polyMeshGenFaceZones(const faceListPMG&);

    // Destructor
        ~polyMeshGenFaceZones();

    // Member functions
        //- add face zone
        inline label addFaceZone(const word& name);

        //- remove face zone
        inline label removeFaceZone(const word& name);

        //- available face zones
        inline void faceZoneIndices(DynList<label>&) const;
        inline void faceZoneNames(DynList<word>&) const;

        //- add face to zone
        inline void addFaceToZone(const label zoneId, const label faceI);

        //- remove face from zone
        inline void removeFaceFromZone(const label faceI);

        //- return the index of a face zone with a given name
        //- returns -1 if the zones does not exist
        inline label faceZoneIndex(const word&) const;

        //- return the name of a face zone with a given index
        inline word faceZoneName(const label) const;

        //- check if the face is in a zone
        inline bool isFaceInZone(const word& name, const label faceI) const;
        inline bool isFaceInZone(const label zoneId, const label faceI) const;

        //- return a zone index of a face
        inline label faceZone(const label faceI) const;

        //- return whether a face is flipped or not
        inline bool isFaceFlipped(const label faceI) const;

        //- set flip for face
        inline void setFlip(const label faceI, const bool);

        //- update face zones from a graph-type information
        void updateFaceZones(const VRWGraph&);

        //- update face zones from a list-type information
        template<class ListType>
        void updateFaceZones(const ListType&);

    // Read zones
        void readFaceZones
        (
            const Time&,
            const fileName& instance,
            const fileName& meshDir
        );

    // Write zones
        void writeFaceZones
        (
            const Time&,
            const fileName& instance,
            const fileName& meshDir
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenFaceZonesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
