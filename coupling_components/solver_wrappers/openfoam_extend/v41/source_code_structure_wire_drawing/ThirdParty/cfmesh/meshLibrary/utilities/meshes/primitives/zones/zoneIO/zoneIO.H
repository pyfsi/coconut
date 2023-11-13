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

#ifndef zoneIO_H
#define zoneIO_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Istream.H"
#include "Ostream.H"
#include "IOdictionary.H"


namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class zoneIO
{
    // Private data
        //- name of the zone
        word name_;

        //- dictionary containing data
        dictionary dict_;

 public:

    // Constructors
        //- Null construct
        zoneIO();

        //- Construct from components
        zoneIO(const word& name, const dictionary& dict);

        //- Construct from Istream
        zoneIO(Istream&);

    // Destructor
        ~zoneIO();

    // Public member functions
        //- return zone name
        const word& name() const;
        word& name();

        //- set the name
        void setName(const word& name);

        //- return the dictionary
        const dictionary& dict() const;
        dictionary& dict();

        //- set the dictionary
        void setDictionary(const dictionary& d);

    // Member operators
        bool operator!=(const zoneIO& zIO) const;

    // Friend operators
        friend Ostream& operator<<(Ostream&, const zoneIO&);
        friend Istream& operator>>(Istream&, zoneIO&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
