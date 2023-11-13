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

Description

\*---------------------------------------------------------------------------*/

#include "zoneIO.H"
#include "labelList.H"
#include "boolList.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

zoneIO::zoneIO()
:
    name_(),
    dict_()
{}

zoneIO::zoneIO(const word& name, const dictionary& dict)
:
    name_(name),
    dict_(dict)
{}

zoneIO::zoneIO(Istream& is)
:
    name_(is),
    dict_(is)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

zoneIO::~zoneIO()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const word& zoneIO::name() const
{
    return name_;
}

word& zoneIO::name()
{
    return name_;
}

void zoneIO::setName(const word& name)
{
    name_ = name;
}

const dictionary& zoneIO::dict() const
{
    return dict_;
}

dictionary& zoneIO::dict()
{
    return dict_;
}

void zoneIO::setDictionary(const dictionary& d)
{
    dict_ = d;
}

bool zoneIO::operator!=(const zoneIO& zIO) const
{
    return (name_ != zIO.name_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const zoneIO& zIO)
{
    os.format(IOstream::ASCII);

    os << zIO.name() << nl << token::BEGIN_BLOCK << nl;

    const word zoneType(zIO.dict().lookup("type"));

    if( zoneType == "cellZone" )
    {
        const labelList cellsInZone(zIO.dict().lookup("cellLabels"));

        os << "type " << "cellZone;" << nl;

        cellsInZone.writeEntry("cellLabels", os);
    }
    else if( zoneType == "faceZone" )
    {
        const labelList facesInZone(zIO.dict().lookup("faceLabels"));

        os << "type " << "faceZone;" << nl;

        facesInZone.writeEntry("faceLabels", os);

        const boolList flipMap(zIO.dict().lookup("flipMap"));

        flipMap.writeEntry("flipMap", os);
    }
    else if( zoneType == "pointZone" )
    {
        const labelList pointsInZone(zIO.dict().lookup("pointLabels"));

        os << "type " << "pointZone;" << nl;

        pointsInZone.writeEntry("pointLabels", os);
    }

    os << nl << token::END_BLOCK << nl;

    return os;
}

Istream& operator>>(Istream& is, zoneIO& zIO)
{
    is >> zIO.name();
    is >> zIO.dict();

    return is;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
