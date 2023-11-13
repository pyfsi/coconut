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

#include "processorBoundaryPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "Ostream.H"
#include "Istream.H"
#include "token.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorBoundaryPatch, 0);
addToRunTimeSelectionTable(boundaryPatchBase, processorBoundaryPatch, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

processorBoundaryPatch::processorBoundaryPatch
(
    const word& name,
    const word& type,
    const label nFaces,
    const label startFace,
    const label myProcNo,
    const label neighbProcNo
)
:
    boundaryPatchBase(name, type, nFaces, startFace),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo)
{}

processorBoundaryPatch::processorBoundaryPatch
(
    const word& name,
    const dictionary& dict
)
:
    boundaryPatchBase(name, dict),
    myProcNo_(readLabel(dict.lookup("myProcNo"))),
    neighbProcNo_(readLabel(dict.lookup("neighbProcNo")))
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const dictionary& processorBoundaryPatch::dict() const
{
    dict_.add("type", type_, true);

    dict_.add("nFaces", nFaces_, true);
    dict_.add("startFace", startFace_, true);
    dict_.add("myProcNo", myProcNo_, true);
    dict_.add("neighbProcNo", neighbProcNo_, true);

    return dict_;
}

dictionary& processorBoundaryPatch::dict()
{
    dict_.add("type", type_, true);

    dict_.add("nFaces", nFaces_, true);
    dict_.add("startFace", startFace_, true);
    dict_.add("myProcNo", myProcNo_, true);
    dict_.add("neighbProcNo", neighbProcNo_, true);

    return dict_;
}

void processorBoundaryPatch::write(Ostream& os) const
{
    this->operator<<(os);
}

void processorBoundaryPatch::writeDict(Ostream& /*os*/) const
{

}

Ostream& processorBoundaryPatch::operator<<(Ostream& os) const
{
    os << patchName() << nl << dict() << nl;

    return os;
}

Istream& processorBoundaryPatch::operator>>(Istream& is)
{
    is >> name_ >> dict_;

	type_ = word(dict_.lookup("type"));
	nFaces_ = readLabel(dict_.lookup("nFaces"));
	startFace_ = readLabel(dict_.lookup("startFace"));
	myProcNo_ = readLabel(dict_.lookup("myProcNo"));
    neighbProcNo_ = readLabel(dict_.lookup("neighbProcNo"));

    return is;
}

void processorBoundaryPatch::operator=(const processorBoundaryPatch& wp)
{
    name_ = wp.name_;
    type_ = wp.type_;
    nFaces_ = wp.nFaces_;
    startFace_ = wp.startFace_;
    myProcNo_ = wp.myProcNo_;
    neighbProcNo_ = wp.neighbProcNo_;
    dict_ = wp.dict_;
}

bool processorBoundaryPatch::operator!=(const processorBoundaryPatch& wp) const
{
    if( name_ != wp.name_ )
    {
        return true;
    }
    else if( type_ != wp.name_ )
    {
        return true;
    }
    else if( (nFaces_ != wp.nFaces_) || (startFace_ != wp.startFace_) )
    {
        return true;
    }
    else if(
        (myProcNo_ != wp.myProcNo_) || (neighbProcNo_ != wp.neighbProcNo_)
    )
    {
        return true;
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
