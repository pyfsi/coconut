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

\*---------------------------------------------------------------------------*/

#include "taperedBoxRefinement.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(taperedBoxRefinement, 0);
addToRunTimeSelectionTable(objectRefinement, taperedBoxRefinement, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void taperedBoxRefinement::calculateHelperData()
{
    sinXY_ = Foam::sin(angleXY_ * M_PI / 180.0);
    sinYZ_ = Foam::sin(angleXZ_ * M_PI / 180.0);
    sinYX_ = Foam::sin(angleYX_ * M_PI / 180.0);
    sinYZ_ = Foam::sin(angleYZ_ * M_PI / 180.0);
    sinZX_ = Foam::sin(angleZX_ * M_PI / 180.0);
    sinZY_ = Foam::sin(angleZY_ * M_PI / 180.0);

    fCentres_[0] = fCentres_[1] = centre_.x();
    fCentres_[0] -= 0.5 * lengthX_;
    fCentres_[1] += 0.5 * lengthX_;

    fCentres_[2] = fCentres_[3] = centre_.y();
    fCentres_[2] -= 0.5 * lengthY_;
    fCentres_[3] += 0.5 * lengthY_;

    fCentres_[4] = fCentres_[5] = centre_.z();
    fCentres_[4] -= 0.5 * lengthZ_;
    fCentres_[5] += 0.5 * lengthZ_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

taperedBoxRefinement::taperedBoxRefinement()
:
    objectRefinement(),
    centre_(),
    lengthX_(-1.0),
    lengthY_(-1.0),
    lengthZ_(-1.0),
    angleXY_(0.0),
    angleXZ_(0.0),
    angleYX_(0.0),
    angleYZ_(0.0),
    angleZX_(0.0),
    angleZY_(0.0),
    sinXY_(0.0),
    sinXZ_(0.0),
    sinYX_(0.0),
    sinYZ_(0.0),
    sinZX_(0.0),
    sinZY_(0.0),
    fCentres_()
{}

taperedBoxRefinement::taperedBoxRefinement
(
    const word& name,
    const scalar cellSize,
    const direction additionalRefLevels,
    const point& centre,
    const scalar lengthX,
    const scalar lengthY,
    const scalar lengthZ,
    const scalar angleXY,
    const scalar angleXZ,
    const scalar angleYX,
    const scalar angleYZ,
    const scalar angleZX,
    const scalar angleZY
)
:
    objectRefinement(),
    centre_(centre),
    lengthX_(lengthX),
    lengthY_(lengthY),
    lengthZ_(lengthZ),
    angleXY_(angleXY),
    angleXZ_(angleXZ),
    angleYX_(angleYX),
    angleYZ_(angleYZ),
    angleZX_(angleZX),
    angleZY_(angleZY),
    sinXY_(0.0),
    sinXZ_(0.0),
    sinYX_(0.0),
    sinYZ_(0.0),
    sinZX_(0.0),
    sinZY_(0.0),
    fCentres_()
{
    calculateHelperData();

    setName(name);
    setCellSize(cellSize);
    setAdditionalRefinementLevels(additionalRefLevels);
}

taperedBoxRefinement::taperedBoxRefinement
(
    const word& name,
    const dictionary& dict
)
:
    objectRefinement(name, dict)
{
    this->operator=(dict);

    calculateHelperData();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool taperedBoxRefinement::intersectsObject(const boundBox& bb) const
{
    const point c = bb.midpoint();

    const vector d = c - centre_;

    //- x-min direction
    if( c.x() < (fCentres_[0] - d.y() * sinXY_ - d.z() * sinXZ_) )
        return false;

    //- x-max direction
    if( c.x() > (fCentres_[1] + d.y() * sinXY_ + d.z() * sinXZ_) )
        return false;

    //- y-min direction
    if( c.y() < (fCentres_[2] - d.x() * sinYX_ - d.z() * sinYZ_) )
        return false;

    //- y-max direction
    if( c.y() > (fCentres_[3] + d.x() * sinYX_ + d.z() * sinYZ_) )
        return false;

    //- z-min direction
    if( c.z() < (fCentres_[4] - d.x() * sinZX_ - d.y() * sinZY_) )
        return false;

    //- z-max direction
    if( c.z() > (fCentres_[5] + d.x() * sinZX_ + d.y() * sinZY_) )
        return false;

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dictionary taperedBoxRefinement::dict(bool /*ignoreType*/) const
{
    dictionary dict;

    if( additionalRefinementLevels() == 0 && cellSize() >= 0.0 )
    {
        dict.add("cellSize", cellSize());
    }
    else
    {
        dict.add("additionalRefinementLevels", additionalRefinementLevels());
    }

    dict.add("type", type());

    dict.add("centre", centre_);
    dict.add("lengthX", lengthX_);
    dict.add("lengthY", lengthY_);
    dict.add("lengthZ", lengthZ_);
    dict.add("angleXY", angleXY_);
    dict.add("angleXZ", angleXZ_);
    dict.add("angleYX", angleYX_);
    dict.add("angleYZ", angleYZ_);
    dict.add("angleZX", angleZX_);
    dict.add("angleZY", angleZY_);

    return dict;
}

void taperedBoxRefinement::write(Ostream& os) const
{
    os  << " type:   " << type()
        << " centre: " << centre_
        << " lengthX: " << lengthX_
        << " lengthY: " << lengthY_
        << " lengthZ: " << lengthZ_
        << " angleXY: " << angleXY_
        << " angleXZ: " << angleXZ_
        << " angleYX: " << angleYX_
        << " angleYZ: " << angleYZ_
        << " angleZX: " << angleZX_
        << " angleZY: " << angleZY_;
}

void taperedBoxRefinement::writeDict(Ostream& os, bool subDict) const
{
    if( subDict )
    {
        os << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    if( additionalRefinementLevels() == 0 && cellSize() >= 0.0 )
    {
        os.writeKeyword("cellSize") << cellSize() << token::END_STATEMENT << nl;
    }
    else
    {
        os.writeKeyword("additionalRefinementLevels")
                << additionalRefinementLevels()
                << token::END_STATEMENT << nl;
    }

    // only write type for derived types
    if( type() != typeName_() )
    {
        os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    }

    os.writeKeyword("centre") << centre_ << token::END_STATEMENT << nl;
    os.writeKeyword("lengthX") << lengthX_ << token::END_STATEMENT << nl;
    os.writeKeyword("lengthY") << lengthY_ << token::END_STATEMENT << nl;
    os.writeKeyword("lengthZ") << lengthZ_ << token::END_STATEMENT << nl;
    os.writeKeyword("angleXY") << angleXY_ << token::END_STATEMENT << nl;
    os.writeKeyword("angleXZ") << angleXZ_ << token::END_STATEMENT << nl;
    os.writeKeyword("angleYX") << angleYX_ << token::END_STATEMENT << nl;
    os.writeKeyword("angleYZ") << angleYZ_ << token::END_STATEMENT << nl;
    os.writeKeyword("angleZX") << angleZX_ << token::END_STATEMENT << nl;
    os.writeKeyword("angleZY") << angleZY_ << token::END_STATEMENT << nl;

    if( subDict )
    {
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
}

void taperedBoxRefinement::operator=(const dictionary& d)
{
    // allow as embedded sub-dictionary "coordinateSystem"
    const dictionary& dict =
    (
        d.found(typeName_())
      ? d.subDict(typeName_())
      : d
    );

    // unspecified centre is (0 0 0)
    if( dict.found("centre") )
    {
        dict.lookup("centre") >> centre_;
    }
    else
    {
        FatalErrorIn
        (
            "void taperedBoxRefinement::operator=(const dictionary& d)"
        ) << "Entry centre is not specified!" << exit(FatalError);
        centre_ = vector::zero;
    }

    // specify lengthX
    if( dict.found("lengthX") )
    {
        lengthX_ = readScalar(dict.lookup("lengthX"));
    }
    else
    {
        FatalErrorIn
        (
            "void taperedBoxRefinement::operator=(const dictionary& d)"
        ) << "Entry lengthX is not specified!" << exit(FatalError);
        lengthX_ = -1.0;
    }

    // specify lengthY
    if( dict.found("lengthY") )
    {
        lengthY_ = readScalar(dict.lookup("lengthY"));
    }
    else
    {
        FatalErrorIn
        (
            "void taperedBoxRefinement::operator=(const dictionary& d)"
        ) << "Entry lengthY is not specified!" << exit(FatalError);
        lengthY_ = -1.0;
    }

    // specify lengthZ
    if( dict.found("lengthZ") )
    {
        lengthZ_ = readScalar(dict.lookup("lengthZ"));
    }
    else
    {
        FatalErrorIn
        (
            "void taperedBoxRefinement::operator=(const dictionary& d)"
        ) << "Entry lengthZ is not specified!" << exit(FatalError);
        lengthZ_ = -1.0;
    }

    // specify angleXY
    if( dict.found("angleXY") )
    {
        angleXY_ = readScalar(dict.lookup("angleXY"));
    }
    else
    {
        FatalErrorIn
        (
            "void taperedBoxRefinement::operator=(const dictionary& d)"
        ) << "Entry angleXY is not specified!" << exit(FatalError);
        angleXY_ = 0.0;
    }

    // specify angleXZ
    if( dict.found("angleXZ") )
    {
        angleXZ_ = readScalar(dict.lookup("angleXZ"));
    }
    else
    {
        FatalErrorIn
        (
            "void taperedBoxRefinement::operator=(const dictionary& d)"
        ) << "Entry angleXZ is not specified!" << exit(FatalError);
        angleXZ_ = 0.0;
    }

    // specify angleYX
    if( dict.found("angleYX") )
    {
        angleYX_ = readScalar(dict.lookup("angleYX"));
    }
    else
    {
        FatalErrorIn
        (
            "void taperedBoxRefinement::operator=(const dictionary& d)"
        ) << "Entry angleYX is not specified!" << exit(FatalError);
        angleYX_ = 0.0;
    }

    // specify angleYZ
    if( dict.found("angleYZ") )
    {
        angleYZ_ = readScalar(dict.lookup("angleYZ"));
    }
    else
    {
        FatalErrorIn
        (
            "void taperedBoxRefinement::operator=(const dictionary& d)"
        ) << "Entry angleYZ is not specified!" << exit(FatalError);
        angleYZ_ = 0.0;
    }

    // specify angleZX
    if( dict.found("angleZX") )
    {
        angleZX_ = readScalar(dict.lookup("angleZX"));
    }
    else
    {
        FatalErrorIn
        (
            "void taperedBoxRefinement::operator=(const dictionary& d)"
        ) << "Entry angleZX is not specified!" << exit(FatalError);
        angleZX_ = 0.0;
    }

    // specify angleZY
    if( dict.found("angleZY") )
    {
        angleZY_ = readScalar(dict.lookup("angleZY"));
    }
    else
    {
        FatalErrorIn
        (
            "void taperedBoxRefinement::operator=(const dictionary& d)"
        ) << "Entry angleZY is not specified!" << exit(FatalError);
        angleZY_ = 0.0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& taperedBoxRefinement::operator<<(Ostream& os) const
{
    os << "name " << name() << nl;
    os << "cell size " << cellSize() << nl;
    os << "additionalRefinementLevels " << additionalRefinementLevels() << endl;

    write(os);
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
