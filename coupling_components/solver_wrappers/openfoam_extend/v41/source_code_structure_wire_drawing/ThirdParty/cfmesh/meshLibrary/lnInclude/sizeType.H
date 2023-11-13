/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Primitive
    long long

Description
    A long long (64 bits on many systems)

SourceFiles
    longLongIO.C

\*---------------------------------------------------------------------------*/

#ifndef sizeType_H
#define sizeType_H

#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef long long sizeType;

# define FOAM_SIZETYPE_MAX LLONG_MAX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return a string representation of a sizeType
word name(const sizeType);

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

sizeType readSizeType(Istream&);
bool readSizeType(const char*, sizeType&);

Istream& operator>>(Istream&, sizeType&);
Ostream& operator<<(Ostream&, const sizeType);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# ifdef NoSizeType
#define MAXMIN(retType, type1, type2)              \
                                                   \
inline retType max(const type1 s1, const type2 s2) \
{                                                  \
    return (s1 > s2)? s1: s2;                      \
}                                                  \
                                                   \
inline retType min(const type1 s1, const type2 s2) \
{                                                  \
    return (s1 < s2)? s1: s2;                      \
}

MAXMIN(sizeType, sizeType, sizeType)
MAXMIN(sizeType, label, sizeType)
MAXMIN(sizeType, sizeType, label)

# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
