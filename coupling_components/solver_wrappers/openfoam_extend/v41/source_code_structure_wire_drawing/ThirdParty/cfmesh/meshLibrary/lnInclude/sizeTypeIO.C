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

Description
    Reads a long long from an input stream, for a given version
    number and File format. If an ascii File is being read, then the line
    numbers are counted and an erroneous read ised.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "sizeType.H"
#include "IOstreams.H"

#include <sstream>

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

word name(const sizeType val)
{
    std::ostringstream buf;
    buf << val;
    return buf.str();
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Istream& operator>>(Istream& is, sizeType& l)
{
    l = readSizeType(is);

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, sizeType&)");

    return is;
}


sizeType readSizeType(Istream& is)
{
    sizeType result = 0;

    //- read the digits and store them to a string
    std::string s;
    s.reserve(20);
    char c;
    bool foundDigit(false);
    while( is.read(c) )
    {
        if( !isdigit(c) )
        {
            //- character is not digit
            //- this can happen only before or after the number
            if( foundDigit )
            {
                break;
            }
            else
            {
                continue;
            }
        }

        s += c;
        foundDigit = true;
    }

    //- calculate the maximum power of 10
    const label maxPow = s.size() - 1;
    sizeType pow(1);
    for(label currPow=0;currPow<maxPow;++currPow)
    {
        const sizeType prevPow = pow;
        pow *= 10;

        if( pow < prevPow )
        {
            FatalErrorIn("sizeType readSizeType(Istream& is)")
                << "Read value " << s
                << " cannot be stored in sizeType!" << exit(FatalError);
        }
    }

    //- go digit by digit and sum it up to the result
    for(size_t i=0;i<s.size();++i)
    {
        if( !isdigit(s[i]) )
        {
            FatalErrorIn("sizeType readSizeType(Istream& is)")
                << "Read value " << s
                << " is not an integer" << exit(FatalError);
        }

        const sizeType digit = label(s[i] - '0');
        if( (digit < 0) || (digit > 9) )
        {
            FatalErrorIn("sizeType readSizeType(Istream& is)")
                << "Digit at position " << label(i)
                << " is not converted to an integer" << exit(FatalError);
        }

        result += digit * pow;
        pow /= 10;
    }

    return result;
}


bool readSizeType(const char* buf, sizeType& s)
{
    char *endptr = NULL;
    s = strtoll(buf, &endptr, 10);
    return (*endptr == 0);
}


Ostream& operator<<(Ostream& os, const sizeType l)
{
    const word val = name(l);

    for(size_t i=0;i<val.size();++i)
        os << val[i];

    os.check("Ostream& operator<<(Ostream&, const sizeType)");
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
