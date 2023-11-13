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
    refLabelledScalar

Description
    A class containing point label and its coordinates. It is used for
    exchanging data over processors

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef refLabelledScalar_H
#define refLabelledScalar_H

#include "labelledScalar.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class refLabelledScalar Declaration
\*---------------------------------------------------------------------------*/

class refLabelledScalar
{
    // Private data
        //- label
        label objectLabel_;

        //- value
        labelledScalar s_;

    public:

    // Constructors
        //- Null construct
        refLabelledScalar()
        :
            objectLabel_(),
            s_()
        {}

        //- Construct from components
        refLabelledScalar(const label sl, const labelledScalar& s)
        :
            objectLabel_(sl),
            s_(s)
        {}

    // Destructor
        ~refLabelledScalar()
        {}

    // Member functions
        //- return scalar label
        inline label objectLabel() const
        {
            return objectLabel_;
        }

        //- return the value
        inline const labelledScalar& lScalar() const
        {
            return s_;
        }

    // Member operators

        inline void operator=(const refLabelledScalar& ls)
        {
            objectLabel_ = ls.objectLabel_;
            s_ = ls.s_;
        }

        inline bool operator==(const refLabelledScalar& ls) const
        {
            if( objectLabel_ == ls.objectLabel_ )
                return true;

            return false;
        }

        inline bool operator<(const refLabelledScalar& ls) const
        {
            if( objectLabel_ < ls.objectLabel_ )
                return true;

            return false;
        }

        inline bool operator!=(const refLabelledScalar& ls) const
        {
            return !this->operator==(ls);
        }

    // Friend operators
        friend Ostream& operator<<(Ostream& os, const refLabelledScalar& ls)
        {
            os << token::BEGIN_LIST;
            os << ls.objectLabel_ << token::SPACE;
            os << ls.s_ << token::END_LIST;

            // Check state of Ostream
            os.check("operator<<(Ostream&, const refLabelledScalar&");

            return os;
        }

        friend Istream& operator>>(Istream& is, refLabelledScalar& ls)
        {
            // Read beginning of refLabelledScalar
            is.readBegin("refLabelledScalar");

            is >> ls.objectLabel_;
            is >> ls.s_;

            // Read end of refLabelledScalar
            is.readEnd("refLabelledScalar");

            // Check state of Istream
            is.check("operator>>(Istream&, refLabelledScalar");

            return is;
        }
};

//- Specify data associated with refLabelledScalar type is contiguous
template<>
inline bool contiguous<refLabelledScalar>() {return true;}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
