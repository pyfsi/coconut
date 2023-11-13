/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::DevSymmTensor4thOrder

Description

Note:
    Order is currently hacked due to missing operations on symmetry planes
    HJ, 29/Oct/2013

SourceFiles
    DevSymmTensor4thOrderTemplateI.H

Author
    Michael Clancy UCD

\*---------------------------------------------------------------------------*/

#ifndef DevSymmTensor4thOrderTemplate_H
#define DevSymmTensor4thOrderTemplate_H

#include "VectorSpace.H"
#include "SymmTensorTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DevSymmTensor4thOrder Declaration
\*---------------------------------------------------------------------------*/

template <class Cmpt>
class DevSymmTensor4thOrder
:
    public VectorSpace<DevSymmTensor4thOrder<Cmpt>, Cmpt, 15>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef DevSymmTensor4thOrder<label> labelType;

    // Member constants

        enum
        {
            rank = 2 // Hacked by PC: Actual rank of DevSymmTensor4thOrder is 4
        };

        //HJ, reconsider and fix.  HJ, 23/Feb/2014


    // Static data members

        static const char* const typeName;
        static const char* componentNames[];

    static const DevSymmTensor4thOrder zero;
    static const DevSymmTensor4thOrder one;
    static const DevSymmTensor4thOrder max;
    static const DevSymmTensor4thOrder min;




        //- Component labeling enumeration
        enum components
        {
            XXXY,
            YZYZ,
            XXYY,
            YYYZ,
            XZYY,
            XYXY,
            XYYY,
            XXXZ,
            XXYZ,
            XZYZ,
            YYYY,
            XXXX,
            XYYZ,
            XZXZ,
            XYXZ
        };


    // Constructors

        //- Construct null
        inline DevSymmTensor4thOrder();

        //- Construct given VectorSpace
        inline DevSymmTensor4thOrder
        (
            const VectorSpace<DevSymmTensor4thOrder<Cmpt>, Cmpt, 15>&
        );

        //- Construct given SymmTensor (for compilation)
        inline DevSymmTensor4thOrder(const SymmTensor<Cmpt>&);

        //- Construct given the 15 components
        inline DevSymmTensor4thOrder
        (
            // TODO: tidy these files
      const Cmpt txxxy,
      const Cmpt tyzyz,
      const Cmpt txxyy,
      const Cmpt tyyyz,
      const Cmpt txzyy,
      const Cmpt txyxy,
      const Cmpt txyyy,
      const Cmpt txxxz,
      const Cmpt txxyz,
      const Cmpt txzyz,
      const Cmpt tyyyy,
      const Cmpt txxxx,
      const Cmpt txyyz,
      const Cmpt txzxz,
      const Cmpt txyxz
        );


        //- Construct from Istream
        DevSymmTensor4thOrder(Istream&);

    // Member Functions

        // Access

    inline const Cmpt& xxxy() const;
    inline const Cmpt& yzyz() const;
    inline const Cmpt& xxyy() const;
    inline const Cmpt& yyyz() const;
    inline const Cmpt& xzyy() const;
    inline const Cmpt& xyxy() const;
    inline const Cmpt& xyyy() const;
    inline const Cmpt& xxxz() const;
    inline const Cmpt& xxyz() const;
    inline const Cmpt& xzyz() const;
    inline const Cmpt& yyyy() const;
    inline const Cmpt& xxxx() const;
    inline const Cmpt& xyyz() const;
    inline const Cmpt& xzxz() const;
    inline const Cmpt& xyxz() const;

    inline Cmpt& xxxy();
    inline Cmpt& yzyz();
    inline Cmpt& xxyy();
    inline Cmpt& yyyz();
    inline Cmpt& xzyy();
    inline Cmpt& xyxy();
    inline Cmpt& xyyy();
    inline Cmpt& xxxz();
    inline Cmpt& xxyz();
    inline Cmpt& xzyz();
    inline Cmpt& yyyy();
    inline Cmpt& xxxx();
    inline Cmpt& xyyz();
    inline Cmpt& xzxz();
    inline Cmpt& xyxz();

        //- Assign to a SymmTensor
        //- maybe I need this to compile
        inline void operator=(const SymmTensor<Cmpt>&);
};


// template<class Cmpt>
// class typeOfRank<Cmpt, 4>
// {
// public:

//     typedef DevSymmTensor4thOrder<Cmpt> type;
// };


// template<class Cmpt>
// class symmTypeOfRank<Cmpt, 4>
// {
// public:

//     typedef DevSymmTensor4thOrder<Cmpt> type;
// };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "DevSymmTensor4thOrderTemplateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
