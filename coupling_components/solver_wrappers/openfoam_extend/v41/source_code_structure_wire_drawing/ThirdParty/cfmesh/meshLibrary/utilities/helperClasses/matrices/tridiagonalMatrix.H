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
    tridiagonalMatrix

Description
    Implementation of a tridiagonal matrix and the TDMA solver

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef tridiagonalMatrix_H
#define tridiagonalMatrix_H

#include "scalar.H"
#include "FixedList.H"
#include "scalarLongList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class tridiagonalMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class tridiagonalMatrix
{
    // Private members
        //- matrix coefficients
        scalarLongList a_;
        scalarLongList b_;
        scalarLongList c_;

        //- source term
        LongList<T> source_;

    // Private member functions
        //- update coefficients by forward elimination
        inline void updateCoefficients();

        //- backward substitution
        inline void backwardSubstitution(LongList<T>& solution);

public:

    // Constructors
        //- Null constructor
        explicit inline tridiagonalMatrix();

        //- Construct from the number of rows
        inline tridiagonalMatrix(const label);

        //- Copy constructor
        inline tridiagonalMatrix(const tridiagonalMatrix&);

    // Destructor
        inline ~tridiagonalMatrix();

    // Member functions
        //- get the size of the matrix
        inline label size() const;

        //- set the number of rows of the matrix
        inline void setSize(const label);

        //- get the coefficient of the matrix
        inline scalar operator()(const label, const label) const;

        //- set the coefficient of the matrix
        inline void setCoeff(const label, const label, const scalar&);

        //- const access to the source
        inline const LongList<T>& source() const;

        //- set the source term at a given position
        inline void setSource(const label, const T&);

        //- find the solution of the system with the given rhs
        inline LongList<T> solve();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tridiagonalMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
