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
    Implementation od 3 x 3 matrix

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
inline void tridiagonalMatrix<T>::updateCoefficients()
{
    if( b_.size() == 0 )
    {
        FatalErrorIn
        (
            "template<class T> "
            "inline void tridiagonalMatrix<T>::updateCoefficients()"
        ) << "Matrix has zero size" << abort(FatalError);
    }

    //- update the first element
    c_[0] /= b_[0];
    source_[0] /= b_[0];

    //- update remaining elements
    for(label i=1;i<b_.size();++i)
    {
        c_[i] /= b_[i] - a_[i] * c_[i-1];

        source_[i] =
            (source_[i] - a_[i] * source_[i-1]) / (b_[i] - a_[i] * c_[i-1]);
    }
}

template<class T>
inline void tridiagonalMatrix<T>::backwardSubstitution(LongList<T>& solution)
{
    const label maxRow = source_.size() - 1;

    solution.setSize(source_.size());

    solution[maxRow] = source_[maxRow];

    for(label i=maxRow-1;i>=0;--i)
    {
        solution[i] = source_[i] - c_[i] * solution[i+1];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
inline tridiagonalMatrix<T>::tridiagonalMatrix()
:
    a_(),
    b_(),
    c_(),
    source_()
{}

template<class T>
inline tridiagonalMatrix<T>::tridiagonalMatrix(const label s)
:
    a_(s, 0.0),
    b_(s, 0.0),
    c_(s, 0.0),
    source_(s, T::zero)
{}

template<class T>
inline tridiagonalMatrix<T>::tridiagonalMatrix(const tridiagonalMatrix<T>& tdm)
:
    a_(tdm.a_),
    b_(tdm.b_),
    c_(tdm.c_),
    source_(tdm.source_)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
inline tridiagonalMatrix<T>::~tridiagonalMatrix()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
inline label tridiagonalMatrix<T>::size() const
{
    return b_.size();
}

template<class T>
inline void tridiagonalMatrix<T>::setSize(const label s)
{
    //- set all matrix coefficients to zero
    a_.setSize(s);
    a_ = 0.0;
    b_.setSize(s);
    b_ = 0.0;
    c_.setSize(s);
    c_ = 0.0;

    source_.setSize(s);
    source_ = T::zero;
}

template<class T>
inline scalar tridiagonalMatrix<T>::operator()
(
    const label i,
    const label j
) const
{
    if( i == j )
    {
        return b_[i];
    }
    else if( j == (i+1) )
    {
        return c_[i];
    }
    else if( j == (i-1) )
    {
        return a_[i];
    }
    else
    {
        FatalErrorIn
        (
            "template<class T> inline scalar tridiagonalMatrix<T>::operator()"
            "(const label, const label) const"
        ) << "Indices (" << i << ", " << j << ") are not valid"
          << abort(FatalError);
    }
}

template<class T>
inline void tridiagonalMatrix<T>::setCoeff
(
    const label i,
    const label j,
    const scalar& coeff
)
{
    if( i == j )
    {
        b_[i] = coeff;
    }
    else if( j == (i+1) )
    {
        c_[i] = coeff;
    }
    else if( j == (i-1) )
    {
        a_[i] = coeff;
    }
    else
    {
        FatalErrorIn
        (
            "template<class T> inline scalar tridiagonalMatrix<T>::operator()"
            "(const label, const label) const"
        ) << "Indices (" << i << ", " << j << ") are not valid"
          << abort(FatalError);
    }
}

template<class T>
inline const LongList<T>& tridiagonalMatrix<T>::source() const
{
    return source_;
}

template<class T>
inline void tridiagonalMatrix<T>::setSource(const label i, const T& s)
{
    source_[i] = s;
}

template<class T>
inline LongList<T> tridiagonalMatrix<T>::solve()
{
    //- update coefficients of the matrix
    updateCoefficients();

    LongList<T> ret(source_.size());

    backwardSubstitution(ret);

    return ret;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
