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

Typedef
    Foam::adTensor

Description
    AdTensor of adScalars.

SourceFiles
    adTensor.C

Author
    Michael Clancy, UCD.

\*---------------------------------------------------------------------------*/

#ifndef fadbadTensor7_H
#define fadbadTensor7_H

#include "fadiff.h"
#include "badiff.h"
#include "TensorTemplate.H"
#include "contiguous.H"
#include "SphericalTensorTemplate.H"
#include "SymmTensorTemplate.H"
#include "symmTensor.H"
#include "DevSymmTensor4thOrderTemplate.H"
#include "devSymmTensor4thOrder.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    typedef fadbad::F<scalar, 7> adScalar7;
    typedef Tensor<adScalar7> adTensor7;
    typedef SymmTensor<adScalar7> adSymmTensor7;
    typedef Vector<adScalar7> adVector7;
    typedef DevSymmTensor4thOrder<adScalar7> adDevSymmTensor4thOrder7;




    // template class Tensor<fadbad::F<scalar, 7> >;
    // template class SymmTensor<fadbad::F<scalar, 7> >;
    // template class Vector<fadbad::F<scalar, 7> >;


    //adScalar7 sgn(adScalar7 val);
    //scalar sgn(scalar val);


    inline Tensor<adScalar7 >
    operator*(const Tensor<adScalar7 >& T, const scalar a)
    {
        return Tensor<adScalar7 >
            (
                T.xx()*a, T.xy()*a, T.xz()*a,
                T.yx()*a, T.yy()*a, T.yz()*a,
                T.zx()*a, T.zy()*a, T.zz()*a
            );
    }

    inline Tensor<adScalar7 >
    operator*(const scalar a, const Tensor<adScalar7 >& T)
    {
        return T*a;
    }


    inline adVector7
    operator*(const adVector7& T, const adScalar7 a)
    {
        return adVector7
            (
                T.x()*a, T.y()*a, T.z()*a
            );
    }
    inline adVector7
    operator*(const adScalar7 a, const adVector7& T)
    {
        return T*a;
    }

    inline adVector7
    operator/(const adVector7& T, const adScalar7 a)
    {
        return adVector7
            (
                T.x()/a, T.y()/a, T.z()/a
            );
    }



    inline adSymmTensor7
    operator*(const adSymmTensor7& T, const scalar a)
    {
        return adSymmTensor7
            (
                T.xx()*a, T.xy()*a, T.xz()*a,
                T.yy()*a, T.yz()*a,
                T.zz()*a
            );
    }

    inline adSymmTensor7
    operator*(const scalar a, const adSymmTensor7& T)
    {
        return T*a;
    }



    inline adSymmTensor7
    operator/(const adSymmTensor7& T, const scalar a)
    {
        return adSymmTensor7
            (
                T.xx()/a, T.xy()/a, T.xz()/a,
                T.yy()/a, T.yz()/a,
                T.zz()/a
            );
    }


    inline adVector7
    operator/(const adVector7& T, const scalar a)
    {
        return adVector7
            (
                T.x()/a, T.y()/a, T.z()/a
            );
    }


    inline adSymmTensor7
    operator*(const adSymmTensor7& T, const adScalar7& a)
    {
        return adSymmTensor7
            (
                T.xx()*a, T.xy()*a, T.xz()*a,
                T.yy()*a, T.yz()*a,
                T.zz()*a
            );
    }

    inline adSymmTensor7
    operator*(const adScalar7& a, const adSymmTensor7& T)
    {
        return T*a;
    }



    inline adSymmTensor7
    operator/(const adSymmTensor7& T, const adScalar7 a)
    {
        return adSymmTensor7
            (
                T.xx()/a, T.xy()/a, T.xz()/a,
                T.yy()/a, T.yz()/a,
                T.zz()/a
            );
    }



    inline Tensor<adScalar7 >
    operator/(const Tensor<adScalar7 >& T, const adScalar7 a)
    {
        return Tensor<adScalar7 >
            (
                T.xx()/a, T.xy()/a, T.xz()/a,
                T.yx()/a, T.yy()/a, T.yz()/a,
                T.zx()/a, T.zy()/a, T.zz()/a
            );
    }







    scalar expAD(scalar a);
    scalar powAD(scalar a, scalar b);
    scalar sqrtAD(scalar a);
    scalar logAD(scalar a);
    scalar magAD(symmTensor T);
    scalar magAD(vector T);



    adScalar7 expAD(adScalar7 a);
    adScalar7 powAD(adScalar7 a, adScalar7 b);
    adScalar7 sqrtAD(adScalar7 a);
    adScalar7 logAD(adScalar7 a);
    adScalar7 magAD(adSymmTensor7 T);
    adScalar7 magAD(Vector<adScalar7> T);


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // these should eventually be replaced with implicit conversions
    // convert to fadbad objects

    adSymmTensor7
    fadbadConvert7( symmTensor T);

    adTensor7
    fadbadConvert7( Tensor<double> T);

    adDevSymmTensor4thOrder7
    fadbadConvert7( devSymmTensor4thOrder T);

    // convert fadbad objects back to scalar objects
    symmTensor
    fadbadConvert7( adSymmTensor7 T);

    Tensor<double>
    fadbadConvert7( Tensor<adScalar7> T);

    devSymmTensor4thOrder
    fadbadConvert7( adDevSymmTensor4thOrder7 T);


    Istream& operator<<(Istream& o, adScalar7 a);

    Istream& operator>>(Istream& o, adScalar7 a);

    Ostream& operator<<(Ostream& o, adScalar7 a);


    inline SymmTensor<adScalar7> operator+(const SphericalTensor<double>& A, const SymmTensor<adScalar7>& B)
    {
        return SymmTensor<adScalar7>(A.ii()+B.xx(), B.xy(), B.xz(), A.ii()+B.yy(), B.yz(), A.ii()+B.zz());
    }

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
