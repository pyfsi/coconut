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

#ifndef fadbadTensor14_H
#define fadbadTensor14_H

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
    typedef fadbad::F<scalar, 14> adScalar14;
    typedef Tensor<adScalar14> adTensor14;
    typedef SymmTensor<adScalar14> adSymmTensor14;
    typedef Vector<adScalar14> adVector14;
    typedef DevSymmTensor4thOrder<adScalar14> adDevSymmTensor4thOrder14;




    // template class Tensor<fadbad::F<scalar, 14> >;
    // template class SymmTensor<fadbad::F<scalar, 14> >;
    // template class Vector<fadbad::F<scalar, 14> >;


    //adScalar14 sgn(adScalar14 val);
    //scalar sgn(scalar val);


    inline Tensor<adScalar14 >
    operator*(const Tensor<adScalar14 >& T, const scalar a)
    {
        return Tensor<adScalar14 >
            (
                T.xx()*a, T.xy()*a, T.xz()*a,
                T.yx()*a, T.yy()*a, T.yz()*a,
                T.zx()*a, T.zy()*a, T.zz()*a
            );
    }

    inline Tensor<adScalar14 >
    operator*(const scalar a, const Tensor<adScalar14 >& T)
    {
        return T*a;
    }


    inline adVector14
    operator*(const adVector14& T, const adScalar14 a)
    {
        return adVector14
            (
                T.x()*a, T.y()*a, T.z()*a
            );
    }
    inline adVector14
    operator*(const adScalar14 a, const adVector14& T)
    {
        return T*a;
    }

    inline adVector14
    operator/(const adVector14& T, const adScalar14 a)
    {
        return adVector14
            (
                T.x()/a, T.y()/a, T.z()/a
            );
    }



    inline adSymmTensor14
    operator*(const adSymmTensor14& T, const scalar a)
    {
        return adSymmTensor14
            (
                T.xx()*a, T.xy()*a, T.xz()*a,
                T.yy()*a, T.yz()*a,
                T.zz()*a
            );
    }

    inline adSymmTensor14
    operator*(const scalar a, const adSymmTensor14& T)
    {
        return T*a;
    }



    inline adSymmTensor14
    operator/(const adSymmTensor14& T, const scalar a)
    {
        return adSymmTensor14
            (
                T.xx()/a, T.xy()/a, T.xz()/a,
                T.yy()/a, T.yz()/a,
                T.zz()/a
            );
    }


    inline adVector14
    operator/(const adVector14& T, const scalar a)
    {
        return adVector14
            (
                T.x()/a, T.y()/a, T.z()/a
            );
    }


    inline adSymmTensor14
    operator*(const adSymmTensor14& T, const adScalar14& a)
    {
        return adSymmTensor14
            (
                T.xx()*a, T.xy()*a, T.xz()*a,
                T.yy()*a, T.yz()*a,
                T.zz()*a
            );
    }

    inline adSymmTensor14
    operator*(const adScalar14& a, const adSymmTensor14& T)
    {
        return T*a;
    }



    inline adSymmTensor14
    operator/(const adSymmTensor14& T, const adScalar14 a)
    {
        return adSymmTensor14
            (
                T.xx()/a, T.xy()/a, T.xz()/a,
                T.yy()/a, T.yz()/a,
                T.zz()/a
            );
    }



    inline Tensor<adScalar14 >
    operator/(const Tensor<adScalar14 >& T, const adScalar14 a)
    {
        return Tensor<adScalar14 >
            (
                T.xx()/a, T.xy()/a, T.xz()/a,
                T.yx()/a, T.yy()/a, T.yz()/a,
                T.zx()/a, T.zy()/a, T.zz()/a
            );
    }






    adScalar14 expAD(adScalar14 a);
    adScalar14 powAD(adScalar14 a, adScalar14 b);
    adScalar14 sqrtAD(adScalar14 a);
    adScalar14 logAD(adScalar14 a);
    adScalar14 magAD(adSymmTensor14 T);
    adScalar14 magAD(Vector<adScalar14> T);


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // these should eventually be replaced with implicit conversions
    // convert to fadbad objects

    adSymmTensor14
    fadbadConvert14( symmTensor T);

    adTensor14
    fadbadConvert14( Tensor<double> T);

    adDevSymmTensor4thOrder14
    fadbadConvert14( devSymmTensor4thOrder T);

    // convert fadbad objects back to scalar objects
    symmTensor
    fadbadConvert14( adSymmTensor14 T);

    Tensor<double>
    fadbadConvert14( Tensor<adScalar14> T);

    devSymmTensor4thOrder
    fadbadConvert14( adDevSymmTensor4thOrder14 T);


    Istream& operator<<(Istream& o, adScalar14 a);

    Istream& operator>>(Istream& o, adScalar14 a);

    Ostream& operator<<(Ostream& o, adScalar14 a);


    inline SymmTensor<adScalar14> operator+(const SphericalTensor<double>& A, const SymmTensor<adScalar14>& B)
    {
        return SymmTensor<adScalar14>(A.ii()+B.xx(), B.xy(), B.xz(), A.ii()+B.yy(), B.yz(), A.ii()+B.zz());
    }


    inline SphericalTensor<adScalar14> operator*(const SphericalTensor<adScalar14>& A, const adScalar14 B)
    {
        return SphericalTensor<adScalar14>(A.ii()*B);
    }
    inline SphericalTensor<adScalar14> operator*(const adScalar14 B, const SphericalTensor<adScalar14>& A)
    {
        return SphericalTensor<adScalar14>(A.ii()*B);
    }

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
