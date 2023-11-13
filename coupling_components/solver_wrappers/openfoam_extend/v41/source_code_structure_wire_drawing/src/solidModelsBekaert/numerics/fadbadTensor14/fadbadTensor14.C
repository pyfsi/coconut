/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |

     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "fadbadTensor14.H"

namespace Foam{


    // adScalar14 sgn(adScalar14 val) {
    //     return (0 < val) - (val < 0);
    // }
    // scalar sgn(scalar val) {
    //     return (0 < val) - (val < 0);
    // }

    template<>
    const adTensor14 adTensor14::zero
    (
        adScalar14(0.0), adScalar14(0.0), adScalar14(0.0),
        adScalar14(0.0), adScalar14(0.0), adScalar14(0.0),
        adScalar14(0.0), adScalar14(0.0), adScalar14(0.0)
    );

    template<>
    const SymmTensor<adScalar14> SymmTensor<adScalar14>::zero
    (
        adScalar14(0.0), adScalar14(0.0), adScalar14(0.0),
        adScalar14(0.0), adScalar14(0.0),
        adScalar14(0.0)
    );
    template<>
    const Tensor<adScalar14> Tensor<adScalar14>::one
    (
        adScalar14(1.0), adScalar14(0.0), adScalar14(0.0),
        adScalar14(0.0), adScalar14(1.0), adScalar14(0.0),
        adScalar14(0.0), adScalar14(0.0), adScalar14(1.0)
    );

    template<>
    const SymmTensor<adScalar14> SymmTensor<adScalar14>::one
    (
        adScalar14(1.0), adScalar14(0.0), adScalar14(0.0),
        adScalar14(1.0), adScalar14(0.0),
        adScalar14(1.0)
    );

    template<>
    const Vector<adScalar14> Vector<adScalar14>::zero
    (
        adScalar14(0.0), adScalar14(0.0), adScalar14(0.0)
    );



    template class Tensor<fadbad::F<scalar, 14> >;
    template class SymmTensor<fadbad::F<scalar, 14> >;
    template class Vector<fadbad::F<scalar, 14> >;

    adScalar14 expAD(adScalar14 a)
    {
        return fadbad::exp(a);
    }

    adScalar14 powAD(adScalar14 a, adScalar14 b)
    {
        return fadbad::pow(a, b);
    }

    adScalar14 sqrtAD(adScalar14 a)
    {
        return fadbad::sqrt(a);
    }

    adScalar14 logAD(adScalar14 a)
    {
        return fadbad::log(a);
    }

    adScalar14 magAD(adSymmTensor14 T)
    {
        return
            sqrtAD
            (
                T.xx()*T.xx() +T.yy()*T.yy() +T.zz()*T.zz() +
                2.0*T.xy()*T.xy() +
                2.0*T.yz()*T.yz() +
                2.0*T.xz()*T.xz()
            );
    }


    adScalar14 magAD(Vector<adScalar14> T)
    {
        return
            sqrtAD
            (
                T.y()*T.y() +
                T.x()*T.x() +
                T.z()*T.z()
            );
    }




    adSymmTensor14
    fadbadConvert14( symmTensor T)
    {
        return adSymmTensor14
            (
                adScalar14(T.xx()),adScalar14(T.xy()),adScalar14(T.xz()),
                adScalar14(T.yy()),adScalar14(T.yz()),adScalar14(T.zz())
            );
    }


    adTensor14
    fadbadConvert14( Tensor<double> T)
    {
        return Tensor<adScalar14>
            (
                adScalar14(T.xx()), adScalar14(T.xy()), adScalar14(T.xz()),
                adScalar14(T.yx()), adScalar14(T.yy()), adScalar14(T.yz()),
                adScalar14(T.zx()), adScalar14(T.zy()), adScalar14(T.zz())
            );
    }



    adDevSymmTensor4thOrder14
    fadbadConvert14( devSymmTensor4thOrder T)
    {
        return adDevSymmTensor4thOrder14
            (
                adScalar14(T.xxxy()),
                adScalar14(T.yzyz()),
                adScalar14(T.xxyy()),
                adScalar14(T.yyyz()),
                adScalar14(T.xzyy()),
                adScalar14(T.xyxy()),
                adScalar14(T.xyyy()),
                adScalar14(T.xxxz()),
                adScalar14(T.xxyz()),
                adScalar14(T.xzyz()),
                adScalar14(T.yyyy()),
                adScalar14(T.xxxx()),
                adScalar14(T.xyyz()),
                adScalar14(T.xzxz()),
                adScalar14(T.xyxz())
            );
    }


    symmTensor
    fadbadConvert14( adSymmTensor14 T)
    {
        return symmTensor
            (
                (T.xx().x()),(T.xy().x()),(T.xz().x()),
                (T.yy().x()),(T.yz().x()),(T.zz().x())
            );
    }
    Tensor<double>
    fadbadConvert14( Tensor<adScalar14> T)
    {
        return Tensor<double>
            (
                (T.xx().x()),(T.xy().x()),(T.xz().x()),
                (T.yx().x()),(T.yy().x()),(T.yz().x()),
                (T.zx().x()),(T.zy().x()),(T.zz().x())
            );
    }



    devSymmTensor4thOrder
    fadbadConvert14( adDevSymmTensor4thOrder14 T)
    {
        return devSymmTensor4thOrder
            (
                T.xxxy().x(),
                T.yzyz().x(),
                T.xxyy().x(),
                T.yyyz().x(),
                T.xzyy().x(),
                T.xyxy().x(),
                T.xyyy().x(),
                T.xxxz().x(),
                T.xxyz().x(),
                T.xzyz().x(),
                T.yyyy().x(),
                T.xxxx().x(),
                T.xyyz().x(),
                T.xzxz().x(),
                T.xyxz().x()
            );
    }


    Istream& operator<<(Istream& o, adScalar14 a)
    {
        return o;
    }

    Istream& operator>>(Istream& o, adScalar14 a)
    {
        return o;
    }

    Ostream& operator<<(Ostream& o, adScalar14 a)
    {
        o << '[' << a.x() << ',' << a.d(0) << ']';
        return o;
    }




  template<>
    const SphericalTensor<adScalar14> SphericalTensor<adScalar14>::I(adScalar14(1.0));
    template<>
    const SphericalTensor<adScalar14> SphericalTensor<adScalar14>::oneThirdI(adScalar14(1.0/3.0));
    template<>
    const SphericalTensor<adScalar14> SphericalTensor<adScalar14>::twoThirdsI(adScalar14(2.0/3.0));


    //template class SphericalTensor<adScalar14>;

template<>
const adDevSymmTensor4thOrder14 adDevSymmTensor4thOrder14::one
(
          0.0,
      3.0/4.0,
      -0.5,
      0.0,
      0.0,
      3.0/4.0,
      0.0,
      0.0,
      0.0,
      0.0,
      1.0,
      1.0,
      0.0,
      3.0/4.0,
      0.0
);


}
