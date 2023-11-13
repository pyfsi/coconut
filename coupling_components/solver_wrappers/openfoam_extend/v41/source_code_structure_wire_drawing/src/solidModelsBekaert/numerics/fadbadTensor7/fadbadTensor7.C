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

#include "fadbadTensor7.H"

namespace Foam{


    // adScalar7 sgn(adScalar7 val) {
    //     return (0 < val) - (val < 0);
    // }
    // scalar sgn(scalar val) {
    //     return (0 < val) - (val < 0);
    // }

    template<>
    const adTensor7 adTensor7::zero
    (
        adScalar7(0.0), adScalar7(0.0), adScalar7(0.0),
        adScalar7(0.0), adScalar7(0.0), adScalar7(0.0),
        adScalar7(0.0), adScalar7(0.0), adScalar7(0.0)
    );

    template<>
    const SymmTensor<adScalar7> SymmTensor<adScalar7>::zero
    (
        adScalar7(0.0), adScalar7(0.0), adScalar7(0.0),
        adScalar7(0.0), adScalar7(0.0),
        adScalar7(0.0)
    );
    template<>
    const Tensor<adScalar7> Tensor<adScalar7>::one
    (
        adScalar7(1.0), adScalar7(0.0), adScalar7(0.0),
        adScalar7(0.0), adScalar7(1.0), adScalar7(0.0),
        adScalar7(0.0), adScalar7(0.0), adScalar7(1.0)
    );

    template<>
    const SymmTensor<adScalar7> SymmTensor<adScalar7>::one
    (
        adScalar7(1.0), adScalar7(0.0), adScalar7(0.0),
        adScalar7(1.0), adScalar7(0.0),
        adScalar7(1.0)
    );

    template<>
    const Vector<adScalar7> Vector<adScalar7>::zero
    (
        adScalar7(0.0), adScalar7(0.0), adScalar7(0.0)
    );



    template class Tensor<fadbad::F<scalar, 7> >;
    template class SymmTensor<fadbad::F<scalar, 7> >;
    template class Vector<fadbad::F<scalar, 7> >;

    adScalar7 expAD(adScalar7 a)
    {
        return fadbad::exp(a);
    }

    adScalar7 powAD(adScalar7 a, adScalar7 b)
    {
        return fadbad::pow(a, b);
    }

    adScalar7 sqrtAD(adScalar7 a)
    {
        return fadbad::sqrt(a);
    }

    adScalar7 logAD(adScalar7 a)
    {
        return fadbad::log(a);
    }

    adScalar7 magAD(adSymmTensor7 T)
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


    adScalar7 magAD(Vector<adScalar7> T)
    {
        return
            sqrtAD
            (
                T.y()*T.y() +
                T.x()*T.x() +
                T.z()*T.z()
            );
    }



    scalar expAD(scalar a)
    {
        return Foam::exp(a);
    }

    scalar powAD(scalar a, scalar b)
    {
        return Foam::pow(a, b);
    }

    scalar sqrtAD(scalar a)
    {
        return Foam::sqrt(a);
    }

    scalar logAD(scalar a)
    {
        return Foam::log(a);
    }


    scalar magAD(symmTensor T)
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


    scalar magAD(vector T)
    {
        return
            sqrtAD
            (
                T.y()*T.y() +
                T.x()*T.x() +
                T.z()*T.z()
            );
    }




    adSymmTensor7
    fadbadConvert7(symmTensor T)
    {
        return adSymmTensor7
            (
                adScalar7(T.xx()),adScalar7(T.xy()),adScalar7(T.xz()),
                adScalar7(T.yy()),adScalar7(T.yz()),adScalar7(T.zz())
            );
    }


    adTensor7
    fadbadConvert7( Tensor<double> T)
    {
        return Tensor<adScalar7>
            (
                adScalar7(T.xx()), adScalar7(T.xy()), adScalar7(T.xz()),
                adScalar7(T.yx()), adScalar7(T.yy()), adScalar7(T.yz()),
                adScalar7(T.zx()), adScalar7(T.zy()), adScalar7(T.zz())
            );
    }



    adDevSymmTensor4thOrder7
    fadbadConvert7( devSymmTensor4thOrder T)
    {
        return adDevSymmTensor4thOrder7
            (
                adScalar7(T.xxxy()),
                adScalar7(T.yzyz()),
                adScalar7(T.xxyy()),
                adScalar7(T.yyyz()),
                adScalar7(T.xzyy()),
                adScalar7(T.xyxy()),
                adScalar7(T.xyyy()),
                adScalar7(T.xxxz()),
                adScalar7(T.xxyz()),
                adScalar7(T.xzyz()),
                adScalar7(T.yyyy()),
                adScalar7(T.xxxx()),
                adScalar7(T.xyyz()),
                adScalar7(T.xzxz()),
                adScalar7(T.xyxz())
            );
    }


    symmTensor
    fadbadConvert7( adSymmTensor7 T)
    {
        return symmTensor
            (
                (T.xx().x()),(T.xy().x()),(T.xz().x()),
                (T.yy().x()),(T.yz().x()),(T.zz().x())
            );
    }
    Tensor<double>
    fadbadConvert7( Tensor<adScalar7> T)
    {
        return Tensor<double>
            (
                (T.xx().x()),(T.xy().x()),(T.xz().x()),
                (T.yx().x()),(T.yy().x()),(T.yz().x()),
                (T.zx().x()),(T.zy().x()),(T.zz().x())
            );
    }



    devSymmTensor4thOrder
    fadbadConvert7( adDevSymmTensor4thOrder7 T)
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


    Istream& operator<<(Istream& o, adScalar7 a)
    {
        return o;
    }

    Istream& operator>>(Istream& o, adScalar7 a)
    {
        return o;
    }

    Ostream& operator<<(Ostream& o, adScalar7 a)
    {
        o << '[' << a.x() << ',' << a.d(0) << ']';
        return o;
    }

  template<>
    const SphericalTensor<adScalar7> SphericalTensor<adScalar7>::I(adScalar7(1.0));
    template<>
    const SphericalTensor<adScalar7> SphericalTensor<adScalar7>::oneThirdI(adScalar7(1.0/3.0));
    template<>
    const SphericalTensor<adScalar7> SphericalTensor<adScalar7>::twoThirdsI(adScalar7(2.0/3.0));



template<>
const adDevSymmTensor4thOrder7 adDevSymmTensor4thOrder7::one
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
