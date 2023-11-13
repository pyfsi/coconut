/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | Web:         http://www.foam-extend.org
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

#include "SymmTensorTemplate.H"
#include "TensorTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    // Construct null
    template <class Cmpt>
    inline DevSymmTensor4thOrder<Cmpt>::DevSymmTensor4thOrder()
    {}


    // Construct given VectorSpace
    template <class Cmpt>
    inline DevSymmTensor4thOrder<Cmpt>::DevSymmTensor4thOrder
    (
        const VectorSpace<DevSymmTensor4thOrder<Cmpt>, Cmpt, 15>& vs
    )
        :
        VectorSpace<DevSymmTensor4thOrder<Cmpt>, Cmpt, 15>(vs)
    {}


    // Construct given 15 Cmpts
    template <class Cmpt>
    inline DevSymmTensor4thOrder<Cmpt>::DevSymmTensor4thOrder
    (
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
    )
    {
        this->v_[XXXX] = txxxx;
        this->v_[XXXY] = txxxy;
        this->v_[YZYZ] = tyzyz;
        this->v_[XXYY] = txxyy;
        this->v_[YYYZ] = tyyyz;
        this->v_[XZYY] = txzyy;
        this->v_[XYXY] = txyxy;
        this->v_[XYYY] = txyyy;
        this->v_[XXXZ] = txxxz;
        this->v_[XXYZ] = txxyz;
        this->v_[XZYZ] = txzyz;
        this->v_[YYYY] = tyyyy;
        this->v_[XXXX] = txxxx;
        this->v_[XYYZ] = txyyz;
        this->v_[XZXZ] = txzxz;
        this->v_[XYXZ] = txyxz;
    }


    // Construct from Istream
    template <class Cmpt>
    inline DevSymmTensor4thOrder<Cmpt>::DevSymmTensor4thOrder(Istream& is)
        :
        VectorSpace<DevSymmTensor4thOrder<Cmpt>, Cmpt, 15>(is)
    {}


    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::xxxy() const
    {
        return this->v_[XXXY];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::yzyz() const
    {
        return this->v_[YZYZ];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::xxyy() const
    {
        return this->v_[XXYY];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::yyyz() const
    {
        return this->v_[YYYZ];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::xzyy() const
    {
        return this->v_[XZYY];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::xyxy() const
    {
        return this->v_[XYXY];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::xyyy() const
    {
        return this->v_[XYYY];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::xxxz() const
    {
        return this->v_[XXXZ];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::xxyz() const
    {
        return this->v_[XXYZ];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::xzyz() const
    {
        return this->v_[XZYZ];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::yyyy() const
    {
        return this->v_[YYYY];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::xxxx() const
    {
        return this->v_[XXXX];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::xyyz() const
    {
        return this->v_[XYYZ];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::xzxz() const
    {
        return this->v_[XZXZ];
    }
    template <class Cmpt>
    inline const Cmpt& DevSymmTensor4thOrder<Cmpt>::xyxz() const
    {
        return this->v_[XYXZ];
    }

    ///////////////////////////////////////////////////////////

    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::xxxy()
    {
        return this->v_[XXXY];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::yzyz()
    {
        return this->v_[YZYZ];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::xxyy()
    {
        return this->v_[XXYY];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::yyyz()
    {
        return this->v_[YYYZ];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::xzyy()
    {
        return this->v_[XZYY];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::xyxy()
    {
        return this->v_[XYXY];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::xyyy()
    {
        return this->v_[XYYY];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::xxxz()
    {
        return this->v_[XXXZ];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::xxyz()
    {
        return this->v_[XXYZ];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::xzyz()
    {
        return this->v_[XZYZ];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::yyyy()
    {
        return this->v_[YYYY];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::xxxx()
    {
        return this->v_[XXXX];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::xyyz()
    {
        return this->v_[XYYZ];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::xzxz()
    {
        return this->v_[XZXZ];
    }
    template <class Cmpt>
    inline Cmpt& DevSymmTensor4thOrder<Cmpt>::xyxz()
    {
        return this->v_[XYXZ];
    }


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

    // TODO: tidy and indent
    // Put checks and errors if necessary

template <class Cmpt>
inline DevSymmTensor4thOrder<Cmpt>
operator*(const SymmTensor<Cmpt>& a, const SymmTensor<Cmpt>& b)
{
    // TODO put a check in here, the tensors a and b must be equal, best way to
    // do this?

    // s_ij = C_ijkl E_kl
    return DevSymmTensor4thOrder<Cmpt>
    (
        a.xx()*b.xy(),
        a.yz()*b.yz(),
        a.xx()*b.yy(),
        a.yy()*b.yz(),
        a.xz()*b.yy(),
        a.xy()*b.xy(),
        a.xy()*b.yy(),
        a.xx()*b.xz(),
        a.xx()*b.yz(),
        a.xz()*b.yz(),
        a.yy()*b.yy(),
        a.xx()*b.xx(),
        a.xy()*b.yz(),
        a.xz()*b.xz(),
        a.xy()*b.xz()
    );
}


template <class Cmpt>
inline DevSymmTensor4thOrder<Cmpt>
operator*(const DevSymmTensor4thOrder<Cmpt>& a, const Cmpt& b)
{
    //- s_ij = C_ijkl E_kl
    return DevSymmTensor4thOrder<Cmpt>
    (
        b*a.xxxy(),
        b*a.yzyz(),
        b*a.xxyy(),
        b*a.yyyz(),
        b*a.xzyy(),
        b*a.xyxy(),
        b*a.xyyy(),
        b*a.xxxz(),
        b*a.xxyz(),
        b*a.xzyz(),
        b*a.yyyy(),
        b*a.xxxx(),
        b*a.xyyz(),
        b*a.xzxz(),
        b*a.xyxz()
    );
}


template <class Cmpt>
inline DevSymmTensor4thOrder<Cmpt>
operator*(const Cmpt& b, const DevSymmTensor4thOrder<Cmpt>& a)
{
    // s_ij = C_ijkl E_kl
    return Foam::operator*(a,b);
}


    //- Double-dot-product between a fourth order symmetric tensor and
    //  a symmetric tensor
    template <class Cmpt>
    inline SymmTensor<Cmpt>
    operator&&(const DevSymmTensor4thOrder<Cmpt>& t4th, const SymmTensor<Cmpt>& st)
    {
        //- s_ij = C_ijkl E_kl

        return SymmTensor<Cmpt>
            (
                t4th.xxxx()*st.xx() + 2*t4th.xxxy()*st.xy() + 2*t4th.xxxz()*st.xz() + t4th.xxyy()*st.yy() + 2*t4th.xxyz()*st.yz() - (t4th.xxxx() + t4th.xxyy())*st.zz(),
                t4th.xxxy()*st.xx() + 2*t4th.xyxy()*st.xy() + 2*t4th.xyxz()*st.xz() + t4th.xyyy()*st.yy() + 2*t4th.xyyz()*st.yz() - (t4th.xxxy() + t4th.xyyy())*st.zz(),
                t4th.xxxz()*st.xx() + 2*t4th.xyxz()*st.xy() + 2*t4th.xzxz()*st.xz() + t4th.xzyy()*st.yy() + 2*t4th.xzyz()*st.yz() - (t4th.xxxz() + t4th.xzyy())*st.zz(),
                t4th.xxyy()*st.xx() + 2*t4th.xyyy()*st.xy() + 2*t4th.xzyy()*st.xz() + t4th.yyyy()*st.yy() + 2*t4th.yyyz()*st.yz() - (t4th.xxyy() + t4th.yyyy())*st.zz(),
                t4th.xxyz()*st.xx() + 2*t4th.xyyz()*st.xy() + 2*t4th.xzyz()*st.xz() + t4th.yyyz()*st.yy() + 2*t4th.yzyz()*st.yz() - (t4th.xxyz() + t4th.yyyz())*st.zz(),
                -(t4th.xxxx() + t4th.xxyy())*st.xx() - 2*(t4th.xxxy() + t4th.xyyy())*st.xy() - 2*(t4th.xxxz() + t4th.xzyy())*st.xz() -
                (t4th.xxyy() + t4th.yyyy())*st.yy() - 2*(t4th.xxyz() + t4th.yyyz())*st.yz() + (t4th.xxxx() + 2*t4th.xxyy() + t4th.yyyy())*st.zz()
            );
    }


    template <class Cmpt>
    inline Cmpt
    operator&&(const DevSymmTensor4thOrder<Cmpt>& A, const DevSymmTensor4thOrder<Cmpt>& B)
    {
        //- s_ij = C_ijkl E_kl
        return
            (
                A.xxxx()*B.xxxx() +
                A.xxxy()*B.xxxy() +
                A.xxxz()*B.xxxz() +
                A.xxxy()*B.xxxy() +
                A.xxyy()*B.xxyy() +
                A.xxyz()*B.xxyz() +
                A.xxxz()*B.xxxz() +
                A.xxyz()*B.xxyz() +
                (-A.xxxx() - A.xxyy())*(-B.xxxx() - B.xxyy()) +
                A.xxxy()*B.xxxy() +
                A.xyxy()*B.xyxy() +
                A.xyxz()*B.xyxz() +
                A.xyxy()*B.xyxy() +
                A.xyyy()*B.xyyy() +
                A.xyyz()*B.xyyz() +
                A.xyxz()*B.xyxz() +
                A.xyyz()*B.xyyz() +
                (-A.xxxy() - A.xyyy())*(-B.xxxy() - B.xyyy()) +
                A.xxxz()*B.xxxz() +
                A.xyxz()*B.xyxz() +
                A.xzxz()*B.xzxz() +
                A.xyxz()*B.xyxz() +
                A.xzyy()*B.xzyy() +
                A.xzyz()*B.xzyz() +
                A.xzxz()*B.xzxz() +
                A.xzyz()*B.xzyz() +
                (-A.xxxz() - A.xzyy())*(-B.xxxz() - B.xzyy()) +
                A.xxxy()*B.xxxy() +
                A.xyxy()*B.xyxy() +
                A.xyxz()*B.xyxz() +
                A.xyxy()*B.xyxy() +
                A.xyyy()*B.xyyy() +
                A.xyyz()*B.xyyz() +
                A.xyxz()*B.xyxz() +
                A.xyyz()*B.xyyz() +
                (-A.xxxy() - A.xyyy())*(-B.xxxy() - B.xyyy()) +
                A.xxyy()*B.xxyy() +
                A.xyyy()*B.xyyy() +
                A.xzyy()*B.xzyy() +
                A.xyyy()*B.xyyy() +
                A.yyyy()*B.yyyy() +
                A.yyyz()*B.yyyz() +
                A.xzyy()*B.xzyy() +
                A.yyyz()*B.yyyz() +
                (-A.xxyy() - A.yyyy())*(-B.xxyy() - B.yyyy()) +
                A.xxyz()*B.xxyz() +
                A.xyyz()*B.xyyz() +
                A.xzyz()*B.xzyz() +
                A.xyyz()*B.xyyz() +
                A.yyyz()*B.yyyz() +
                A.yzyz()*B.yzyz() +
                A.xzyz()*B.xzyz() +
                A.yzyz()*B.yzyz() +
                (-A.xxyz() - A.yyyz())*(-B.xxyz() - B.yyyz()) +
                A.xxxz()*B.xxxz() +
                A.xyxz()*B.xyxz() +
                A.xzxz()*B.xzxz() +
                A.xyxz()*B.xyxz() +
                A.xzyy()*B.xzyy() +
                A.xzyz()*B.xzyz() +
                A.xzxz()*B.xzxz() +
                A.xzyz()*B.xzyz() +
                (-A.xxxz() - A.xzyy())*(-B.xxxz() - B.xzyy()) +
                A.xxyz()*B.xxyz() +
                A.xyyz()*B.xyyz() +
                A.xzyz()*B.xzyz() +
                A.xyyz()*B.xyyz() +
                A.yyyz()*B.yyyz() +
                A.yzyz()*B.yzyz() +
                A.xzyz()*B.xzyz() +
                A.yzyz()*B.yzyz() +
                (-A.xxyz() - A.yyyz())*(-B.xxyz() - B.yyyz()) +
                (-A.xxxx()-A.xxyy())*(-B.xxxx()-B.xxyy()) +
                (-A.xxxy()-A.xyyy())*(-B.xxxy()-B.xyyy()) +
                (-A.xxxz()-A.xzyy())*(-B.xxxz()-B.xzyy()) +
                (-A.xxxy()-A.xyyy())*(-B.xxxy()-B.xyyy()) +
                (-A.xxyy()-A.yyyy())*(-B.xxyy()-B.yyyy()) +
                (-A.xxyz()-A.yyyz())*(-B.xxyz()-B.yyyz()) +
                (-A.xxxz()-A.xzyy())*(-B.xxxz()-B.xzyy()) +
                (-A.xxyz()-A.yyyz())*(-B.xxyz()-B.yyyz()) +
                (-(-A.xxxx()-A.xxyy()) - (-A.xxyy()-A.yyyy()))*(-(-B.xxxx()-B.xxyy()) - (-B.xxyy()-B.yyyy()))
            );
    }


    //- Double-dot-product between a symmetric tensor and a symmetric
    //  fourth order tensor
    template <class Cmpt>
    inline SymmTensor<Cmpt>
    operator&&(const SymmTensor<Cmpt>& st, const DevSymmTensor4thOrder<Cmpt>& t4th)
    {
        //- s_ij =  E_kl C_klij
        //- as C_ijkl = C_klij
        //- then s_ij = C_ijkl E_kl
        return (t4th && st);
    }

    //- Inner-product between a symmTensor4thOrder and a vector
    //  this is not implemented
    template <class Cmpt>
    inline Vector<Cmpt>
    operator&(const DevSymmTensor4thOrder<Cmpt>& st, const Vector<Cmpt>& v)
    {
        notImplemented
            (
                "operator&(const DevSymmTensor4thOrder<Cmpt>& st, const Vector<Cmpt>& v)"
                " not implemented in DevSymmTensor4thOrderTemplateI.H"
            );

        return Vector<Cmpt>
            (
                v.x(),
                v.y(),
                v.z()
            );
    }


    //- Inner-product between a vector and a symmTensor4thOrder
    //  this is not implemented
    template <class Cmpt>
    inline Vector<Cmpt>
    operator&(const Vector<Cmpt>& v, const DevSymmTensor4thOrder<Cmpt>& st)
    {
        notImplemented
            (
                "operator&(const Vector<Cmpt>& v, const DevSymmTensor4thOrder<Cmpt>& st)"
                " not implemented in DevSymmTensor4thOrderTemplateI.H"
            );

        return Vector<Cmpt>
            (
                v.x(),
                v.y(),
                v.z()
            );
    }


    template<class Cmpt>
    class outerProduct<DevSymmTensor4thOrder<Cmpt>, Cmpt>
    {
    public:

        typedef DevSymmTensor4thOrder<Cmpt> type;
    };

    template<class Cmpt>
    class outerProduct<Cmpt, DevSymmTensor4thOrder<Cmpt> >
    {
    public:

        typedef DevSymmTensor4thOrder<Cmpt> type;
    };

    template<class Cmpt>
    class innerProduct<DevSymmTensor4thOrder<Cmpt>, Vector<Cmpt> >
    {
    public:

        typedef Vector<Cmpt> type;
    };

    template<class Cmpt>
    class innerProduct<Vector<Cmpt>, DevSymmTensor4thOrder<Cmpt> >
    {
    public:

        typedef Vector<Cmpt> type;
    };

    template<class Cmpt>
    class scalarProduct<DevSymmTensor4thOrder<Cmpt>, SymmTensor<Cmpt> >
    {
    public:

        typedef SymmTensor<Cmpt> type;
    };

    template<class Cmpt>
    class scalarProduct<SymmTensor<Cmpt>, DevSymmTensor4thOrder<Cmpt> >
    {
    public:

        typedef SymmTensor<Cmpt> type;
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
