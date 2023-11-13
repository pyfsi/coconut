/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
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

#include "elastoPlasticDistortional.H"

// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //


template<class t>
Foam::SymmTensor<t> Foam::elastoPlasticDistortional::calculateStress
(
    const SymmTensor<t>& C,
    const SymmTensor<t>& CInv,
    const SymmTensor<t>& Cp,
    const SymmTensor<t>& CpInv
)
{
    return
        mu_.value()*(CpInv - CInv)
     + (lambda_.value()/2.0)*(det(C)/det(Cp) - 1.0)*CInv;
}


template<class t>
void Foam::elastoPlasticDistortional::pack_X
(
    List<t>& X,
    const t& deltaGamma,
    const SymmTensor<t>& effectiveStress
)
{
    X.setSize(7);

    X[0] = deltaGamma;
    X[1] = effectiveStress.xx();
    X[2] = effectiveStress.xy();
    X[3] = effectiveStress.xz();
    X[4] = effectiveStress.yy();
    X[5] = effectiveStress.yz();
    X[6] = effectiveStress.zz();
}


template<class t>
void Foam::elastoPlasticDistortional::unpack_X
(
    const List<t>& X,
    t& deltaGamma,
    SymmTensor<t>& effectiveStress
)
{
    if (X.size() != 7)
    {
        FatalErrorIn
        (
            "template<class t>\n"
            "void Foam::elastoPlasticDistortional::unpack_X\n"
            "(\n"
            "    const List<t>& X,\n"
            "    t& deltaGamma,\n"
            "    SymmTensor<t>& effectiveStress\n"
            ")"
        )   << "unpack_X: X not right size" << abort(FatalError);
    }

    deltaGamma = X[0];
    effectiveStress.xx() = X[1];
    effectiveStress.xy() = X[2];
    effectiveStress.xz() = X[3];
    effectiveStress.yy() = X[4];
    effectiveStress.yz() = X[5];
    effectiveStress.zz() = X[6];
}


// ************************************************************************* //
