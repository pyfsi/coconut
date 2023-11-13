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

Class
    Foam::fv::skewCorrectedGaussGrad

Description
    Point Gauss gradient calculation where the cell values are interpolated to
    the points using a least squares method.

SourceFiles
    skewCorrectedGaussGrad.C

Author
    Zeljko Tukovic Zagreb/UCD
    Re-organised by Philip Cardiff UCD/UT

\*---------------------------------------------------------------------------*/

#ifndef skewCorrectedGaussGrad_H
#define skewCorrectedGaussGrad_H

#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class skewCorrectedGaussGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class skewCorrectedGaussGrad
:
    public fv::gaussGrad<Type>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        skewCorrectedGaussGrad(const skewCorrectedGaussGrad&);

        //- Disallow default bitwise assignment
        void operator=(const skewCorrectedGaussGrad&);


public:

    //- Runtime type information
    TypeName("skewCorrectedGauss");


    // Constructors

        //- Construct from mesh
        skewCorrectedGaussGrad(const fvMesh& mesh)
        :
            gaussGrad<Type>(mesh)
        {}

        //- Construct from Istream
        skewCorrectedGaussGrad(const fvMesh& mesh, Istream&)
        :
            gaussGrad<Type>(mesh)
        {}


    // Member Functions

        //- Calculate and return the grad of the given field.
        //  Used by grad either to recalculate the cached gradient when it is
        //  out of date with respect to the field or when it is not cached.
        //- Volume gradient
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const word& name
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "skewCorrectedGaussGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
