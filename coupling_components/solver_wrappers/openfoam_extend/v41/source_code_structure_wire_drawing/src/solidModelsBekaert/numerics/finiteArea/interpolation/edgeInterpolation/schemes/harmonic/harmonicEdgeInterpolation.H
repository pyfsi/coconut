/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    harmonicEdgeInterpolation

Description
    Harmonic-mean differencing scheme class.

    This scheme interpolates 1/field using a scheme specified at run-time
    and return the reciprocal of the interpolate.

    The scheme is based on the harmonic-mean differencing scheme class from the
    finiteVolume.

Authors
    Vanja Skuric, FMENA Zagreb. All rights reserved.

SourceFiles
    harmonicEdgeInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef harmonicEdgeInterpolation_H
#define harmonicEdgeInterpolation_H

#include "edgeInterpolationScheme.H"
#include "areaFields.H"
#include "magLongDeltaFa.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class harmonicEdgeInterpolation
:
    virtual public edgeInterpolationScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const harmonicEdgeInterpolation&);

public:

    //- Runtime type information
    TypeName("harmonic");


    // Constructors

        //- Construct from mesh
        harmonicEdgeInterpolation(const faMesh& mesh)
        :
            edgeInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        harmonicEdgeInterpolation(const faMesh& mesh, Istream&)
        :
            edgeInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        harmonicEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField&,
            Istream&
        )
        :
            edgeInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<edgeScalarField> weights
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        tmp<scalarField> weights
        (
            const Field<Type>& fOwn,
            const Field<Type>& fNei,
            const faPatch& patch
        ) const;
};


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
harmonicEdgeInterpolate(const GeometricField<Type, faPatchField, areaMesh>& vf)
{
    return edgeInterpolationScheme<Type>::interpolate
    (
        vf,
        vf.mesh().edgeInterpolation::weights()
    );
}


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
harmonicEdgeInterpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh> > tinterp =
        harmonicEdgeInterpolate(tvf());
    tvf.clear();
    return tinterp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "harmonicEdgeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
