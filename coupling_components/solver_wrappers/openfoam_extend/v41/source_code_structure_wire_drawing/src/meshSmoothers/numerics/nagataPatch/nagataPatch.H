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
    Foam::nagataPatch

Description
    Construct a Nagata patch on a triSurface e.g. as described in:
    D.M. Neto, M.C. Oliveira, L.F. Menezes, J.L. Alves. Improving Nagata patch
    interpolation applied for tool surface description in sheet metal forming
    simulation. Computer-Aided Design 45 (2013) 639–656, and in the original
    paper by Nagata.

    A Nagata patch is a local parabolic surface fit to each triangle. We can
    then use this to approximate the original smooth geometry. It can be useful
    for projecting refined surface meshes or smoothly calculating distances to
    the surface.

Author
    Philip Cardiff, UCD.  All rights reserved

SourceFiles
    nagataPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nagataPatch_H
#define nagataPatch_H

#include "typeInfo.H"
#include "triSurface.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class nagataPatch Declaration
\*---------------------------------------------------------------------------*/

class nagataPatch
{
    // Private data

        //- Reference to the triSurface
        const triSurface& triSurf_;

        //- Coefficient for switching to linear interpolation instead of
        //  quadratic to control unrealistic curvatures
        //  For no limiting, set epsilon1 to zero and epsilon2 to a large number
        const scalar epsilon1_;
        const scalar epsilon2_;

        //- Correct point normals on patch edges
        const bool correctPointNormalsOnPatchEdges_;

        //- 6 interpolation coefficients for each face
        mutable autoPtr< List< List<vector> > > faceInterpCoeffsPtr_;

        //- 3 curvature vectors for each face
        mutable autoPtr< List< List<vector> > > faceCurvatureCoeffsPtr_;

        //- faceOnEdgeOfPatch list: it is true for faces at the
        //  edge of a patch and zero otherwise
        mutable autoPtr<boolList> faceOnEdgeOfPatchPtr_;


    // Private Member Functions

        //- Calculate the face interpolation coefficients
        void calcFaceInterpCoeffs() const;

        //- Return a reference to the face interpolation coefficients
        const List< List<vector> >&  faceInterpCoeffs() const;

        //- Calculate the face curvative coefficients
        void calcFaceCurvatureCoeffs() const;

        //- Return a reference to the face curvative coefficients
        const List< List<vector> >& faceCurvatureCoeffs() const;

        //- Determine local coordinates of the given point within the face
        void localCoordinates
        (
            scalar& eta, scalar& xi, const point& pt, const label triFaceID
        );

        //- Calculate curvature coefficient
        void curvatureCoeff
        (
            vector& coeff,
            const vector& x0,
            const vector& x1,
            const vector& n0,
            const vector& n1
        ) const;

        //- Multiply rectangular matrices
        //  Result returned as a copy
        RectangularMatrix<scalar> multiply
        (
            const RectangularMatrix<scalar>& a,
            const RectangularMatrix<scalar>& b
        ) const;

        //- Calculate point normal using neighbour faces in the given region
        void updatePointNormalUsingFacesNormalInRegion
        (
            vector& n,
            const label pointID,
            const label curRegionID,
            const triSurface& triSurf
        ) const;

        //- Return faceOnEdgeOfPatch bool list
        const boolList& faceOnEdgeOfPatch() const;

        //- Calculate the faceOnEdgeOfPatch list
        void calcFaceOnEdgeOfPatch() const;

        //- Disallow default bitwise copy construct
        nagataPatch(const nagataPatch&);

        //- Disallow default bitwise assignment
        void operator=(const nagataPatch&);


public:

    //- Runtime type information
    TypeName("nagataPatch");


    // Constructors

        //- Construct from triSurface
        nagataPatch
        (
            const triSurface& triSurf,
            const scalar epsilon1 = 0.00625, //0.0125,
            const scalar epsilon2 = 0.2, //0.1,
            const bool correctPointNormalsOnPatchEdges = true
        );


    //- Destructor
    virtual ~nagataPatch();


    // Member Functions

        //- Project point to Nagata patch surface, given the face Nagata
        //  coefficients of the triangle we are projecting on
        void projectPoint
        (
            point& projectedPoint,
            const List<vector>& faceCoeffs,
            const scalar eta,
            const scalar xi
        );

        //- Project points to Nagata patch surface
        tmp<pointField> projectPoints
        (
            const pointField& origPoints,
            const labelList& triFaceIDs
        );

        //- Refine and project patches
        triSurface refineAndProjectPatches
        (
            const wordList& patchesToProject, const int nRefinementLevels
        );

        //- Remove empty patches from the given triSurf
        //  We pass the triSurf as an argument so it is clear that we are not
        //  changing the triSurf stored as a reference in any other function
        //  It might make more sense to defined this function outside of the
        //  class
        void removeEmptyPacthes(triSurface& triSurf) const;

        //- Clearout demand driven data
        void clearOut();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// #ifdef NoRepository
// #   include "nagataPatchTemplates.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
