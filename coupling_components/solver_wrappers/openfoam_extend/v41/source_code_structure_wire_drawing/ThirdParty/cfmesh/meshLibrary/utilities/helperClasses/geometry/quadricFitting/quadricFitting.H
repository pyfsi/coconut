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
    quadricFitting

Description
    A class for fitting a quadric surface through the set of points. It is
    useful for curvature estimation and for calculating normals.

SourceFiles
    quadricFittingI.H

\*---------------------------------------------------------------------------*/

#ifndef quadricFitting_H
#define quadricFitting_H

#include "DynList.H"
#include "FixedList.H"
#include "plane.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        class quadricFitting Declaration
\*---------------------------------------------------------------------------*/

class quadricFitting
{
    // Private data
        //- seed point for fitting the quadric surface
        const point origin_;

        //- normal vector at the origin
        vector normal_;

        //- x-coordinate vector for the transformed coordinate system
        vector vecX_;

        //- y-coordinate vector for the transformed coordinate system
        vector vecY_;

        //- other points
        DynList<point, 32> otherPoints_;

        //- transformed coordinates of other points
        DynList<point, 32> transformedPoints_;

        //- coefficients of the quadric surface
        FixedList<scalar, 5> coefficients_;

        //- sum of approximation error
        scalar cumulativeError_;

        //- it is a valid quadric patch
        bool isValid_;

    // Private member functions
        //- initial estimate of the normal vector
        inline void calculateNormalVector();

        //- calculate transformed coordinate system
        inline void calculateCoordinateSystem();

        //- calculate coefficients of the quadric surface fit
        inline void calculateQuadricCoeffs();

        //- iterative procedure for updating normals
        //- until the best fit is achieved
        inline void calculateBestFit();

        //- Disallow default bitwise copy construct
        quadricFitting(const quadricFitting&);

        //- Disallow default bitwise assignment
        void operator=(const quadricFitting&);

public:

    // Constructors

        //- Construct from point, normal, and neighbouring points
        template<class ListType>
        inline quadricFitting(const point&, const vector, const ListType&);

        //- Construct from point and other points
        template<class ListType>
        inline quadricFitting(const point&, const ListType&);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Destructor
        inline ~quadricFitting();

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Member Functions

        // Access
            //- provides information if it is a valid quadric patch
            inline bool isValid() const;

            //- return the error as the sum of squared distances between the
            //- original point coordinates and the quadric patch
            inline scalar cumulativeError() const;

            //- Return surface normal
            inline const vector& normal() const;

            //- Return min curvature
            inline scalar minCurvature() const;

            //- Return max curvature
            inline scalar maxCurvature() const;

            //- Return mean curvature
            inline scalar meanCurvature() const;

            //- Return Gaussian curvature
            inline scalar gaussianCurvature() const;

            //- Return min curvature vector
            inline vector minCurvatureVector() const;

            //- Return max curvature vector
            inline vector maxCurvatureVector() const;

            //- return the transformed coordinates in the 2D tangent plane
            inline const DynList<point, 32>& transformedPoints() const;

            //- calculate point at the quadric patch with the given
            //- parametric coordinates
            inline point pointCoordinates(const scalar u, const scalar v) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "quadricFittingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
