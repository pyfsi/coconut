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

Description

\*---------------------------------------------------------------------------*/

#include "quadricFitting.H"

//#define DEBUGQuadric

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * * * //

inline void quadricFitting::calculateNormalVector()
{
    if( otherPoints_.size() == 0 )
    {
        isValid_ = false;
        return;
    }

    //- assemble a covariance matrix
    symmTensor mat(symmTensor::zero);

    forAll(otherPoints_, i)
    {
        const vector d = (otherPoints_[i] - origin_);

        mat += symm(d * d);
    }

    mat /= otherPoints_.size();

    //- small determinants are also invalid
    if( mag(det(mat)) < SMALL )
    {
        normal_ = vector::zero;
        isValid_ = false;
        return;
    }

    //- find the eigenvalues of the tensor
    try
    {
        const vector ev = eigenValues(mat);

        //- estimate the normal as the eigenvector associated
        //- to the smallest eigenvalue
        normal_ = eigenVector(mat, ev[0]);

        # ifdef DEBUGQuadric
        Info << "Estimated normal of the quadric patch is " << normal_ << endl;
        Info << "ev " << ev << endl;
        Info << "Middle eigenvector " << eigenVector(mat, ev[1]) << endl;
        Info << "Max eigenvector " << eigenVector(mat, ev[2]) << endl;
        # endif
    }
    catch(const char*)
    {
        normal_ = vector::zero;
        isValid_ = false;
    }
}

inline void quadricFitting::calculateCoordinateSystem()
{
    if( mag(normal_) < VSMALL )
    {
        vecX_ = vecY_ = vector::zero;
        forAll(transformedPoints_, i)
            transformedPoints_[i] = vector::zero;

        isValid_ = false;
        return;
    }

    const plane pl(origin_, normal_);

    bool valid(false);

    label pI(0);
    while( !valid )
    {
        const point vp = pl.nearestPoint(otherPoints_[pI]);
        vecX_ = vp - origin_;

        if( magSqr(vecX_) < VSMALL )
        {
            ++pI;
            continue;
        }

        vecX_ /= mag(vecX_);
        vecY_ = normal_ ^ vecX_;
        valid = true;
    }

    if( !valid )
    {
        normal_ = vecX_ = vecY_ = vector::zero;
        isValid_ = false;
    }

    # ifdef DEBUGQuadric
    Info << "vecX " << vecX_ << endl;
    Info << "vecY " << vecY_ << endl;
    Info << "normal " << normal_ << endl;
    # endif

    //- transform other points to this new coordinate system
    transformedPoints_.setSize(otherPoints_.size());

    forAll(otherPoints_, i)
    {
        const vector delta = (otherPoints_[i] - origin_);

        # ifdef DEBUGQuadric
        Info << "Other point " << otherPoints_[i] << endl;
        Info << "Delta " << delta << endl;
        # endif

        point& p = transformedPoints_[i];

        p.x() = delta & vecX_;
        p.y() = delta & vecY_;
        p.z() = delta & normal_;
    }
}

inline void quadricFitting::calculateQuadricCoeffs()
{
    if( mag(normal_) < VSMALL )
    {
        coefficients_.setSize(5);
        coefficients_ = 0.0;
        isValid_ = false;

        return;
    }

    simpleMatrix<scalar> mat(5, 0.0, 0.0);
    DynList<scalar> helper;
    helper.setSize(5);

    forAll(transformedPoints_, i)
    {
        const point& p = transformedPoints_[i];

        helper[0] = sqr(p.x());
        helper[1] = sqr(p.y());
        helper[2] = p.x() * p.y();
        helper[3] = p.x();
        helper[4] = p.y();

        for(label rowI=0;rowI<5;++rowI)
        {
            for(label colI=rowI;colI<5;++colI)
                mat[rowI][colI] += helper[rowI] * helper[colI];

            mat.source()[rowI] += helper[rowI] * p.z();
        }
    }

    for(label rowI=0;rowI<5;++rowI)
    {
        for(label colI=0;colI<rowI;++colI)
            mat[rowI][colI] = mat[colI][rowI];

        if( mag(mat[rowI][rowI]) < SMALL )
            mat[rowI][rowI] = SMALL;
    }

    # ifdef DEBUGQuadric
    Info << "Matrix " << mat << endl;
    # endif

    coefficients_ = mat.LUsolve();

    cumulativeError_ = 0.0;
    forAll(transformedPoints_, pI)
    {
        const point& tp = transformedPoints_[pI];
        const point onQuadric = pointCoordinates(tp.x(), tp.y());

        cumulativeError_ += magSqr(onQuadric - tp);
    }

    # ifdef DEBUGQuadric
    Info << "Quadric error " << cumulativeError_ << endl;
    # endif
}

inline void quadricFitting::calculateBestFit()
{
    label iteration(0);

    while( iteration++ < 10 )
    {
        calculateCoordinateSystem();

        calculateQuadricCoeffs();

        if( (mag(coefficients_[3]) > SMALL) || (mag(coefficients_[4]) > SMALL) )
        {
            //- correct the normal
            const scalar d =
                Foam::sqrt(1.0 + sqr(coefficients_[3]) + sqr(coefficients_[4]));

            const vector newNormal
            (
                normal_ / d -
                coefficients_[3] * vecX_ / d -
                coefficients_[4] * vecY_ / d
            );

            normal_ = newNormal;

            # ifdef DEBUGQuadric
            Info << "Updated normal " << normal_ << endl;
            # endif
        }
        else
        {
            break;
        }
    }

    calculateCoordinateSystem();

    # ifdef DEBUGQuadric
    Info << "vecX " << vecX_ << endl;
    Info << "vecY " << vecY_ << endl;
    Info << "normal " << normal_ << endl;
    # endif

    calculateQuadricCoeffs();

    # ifdef DEBUGQuadric
    Info << "Other points " << otherPoints_ << endl;
    Info << "Transformed points " << transformedPoints_ << endl;
    Info << "normal " << normal_ << endl;
    Info << "Coefficients " << coefficients_ << endl;
    Info << "isValid " << isValid_ << endl;
    # endif
}

// * * * * * * * * * * * * * * * *  Constructors * * * * * * * * * * * * * * //

template<class ListType>
inline quadricFitting::quadricFitting
(
    const point& origin,
    const vector normal,
    const ListType& otherPoints
)
:
    origin_(origin),
    normal_(normal),
    vecX_(vector::zero),
    vecY_(vector::zero),
    otherPoints_(),
    transformedPoints_(),
    coefficients_(),
    cumulativeError_(),
    isValid_(true)
{
    otherPoints_.setSize(otherPoints.size());
    forAll(otherPoints, i)
        otherPoints_[i] = otherPoints[i];

    const scalar magN = mag(normal_);
    if( magN < VSMALL )
    {
        normal_ = vector::zero;

        calculateNormalVector();
    }
    else
    {
        normal_ /= (magN + VSMALL);
    }

    # ifdef DEBUGQuadric
    Info << "Constructing quadric patch from points " << origin
         << " normal " << normal_ << " other points " << otherPoints_ << endl;
    # endif

    calculateBestFit();
}

//- Construct from point and other points
template<class ListType>
inline quadricFitting::quadricFitting
(
    const point& origin,
    const ListType& otherPoints
)
:
    origin_(origin),
    normal_(),
    vecX_(),
    vecY_(),
    otherPoints_(),
    transformedPoints_(),
    coefficients_(),
    cumulativeError_(),
    isValid_(true)
{
    otherPoints_.setSize(otherPoints.size());
    forAll(otherPoints, i)
        otherPoints_[i] = otherPoints[i];

    calculateNormalVector();

    calculateBestFit();
}

inline quadricFitting::~quadricFitting()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool quadricFitting::isValid() const
{
    return isValid_;
}

inline scalar quadricFitting::cumulativeError() const
{
    return cumulativeError_;
}

inline const vector& quadricFitting::normal() const
{
    return normal_;
}

//- Return min curvature
inline scalar quadricFitting::minCurvature() const
{
    return
    (
        coefficients_[0] + coefficients_[1] -
        sqrt
        (
            max
            (
                sqr(coefficients_[0] - coefficients_[1]) +
                sqr(coefficients_[2]),
                VSMALL
            )
        )
    );
}

//- Return max curvature
inline scalar quadricFitting::maxCurvature() const
{
    return
    (
        coefficients_[0] + coefficients_[1] +
        sqrt
        (
            max
            (
                sqr(coefficients_[0] - coefficients_[1]) +
                sqr(coefficients_[2]),
                VSMALL
            )
        )
    );
}

//- Return mean curvature
inline scalar quadricFitting::meanCurvature() const
{
    return 0.5 * (minCurvature() + maxCurvature());
}

//- Return Gaussian curvature
inline scalar quadricFitting::gaussianCurvature() const
{
    return minCurvature() * maxCurvature();
}

//- Return min curvature vector
inline vector quadricFitting::minCurvatureVector() const
{
    const scalar theta =
    0.5 * Foam::atan
    (
        coefficients_[2] /
        stabilise(coefficients_[0] - coefficients_[1], VSMALL)
    );

    return Foam::cos(theta) * vecX_ + Foam::sin(theta) * vecY_;
}

//- Return max curvature vector
inline vector quadricFitting::maxCurvatureVector() const
{
    const scalar theta =
        0.5 * Foam::atan
        (
            coefficients_[2] /
            stabilise(coefficients_[0] - coefficients_[1], VSMALL)
        );

    return Foam::sin(theta) * vecX_ + Foam::cos(theta) * vecY_;
}

inline const DynList<point, 32>& quadricFitting::transformedPoints() const
{
    return transformedPoints_;
}

inline point quadricFitting::pointCoordinates
(
    const scalar u,
    const scalar v
) const
{
    point newP = origin_;
    newP += u * vecX_;
    newP += v * vecY_;
    newP +=
        (
            sqr(u) * coefficients_[0] +
            sqr(v) * coefficients_[1] +
            u * v * coefficients_[2] +
            u * coefficients_[3] +
            v * coefficients_[4]
        ) * normal_;


    return newP;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
