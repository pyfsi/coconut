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

#include "demandDrivenData.H"
#include "surfaceOptimizerHeight.H"
#include "matrix2D.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const vector surfaceOptimizerHeight::dirVecs[4] =
    {
        vector(-1.0, -1.0, 0.0),
        vector(1.0, -1.0, 0.0),
        vector(-1.0, 1.0, 0.0),
        vector(1.0, 1.0, 0.0)
    };

const tensor surfaceOptimizerHeight::gradGradLsq_
(
    4. ,0. ,0.,
    0. ,4. ,0.,
    0. ,0. ,0.
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar surfaceOptimizerHeight::evaluateStabilisationFactor() const
{
    //- find the maximum Lsq
    scalar LsqMax(0.0);
    forAll(trias_, triI)
    {
        const point& p0 = pts_[trias_[triI][0]];
        const point& p1 = pts_[trias_[triI][1]];
        const point& p2 = pts_[trias_[triI][2]];

        # ifdef DEBUGSmooth
        Info << "Triangle " << triI << " area " << Atri << endl;
        # endif

        const scalar LSqrTri
        (
            magSqr(p0 - p1) +
            magSqr(p2 - p0)
        );

        LsqMax = Foam::max(LsqMax, LSqrTri);
    }

    # ifdef DEBUGSmooth
    Info << "Amin " << Amin << endl;
    Info << "LsqMax " << LsqMax << endl;
    # endif

    //- K is greater than zero in case the stabilisation is needed
    scalar K = SMALL * LsqMax;

    return K;
}

scalar surfaceOptimizerHeight::evaluateFunc(const scalar& K) const
{
    scalar func(0.0);
    scalar sumW(0.0);

    forAll(trias_, triI)
    {
        const triFace& tri = trias_[triI];

        const point& p0 = pts_[tri[0]];
        const point& p1 = pts_[tri[1]];
        const point& p2 = pts_[tri[2]];

        const scalar w = weights_[triI];

        const vector v = p2 - p1;
        vector n(-v.y(), v.x(), 0.0);
        n /= (sqrt(sqr(v.x()) + sqr(v.y())) + VSMALL);

        const scalar h = (p0 - p1) & n;

        const scalar stab = sqrt(sqr(h) + K);

        # ifdef DEBUGSmooth
        Info << "Triangle " << triI << " height " << h << endl;
        # endif

        const scalar LSqrTri
        (
            magSqr(p0 - p1) +
            magSqr(p2 - p0)
        );

        const scalar Hstab = Foam::max(VSMALL, 0.5 * (h + stab));
        func += w * LSqrTri / Hstab;
        sumW += w;
    }

    return func / (sumW + VSMALL);
}

//- evaluate gradients needed for optimisation
void surfaceOptimizerHeight::evaluateGradients
(
    const scalar& K,
    vector& gradF,
    tensor& gradGradF
) const
{
    gradF = vector::zero;
    gradGradF = tensor::zero;

    scalar sumW(0.0);

    forAll(trias_, triI)
    {
        const triFace& tri = trias_[triI];

        const point& p0 = pts_[tri[0]];
        const point& p1 = pts_[tri[1]];
        const point& p2 = pts_[tri[2]];

        if( magSqr(p1 - p2) < VSMALL )
            continue;

        const scalar w = weights_[triI];

        const scalar LSqrTri
        (
            magSqr(p0 - p1) +
            magSqr(p2 - p0)
        );

        const vector v = p2 - p1;
        vector n(-v.y(), v.x(), 0.0);
        n /= (sqrt(sqr(v.x()) + sqr(v.y())) + VSMALL);

        const scalar Htri = (p0 - p1) & n;

        const scalar stab = sqrt(sqr(Htri) + K);

        const scalar Hstab = Foam::max(VSMALL, 0.5 * (Htri + stab));

        const vector gradHstab = 0.5 * (n + Htri * n / stab);

        const tensor gradGradHstab =
            0.5 *(1.0 - sqr(Htri/stab)) * (n * n) / stab;

        const vector gradLt(4.0 * p0 - 2.0 * p1 - 2.0 * p2);

        //- calculate the gradient
        const scalar sqrHstab = sqr(Hstab);
        gradF += w * gradLt / Hstab - (LSqrTri * gradHstab) / sqrHstab;

        //- calculate the second gradient
        gradGradF +=
            w *
            (
                gradGradLsq_ / Hstab -
                twoSymm(gradLt * gradHstab) / sqrHstab -
                gradGradHstab * LSqrTri / sqrHstab +
                2.0 * LSqrTri * (gradHstab * gradHstab) / (sqrHstab * Hstab)
            );

        sumW += w;
    }

    gradF /= (sumW + VSMALL);
    gradGradF /= (sumW + VSMALL);
}

scalar surfaceOptimizerHeight::optimiseDivideAndConquer(const scalar tol)
{
    point& pOpt = pts_[trias_[0][0]];

    pOpt = 0.5 * (pMax_ + pMin_);
    point currCentre = pOpt;
    scalar dx = (pMax_.x() - pMin_.x()) / 2.0;
    scalar dy = (pMax_.y() - pMin_.y()) / 2.0;

    label iter(0);

    //- find the value of the functional in the centre of the bnd box
    scalar K = evaluateStabilisationFactor();
    scalar funcBefore, funcAfter(evaluateFunc(K));

    do
    {
        funcBefore = funcAfter;

        funcAfter = VGREAT;
        point minCentre(vector::zero);

        for(label i=0;i<4;++i)
        {
            pOpt.x() = currCentre.x() + 0.5 * dirVecs[i].x() * dx;
            pOpt.y() = currCentre.y() + 0.5 * dirVecs[i].y() * dy;

            K = evaluateStabilisationFactor();
            const scalar func = evaluateFunc(K);

            if( func < funcAfter )
            {
                minCentre = pOpt;
                funcAfter = func;
            }
        }

        //- set the centre with the minimum value
        //- as the centre for future search
        currCentre = minCentre;
        pOpt = minCentre;

        //- halve the search range
        dx *= 0.5;
        dy *= 0.5;

        //- calculate the tolerence
        const scalar t = mag(funcAfter - funcBefore) / funcAfter;

        # ifdef DEBUGSmooth
        Info << "Point position " << pOpt << endl;
        Info << "Func before " << funcBefore << endl;
        Info << "Func after " << funcAfter << endl;
        Info << "Normalised difference " << t << endl;
        # endif

        if( t < tol )
            break;
    } while( ++iter < 100 );

    return funcAfter;
}

scalar surfaceOptimizerHeight::optimiseSteepestDescent(const scalar tol)
{
    point& pOpt = pts_[trias_[0][0]];

    //- find the bounding box
    const scalar avgEdge = Foam::mag(pMax_ - pMin_);

    //- find the minimum value on the 5 x 5 raster
    scalar K = evaluateStabilisationFactor();
    scalar funcBefore, funcAfter(evaluateFunc(K));

    //- start with steepest descent optimisation
    vector gradF;
    tensor gradGradF;
    vector disp;
    disp.z() = 0.0;

    direction nIterations(0);
    do
    {
        funcBefore = funcAfter;

        evaluateGradients(K, gradF, gradGradF);

        //- store data into a matrix
        matrix2D mat;
        mat[0][0] = gradGradF.xx();
        mat[0][1] = gradGradF.xy();
        mat[1][0] = gradGradF.yx();
        mat[1][1] = gradGradF.yy();
        FixedList<scalar, 2> source;
        source[0] = gradF.x();
        source[1] = gradF.y();

        //- calculate the determinant
        const scalar det = mat.determinant();

        if( mag(det) < VSMALL )
        {
            disp = vector::zero;
        }
        else
        {
            disp.x() = mat.solveFirst(source);
            disp.y() = mat.solveSecond(source);

            if( mag(disp) > 0.2 * avgEdge )
            {
                vector dir = disp / mag(disp);

                disp = dir * 0.2 * avgEdge;
            }
        }

        # ifdef DEBUGSmooth
        Info << "Second gradient " << gradGradF << endl;
        Info << "Gradient " << gradF << endl;
        Info << "Displacement " << disp << endl;
        Info << "K = " << K << endl;
        # endif

        pOpt -= disp;

        K = evaluateStabilisationFactor();
        funcAfter = evaluateFunc(K);

        if( mag(funcAfter - funcBefore) / funcBefore < tol )
            break;

        #ifdef DEBUGSmooth
        Info << "New coordinates " << pOpt << endl;
        # endif

    } while( ++nIterations < 100 );

    return funcAfter;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

surfaceOptimizerHeight::surfaceOptimizerHeight
(
    DynList<point>& pts,
    const DynList<triFace>& trias,
    const DynList<scalar>& weights
)
:
    pts_(pts),
    trias_(trias),
    weights_(trias.size(), 1.0),
    pMin_(),
    pMax_()
{
    if( weights.size() )
    {
        if( weights.size() == trias.size() )
        {
            weights_ = weights;
        }
        else
        {
            FatalError << "The number of weights " << weights.size()
                << " is not equal to the number of triangles " << trias.size()
                << abort(FatalError);
        }
    }

    pMin_ = pts_[trias_[0][1]];
    pMax_ = pMin_;

    forAll(trias_, triI)
    {
        const triFace& tf = trias_[triI];

        for(label i=1;i<3;++i)
        {
            pMin_ = Foam::min(pMin_, pts_[tf[i]]);
            pMax_ = Foam::max(pMax_, pts_[tf[i]]);
        }
    }
}

surfaceOptimizerHeight::~surfaceOptimizerHeight()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

point surfaceOptimizerHeight::optimizePoint(const scalar tol)
{
    const scalar scale = mag(pMax_ - pMin_);
    forAll(pts_, i)
        pts_[i] /= scale;
    pMin_ /= scale;
    pMax_ /= scale;

    point& pOpt = pts_[trias_[0][0]];

    const point pOrig = pOpt;

    try
    {
        const scalar funcDivide = optimiseDivideAndConquer(tol);
        const point newPoint = pOpt;

        const scalar funcSteepest = optimiseSteepestDescent(tol);

        if( funcSteepest > funcDivide )
            pOpt = newPoint;
    }
    catch(...)
    {
        pOpt = pOrig;
    }

    forAll(pts_, i)
        pts_[i] *= scale;
    pMin_ *= scale;
    pMax_ *= scale;

    return pOpt;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
