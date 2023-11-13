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
#include "surfaceQuadricMetric.H"
#include "partTetMeshSimplex.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Private member functions

scalar surfaceQuadricMetric::evaluateMetric() const
{
    scalar val(0.0);

    forAll(normals_, nI)
    {
        const scalar d = normals_[nI] & (p_ - centres_[nI]);
        const scalar stab = sqrt(sqr(d) + Ksq_);

        val += 4.0 / sqr(stab + d);
    }

    return val;
}

void surfaceQuadricMetric::evaluateGradients
(
    vector& grad,
    tensor& gradGrad
) const
{
    grad = vector::zero;
    gradGrad = tensor::zero;

    forAll(normals_, nI)
    {
        const vector& n = normals_[nI];

        const scalar d = (n & (p_ - centres_[nI]));
        const scalar stab = sqrt(sqr(d) + Ksq_);

        grad += -8.0 * (((d * n) / stab) + n) / pow(stab + d, 3);
        gradGrad +=
            24.0 * (d * n / stab + n) * (d * n / stab + n) / pow(d + stab, 4)
            -8.0 * (n * n / stab - sqr(d) * n * n / pow(stab, 3)) /
            pow(d + stab, 3);
    }

    gradGrad.zz() = SMALL * (gradGrad.xx() + gradGrad.yy());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

surfaceQuadricMetric::surfaceQuadricMetric
(
    DynList<point>& pts,
    const DynList<triFace>& trias
)
:
    pts_(pts),
    trias_(trias),
    p_(pts[trias[0][0]]),
    normals_(),
    centres_(),
    bb_(),
    Ksq_(SMALL)
{
    bb_.min() = bb_.max() = pts_[trias[0][1]];

    forAll(trias_, tI)
    {
        const triFace& pt = trias_[tI];

        for(direction dir=0;dir<vector::nComponents;++dir)
        {
            bb_.min()[dir] = std::min(bb_.min()[dir], pts_[pt[1]][dir]);
            bb_.max()[dir] = std::max(bb_.max()[dir], pts_[pt[1]][dir]);
            bb_.min()[dir] = std::min(bb_.min()[dir], pts_[pt[2]][dir]);
            bb_.max()[dir] = std::max(bb_.max()[dir], pts_[pt[2]][dir]);
        }

        vector n = vector(0., 0., 1.0) ^ (pts_[pt[2]] - pts_[pt[1]]);
        n /= (mag(n) + VSMALL);

        centres_.append(0.5 * (pts_[pt[2]] + pts_[pt[1]]));
        normals_.append(n);
    }

    Ksq_ = 1e-6 * bb_.mag();
}


surfaceQuadricMetric::~surfaceQuadricMetric()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

point surfaceQuadricMetric::optimizePoint(const scalar)
{
    const scalar scale = 1.0 / bb_.mag();
    forAll(pts_, pI)
        pts_[pI] *= scale;
    forAll(centres_, cI)
        centres_[cI] *= scale;
    bb_.min() *= scale;
    bb_.max() *= scale;

    if( !bb_.contains(p_) )
        p_ = 0.5 * (bb_.min() + bb_.max());

    const scalar tol = sqr(2.0 * SMALL) * magSqr(bb_.min() - bb_.max());

    label iterI, outerIter(0);

    vector gradF, disp;
    tensor gradGradF;

    scalar func, lastFunc;

    # ifdef DEBUGSmooth
    forAll(normals_, nI)
    {
        const scalar fx = normals_[nI] & (p_ - centres_[nI]);
        Info << "Tet " << nI << " has distance " << fx << " func "
            << CfMesh::sqr(mag(fx) - fx) << endl;
    }
    Info << "BoundBox size " << (bb_.max() - bb_.min()) << endl;
    Info << "Tolerance " << tol << endl;
    # endif

    bool finished;
    do
    {
        finished = true;

        lastFunc = evaluateMetric();

        iterI = 0;
        do
        {
            # ifdef DEBUGSmooth
            Info << "Iteration " << iterI << endl;
            Info << "Initial metric value " << lastFunc << endl;
            # endif

            //- store previous value
            const point pOrig = p_;

            //- evaluate gradients
            evaluateGradients(gradF, gradGradF);

            //- calculate displacement
            const scalar determinant = det(gradGradF);
            if( mag(determinant) > SMALL )
            {
                disp = (inv(gradGradF, determinant) & gradF);

                for(direction i=0;i<vector::nComponents;++i)
                {
                    const scalar& val = disp[i];
                    if( (val != val) || ((val - val) != (val - val)) )
                    {
                        disp = vector::zero;
                        break;
                    }
                }

                p_ -= disp;

                func = evaluateMetric();

                # ifdef DEBUGSmooth
                Info << "Second grad " << gradGradF << endl;
                Info << "inv(gradGradF, determinant) "
                    << inv(gradGradF, determinant) << endl;
                Info << "Gradient " << gradF << endl;
                Info << "Determinant " << determinant << endl;
                Info << "Displacement " << disp << endl;
                Info << "New metric value " << func << endl;
                # endif

                scalar relax(0.8);
                label nLoops(0);
                while( func > lastFunc )
                {
                    p_ = pOrig - relax * disp;
                    relax *= 0.5;
                    func = evaluateMetric();

                    if( func < lastFunc )
                        continue;

                    //- it seems that this direction is wrong
                    if( ++nLoops == 5 )
                    {
                        p_ = pOrig;
                        disp = vector::zero;
                        func = 0.0;
                    }
                }

                lastFunc = func;
            }
            else
            {
                Info << "Tu sam. Determinant " << determinant << " gradGradF " << gradGradF << endl;
                disp = vector::zero;
            }
        } while( (magSqr(disp) > tol) && (++iterI < 10) );

        if( lastFunc < VSMALL )
            finished = false;
    } while( !finished && (++outerIter < 5) );

    //- scale back to the original size
    forAll(pts_, pI)
        pts_[pI] /= scale;
    forAll(centres_, cI)
        centres_[cI] /= scale;
    bb_.min() /= scale;
    bb_.max() /= scale;

    # ifdef DEBUGSmooth
    Info << "Last value " << lastFunc << endl;
    Info << "Metric value " << evaluateMetric() << endl;
    forAll(normals_, nI)
    {
        const scalar fx = normals_[nI] & (p_ - centres_[nI]);
        Info << "Tet " << nI << " has distance " << fx << endl;
    }
    # endif

    return p_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
