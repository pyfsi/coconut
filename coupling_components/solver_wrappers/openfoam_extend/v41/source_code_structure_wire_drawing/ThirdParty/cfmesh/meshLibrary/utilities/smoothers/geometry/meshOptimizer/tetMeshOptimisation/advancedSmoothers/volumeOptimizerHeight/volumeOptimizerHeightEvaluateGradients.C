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
#include "volumeOptimizerHeight.H"
#include "tetrahedron.H"
#include "helperFunctions.H"

#include <stdexcept>

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar volumeOptimizerHeight::evaluateFunc() const
{
    const scalar K = evaluateStabilisationFactor();

    scalar func(0.0);

    forAll(tets_, tetI)
    {
        const partTet& pt = tets_[tetI];
        const tetrahedron<point, point> tet
        (
            points_[pt.a()],
            points_[pt.b()],
            points_[pt.c()],
            points_[pt.d()]
        );

        const scalar LSqrTri
        (
            magSqr(tet.d() - tet.a()) +
            magSqr(tet.d() - tet.b()) +
            magSqr(tet.d() - tet.c())
        );

        vector nTri = (tet.b() - tet.a()) ^ (tet.c() - tet.a());
        nTri /= (mag(nTri) + VSMALL);

        const scalar hTet = (tet.d() - tet.a()) & nTri;

        const scalar hStab =
            Foam::max(0.5 * (hTet + sqrt(sqr(hTet) + K)), SMALL * LSqrTri);

        func += LSqrTri / hStab;
    }

    return func;
}

scalar volumeOptimizerHeight::evaluateStabilisationFactor() const
{
    scalar K = 0.0;

    scalar Vmin(VGREAT), LSqMax(0.0);

    forAll(tets_, tetI)
    {
        const partTet& pt = tets_[tetI];
        const tetrahedron<point, point> tet
        (
            points_[pt.a()],
            points_[pt.b()],
            points_[pt.c()],
            points_[pt.d()]
        );

        const scalar Vtri = tet.mag();

        Vmin = Foam::min(Vmin, Vtri);

        const scalar LSqrTri
        (
            magSqr(tet.d() - tet.a()) +
            magSqr(tet.d() - tet.b()) +
            magSqr(tet.d() - tet.c())
        );

        LSqMax = Foam::max(LSqMax, LSqrTri);
    }

    if( Vmin < SMALL * LSqMax )
        K = max(ROOTVSMALL, SMALL * LSqMax);
    if( Vmin < VSMALL )
    {
        K = max(K, -10.0 * SMALL * Vmin);
        K = max(ROOTVSMALL, K);
    }

    return K;
}

void volumeOptimizerHeight::evaluateGradientsExact
(
    vector& gradF,
    tensor &gradGradF
) const
{
    gradF = vector::zero;
    gradGradF = tensor::zero;

    const scalar K = evaluateStabilisationFactor();

    tensor gradGradLsq(tensor::zero);
    gradGradLsq.xx() = 6.0;
    gradGradLsq.yy() = 6.0;
    gradGradLsq.zz() = 6.0;

    const point& p = points_[pointI_];

    forAll(tets_, tetI)
    {
        const partTet& pt = tets_[tetI];
        const tetrahedron<point, point> tet
        (
            points_[pt.a()],
            points_[pt.b()],
            points_[pt.c()],
            points_[pt.d()]
        );

        //- calculate the height and its gradients
        vector n = 0.5 * ((tet.b() - tet.a()) ^(tet.c() - tet.a()));
        n /= (mag(n) + VSMALL);

        const scalar h = (tet.d() - tet.a()) & n;

        const scalar stab = sqrt(sqr(h) + K);

        const scalar hStab = 0.5 * (h + stab);
        const scalar sqHstab = sqr(hStab);

        const vector gradStab = h * n / stab;

        const tensor gradGradHstab =
            (gradStab * gradStab / stab) * (1.0 - sqr(h / stab));

        const vector gradHstab = 0.5 * (n + gradStab);

        if( hStab < VSMALL )
        {
            //- do not consider this tet
            continue;
        }

        //- calculate the sum of edge lengths
        const scalar LSqrTri
        (
            magSqr(tet.d() - tet.a()) +
            magSqr(tet.d() - tet.b()) +
            magSqr(tet.d() - tet.c())
        );

        //- calculate the gradient of the Frobenius norm
        const vector gradLsq = 2. * (3. * p - tet.a() - tet.b() - tet.c());

        //- calculate the gradient of the functional
        gradF += gradLsq / hStab - LSqrTri * gradHstab / sqHstab;

        //- calculate the second gradient
        gradGradF +=
            gradGradLsq / hStab +
            2.0 * LSqrTri * (gradHstab * gradHstab) / (sqHstab * hStab) -
            (
                LSqrTri * gradGradHstab +
                gradLsq * gradHstab +
                gradHstab * gradLsq
            ) / sqHstab;
    }

    if( help::isnan(gradF) || help::isinf(gradF) )
        throw std::range_error("gradF is out of range");
    if( help::isnan(gradGradF) || help::isinf(gradGradF) )
        throw std::range_error("gradGradF is not within a valid range");
}

scalar volumeOptimizerHeight::optimiseDivideAndConquer(const scalar tol)
{
    point& pOpt = points_[pointI_];

    pOpt = 0.5 * (bb_.max() + bb_.min());
    point currCentre = pOpt;
    scalar dx = (bb_.max().x() - bb_.min().x()) / 2.0;
    scalar dy = (bb_.max().y() - bb_.min().y()) / 2.0;
    scalar dz = (bb_.max().z() - bb_.min().z()) / 2.0;

    label iter(0);

    //- find the value of the functional in the centre of the bnd box
    scalar funcBefore, funcAfter(evaluateFunc());

    do
    {
        funcBefore = funcAfter;

        funcAfter = VGREAT;
        point minCentre(vector::zero);

        for(label i=0;i<8;++i)
        {
            pOpt.x() = currCentre.x() + 0.5 * dirVecs[i].x() * dx;
            pOpt.y() = currCentre.y() + 0.5 * dirVecs[i].y() * dy;
            pOpt.z() = currCentre.z() + 0.5 * dirVecs[i].z() * dz;

            const scalar func = evaluateFunc();

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
        dz *= 0.5;

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

scalar volumeOptimizerHeight::optimiseSteepestDescent(const scalar tol)
{
    label iter(0);

    point& p = points_[pointI_];

    # ifdef DEBUGSmooth
    Info << nl << "Smoothing point " << pointI_
         << " with coordinates " << p << endl;
    scalar Vmina(VGREAT);
    forAll(tets_, tetI)
    Vmina = Foam::min(Vmina, tets_[tetI].mag(points_));
    Info << "Vmin before " << Vmina << endl;
    # endif

    vector gradF;
    vector disp(vector::zero);
    tensor gradGradF;
    point pOrig;

    scalar funcBefore, funcAfter(evaluateFunc());

    bool finished;
    do
    {
        finished = false;
        pOrig = p;
        funcBefore = funcAfter;

        evaluateGradientsExact(gradF, gradGradF);

        const scalar determinant = Foam::det(gradGradF);
        if( determinant > SMALL )
        {
            disp = (inv(gradGradF, determinant) & gradF);

            p -= disp;

            funcAfter = evaluateFunc();

            # ifdef DEBUGSmooth
            Info << nl << "gradF " << gradF << endl;
            Info << "gradGradF " << gradGradF << endl;
            Info << "det(gradGradF) " << determinant << endl;
            Info << "disp " << disp << endl;
            Info << "Func before " << funcBefore << endl;
            Info << "Func after " << funcAfter << endl;
            # endif

            scalar relax(0.8);
            label nLoops(0);

            while( funcAfter > funcBefore )
            {
                p = pOrig - relax * disp;
                relax *= 0.5;
                funcAfter = evaluateFunc();

                if( funcAfter < funcBefore )
                    continue;

                if( ++nLoops == 5 )
                {
                    //- it seems that this direction is wrong, stop the loop
                    p = pOrig;
                    disp = vector::zero;
                    finished = true;
                    funcAfter = funcBefore;
                }
            }

            if( mag(funcBefore - funcAfter) / funcBefore < tol )
                finished = true;
        }
        else
        {
            //- move in random direction
            //- this is usually needed to move the point of the zero volume
            disp = vector::zero;
            forAll(tets_, tetI)
            {
                const partTet& tet = tets_[tetI];
                const scalar Vtri = tet.mag(points_);

                if( Vtri < SMALL )
                {
                    triangle<point, point> tri
                    (
                        points_[tet.a()],
                        points_[tet.b()],
                        points_[tet.c()]
                    );

                    vector n = tri.normal();
                    const scalar d = mag(n);

                    if( d > VSMALL )
                        disp += 0.01 * (n / d);
                }
            }

            p += disp;
            funcAfter = evaluateFunc();
        }
    } while( (++iter < 100) && !finished );

    # ifdef DEBUGSmooth
    scalar Vmin(VGREAT);
    forAll(tets_, tetI)
        Vmin = Foam::min(Vmin, tets_[tetI].mag(points_));

    Info << nl << "New coordinates for point "
        << pointI_ << " are " << p << endl;
    Info << "Num iterations " << iter << " gradient " << gradF << endl;
    Info << "Vmin " << Vmin << endl;
    # endif

    return funcAfter;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
