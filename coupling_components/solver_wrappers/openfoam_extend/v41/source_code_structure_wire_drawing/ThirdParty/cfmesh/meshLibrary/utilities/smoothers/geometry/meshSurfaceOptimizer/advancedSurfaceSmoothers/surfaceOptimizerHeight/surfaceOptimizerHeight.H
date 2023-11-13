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
    surfaceOptimizerHeight

Description
    Performs optimisation of a central vertex in a simplex

SourceFiles
    surfaceOptimizerHeight.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceOptimizerHeight_H
#define surfaceOptimizerHeight_H

#include "point.H"
#include "triFace.H"
#include "DynList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class surfaceOptimizerHeight Declaration
\*---------------------------------------------------------------------------*/

class surfaceOptimizerHeight
{
    // Private static data
        //- direction vectors for divide and conquer algorithm
        static const vector dirVecs[4];

        //- second gradient of the sum of squares of edge length
        static const tensor gradGradLsq_;

    // Private data
        //- reference to the simplex points
        DynList<point>& pts_;

        //- reference to the triangles forming a simplex
        const DynList<triFace>& trias_;

        //- weights for each triangle
        DynList<scalar> weights_;

        //- min position of the bnd box
        point pMin_;

        //- max position of the bnd box
        point pMax_;

    // Private member functions
        //- evaluate stabilisation factor
        scalar evaluateStabilisationFactor() const;

        //- evaluate the functional
        scalar evaluateFunc(const scalar& K) const;

        //- evaluate gradients needed for optimisation
        void evaluateGradients(const scalar&, vector&, tensor&) const;

        //- optimise point position using the divide and conquer technique
        scalar optimiseDivideAndConquer(const scalar tol);

        //- optimise point position via the steepest descent method
        scalar optimiseSteepestDescent(const scalar tol);

        //- Disallow default bitwise copy construct
        surfaceOptimizerHeight(const surfaceOptimizerHeight&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceOptimizerHeight&);

public:

    // Constructors

        //- Construct from transformed points and triangles forming a simplex
        surfaceOptimizerHeight
        (
            DynList<point>& pts,
            const DynList<triFace>& trias,
            const DynList<scalar>& weights = DynList<scalar>()
        );


    // Destructor

        ~surfaceOptimizerHeight();

    // Member Functions

        //- optimizes position of a central point in the simplex
        point optimizePoint(const scalar tol = 0.1);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
