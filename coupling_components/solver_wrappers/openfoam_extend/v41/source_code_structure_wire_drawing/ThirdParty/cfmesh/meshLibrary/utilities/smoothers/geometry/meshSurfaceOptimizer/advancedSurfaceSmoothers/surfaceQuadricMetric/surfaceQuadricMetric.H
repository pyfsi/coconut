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
    surfaceQuadricMetric

Description
    Performs optimisation of a central vertex in a simplex

SourceFiles
    surfaceQuadricMetric.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceQuadricMetric_H
#define surfaceQuadricMetric_H

#include "DynList.H"
#include "boundBox.H"
#include "point.H"
#include "triFace.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class partTetMeshSimplex;

/*---------------------------------------------------------------------------*\
                    Class surfaceQuadricMetric Declaration
\*---------------------------------------------------------------------------*/

class surfaceQuadricMetric
{
    // Private data
        //- reference to the simplex points
        DynList<point>& pts_;

        //- reference to the triangles forming a simplex
        const DynList<triFace>& trias_;

        //- free vertex
        point& p_;

        //- normals of edges forming the outer hull
        DynList<vector> normals_;

        //- centres of edges forming the outer hull
        DynList<point> centres_;

        //- bounding box
        boundBox bb_;

        //- stability constant
        scalar Ksq_;

    // Private member functions
        //- evaluate the value of the metric
        scalar evaluateMetric() const;

        //- evaluate metric gradients
        void evaluateGradients(vector& grad, tensor& gradGrad) const;

    public:

    // Constructor
        surfaceQuadricMetric
        (
            DynList<point>& pts,
            const DynList<triFace>& trias
        );

    // Destructor
        ~surfaceQuadricMetric();

    // Member functions
        //- find the best position such that all tets making a simplex
        //- have a positive volume
        point optimizePoint(const scalar tol = 0.01);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
