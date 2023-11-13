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
    boundaryPatch

Description
    Implementation of poly line

SourceFiles
    boundaryPatch.C

\*---------------------------------------------------------------------------*/

#ifndef polyLine_H
#define polyLine_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "splineBase.H"
#include "typeInfo.H"
#include "scalarField.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyLine
:   public splineBase
{
    // Private data
        //- interval coefficient
        scalarField startingCoeffs_;

    // Private member functions
        //- calculate coefficients
        void calculateCoefficients();
public:

    //- Runtime type information
    TypeName("polyLine");

    // Constructors
        //- Null construct
        polyLine();

        //- construct from control points and the type
        polyLine(const LongList<point>& points, const word& type);

    // Return clone of the object
        virtual autoPtr<splineBase> clone() const
        {
            return autoPtr<splineBase>
            (
                new polyLine
                (
                    this->points(),
                    this->type()
                )
            );
        }

    // Public member functions
        //- evaluate a point on the spline for the given parameter
        virtual point evaluate(const scalar t) const;

        //- create a polyLine that differs from the spline less
        //- than the given tolerance
        virtual void createPolyLine
        (
            const scalar tol,
            LongList<point>&
        ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
