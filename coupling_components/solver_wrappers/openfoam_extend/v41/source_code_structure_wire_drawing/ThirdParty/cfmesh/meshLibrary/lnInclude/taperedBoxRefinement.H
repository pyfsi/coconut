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
    taperedBoxRefinement

Description
    Checks if a box is contained inside the box object

SourceFiles
    taperedBoxRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef taperedBoxRefinement_H
#define taperedBoxRefinement_H

#include "objectRefinement.H"
#include "point.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class taperedBoxRefinement Declaration
\*---------------------------------------------------------------------------*/

class taperedBoxRefinement
:
    public objectRefinement
{
    // Private data
        //- centre of the box
        point centre_;

        //- length of box sides
        scalar lengthX_;
        scalar lengthY_;
        scalar lengthZ_;

        //- tapper angles in each direction
        scalar angleXY_;
        scalar angleXZ_;
        scalar angleYX_;
        scalar angleYZ_;
        scalar angleZX_;
        scalar angleZY_;

        //- helper data to speed-up inside evaluations
        scalar sinXY_;
        scalar sinXZ_;
        scalar sinYX_;
        scalar sinYZ_;
        scalar sinZX_;
        scalar sinZY_;

        FixedList<scalar, 6> fCentres_;

    // Private member functions
        //- calculate helper data
        void calculateHelperData();

public:

    //- Runtime type information
    TypeName("taperedBox");


    // Constructors

        //- Null construct
        taperedBoxRefinement();

        //- Construct from name, cell size, centre and sizes
        taperedBoxRefinement
        (
            const word& name,
            const scalar cellSize,
            const direction additionalRefLevels,
            const point& centre,
            const scalar lengthX,
            const scalar lengthY,
            const scalar lengthZ,
            const scalar angleXY,
            const scalar angleXZ,
            const scalar angleYX,
            const scalar angleYZ,
            const scalar angleZX,
            const scalar angleZY
        );

        //- Construct from dictionary
        taperedBoxRefinement(const word& name, const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<objectRefinement> clone
        (
            const taperedBoxRefinement& sr
        ) const
        {
            return autoPtr<objectRefinement>
            (
                new taperedBoxRefinement
                (
                    sr.name(),
                    sr.cellSize(),
                    sr.additionalRefinementLevels(),
                    sr.centre_,
                    sr.lengthX_,
                    sr.lengthY_,
                    sr.lengthZ_,
                    sr.angleXY_,
                    sr.angleXZ_,
                    sr.angleYX_,
                    sr.angleYZ_,
                    sr.angleZX_,
                    sr.angleZY_
                )
            );
        }

    // Member Functions

        //- check if a boundBox intersects or is inside the object
        bool intersectsObject(const boundBox&) const;

    //- Return as dictionary of entries
        dictionary dict(bool ignoreType = false) const;

      // Write

        //- Write
        void write(Ostream&) const;

        //- Write dictionary
        void writeDict(Ostream&, bool subDict = true) const;

    // Member Operators

        //- assign from dictionary
        void operator=(const dictionary&);

    // IOstream Operators

        Ostream& operator<<(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
