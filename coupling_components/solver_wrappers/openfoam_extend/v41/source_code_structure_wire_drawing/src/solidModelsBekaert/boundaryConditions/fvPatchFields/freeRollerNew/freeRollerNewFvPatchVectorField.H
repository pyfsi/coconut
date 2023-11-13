/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    freeRollerNewFvPatchVectorField

Description
    A free roller boundary condition.

    A rotational source term is added to the roller cells to artificially
    accelerate the roller to the approximate RPM; this allows us to more quickly
    jump to steady state rolling.

SourceFiles
    freeRollerNewFvPatchVectorField.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef freeRollerNewFvPatchVectorField_H
#define freeRollerNewFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class freeRollerNewFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{
    // Private data

        //- Name of the displacement field
        const word fieldName_;

        //- Displacement to compress wire
        // This is increased each pass
        const vector displacement_;

        //- Displacement is linearly applied over dispRampTime_ time
        scalar dispRampTime_;

        //- Rotation axis
        //  Displacement component in this direction will be zero
        vector rotationAxis_;

        //- Approximate time to apply acceleration
        const scalar accelerationTime_;

        //- Initial origin of rotation
        vector initialRotationOrigin_;

        //- The current axis displacement
        //  This may be looked up by e.g. function objects
        vector currentAxisDisplacement_;

        //- The total axis translation displacement has been reached
        bool totalAxisDisplacementReached_;

        //- Name of the cell zone to be accelerated
        const word zoneName_;

        //- Approximate RPM to be set
        const scalar rpm_;

        //- Current time index
        label curTimeIndex_;

    // Private member functions

        //- Return the cell zone index
        label zoneID() const;

public:

    //- Runtime type information
    TypeName("freeRollerNew");


    // Constructors

        //- Construct from patch and internal field
        freeRollerNewFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        freeRollerNewFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // freeRollerNewFvPatchVectorField
        //  onto a new patch
        freeRollerNewFvPatchVectorField
        (
            const freeRollerNewFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new freeRollerNewFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        freeRollerNewFvPatchVectorField
        (
            const freeRollerNewFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new freeRollerNewFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~freeRollerNewFvPatchVectorField();


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
