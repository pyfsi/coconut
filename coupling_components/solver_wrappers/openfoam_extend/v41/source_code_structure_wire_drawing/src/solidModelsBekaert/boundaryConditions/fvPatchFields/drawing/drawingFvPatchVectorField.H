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
    drawingFvPatchVectorField

Description
    The drawing component is a special case of the fixedDisplacementZeroShear
    condition, where the normal component of displacement is specified and the
    force/stress in the tangential/shear direction is zero.

    The velocity is ramped up over the velocity ramp time, where the velocity
    ramp time is optional and defaults to 0.

SourceFiles
    drawingFvPatchVectorField.C

Author
    Philip Cardiff, UCD, based on fixedDisplacementZeroShear

\*---------------------------------------------------------------------------*/

#ifndef drawingFvPatchVectorField_H
#define drawingFvPatchVectorField_H

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

class drawingFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{
    // Private data

        //- Velocity in patch normal direction
        const scalar velocity_;

        //- Velocity ramp time (defaults to 0)
        const scalar velocityRampTime_;

    // Private functions

        //- Set face displacement (or displacement increment) field
        void setFaceDisplacementField(vectorField& disp) const;

public:

    //- Runtime type information
    TypeName("drawing");


    // Constructors

        //- Construct from patch and internal field
        drawingFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        drawingFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // drawingFvPatchVectorField
        //  onto a new patch
        drawingFvPatchVectorField
        (
            const drawingFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new drawingFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        drawingFvPatchVectorField
        (
            const drawingFvPatchVectorField&,
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
                new drawingFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~drawingFvPatchVectorField()
        {}


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

            //- Return the average velocity
            scalar averageVelocity() const;


        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
