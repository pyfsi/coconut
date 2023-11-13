/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::timeVaryingVelocityFvPatchVectorField

Description
    Time varying velocity boundary condition with boundary
    non-orthogonal correction for stress solvers.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            timeVaryingVelocity;
        fileName        "$FOAM_CASE/time-series";
        outOfBounds     clamp;           // (error|warn|clamp|repeat)
    }
    @endverbatim

SourceFiles
    timeVaryingVelocityFvPatchVectorField.C

Authro
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingVelocityFvPatchVectorField_H
#define timeVaryingVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class timeVaryingVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingVelocityFvPatchVectorField
:
    public fixedDisplacementFvPatchVectorField
{

    // Private Data

        //- The time series being used, including the bounding treatment
        interpolationTable<vector> timeSeries_;

public:

    //- Runtime type information
    TypeName("timeVaryingVelocity");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // timeVaryingVelocityFvPatchVectorField
        //  onto a new patch
        timeVaryingVelocityFvPatchVectorField
        (
            const timeVaryingVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingVelocityFvPatchVectorField
        (
            const timeVaryingVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingVelocityFvPatchVectorField
        (
            const timeVaryingVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingVelocityFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return the time series used
            const interpolationTable<vector>& timeSeries() const
            {
                return timeSeries_;
            }

        // Evaluation functions

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
