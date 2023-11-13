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
    fixedTangentialNormalPressureFvPatchVectorField

Description
    The component of the displacement tangential to the patch is fixed, and the
    normal pressure is specified as fixed or, optionally, the as a time series.
    Patch contains non-orthogonal correction.

SourceFiles
    fixedTangentialNormalPressureFvPatchVectorField.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef fixedTangentialNormalPressureFvPatchVectorField_H
#define fixedTangentialNormalPressureFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class fixedTangentialNormalPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedTangentialNormalPressureFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{
    // Private data

        //- Fixed pressure
        //  Note: if pressure is given as a spatially varying field, it is not
        //  mapped/advected when there is ALE/mesh smoothing
        scalarField pressure_;

        //- Normal direction
        vector normal_;

        //- Or optionally pressure specified as a time series
        interpolationTable<scalar> pressureTimeSeries_;

        //- Update average tangential displacement to reduce the average shear
        //  stress on the patch
        Switch updateTangentialDisplacement_;

        //- The updateTangentialDisplacement is set to false once the specified
        //  time is reached
        const scalar updateTangentialDisplacementEndTime_;

        //- Under-relaxation factor for apply updateTangentialDisplacement
        const scalar updateTangentialDisplacementRelaxFactor_;

public:

    //- Runtime type information
    TypeName("fixedTangentialNormalPressure");


    // Constructors

        //- Construct from patch and internal field
        fixedTangentialNormalPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedTangentialNormalPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // fixedTangentialNormalPressureFvPatchVectorField
        //  onto a new patch
        fixedTangentialNormalPressureFvPatchVectorField
        (
            const fixedTangentialNormalPressureFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new fixedTangentialNormalPressureFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedTangentialNormalPressureFvPatchVectorField
        (
            const fixedTangentialNormalPressureFvPatchVectorField&,
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
                new fixedTangentialNormalPressureFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~fixedTangentialNormalPressureFvPatchVectorField()
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


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
