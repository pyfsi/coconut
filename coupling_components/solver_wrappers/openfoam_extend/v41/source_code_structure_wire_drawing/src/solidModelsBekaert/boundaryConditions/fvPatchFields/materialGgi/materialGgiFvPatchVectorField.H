/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    materialGgiFvPatchVectorField

Description
    Generalized grid interface patch field.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

SourceFiles
    materialGgiFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef materialGgiFvPatchVectorField_H
#define materialGgiFvPatchVectorField_H

#include "ggiFvPatchField.H"
#include "simpleCohesiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class materialGgiFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class materialGgiFvPatchVectorField
:
    public ggiFvPatchField<vector>
{

    // Private data

    // Private functions

public:

    //- Runtime type information
    TypeName("materialGgi");


    // Constructors

        //- Construct from patch and internal field
        materialGgiFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        materialGgiFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given materialGgiFvPatchVectorField 
        //  onto a new patch
        materialGgiFvPatchVectorField
        (
            const materialGgiFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        materialGgiFvPatchVectorField
        (
            const materialGgiFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );


        //- Construct as copy
        materialGgiFvPatchVectorField
        (
            const materialGgiFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new materialGgiFvPatchVectorField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new materialGgiFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return neighbour field given internal cell data
            template<class Type>
            tmp<Field<Type> > patchNeighbourField
            (
                const Field<Type>& iField
            ) const;

            //- Return neighbour field given internal cell data
            template<class Type>
            tmp<Field<Type> > patchNeighbourField
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& sf
            ) const;

        // Evaluation functions

            //- Return patch-normal gradient
            virtual tmp<Field<vector> > snGrad() const;

            //- Return local patch-normal gradient
            tmp<Field<tensor> > localGrad() const;

            //- Return local patch-normal gradient
            tmp<Field<vector> > localSnGrad() const;

            //- Return local patch-normal gradient
            tmp<Field<vector> > modifiedTraction() const;

            //- Return effective material propertie at the interface
            tmp<scalarField> iK() const;

            //- Return neighbour field given internal cell data
            virtual tmp<Field<vector> > patchNeighbourField() const;

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Calculate patch face interpolate given weighting factors
            virtual void patchInterpolate
            (
                GeometricField<vector, fvsPatchField, surfaceMesh>& fField,
                const scalarField& pL
            ) const;

            //- Calculate patch face interpolate given two weighting factors
            virtual void patchInterpolate
            (
                GeometricField<vector, fvsPatchField, surfaceMesh>& fField,
                const scalarField& pL,
                const scalarField& pY
            ) const;

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "materialGgiFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
