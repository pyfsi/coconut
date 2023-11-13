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
    materialInterface

Description
    Generalized grid interface patch field.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

    Adapted by Philip Cardiff to remove coheisve law

SourceFiles
    materialInterface.C

\*---------------------------------------------------------------------------*/

#ifndef materialInterface_H
#define materialInterface_H

#include "ggiFvPatchField.H"
//#include "simpleCohesiveLaw.H"
// #include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class materialInterface Declaration
\*---------------------------------------------------------------------------*/

class materialInterfaceFvPatchVectorField
:
    public ggiFvPatchField<vector>
{

    // Private data

        //- Cohesive law
        //simpleCohesiveLaw* cohesiveLawPtr_;

        //- Whether this patch is master
        //scalarField debonded_;

        //- Displacement jump
        //vectorField displacementJump_;

        //- Traction
        //vectorField traction_;

        //- Normal traction
        //scalarField normalTraction_;

        //- Current traction
        //vectorField initiationTraction_;

        //- Current normal gradient
        vectorField grad_;

        //- Current time index
        label curTimeIndex_;

        //- Iteration counter
        label iCorr_;

        //- Iteration counter
        //bool jumpConverged_;

        //- Number of linear solver iterations
        mutable labelList nIter_;

    // Private functions

public:

    //- Runtime type information
    TypeName("materialInterface");


    // Constructors

        //- Construct from patch and internal field
        materialInterfaceFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        materialInterfaceFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given materialInterface 
        //  onto a new patch
        materialInterfaceFvPatchVectorField
        (
            const materialInterfaceFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        materialInterfaceFvPatchVectorField
        (
            const materialInterfaceFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );


        //- Construct as copy
        materialInterfaceFvPatchVectorField
        (
            const materialInterfaceFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new materialInterfaceFvPatchVectorField(*this)
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
                new materialInterfaceFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return cohesive law
            //const simpleCohesiveLaw& law() const;

            //- Return true if face debonded
            // const scalarField& debonded() const
            // {
            //     return debonded_;
            // }

            //- Return reference to traction at the interface
            // const vectorField& initiationTraction() const
            // {
            //     return initiationTraction_;
            // }

            //- Return reference to traction at the interface
            // const vectorField& traction() const
            // {
            //     return traction_;
            // }

            //- Return reference to traction at the interface
            // const scalarField& normalTraction() const
            // {
            //     return normalTraction_;
            // }

            //- Return reference to displacement jump at the interface
            // const vectorField& displacementJump() const
            // {
            //     return displacementJump_;
            // }

            //- Return reference to displacement jump at the interface
            // vectorField& displacementJump()
            // {
            //     return displacementJump_;
            // }

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
            tmp<Field<vector> > localSnGrad() const;

//             //- Return neighbour field given internal cell data
//             tmp<Field<vector> > patchInternalField
//             (
//                 const vectorField& internalField
//             ) const;

            //- Return neighbour field given internal cell data
            virtual tmp<Field<vector> > patchNeighbourField() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

//             //- Evaluate the patch field
//             virtual void evaluate
//             (
//                 const Pstream::commsTypes commsType
//             );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector> > gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<vector> > gradientBoundaryCoeffs() const;

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

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<vector>& matrix);

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "materialInterfaceFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
