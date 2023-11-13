/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    coupledSolidTractionFvPatchVectorField

Description
    Region couple patch field

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    coupledSolidTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledSolidTractionFvPatchVectorField_H
#define coupledSolidTractionFvPatchVectorField_H

#include "coupledFvPatchField.H"
#include "ggiLduInterfaceField.H"
#include "regionCoupleFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class coupledSolidTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class coupledSolidTractionFvPatchVectorField
:
    public ggiLduInterfaceField,
    public coupledFvPatchField<vector>
{
    // Private data

        //- Local reference cast into the regionCouple patch
        //const regionCoupleFvPatch& regionCouplePatch_;

        //- Name of remote field to couple to
        //word remoteFieldName_;

        //- Matrix update buffer
        mutable scalarField matrixUpdateBuffer_;

        //- Original patch field.  Required for correct evaluation
        //  in harmonic averaging
        //mutable Field<Type> originalPatchField_;

        //- Current time index used to store originalPatchField_
        //mutable label curTimeIndex_;


// protected:

//         //- Set remote field name
//         void setRemoteFieldName(const word& name)
//         {
//             remoteFieldName_ = name;
//         }

//         //- Return non-const access to matrixUpdateBuffer
//         void setMatrixUpdateBuffer(const scalarField& field) const
//         {
//             matrixUpdateBuffer_ = field;
//         }

//         //- Return original patch field
//         const Field<Type>& originalPatchField() const;


public:

    //- Runtime type information
    TypeName("coupledSolidTraction");


    // Constructors

        //- Construct from patch and internal field
        coupledSolidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coupledSolidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given coupledSolidTractionFvPatchVectorField
        //  onto a new patch
        coupledSolidTractionFvPatchVectorField
        (
            const coupledSolidTractionFvPatchVectorField<vector>&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        coupledSolidTractionFvPatchVectorField
        (
            const coupledSolidTractionFvPatchVectorField<vector>&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new coupledSolidTractionFvPatchVectorField<vector>(*this)
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
                new coupledSolidTractionFvPatchVectorField<vector>(*this, iF)
            );
        }


    // Member functions

        // Access functions

            //- Return remote field name
            // const word& remoteFieldName() const
            // {
            //     return remoteFieldName_;
            // }

            //- Return contents of a matrix update buffer
            const scalarField& matrixUpdateBuffer() const
            {
                return matrixUpdateBuffer_;
            }

            //- Return a named shadow patch field
            // template<class LookupField, class LookupType>
            // const typename LookupField::PatchFieldType&
            // lookupShadowPatchField
            // (
            //     const word& name,
            //     const LookupField* = NULL,
            //     const LookupType* = NULL
            // ) const;

            // //- Return a shadow patch
            // const regionCoupleFvPatch& regionCouplePatch() const
            // {
            //     return regionCouplePatch_;
            // }

            //- Return shadow patch field
            //const coupledSolidTractionFvPatchVectorField<Type>&
            //shadowPatchField() const;


        // Evaluation functions

            //- Return neighbour field given internal cell data
            virtual tmp<Field<Type> > patchNeighbourField() const;

            //- Return named neighbour field given internal cell data
            virtual tmp<Field<Type> > patchNeighbourField
            (
                const word& name
            ) const;

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return patch-normal gradient
            virtual tmp<Field<Type> > snGrad() const;


            //- Transform given patch component field
            virtual void transformCoupleField
            (
                scalarField& f,
                const direction cmpt
            ) const
            {
                ggiLduInterfaceField::transformCoupleField(f, cmpt);
            }

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

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;


        // GGI coupled interface functions

            //- Does the patch field perform the transfromation
            virtual bool doTransform() const
            {
                // return
                //     !(
                //         regionCouplePatch_.parallel()
                //      || pTraits<Type>::rank == 0
                //     );

                notImplemented("doTransform");

                return false;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                //return regionCouplePatch_.forwardT();

                notImplemented("forwardT");

                return tensorField();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                //return regionCouplePatch_.reverseT();

                notImplemented("reverseT");

                return tensorField();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return pTraits<vector>::rank;
            }


        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "coupledSolidTractionFvPatchVectorField.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
