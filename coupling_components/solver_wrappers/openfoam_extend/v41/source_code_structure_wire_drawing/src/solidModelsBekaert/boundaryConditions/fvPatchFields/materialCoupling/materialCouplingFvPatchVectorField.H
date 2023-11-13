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
    materialCouplingFvPatchVectorField

Description
    Region couple patch field

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

SourceFiles
    materialCouplingFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef materialCouplingFvPatchVectorField_H
#define materialCouplingFvPatchVectorField_H

#include "coupledFvPatchField.H"
#include "ggiLduInterfaceField.H"
#include "regionCoupleFvPatch.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class materialCouplingFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class materialCouplingFvPatchVectorField
:
    public ggiLduInterfaceField,
    public coupledFvPatchField<vector>
{
    // Private data

        //- Local reference cast into the regionCouple patch
        const regionCoupleFvPatch& regionCouplePatch_;

        //- Name of remote field to couple to
        word remoteFieldName_;

        //- Matrix update buffer
        mutable scalarField matrixUpdateBuffer_;

        //- Original patch field.  Required for correct evaluation
        //  in harmonic averaging
        mutable Field<vector> originalPatchField_;

        //- Current time index used to store originalPatchField_
        mutable label curTimeIndex_;


protected:

        //- Set remote field name
        void setRemoteFieldName(const word& name)
        {
            remoteFieldName_ = name;
        }

        //- Return non-const access to matrixUpdateBuffer
        void setMatrixUpdateBuffer(const scalarField& field) const
        {
            matrixUpdateBuffer_ = field;
        }

        //- Return original patch field
        const Field<vector>& originalPatchField() const;

        //- Return local patch-normal gradient
        tmp<Field<vector> > modifiedTraction() const;

        //- Return effective material propertie at the interface
        tmp<scalarField> iK() const;

public:

    //- Runtime type information
    TypeName("materialCoupling");


    // Constructors

        //- Construct from patch and internal field
        materialCouplingFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        materialCouplingFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given materialCouplingFvPatchVectorField
        //  onto a new patch
        materialCouplingFvPatchVectorField
        (
            const materialCouplingFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        materialCouplingFvPatchVectorField
        (
            const materialCouplingFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new materialCouplingFvPatchVectorField(*this)
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
                new materialCouplingFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access functions

            //- Return true if this patch field is derived from
            //  coupledFvPatchField<Type>.
            virtual bool coupled() const
            {
                return regionCouplePatch_.coupled();
            }
  
            //- Return remote field name
            const word& remoteFieldName() const
            {
                return remoteFieldName_;
            }

            //- Return contents of a matrix update buffer
            const scalarField& matrixUpdateBuffer() const
            {
                return matrixUpdateBuffer_;
            }

            //- Return a named shadow patch field
            template<class LookupField, class LookupType>
            const typename LookupField::PatchFieldType&
            lookupShadowPatchField
            (
                const word& name,
                const LookupField* = NULL,
                const LookupType* = NULL
            ) const;

            //- Return a shadow patch
            const regionCoupleFvPatch& regionCouplePatch() const
            {
                return regionCouplePatch_;
            }

            //- Return shadow patch field
            const materialCouplingFvPatchVectorField& shadowPatchField() const;


        // Evaluation functions

            //- Return neighbour field given internal cell data
            virtual tmp<Field<vector> > patchNeighbourField() const;

            //- Return named neighbour field given internal cell data
            virtual tmp<Field<vector> > patchNeighbourField
            (
                const word& name
            ) const;

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Evaluate the patch field
            virtual void evaluate(const Pstream::commsTypes commsType);

            // //- Update the coefficients associated with the patch field
            // virtual void updateCoeffs();

            //- Return patch-normal gradient
            virtual tmp<Field<vector> > snGrad() const;

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
                return
                    !(
                        regionCouplePatch_.parallel()
                     || pTraits<vector>::rank == 0
                    );
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return regionCouplePatch_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return regionCouplePatch_.reverseT();
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

#ifdef NoRepository
#   include "materialCouplingFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
