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
    Foam::thermalContactFvPatchScalarField

Description
    Foam::thermalContactFvPatchScalarField

SourceFiles
    thermalContactFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef thermalContactFvPatchScalarField_H
#define thermalContactFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "solidContactFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class thermalContactFvPatch Declaration
\*---------------------------------------------------------------------------*/

class thermalContactFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Is this patch the master
        Switch master_;

        // Copy of contact dictionary
        const dictionary dict_;

        //- Under-relaxation factor
        scalar underRelaxation_;

        //- Heat transfer coefficient (outside the contact region)
        scalarField alpha_;

        //- Surrounding temperature
        scalar Tinf_;

        //- Reference contact resistance
        scalar Rc_;

        //- Pressure sensitive factor for contact conductance
        scalar beta_;

        //- Approximate Ultimate Tensile Strength of softer material in contact
        // This is used in the calculation of the pressure dependent contact
        // conductance
        scalar UTS_;

        //- Displacement field name
        word DUName_;

        //- Current time index
        label curTimeIndex_;

        //- Zone-to-zone interpolation
        // mutable PtrList<newGgiZoneInterpolation> zoneToZones_;


    // Private member functions

        //- Check that thermalContact and solidContact both use the same master
        //  patch
        void checkConsistentMaster() const;

        //- Return the shadow patch
        //  This function is not allowed if shadowPatchIndices.size() > 1
        //  i.e. it is always OK for the slave(s) to call this function but not
        //  for the master when there are multiple slaves
        const thermalContactFvPatchScalarField& shadowPatchField() const;

        //- Return the shadow patch for shadowI
        //  All slaves will only have one shadow whereas the master can have
        //  many
        const thermalContactFvPatchScalarField& shadowPatchField
        (
            const label shadI
        ) const;

        //- Return the solidContact patch field on this patch
        const solidContactFvPatchVectorField& solidContactPatch() const;


public:

    //- Runtime type information
    TypeName("thermalContact");


    // Constructors

        //- Construct from patch and internal field
        thermalContactFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermalContactFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given thermalContactFvPatchScalarField
        //  onto a new patch
        thermalContactFvPatchScalarField
        (
            const thermalContactFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        thermalContactFvPatchScalarField
        (
            const thermalContactFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new thermalContactFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        thermalContactFvPatchScalarField
        (
            const thermalContactFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new thermalContactFvPatchScalarField(*this, iF)
            );
        }


    // Destructor

        virtual ~thermalContactFvPatchScalarField();

    // Member functions

        // Access

            //- Return if this patch is the master
            Switch master() const
            {
                return master_;
            }

            //- Return the shadow patch names
            const wordList& shadowPatchNames() const;

            //- Return the shadow patch indices
            const labelList& shadowPatchIndices() const;

            //- Return zone index
            // label zoneIndex() const;

            //- Return index of the shadow zones
            // const labelList& shadowZoneIndices() const;

            //- Return under-relaxation factor
            scalar underRelaxation() const;

            //- Return contact reference resistance
            scalar Rc() const;

            //- Return pressure sensitivity parameter for contact conductance
            scalar beta() const;

            //- Return Ultimate Tensile Strength
            scalar UTS() const;

            //- Return heat transfer coefficient
            const scalarField& alpha() const
            {
                return alpha_;
            }

            //- Return surrounding temperature
            scalar Tinf() const
            {
                return Tinf_;
            }

            //- Return pressure sensitive contact conductance
            tmp<scalarField> Hc() const;

            //- Return friction flux for this patch
            //  This partitions the frictionFluxRate from solidContact into the
            //  master contribution and slave contribution
            //  Note: this function was formerly known as Qc() and we renamed it
            //  to be more clear
            tmp<scalarField> frictionFluxRateForThisPatch() const;

            //- Return const reference to zone-to-zone interpolation
            // const PtrList<newGgiZoneInterpolation>& zoneToZones() const;

            //- Return reference to zone-to-zone interpolation
            // PtrList<newGgiZoneInterpolation>& zoneToZones();

            //- Return const reference to zone-to-zone interpolation for the
            //  current slave patch; an error is thrown if the master calls this
            //  function
            // const newGgiZoneInterpolation& zoneToZoneForThisSlave() const;

            //- Return contact indicator
            const scalarField& contact() const;

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Return gradient at boundary
            //virtual tmp<Field<scalar> > snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "thermalContactFvPatchScalarFieldTemplates.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
