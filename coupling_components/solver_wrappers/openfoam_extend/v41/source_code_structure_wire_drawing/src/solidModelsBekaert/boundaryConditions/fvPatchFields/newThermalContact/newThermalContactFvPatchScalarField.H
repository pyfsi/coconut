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
    Foam::newThermalContactFvPatchScalarField

Description
    Foam::newThermalContactFvPatchScalarField

SourceFiles
    newThermalContactFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef newThermalContactFvPatchScalarField_H
#define newThermalContactFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "newGgiInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class newThermalContactFvPatch Declaration
\*---------------------------------------------------------------------------*/

class newThermalContactFvPatchScalarField
:
    // public mixedFvPatchScalarField
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Is this patch the master
        Switch master_;

        //- Shadow patch name
        word shadowPatchName_;

        //- Shadow patch index
        mutable label shadowPatchIndex_;

        //- Zone index
        mutable label zoneIndex_;

        //- Shadow zone index
        mutable label shadowZoneIndex_;

        //- Under-relaxation factor
        scalar underRelaxation_;

        //- Heat transfer coefficient for uncovered parts
        scalar alpha_;

        //- Surrounding temperature
        scalar Tinf_;

        //- Pressure sensitive factor for contact conductance
        scalar beta_;

        //- Reference contact resistance
        scalar Rc_;

        //- Approximate Ultimate Tensile Strength of softer material in contact
        // This is used in the calculation of the pressure dependent contact
        // conductance
        scalar UTS_;

        //- Heat generation at the contact surface
        scalarField Qc_;

        //- If true then Qc is calculated from dissipated frictional
        //  energy from mechanical contact
        bool useFrictionQc_;

        //- Displacement field name
        word DUName_;

        //- Current time index
        label curTimeIndex_;

        //- Zone-to-zone interpolation
        mutable newGgiZoneInterpolation* zoneToZonePtr_;

        //- Pressure sensitive contact conductivity
        mutable scalarField* HcPtr_;

        //- Contact area (0 - no contact, 1 - full contact, > 0 & < 1 - partial
        //  contact)
        mutable scalarField* contactPtr_;


    // Private member functions

        //- Return the shadow patch index
        label shadowPatchIndex() const;

        //- Set the shadow patch index
        void setShadowPatchIndex() const;

        //- Return zone index
        label zoneIndex() const;

        //- Set the zone index
        void setZoneIndex() const;

        //- Return shadow zone index
        label shadowZoneIndex() const;

        //- Set the shadow zone index
        void setShadowZoneIndex() const;

        //- Calculate interpolation
        void calcZoneToZone() const;

        //- Calculate contact indicator
        void calcContact() const;

        //- Calculate total contact conductivity
        void calcHc() const;

        //- Return thermal conductivity of a current patch
        const scalarField& curPatchK() const;

        //- Interpolate thermal conductivity of a shadow patch onto current
        //  patch
        const tmp<scalarField> shadowPatchKOnCurPatch() const;

        //- Return rho*c of a current patch
        const scalarField& curPatchRhoC() const;

        //- Interpolate rho*C of a shadow patch onto current  patch
        const tmp<scalarField> shadowPatchRhoCOnCurPatch() const;

        //- Interpolate temperature field of a shadow patch onto current
        //  patch
        const tmp<scalarField> shadowPatchTempOnCurPatch() const;

        //- Get zone field from patch field
        template<class Type>
        tmp<Field<Type> > zoneField
        (
            const label zoneIndex,
            const label patchIndex,
            const Field<Type>& patchField
        ) const;

        //- Get patch field from zone field
        template<class Type>
        tmp<Field<Type> > patchField
        (
            const label patchIndex,
            const label zoneIndex,
            const Field<Type>& zoneField
        ) const;

public:

    //- Runtime type information
    TypeName("newThermalContact");


    // Constructors

        //- Construct from patch and internal field
        newThermalContactFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        newThermalContactFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  newThermalContactFvPatchScalarField onto a new patch
        newThermalContactFvPatchScalarField
        (
            const newThermalContactFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        newThermalContactFvPatchScalarField
        (
            const newThermalContactFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new newThermalContactFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        newThermalContactFvPatchScalarField
        (
            const newThermalContactFvPatchScalarField&,
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
                new newThermalContactFvPatchScalarField(*this, iF)
            );
        }


    // Destructor

        virtual ~newThermalContactFvPatchScalarField();


    // Member functions

        // Access

            //- Return if this patch is the master
            Switch master() const
            {
                return master_;
            }

            //- Return the name of the shadow patch
            const word& shadowPatchName() const
            {
                return shadowPatchName_;
            }

            //- Return under-relaxation factor
            scalar underRelaxation() const;

            //- Return reference contact resistance
            scalar Rc() const;

            //- Return pressure sensitivity parameter for contact conductance
            scalar beta() const;

            //- Return Ultimate Tensile Strength
            scalar UTS() const;

            //- Return heat transfer coefficient
            scalar alpha() const
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

            //- Return surface friction generated heat
            tmp<scalarField> Qc() const;

            //- Return contact area
            const scalarField& contact() const;

            //- Return useFrictionQc_
            bool useFrictionQc() const;

            //- Return reference to zone-to-zone interpolation
            const newGgiZoneInterpolation& zoneToZone() const;


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

            //- Clear out demand driven data
            virtual void clearOut();


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

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

#ifdef NoRepository
#   include "newThermalContactFvPatchScalarFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
