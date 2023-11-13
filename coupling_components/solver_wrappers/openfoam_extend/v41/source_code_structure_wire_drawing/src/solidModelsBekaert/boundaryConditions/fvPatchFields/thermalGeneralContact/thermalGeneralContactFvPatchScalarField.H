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
    thermalGeneralContactFvPatchScalarField

Description
    General thermal contact boundary conditions for finite volume solidMechanics
    solvers, for heat conduction between contacting solids.

    Multiple patches may be defined as thermalGeneralContact and master-slave
    contact pairs will be created between each.

    This boundary condition mimics the solidGeneralContact

SourceFiles
    thermalGeneralContactFvPatchScalarField.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef thermalGeneralContactFvPatchScalarField_H
#define thermalGeneralContactFvPatchScalarField_H

// #include "mixedFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
//#include "newGgiInterpolation.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class thermalGeneralContactFvPatchScalarField
:
    // public mixedFvPatchScalarField
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the displacement field
        const word fieldName_;

        //- Is this the global master patch: the global master patches keeps a
        // copy of mesh.allPoints() and updates them to the deformed position
        // each time before the contact is corrected
        mutable bool* globalMasterPtr_;

        //- Patch index of global master patch
        mutable label* globalMasterIndexPtr_;

        //- If the current patch is the local slave of the contact pair
        // The local slave is the patch with the lower patch index
        mutable boolList* localSlavePtr_;

        //- List of shadow patch names
        mutable List<word>* shadowPatchNamesPtr_;

        //- List of shadow patch indices
        mutable List<label>* shadowPatchIndicesPtr_;

        //- Zone index i.e. index of global face zone corresponding to the
        // current patch
        label zoneIndex_;

        //- List of shadow zone names
        mutable List<word>* shadowZoneNamesPtr_;

        //- List of shadow zone indices
        mutable List<label>* shadowZoneIndicesPtr_;

        //- Under-relaxation factor
        scalar underRelaxation_;

        //- Thermal conductivity field name
        word thermalConductivityName_;

        //- Heat transfer coefficient for uncovered parts
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

        //- Heat generation at the contact surface
        scalarField Qc_;

        //- If true then Qc is calculated from dissipated frictional
        //- energy from mechanical contact
        bool useFrictionQc_;

        //- Displacement increment field name
        word DUName_;

        //- Current time index
        label curTimeIndex_;

        //- Pressure sensitivity contact conductivity
        mutable scalarField* HcPtr_;

        //- Contact indicator (0 - no contact, 1 - contact)
        mutable scalarField* contactPtr_;

        //- The contact indicator for each individual contact pair
        mutable List<scalarField>* contactsPtr_;


    // Private Member Functions

            //- Check global master patch
            void calcGlobalMaster() const;

            //- Check global master patch index
            void calcGlobalMasterIndex() const;

            //- Calculate local slaves
            void calcLocalSlave() const;

            //- Return local slaves
            const boolList& localSlave() const;

            //- Calculate shadow patch names
            void calcShadowPatchNames() const;

            //- Calculate shadow zone names
            void calcShadowZoneNames() const;

            //- Calculate interpolation
            void calcZoneToZones() const;

            //- Calculate contact indicator
            void calcContact() const;

            //- Calculate contact indicator for each individual pair
            void calcContacts() const;

            //- Calculate pressure sensitive contact conductivity
            void calcHc() const;

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
    TypeName("thermalGeneralContact");


    // Constructors

        //- Construct from patch and internal field
        thermalGeneralContactFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermalGeneralContactFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given thermalGeneralContactFvPatchScalarField
        //  onto a new patch
        thermalGeneralContactFvPatchScalarField
        (
            const thermalGeneralContactFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        thermalGeneralContactFvPatchScalarField
        (
            const thermalGeneralContactFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new thermalGeneralContactFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        thermalGeneralContactFvPatchScalarField
        (
            const thermalGeneralContactFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new thermalGeneralContactFvPatchScalarField(*this, iF)
            );
        }


    // Destructor

        virtual ~thermalGeneralContactFvPatchScalarField();


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return true;
            }

            //- Return true if this patch is the global master
            bool globalMaster() const;

            //- Return index of global master patch
            label globalMasterIndex() const;

            //- Return the list of shadow patch names
            const List<word>& shadowPatchNames() const;

            //- Return the list of shadow patch indices
            const List<label>& shadowPatchIndices() const;

            //- Return index of the zone
            label zoneIndex() const
            {
                return zoneIndex_;
            }

            //- Return the list of shadow zone names
            const List<word>& shadowZoneNames() const;

            //- Return the list of shadow zone indices
            const List<label>& shadowZoneIndices() const;

            //- Return under-relaxation factor
            scalar underRelaxation() const
            {
                return underRelaxation_;
            }

            //- Return the name of the thermal conductivity field
            const word& thermalConductivityName() const
            {
                return thermalConductivityName_;
            }

            //- Return contact reference resistance
            scalar Rc() const
            {
                return Rc_;
            }

            //- Return pressure sensitivity parameter for contact conductance
            scalar beta() const
            {
                return beta_;
            }

            //- Return Ultimate Tensile Strength
            scalar UTS() const
            {
                return UTS_;
            }

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
            const scalarField& Hc() const;

            //- Return contact resistance
            tmp<scalarField> Qc() const;

            //- Return const reference to patch-to-patch interpolation
            // const extendedGgiZoneInterpolation& zoneToZone
            // (
            //     const label shadowI
            // ) const;

            //- Return reference to patch-to-patch interpolation
            //extendedGgiZoneInterpolation& zoneToZone(const label shadowI);

            //- Return contact indicator
            const scalarField& contact() const;

            //- Return contact indicator for shadowI
            //tmp<scalarField> contact(const label shadowI) const;
            const scalarField& contact(const label shadowI) const;

            //- Given the patch index of a shadow, return the position in the
            // shadow list of the current patch
            label findShadowID(const label patchID) const;

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
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

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            //virtual tmp<Field<scalar> > gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            //virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "thermalGeneralContactFvPatchScalarFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
