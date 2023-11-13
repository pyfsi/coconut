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
    solidGeneralContactFvPatchVectorField

Description
    General contact boundary conditions for finite volume solidMechanics
    solvers.

    Multiple patches may be defined as solidGeneralContact and master-slave
    contact pairs will be created between each.

    As calculation of point distances is expensive, an initial quick check will
    be performed to check if a contact pair is active.

SourceFiles
    solidGeneralContactFvPatchVectorField.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef solidGeneralContactFvPatchVectorField_H
#define solidGeneralContactFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"
#include "normalContactModel.H"
#include "frictionContactModel.H"
#include "newGgiInterpolation.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class solidGeneralContactFvPatchVectorField
:
        public solidTractionFvPatchVectorField
{
    // Private data

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
        mutable label zoneIndex_;

        //- List of shadow zone names
        mutable List<word>* shadowZoneNamesPtr_;

        //- List of shadow zone indices
        mutable List<label>* shadowZoneIndicesPtr_;

        //- Is the current patch considered rigid in contact
        Switch rigidMaster_;

        // Copy of contact dictionary
        const dictionary dict_;

        //- List of normal contact models
        mutable PtrList<normalContactModel> normalModels_;

        //- List of friction contact models
        mutable PtrList<frictionContactModel> frictionModels_;

        // Master face zone kept in deformed configuration
        mutable standAlonePatch* zonePtr_;

        //- List of zone-to-zone interpolators
        mutable PtrList<newGgiStandAlonePatchInterpolation> zoneToZones_;

        //- Patch to patch algorithm
        const intersection::algorithm alg_;

        //- Patch to patch direction
        const intersection::direction dir_;

        //- Current time index
        label curTimeIndex_;

        //- Global face zones in the mesh
        mutable labelList* globalFaceZonesPtr_;

        // Point map from global to local for global face zones
        mutable IOList<labelList>* globalToLocalFaceZonePointMapPtr_;

        //- The current patch tractions for each contatc pair
        mutable List<vectorField>* curPatchTractionPtr_;

        // Rate of friction energy dissipation
        mutable scalarField* QcPtr_;

        // Rate of friction energy dissipation for each contact pair
        mutable List<scalarField>* QcsPtr_;

        //- Offset used if bounding box of patch has a zero dimension
        mutable scalar bbOffset_;

    // Private Member Functions

        //- Is a moving mesh (e.g. updated Lagrangian) approach be used
        bool movingMesh() const;

        //- Move the master and slave zones to the deformed configuration
        void moveZonesToDeformedConfiguration();

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

        //- Create the normal contact model
        void calcNormalModels() const;

        //- Create the friction contact model
        void calcFrictionModels() const;

        //- Calculate zone index
        void calcZoneIndex() const;

        //- Calculate zone
        void calcZone() const;

        //- Calculate zone old points
        void calcOldZonePoints() const;

        //- Calculate interpolation
        void calcZoneToZones() const;

        //- Make the current patch traction list
        void makeCurPatchTractions() const;

        //- Calculate rate of dissipated friction energy
        void calcQc() const;

        //- Calculate rate of dissipated friction energy for each contact
        //  pair
        void calcQcs() const;

        //- Calculate the bounding box offset
        void calcBbOffset() const;

        //- Return the bounding box offset
        scalar bbOffset() const;

        // Contact dictionary
        const dictionary& dict() const
        {
            return dict_;
        };

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

        //- Interpolate zone face values to points
        template<class Type>
        tmp<Field<Type> > zoneFaceToPointInterpolate
        (
            const label zoneIndex,
            const Field<Type>& zoneField,
            const label shadowI
        ) const;

        //- Interpolate zone point values to faces
        template<class Type>
        tmp<Field<Type> > zonePointToFaceInterpolate
        (
            const label zoneIndex,
            const Field<Type>& zonePointField,
            const label shadowI
        ) const;

public:

    //- Runtime type information
    TypeName("solidGeneralContact");


    // Constructors

        //- Construct from patch and internal field
        solidGeneralContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidGeneralContactFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidGeneralContactFvPatchVectorField
        // onto a new patch
        solidGeneralContactFvPatchVectorField
        (
            const solidGeneralContactFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidGeneralContactFvPatchVectorField
        (
            const solidGeneralContactFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<vector> > clone() const
        {
            return tmp<fvPatchField<vector> >
            (
                new solidGeneralContactFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidGeneralContactFvPatchVectorField
        (
            const solidGeneralContactFvPatchVectorField&,
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
                new solidGeneralContactFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~solidGeneralContactFvPatchVectorField();


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
            label zoneIndex() const;

            //- Return the list of shadow zone names
            const List<word>& shadowZoneNames() const;

            //- Return the list of shadow zone indices
            const List<label>& shadowZoneIndices() const;

            //- Return reference to the normal contact model
            normalContactModel& normalModel(const label shadowI);

            //- Return const reference to the normal contact model
            const normalContactModel& normalModel(const label shadowI) const;

            //- Return reference to the friction contact model
            frictionContactModel& frictionModel(const label shadowI);

            //- Return const reference to the friction contact model
            const frictionContactModel& frictionModel
            (
                const label shadowI
            ) const;

            //- Return const reference to deformed master zone
            const standAlonePatch& zone() const;

            //- Return reference to deformed master zone
            standAlonePatch& zone();

            //- Return const reference to deformed master zone
            const standAlonePatch& shadowZone(const label shadowI) const;

            //- Return reference to deformed master zone
            standAlonePatch& shadowZone(const label shadowI);

            //- Return const reference to patch-to-patch interpolation
            const newGgiStandAlonePatchInterpolation& zoneToZone
            (
                const label shadowI
            ) const;

            //- Return reference to patch-to-patch interpolation
            newGgiStandAlonePatchInterpolation& zoneToZone
            (
                const label shadowI
            );

            //- Given the patch index of a shadow, return the position in the
            // shadow list of the current patch
            label findShadowID(const label patchID) const;

            //- Return the current patch traction for shadowI
            const vectorField& curPatchTractions(const label shadowI) const;

            //- Return the current patch traction for shadowI
            vectorField& curPatchTractions(const label shadowI);

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


        // Member functions


            //- Dissipated energy due to friction
            virtual const scalarField& Qc() const;

            //- Dissipated energy due to friction from contact with shadowI
            virtual const scalarField& Qc(const label shadowI) const;

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "solidGeneralContactFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
