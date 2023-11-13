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
    solidCohesiveFvPatchVectorField

Description
    Cohesive law fv patch field for arbitrary crack procedure,
    where the mode mixity is allowed to vary depending on how
    the crack opens.

    This procedure allows sigmaMax, tauMax, GIc, and GIIc
    to be independently specified.

    The cohesive law surface field is looked up from mechanical.

    The dissipated fracture energy is constantly monitored, by
    integrating in time.

Author
    Philip Cardiff UCD
    based on original Tukovic cohesive law boundary condition.

SourceFiles
    solidCohesiveFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidCohesiveFvPatchVectorField_H
#define solidCohesiveFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"
#include "cohesiveFvPatch.H"
#include "Switch.H"
#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class solidCohesiveFvPatch Declaration
\*---------------------------------------------------------------------------*/

class solidCohesiveFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{

    // Private Data

        //- Name of the displacement field
        const word fieldName_;

        //- Cohesive law - now part of mechanicalModel
        //cohesiveZoneModel* cohesiveZoneModelPtr_;

        //- Under-relaxation factor
        scalar relaxationFactor_;

        //- Current traction
        vectorField traction_;

        //- This field indicates separate cracks
        scalarField identity_;

        //- Crack pressure field which can be used to present hydraulic pressure
        scalarField pressure_;

        //- Current deltas: displacement difference between a crakc face and its
        //  shadow face
        vectorField delta_;

        //- Minimal unloading relative separation distance
        //scalar minUnloadingSeparationDistance_;

        //- If yes, contact procedure will be used to calculate
        //  cohesive traction for the negative normal separation distance
        Switch contact_;

        //- if a face has cracked
        //Field<bool> cracked_;
        scalarField cracked_;

        //- Current normal traction
        scalarField curTractionN_;

        //- Old normal traction
        scalarField oldTractionN_;

        //- Current shear traction
        scalarField curTractionS_;

        //- Old shear traction
        scalarField oldTractionS_;

        //- Current normal delta
        scalarField deltaN_;

        //- Old normal delta
        scalarField oldDeltaN_;

        //- Current shear delta
        scalarField deltaS_;

        //- Old shear delta
        scalarField oldDeltaS_;

        //- Unloading effective delta: this is the maximum delta effective that
        //  has been reached by the face, so we can monitor if a face starts to
        //  unload
        scalarField unloadingDeltaEff_;

        //- Current energy dissipated in mode-I
        scalarField currentGI_;

        //- Old energy dissipated in mode-I
        scalarField oldGI_;

        //- Current energy dissipated in mode-II
        scalarField currentGII_;

        //- Old energy dissipated in mode-II
        scalarField oldGII_;

        //- Current time index
        label curTimeIndex_;

        //- Penalty factor and returning penalty factor
        //  for penalty contact
        scalar* penaltyFactorPtr_;

        //- Penalty factor scale: this is a user option that scales the penalty
        //  factor
        scalar penaltyScale_;

        //- Coulomb friction coefficient for contact
        scalar frictionCoeff_;

        //- explicitSeparationDistance is a switch to use old separation
        //  distances, in an attempt to provide better convergence
        Switch explicitSeparationDistance_;

        //- curDeltaN/curDeltaS is always latest delta
        //  whereas deltaN/deltaS is only latest when explicitDist is false
        scalarField curDeltaN_;
        scalarField curDeltaS_;

        //- Update global patch material field needed for multimaterial cohesive
        //  laws
        bool updateGlobalPatchMaterials_;

        //- Global patch material field needed for multimaterial cohesive laws
        scalarField globalPatchMaterials_;

        //- Disable directAddresing in the case of mesh topological changes
        label nOldCells_;

        //- Offset with the contact law: defaults to 0
        scalar offset_;

        //- Histeretic offset with the contact law: defaults to 0: may not be
        //  needed?
        scalar histereticOffset_;


     // Private member functions

        //- penalty factor
        void calcPenaltyFactor();
        virtual scalar penaltyFactor()
        {
            if (!penaltyFactorPtr_)
            {
                calcPenaltyFactor();
            }

            return *penaltyFactorPtr_;
        }


public:

    //- Runtime type information
    TypeName("solidCohesive");


    // Constructors

        //- Construct from patch and internal field
        solidCohesiveFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidCohesiveFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidCohesiveFvPatchVectorField onto a new patch
        solidCohesiveFvPatchVectorField
        (
            const solidCohesiveFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidCohesiveFvPatchVectorField
        (
            const solidCohesiveFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new solidCohesiveFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        solidCohesiveFvPatchVectorField
        (
            const solidCohesiveFvPatchVectorField&,
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
                new solidCohesiveFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return cohesive law
            //const cohesiveZoneModel& law() const;

            //- Return current traction
            // const vectorField& traction() const
            vectorField& traction()
            {
                return traction_;
            }

            // Eric added:
            //- Return pressure etc.
            scalarField& identity()
            {
                return identity_;
            }

            scalarField& pressure()
            {
                return pressure_;
            }

            const scalarField& pressure() const
            {
                return pressure_;
            }

            // Return current deltas
            vectorField& delta()
            {
                return delta_;
            }

            //- crack and damage field for visualisation
            //  cracked faces set to 2, damaged but not cracked set to 1
            tmp<scalarField> crackingAndDamage() const;

            //- GI field for visualisation
            tmp<scalarField> GI() const;

            //- GII field for visualisation
            tmp<scalarField> GII() const;

            //- Print number of damaged and cracked faces
            bool cracking();

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchVectorField&) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
