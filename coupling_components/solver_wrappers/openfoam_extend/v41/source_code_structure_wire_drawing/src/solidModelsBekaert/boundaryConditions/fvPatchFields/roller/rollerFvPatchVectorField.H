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
    Foam::rollerFvPatchVectorField

Description
    Boundary condition for a roller,

    A rolling pass consists of six phases:
        1) displacement the roller to the required reduction;
        2) angular acceleration of the roller;
        3) rolling at a constant RPM;
        4) angular decceleration of the roller;
        5) roller is retracted from the wire;
        6) translation the roller into position for the start of the next
           pass, and angle is linearly reset to zero.

    Example of the boundary condition specification:
    @verbatim
    rollerAxis
    {
        type                   roller;
        initialRotationOrigin    ( 0.04 0.090355 0 );
        rotationAxis             ( 0 0 1 );
        displacement             ( 0 -0.000271 0 );
        dispRampTime             0.005;
        rpm                      50;
        rotationRampTime         0.005;
        rotationConstantTime     0.03;
        translation              ( 0.01 0 0 );
        translationTime          0.01;
        value                  uniform (0 0 0);
    }
    @endverbatim


SourceFiles
    rollerFvPatchVectorField.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef rollerFvPatchVectorField_H
#define rollerFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedDisplacementFvPatchVectorField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class rollerFvPatch Declaration
\*---------------------------------------------------------------------------*/

class rollerFvPatchVectorField
:
    public fixedDisplacementFvPatchVectorField
{

    // Private Data

        //- Displacement to compress wire
        // This is increased each pass
        vector displacement_;

        //- Displacement from first pass is added to the displacement each new
        // pass
        //const vector firstDisplacement_;

        //- Displacement is linearly applied over dispRampTime_ time
        scalar dispRampTime_;

        //- Revolutions-per-minute
        scalar rpm_;

        //- Current revolutions-per-minute
        scalar currentRpm_;

        //- Rotation ramp up time
        // angular velocity is linearly ramped up
        scalar rotRampTime_;

        //- Rotation axis
        const vector rotationAxis_;

        //- Initial origin of rotation
        vector initialRotationOrigin_;

        //- The current axis displacement
        //  This may be looked up by e.g. function objects
        vector currentAxisDisplacement_;

        //- Reference face centres
        mutable vectorField* origFaceCentresPtr_;

        //- Current angle - relative to initial position
        scalar curAngle_;

        //- Current time index
        label curTimeIndex_;

    // Private Member Functions

        // Make reference face centres
        void makeReferenceCentre() const;

public:

    //- Runtime type information
    TypeName("roller");


    // Constructors

        //- Construct from patch and internal field
        rollerFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rollerFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given rollerFvPatchVectorField
        //  onto a new patch
        rollerFvPatchVectorField
        (
            const rollerFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        rollerFvPatchVectorField
        (
            const rollerFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new rollerFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        rollerFvPatchVectorField
        (
            const rollerFvPatchVectorField&,
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
                new rollerFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~rollerFvPatchVectorField();


    // Member functions

        // Access functions

            //- Time at which the roller reaches max rpmp
            scalar maxRpmTime() const
            {
                return dispRampTime_ + rotRampTime_;
            }

            //- Reference to the disp ramp time
            scalar& dispRampTime()
            {
                return dispRampTime_;
            }

            //- Reference to the rot ramp time
            scalar& rotRampTime()
            {
                return rotRampTime_;
            }

            //- Reference to the rpm
            scalar& rpm()
            {
                return rpm_;
            }

            //- Const reference to the rpm
            const scalar& rpm() const
            {
                return rpm_;
            }

            //- Reference to the current rpm
            const scalar& currentRpm() const
            {
                return currentRpm_;
            }

            //- Reference to compression displacement of the roller
            vector& displacement()
            {
                return displacement_;
            }

            //- Const reference to compression displacement of the roller
            const vector& displacement() const
            {
                return displacement_;
            }

            //- Current roller angle
            scalar& curAngle()
            {
                return curAngle_;
            }

            //- Rotation axis
            const vector& rotationAxis() const
            {
                return rotationAxis_;
            }

            //- Point on the initial rotation axis
            const vector& initialRotationOrigin() const
            {
                return initialRotationOrigin_;
            }

            //- Current axis displacement
            const vector& currentAxisDisplacement() const
            {
                return currentAxisDisplacement_;
            }

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


        // Access functions

            //- Reference face centres
            const vectorField& origFaceCentres() const;

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
