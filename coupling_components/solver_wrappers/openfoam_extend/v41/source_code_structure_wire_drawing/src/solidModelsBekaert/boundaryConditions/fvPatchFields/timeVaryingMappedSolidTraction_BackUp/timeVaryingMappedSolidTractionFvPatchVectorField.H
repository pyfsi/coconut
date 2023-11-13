/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::timeVaryingMappedSolidTractionFvPatchVectorField

Description
    Foam::timeVaryingMappedSolidTractionFvPatchVectorField

SourceFiles
    timeVaryingMappedSolidTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMappedSolidTractionFvPatchVectorField_H
#define timeVaryingMappedSolidTractionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "interpolationTable.H"
#include "coordinateSystem.H"
#include "FixedList.H"
#include "instantList.H"
#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
	Class timeVaryingMappedSolidTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingMappedSolidTractionFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{
    // Private data

        //- Name of the field tables
        word pressureFieldTableName_;

        word tractionFieldTableName_;

        //Pressure        
	scalarField pressure_;

	//Traction
	vectorField traction_;

        //- Coordinate system
        autoPtr<coordinateSystem> referenceCS_;

        //- Current interpolation addressing to face centres of underlying
        //  patch
        List<FixedList<label, 3> > nearestVertex_;

        //- Current interpolation factors to face centres of underlying
        //  patch
        List<FixedList<scalar, 3> > nearestVertexWeight_;

        //- List of boundaryData time directories
        instantList sampleTimes_;

        //- Current starting index in sampleTimes
        label startSampleTime_;

        //- Interpolated pressure values from startSampleTime
        Field<scalar> startSampledPressureValues_;

	//- Interpolated traction values from startSampleTime
	Field<vector> startSampledTractionValues_;

        //- Current end index in sampleTimes
        label endSampleTime_;

        //- Interpolated pressure values from endSampleTime
        Field<scalar> endSampledPressureValues_;

	//- Interpolated traction values from endSampleTime 
	Field<vector> endSampledTractionValues_;

	bool findSampleTimeMultiple_;

    // Private Member Functions

        //- Get names of times
        static wordList timeNames(const instantList&);

        //- Find times around current time
        void findTime
        (
            const fileName& instance,
            const fileName& local,
            const scalar timeVal,
            label& lo,
            label& hi
        ) const;

        //- Read boundary points and determine interpolation weights to patch
        //  faceCentres
        void readSamplePoints();

        //- Find boundary data inbetween current time and interpolate
        void checkTable();

        //- Do actual interpolation using current weights
        template<class Type>
        tmp<Field<Type> > interpolate( const Field<Type>&) const;


public:

    //- Runtime type information
    TypeName("timeVaryingMappedSolidTraction");

    // Constructors

        //- Construct from patch and internal field
        timeVaryingMappedSolidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingMappedSolidTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given timeVaryingMappedSolidTractionFvPatchVectorField
        //  onto a new patch
        timeVaryingMappedSolidTractionFvPatchVectorField
        (
            const timeVaryingMappedSolidTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingMappedSolidTractionFvPatchVectorField
        (
            const timeVaryingMappedSolidTractionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingMappedSolidTractionFvPatchVectorField(*this)
            );
        }
        //- Construct as copy setting internal field reference
        timeVaryingMappedSolidTractionFvPatchVectorField
        (
            const timeVaryingMappedSolidTractionFvPatchVectorField&,
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
                new timeVaryingMappedSolidTractionFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the coordinateSystem
            const coordinateSystem& referenceCS() const
            {
                return referenceCS_;
            }
           
            virtual const scalarField& pressure() const
 	    {
		return pressure_;
	    }

	    virtual scalarField& pressure()
            {
		return pressure_;
	    }
	    
            virtual const vectorField& traction() const
	    {
		return traction_;
	    }

	    virtual vectorField& traction()
	    {
		return traction_;
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
                const fvPatchVectorField&,
                const labelList&
            );


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

#endif

// ************************************************************************* //
