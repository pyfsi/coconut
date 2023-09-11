/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::timeVaryingMappedCoupledVelocityFvPatchField

Description
    This boundary conditions interpolates the values from a set of supplied
    points in space and time.

    By default the data files should be provide in
    constant/boundaryData/\<patch name\>/ directory:
      - points             : pointField of locations
      - \<time\>/\<field\> : field of values at time \<time\>

    Alternatively the names and locations of the points and field files may be
    specified explicitly via the optional dictionary entries:
      - dataDir \<optional top-level directory of the points and field data>;
      - points \<optional path including name of points file relative to
                 dataDir\>;
      - sample \<optional name of the sub-directory in the time directories
                containing the fields\>;
    This is particularly useful when mapping data from another case for which
    the \c sample \c functionObject is used to obtain the patch field data for
    mapping.

    The default mode of operation (mapMethod planarInterpolation) is to project
    the points onto a plane (constructed from the first three points) and
    construct a 2D triangulation and finds for the face centres the triangle it
    is in and the weights to the 3 vertices.

    The optional mapMethod nearest will avoid all projection and triangulation
    and just use the value at the nearest vertex.

    Values are interpolated linearly between times.

Usage
    \table
        Property     | Description                | Required    | Default value
        setAverage   | Switch to activate setting of average value | no | false
        perturb      | Perturb points for regular geometries | no | 1e-5
        fieldTable   | Alternative field name to sample | no| this field name
        mapMethod    | Type of mapping            | no | planarInterpolation
        offset       | Offset to mapped values    | no | Zero
        dataDir      | Top-level directory of the points and field data \\
                         | no | constant/boundaryData/\<patch name\>
        points       | Path including name of points file relative to dataDir \\
                         | no | points
        sample       | Name of the sub-directory in the time directories \\
                         containing the fields | no | ""
    \endtable

    \verbatim
    <patch name>
    {
        type            timeVaryingMappedCoupledVelocity;
    }
    \endverbatim

    \verbatim
    <patch name>
    {
        type            timeVaryingMappedCoupledVelocity;
        dataDir         "../<source case name>/postProcessing/sample";
        points          "0/<sample name>/faceCentres";
        sample          <sample name>;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField
    Foam::coconut_timeVaryingMappedFvPatchField
    Foam::Function1s

SourceFiles
    timeVaryingMappedCoupledVelocityFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMappedCoupledVelocityFvPatchField_H
#define timeVaryingMappedCoupledVelocityFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "coconut_timeVaryingMappedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class timeVaryingMappedCoupledVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class timeVaryingMappedCoupledVelocityFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        coconut_timeVaryingMappedFvPatchField<Type> fieldMapper_;


public:

    //- Runtime type information
    TypeName("timeVaryingMappedCoupledVelocity");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingMappedCoupledVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingMappedCoupledVelocityFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given timeVaryingMappedCoupledVelocityFvPatchField
        //  onto a new patch
        timeVaryingMappedCoupledVelocityFvPatchField
        (
            const timeVaryingMappedCoupledVelocityFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        timeVaryingMappedCoupledVelocityFvPatchField
        (
            const timeVaryingMappedCoupledVelocityFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new timeVaryingMappedCoupledVelocityFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        timeVaryingMappedCoupledVelocityFvPatchField
        (
            const timeVaryingMappedCoupledVelocityFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new timeVaryingMappedCoupledVelocityFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchField<Type>&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "timeVaryingMappedCoupledVelocityFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
