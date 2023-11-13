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
    Foam::thermalSymmetryFvPatchScalarField

Description
    Symmetry plane with non-orthogonal correction for scalar field.

SourceFiles
    thermalSymmetryFvPatchScalarField.C

Author
    Philip Cardiff UCD based on symmetryDisplacement

\*---------------------------------------------------------------------------*/

#ifndef thermalSymmetryFvPatchScalarField_H
#define thermalSymmetryFvPatchScalarField_H

#include "symmetryFvPatchFields.H"
#include "symmetryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class thermalSymmetryFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class thermalSymmetryFvPatchScalarField
:
    public symmetryFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("thermalSymmetry");


    // Constructors

        //- Construct from patch and internal field
        thermalSymmetryFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermalSymmetryFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given thermalSymmetryFvPatchScalarField
        //  onto a new patch
        thermalSymmetryFvPatchScalarField
        (
            const thermalSymmetryFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        thermalSymmetryFvPatchScalarField
        (
            const thermalSymmetryFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new thermalSymmetryFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        thermalSymmetryFvPatchScalarField
        (
            const thermalSymmetryFvPatchScalarField&,
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
                new thermalSymmetryFvPatchScalarField(*this, iF)
            );
        }

        tmp<Field<scalar> > snGrad() const;

        void evaluate(const Pstream::commsTypes);

        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
