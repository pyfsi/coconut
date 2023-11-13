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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    hydrostaticPressureFvPatchVectorField

Description
    A pressure boundary condition for solids, where the user inputs the fluid
    height and the boundary condition calculates the local pressure based on the
    coordinates i.e. p = rho*g*h.

    Optionally, the user can specify the reference height (defaults to zero),
    the fluid density (defaults to 1000) and the height direction (defaults to
    vector(0, 1, 0).

Author
   Philip Cardiff, UCD. All rights reserved.

SourceFiles
    hydrostaticPressureFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef hydrostaticPressureFvPatchVectorField_H
#define hydrostaticPressureFvPatchVectorField_H

#include "solidTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class hydrostaticPressureFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class hydrostaticPressureFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{

    // Private Data

        //- Fluid density
        const scalar rho_;

        //- Height direction
        const vector heightDirection_;

        //- Reference position (i.e. zero height is at this position)
        const scalar referencePosition_;


public:

    //- Runtime type information
    TypeName("hydrostaticPressure");


    // Constructors

        //- Construct from patch and internal field
        hydrostaticPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        hydrostaticPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  hydrostaticPressureFvPatchVectorField onto a new patch
        hydrostaticPressureFvPatchVectorField
        (
            const hydrostaticPressureFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        hydrostaticPressureFvPatchVectorField
        (
            const hydrostaticPressureFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new hydrostaticPressureFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        hydrostaticPressureFvPatchVectorField
        (
            const hydrostaticPressureFvPatchVectorField&,
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
                new hydrostaticPressureFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
