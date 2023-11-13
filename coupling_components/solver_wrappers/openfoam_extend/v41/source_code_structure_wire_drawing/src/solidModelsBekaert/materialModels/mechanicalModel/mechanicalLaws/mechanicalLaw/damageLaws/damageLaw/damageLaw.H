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
    damageLaw

Description
    Material mechanical for solids.

SourceFiles
    damageLaw.C
    newDamageLaw.C

\*---------------------------------------------------------------------------*/

#ifndef damageLaw_H
#define damageLaw_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class damageLaw Declaration
\*---------------------------------------------------------------------------*/

class damageLaw
{
    // Private data

        //- Name
        const word name_;

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Copy of the damage law dictionary
        const dictionary dict_;

        //- Is the damage coupled
        const Switch coupled_;

        //- Damage field
        //  undamaged = 1.0
        //  fully damaged = 0.0
        volScalarField* damagePtr_;

        //- Const reference to curMaterial index field: this is 1.0 on the cells
        //  owned by the current material
        const volScalarField& curMaterial_;


    // Private Member Functions

        //- Disallow copy construct
        damageLaw(const damageLaw&);

        //- Disallow default bitwise assignment
        void operator=(const damageLaw&);


public:

    //- Runtime type information
    TypeName("damageLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            damageLaw,
            dictionary,
            (
                const word name,
                const fvMesh& mesh,
                const dictionary& dict,
                const volScalarField& curMaterial
            ),
            (name, mesh, dict, curMaterial)
        );


    // Selectors

        //- Return a reference to the selected mechanical model
        static autoPtr<damageLaw> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const volScalarField& curMaterial
        );


    // Constructors

        //- Construct from dictionary
        damageLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const volScalarField& curMaterial
        );


    // Destructor

        virtual ~damageLaw();


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return coupled switch
        const Switch& coupled() const
        {
            return coupled_;
        }

        //- Return reference to curMaterial
        volScalarField& damage()
        {
            if (!damagePtr_)
            {
                FatalErrorIn("volScalarField& damage()")
                    << "pointer not set!" << abort(FatalError);
            }

            return *damagePtr_;
        }

        //- Return reference to curMaterial
        const volScalarField& curMaterial() const
        {
            return curMaterial_;
        }

        //- Update the damage field
        virtual void updateDamage(const volScalarField DEpsilonPEq) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
