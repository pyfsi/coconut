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
    mechanicalLaw

Description
    Material mechanical for solids.

SourceFiles
    mechanicalLaw.C
    newMechanicalLaw.C
    mechanicalLawTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef mechanicalLaw_H
#define mechanicalLaw_H

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
                         Class mechanicalLaw Declaration
\*---------------------------------------------------------------------------*/

class mechanicalLaw
{
    // Private data

        //- Name
        const word name_;

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Copy of law dictionary
        const dictionary dict_;

        // Current materials index field
        //     1 for the current materials
        //     0 for other materials
        mutable volScalarField* curMaterialPtr_;

        // Current materials index field
        //     1 for the current materials
        //     0 for other materials
        mutable surfaceScalarField* curMaterialfPtr_;

        // Current material index i.e. is the current material number 0, 1,
        // 2, etc.
        mutable label curMatIndex_;


    // Private Member Functions

        // Make current material field
        void calcCurMaterial() const;

        // Make current material field
        void calcCurMaterialf() const;

        //- Disallow copy construct
        mechanicalLaw(const mechanicalLaw&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicalLaw&);


protected:

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Const reference to law dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Set the field pointer by looking up the field from the
        //  objectRegistry; if the field is not found then create it
        //  The initial feidl value is set to defaultValue if the field if not
        //  read from disk
        template<class Type>
        void SetFieldPtr
        (
            GeometricField<Type, fvPatchField, volMesh>*& fieldPtr,
            const word& fieldName,
            const dimensioned<Type>& defaultValue,
            const word& patchTypes = "calculated"
        );

        //- Set the field pointer by looking up the field from the
        //  objectRegistry; if the field is not found then create it
        //  The initial feidl value is set to defaultValue if the field if not
        //  read from disk
        template<class Type>
        void SetFieldPtr
        (
            GeometricField<Type, fvsPatchField, surfaceMesh>*& fieldPtr,
            const word& fieldName,
            const dimensioned<Type>& defaultValue,
            const word& patchTypes = "calculated"
        );

        //- Clean up the pointer
        //  i.e. delete the object if it exists and set the pointer to NULL
        template<class Type>
        void CleanFieldPtr
        (
            GeometricField<Type, fvPatchField, volMesh>*& fieldPtr,
            const word& fieldName
        );

        //- Clean up the pointer
        //  i.e. delete the object if it exists and set the pointer to NULL
        template<class Type>
        void CleanFieldPtr
        (
            GeometricField<Type, fvsPatchField, surfaceMesh>*& fieldPtr,
            const word& fieldName
        );


public:

    //- Runtime type information
    TypeName("mechanicalLaw");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            mechanicalLaw,
            dictionary,
            (
                const word name,
                const fvMesh& mesh,
                const dictionary& dict,
                const label lawIndex
            ),
            (name, mesh, dict, lawIndex)
        );


    // Selectors

        //- Return a reference to the selected mechanical model
        static autoPtr<mechanicalLaw> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const label lawIndex
        );


    // Constructors

        //- Construct from dictionary
        mechanicalLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const label lawIndex
        );


    // Destructor

        virtual ~mechanicalLaw();


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return density
        virtual tmp<volScalarField> rho() const = 0;

        //- Return modulus of elasticity
        virtual tmp<volScalarField> E() const = 0;

        //- Return Poisson's ratio
        virtual tmp<volScalarField> nu() const = 0;

        //- Return bulk modulus
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return Rhie Chow scale factor
        virtual tmp<volScalarField> RhieChowScaleFactor() const;

        //- Update stress based on mechanical law
        virtual void correct(volSymmTensorField& sigma)
        {}

        //- Update stress based on mechanical law
        virtual void correct(volSymmTensorField& sigma, const int flag)
        {
            notImplemented
            (
                "virtual void correct"
                "(volSymmTensorField& sigma, const int flag)"
            );
        }

        //- Update stress based on mechanical law
        virtual void correct(surfaceSymmTensorField& sigma)
        {}

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress
        //  This is called at end of time-step and does not necessarilly have to
        //  do anything with the yield stress: it's name is historical
        virtual void updateYieldStress();

	virtual void resetYieldStress();//MATHIEU: Reset the yield stress in plastic deformation regime when a partitioned FSI-calculation is executed.
					//This is called at the beginning of a FSI-iteration. 

        //- Rate of dissipated energy from plasticity
        virtual tmp<volScalarField> plasticDissipationRate() const
        {
            notImplemented("tmp<volScalarField> plasticDissipationRate()");

            // Keep compiler happy
            tmp<volScalarField>* tmpPtr = NULL;
            return *tmpPtr;
        }

        //- Return the desired new time-step
        virtual scalar newDeltaT() const;

        // Return reference to current material field
        const volScalarField& curMaterial() const;

        // Return reference to current material field
        volScalarField& curMaterial();

        // Return reference to current material field
        const surfaceScalarField& curMaterialf() const;

        // Return reference to current material field
        surfaceScalarField& curMaterialf();

        // Return reference to current material index
        label& curMaterialIndex() const
        {
            return curMatIndex_;
        }

        //- Set current material index; this is only needed for multiMaterial
        //  analyses
        virtual void setMaterialIndex(label curMaterialIndex)
        {
            curMatIndex_ = curMaterialIndex;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "mechanicalLawTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
