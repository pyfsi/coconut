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
    dugdaleModeICohesiveLaw

Description
    DugdaleFixedModeMix cohesive law
    The normal traction stress stays constant during damage at sigmaMax; the
    shear traction is set to zero i.e. mode-I only propagation.

SourceFiles
    dugdaleModeICohesiveLaw.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef dugdaleModeICohesiveLaw_H
#define dugdaleModeICohesiveLaw_H

#include "cohesiveLaw.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "PtrDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dugdaleModeICohesiveLaw Declaration
\*---------------------------------------------------------------------------*/

class dugdaleModeICohesiveLaw
:
    public cohesiveLaw
{

  private:

    // Private data

        //- Criticial fracture energy mode I
        dimensionedScalar GIc_;

        //- Maximum opening cohesive strength
        dimensionedScalar sigmaMax_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dugdaleModeICohesiveLaw(const dugdaleModeICohesiveLaw&);

        //- Disallow default bitwise assignment
        void operator=(const dugdaleModeICohesiveLaw&);

public:

    //- Runtime type information
    TypeName("dugdaleFixedModeMix");

    // Static data members


    // Constructors

        //- Construct from dictionary
        dugdaleModeICohesiveLaw
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~dugdaleModeICohesiveLaw();


    // Member Functions

        //- Return materials;
        virtual tmp<volScalarField> materials() const;

        //- Return opening cohesive strength value
        virtual scalar sigmaMaxValue() const
        {
            return sigmaMax_.value();
        };

        //- Return opening cohesive strength surface field
        virtual tmp<surfaceScalarField> sigmaMax() const;

        //- Return shearing cohesive strength value
        virtual scalar tauMaxValue() const
        {
            notImplemented
            (
                "scalar dugdaleModeICohesiveLaw::tauMaxValue()const"
            );

            return 0.0;
        }

        //- Return shearing cohesive strength surface field
        virtual tmp<surfaceScalarField> tauMax() const;

        //-  Return Mode I fracture energy value
        virtual scalar GIcValue() const
        {
            return GIc_.value();
        };

        //-  Return Mode I fracture energy surface field
        virtual tmp<surfaceScalarField> GIc() const;

        //-  Return Mode II fracture energy value
        virtual scalar GIIcValue() const
        {
            notImplemented
            (
                "scalar dugdaleModeICohesiveLaw::GIIcValue()const"
            );

            return 0.0;
        }

        //-  Return Mode II fracture energy surface field
        virtual tmp<surfaceScalarField> GIIc() const;

        //- Return tractions during damage given current deltas and Gs
        virtual void damageTractions
        (
            scalar& tN, // current normal traction
            scalar& tS, // current shear traction
            const scalar deltaN, // current deltaN
            const scalar deltaS, // current deltaS
            const scalar GI, // current dissipated GI
            const scalar GII, // current dissipated GII
            const label faceID, // needed for multi-mat
            const scalarField& globalPatchMaterials
            ) const;

        //- Return cohsive tractions during damage given current deltas and Gs
        //  for implicit damage procedure
        virtual void damageTractions
        (
            scalar& tN, // current normal traction
            scalar& tS, // current shear traction
            const scalar deltaN, // current deltaN
            const scalar deltaS, // current deltaS
            const scalar GI, // current dissipated GI
            const scalar GII, // current dissipated GII
            const label internalFaceID // needed for multi-mat
        ) const;

        //- Correct the cohesive model
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
