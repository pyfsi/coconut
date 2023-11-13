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
    levanovFrictionModel

Description
    Levanov type friction model, where shear stress is a function of slip,
    pressure and the material yield stress.

SourceFiles
    levanovFrictionModel.C

Author
    Philip Cardiff, UCD. All rights reserved.
    Tian Tang, Bekaert. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef levanovFrictionModel_H
#define levanovFrictionModel_H

#include "frictionContactModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class levanovFrictionModel Declaration
\*---------------------------------------------------------------------------*/

class levanovFrictionModel
:
    public frictionContactModel
{
    // Private data

         //- The contact model dictionary
        dictionary frictionContactModelDict_;

         //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Slave traction
        vectorField slaveTraction_;

        //- Previous slave traction used for under-relaxation
        vectorField prevSlaveTraction_;

        //- Slip
        vectorField slip_;

        //- Under-relaxation factor
        scalar relaxFac_;

        //- Surface constant taking into account the state of the roughness and
        //  the adhesion properties of the contacting surfaces. Have to be
        //  calibrated via measurements. k0 value has to be set between 0 and 1.
        const scalar k0_;

        //-Iinitial flow strength of the steel
        const scalar kf0_;

        //- Model constants with default values
        const scalar C_;
        const scalar a_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        //levanovFrictionModel(const levanovFrictionModel&);

        //- Disallow default bitwise assignment
        void operator=(const levanovFrictionModel&);

public:

    //- Runtime type information
    TypeName("levanov");


    // Constructors

        //- Construct from dictionary
        levanovFrictionModel
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID
        );

        //- Construct as copy
        levanovFrictionModel(const levanovFrictionModel&);

        //- Construct and return a clone
        virtual autoPtr<frictionContactModel> clone() const
        {
            return
                autoPtr<frictionContactModel>
                (
                    new levanovFrictionModel(*this)
                );
        }


    //- Destructor
    virtual ~levanovFrictionModel()
    {}


    // Member Functions

        //- Correct the friction model
        virtual void correct
        (
            const vectorField& slavePressure,
            const vectorField& slaveFaceNormals,
            const scalarField& slavePatchAreaInContact,
            const vectorField& slaveDU,
            const vectorField& masterDUInterpToSlave
        );

        //- Return slave friction traction
        virtual const vectorField& slaveTraction() const
        {
            return slaveTraction_;
        }

        //- Return slip field
        virtual const vectorField& slip() const
        {
            // Note: slip is updated inside the correct function
            return slip_;
        }

        //- Perform any reauired field mapping e.g. after a topo change
        virtual void autoMap(const fvPatchFieldMapper& m);

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
