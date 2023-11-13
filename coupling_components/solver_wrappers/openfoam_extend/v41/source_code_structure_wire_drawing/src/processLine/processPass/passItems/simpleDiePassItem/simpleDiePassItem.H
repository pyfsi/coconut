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
    simpleDiePassItem

Description
    Drawing pass.

SourceFiles
    simpleDiePassItem.C

Author
    Philip Cardiff, UCD.
    Peter De Jaeger, Bekaert.

\*---------------------------------------------------------------------------*/

#ifndef simpleDiePassItem_H
#define simpleDiePassItem_H

#include "passItem.H"
#include "boundaryConditionsContainer.H"
#include "dataHandler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class simpleDiePassItem Declaration
\*---------------------------------------------------------------------------*/

class simpleDiePassItem
:
    public passItem
{

    // Private Data

        //- Switch if DU boundary conditions dictonary is generated
        Switch boundaryConditionsDUGenerated_;

        //- Switch if T boundary conditions dictonary is generated
        Switch boundaryConditionsTGenerated_;

        //- die downstream patch name
        word dieDownstreamPatchName_;

        //- die upstream patch name
        word dieUpstreamPatchName_;

        //- die contact patch name
        word dieContactPatchName_;

        //- die shadow patch name
        word dieContactShadowPatchName_;

        //- die dieToHousing patch name
        word dieToHousingPatchName_;

        //- die XZ symmetry plane patch name
        word dieSymmetryYPatchName_;

        //- die XY symmetry plane patch name
        word dieSymmetryZPatchName_;

        //- die front patch name
        word dieFrontPatchName_;

        //- die back patch name
        word dieBackPatchName_;

        //- die entrance cone patch name
        word dieEntranceConePatchName_;

        //- die exit patch name
        word dieExitConePatchName_;

        //- Friction coefficient
        scalar frictionCoefficientAxial_;

    // Private Member Functions

        //- Write mesh dict
        void writeMeshDict(Time& runTime);

        //- Disallow default bitwise copy construct
        simpleDiePassItem(const simpleDiePassItem&);

        //- Disallow default bitwise assignment
        void operator=(const simpleDiePassItem&);

public:

    //- Runtime type information
    TypeName("simpleDie");

    // Static data members


    // Constructors

        //- Construct from dictionary
        simpleDiePassItem
        (
            const processPass& parent,
            const word& name,
            const dictionary& dict,
            const int& passNr,
            dataHandler* dataContainerPtr
        );


    // Destructor

        virtual ~simpleDiePassItem();


    // Member Functions

        //- Construct the mesh
        virtual void setupMesh(Time& runTime);

        //- Return const reference to the DU boundary conditions subDict
        virtual dictionary& boundaryConditionsDU();

        //- Return const reference to the T boundary conditions subDict
        virtual dictionary& boundaryConditionsT();

        //- Return copy of the contact patch
        virtual word contactPatch();

        //- Return copy of the friction coefficient
        virtual scalar frictionCoefficientAxial();

        //- Return functionObjects for this passItem
        virtual List<dictionary> functionObjects();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
