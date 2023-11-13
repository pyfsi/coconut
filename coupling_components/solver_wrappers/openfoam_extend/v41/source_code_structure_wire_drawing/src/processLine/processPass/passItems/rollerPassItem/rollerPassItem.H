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
    rollerPassItem

Description
    Drawing pass.

SourceFiles
    rollerPassItem.C

Author
    Philip Cardiff, UCD.
    Peter De Jaeger, Bekaert.

\*---------------------------------------------------------------------------*/

#ifndef rollerPassItem_H
#define rollerPassItem_H

#include "passItem.H"
#include "boundaryConditionsContainer.H"
#include "dataHandler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rollerPassItem Declaration
\*---------------------------------------------------------------------------*/

class rollerPassItem
:
    public passItem
{

    // Private Data

        //- Type of roller, i.e., topRoll, bottomRoll, leftRoll or rightRoll
        word rollerType_;

        //- Switch if DU boundary conditions dictonary is generated
        Switch boundaryConditionsDUGenerated_;

        //- Switch if T boundary conditions dictonary is generated
        Switch boundaryConditionsTGenerated_;

        //- roller axis patch name extension
        word rollerAxisPatchName_;

        //- roller contact patch name extension
        word rollerContactPatchName_;

        //- roller shadow patch name extension
        word rollerContactShadowPatchName_;

        //- roller back patch name extension
        word rollerBackPatchName_;

        //- roller front patch name extension
        word rollerFrontPatchName_;

        //- roller XZ symmetry plane patch name
        word rollerSymmetryYPatchName_;

        //- roller XY symmetry plane patch name
        word rollerSymmetryZPatchName_;

        //- roller-to-roller front patch name
        word rollerToRollerFrontPatchName_;

        //- roller-to-roller back patch name
        word rollerToRollerBackPatchName_;

        //- roller-to-air back patch name
        word rollerToAirPatchName_;

        //- Friction coefficient
        scalar frictionCoefficientAxial_;

        //- Roller position after aligning with wire
        vector rollerPosition_;

        //- Motion enabler: either roller (drien) or freeRoller (non-driven)
        word motionEnabler_;

    // Private Member Functions

        //- Write mesh dict
        void writeMeshDict(Time& runTime);

        //- Disallow default bitwise copy construct
        rollerPassItem(const rollerPassItem&);

        //- Disallow default bitwise assignment
        void operator=(const rollerPassItem&);

public:

    //- Runtime type information
    TypeName("roller");

    // Static data members


    // Constructors

        //- Construct from dictionary
        rollerPassItem
        (
            const processPass& parent,
            const word& name,
            const dictionary& dict,
            const int& passNr,
            dataHandler* dataContainerPtr
        );


    // Destructor

        virtual ~rollerPassItem();


    // Member Functions

        //- Translate the passItem mesh by the given vector
        virtual void positionMesh(const Time& runTime, const fvMesh& wireMesh);

        //- Construct the mesh
        virtual void setupMesh(Time& runTime);

        //- Return const reference to the DU boundary conditions subDict
        virtual dictionary& boundaryConditionsDU();

        //- Return const reference to the T boundary conditions subDict
        virtual dictionary& boundaryConditionsT();

        //- Add roller-to-roller boundary conditions
        void rollerToRollerBoundaryConditions();

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
