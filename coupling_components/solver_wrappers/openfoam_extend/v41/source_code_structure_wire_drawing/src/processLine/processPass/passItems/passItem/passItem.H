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
    passItem

Description
    Class holding data and functions of the deformation tools

SourceFiles
    passItem.C
    newPassItem.C

Author
    Philip Cardiff, UCD.
    Peter De Jaeger, Bekaert.

\*---------------------------------------------------------------------------*/

#ifndef passItem_H
#define passItem_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "dataHandler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Class forward declaration
    class processPass;
    class wirePassItem;

/*---------------------------------------------------------------------------* \
                         Class passItem Declaration
\*---------------------------------------------------------------------------*/

class passItem
{
    // Private data

        //- Reference to the parent class (process pass)
        const processPass* parentPtr_;

        //- Name of the tool
        const word name_;

        //- Reference to the pass dictionary
        dictionary dict_;

        //- Reference to the pass number
        const int& passNr_;

        // Mesh for the region
        mutable autoPtr<fvMesh> meshPtr_;

        //- wire length
        scalar wireLength_;

        //- wire axial resolution
        scalar wireAxialResolution_;

        //- Take reference to the process line data container
        dataHandler* dataContainerPtr_;

    // Private Member Functions

        //- Read the mesh from the region and set the mesh pointer
        void readMesh(const Time& runTime) const;

        //- Disallow copy construct
        passItem(const passItem&);

        //- Disallow default bitwise assignment
        void operator=(const passItem&);


public:

    //- Runtime type information
    TypeName("passItem");

        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            passItem,
            dictionary,
            (
                const processPass& parent,
                const word& name,
                const dictionary& dict,
                const int& passNr,
                dataHandler* dataContainerPtr
            ),
            (parent, name, dict, passNr, dataContainerPtr)
        );


    // Selectors

        //- Return a reference to the selected deformation tool
        static autoPtr<passItem> New
        (
            const processPass& parent,
            const word& name,
            const dictionary& dict,
            const int& passNr,
            dataHandler* dataContainerPtr
        );


    // Constructors

        //- Construct from dictionary
        passItem
        (
            const processPass& parent,
            const word& name,
            const dictionary& dict,
            const int& passNr,
            dataHandler* dataContainerPtr
        );


    // Destructor

        virtual ~passItem()
        {}


    // Member Functions

        //- Const reference to the parent class (process pass)
        const processPass& parent() const;

        //- Return reference to the dataContainer
        dataHandler& dataContainer()
        {
            return *passItem::dataContainerPtr_;
        }

        //- Return const reference to the name of the tool
        const word& name() const
        {
            return name_;
        }

        //- Return const reference to passItem dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Return reference to passItem dictionary
        dictionary& dict()
        {
            return dict_;
        }

        //- Return reference to mesh
        fvMesh& mesh(const Time& runTime);

        //- Return const reference to mesh
        const fvMesh& mesh(const Time& runTime) const;

        //- Return reference to mesh input dict
        dictionary& meshInputDict()
        {
            return dict().subDict("mesh");
        }

        //- Return const reference to the phyiscal properties
        const dictionary& physcialProperties() const
        {
            return dict_.subDict("physicalProperties");
        }

        //- Return const reference to the mechanical law subDict
        const dictionary& mechanicalLawDict() const
        {
            return physcialProperties().subDict("mechanical");
        }

        //- Return const reference to the thermal law subDict
        const dictionary& thermalLawDict() const
        {
            return physcialProperties().subDict("thermal");
        }

        //- Return const reference to the boundary conditions subDict
        dictionary& boundaryConditions()
        {
            return dict_.subDict("boundaryConditions");
        }

        //- Return const reference to the DU boundary conditions subDict
        virtual dictionary& boundaryConditionsDU()
        {
            return boundaryConditions().subDict("DU");
        }

        //- Return const reference to the T boundary conditions subDict
        virtual dictionary& boundaryConditionsT()
        {
            return boundaryConditions().subDict("T");
        }

        //- Return const reference to the contact patch
        virtual word contactPatch() = 0;

        //- Return const reference to the contact patch
        virtual scalar frictionCoefficientAxial() = 0;

        //- Translate the passItem mesh by the given vector
        virtual void positionMesh(const Time& runTime, const fvMesh& wireMesh);

        //- Translate the passItem mesh by the given vector
        void positionMesh(const Time& runTime, const vector& translation);

        //- Construct the mesh
        virtual void setupMesh(Time& runTime) = 0;

        //- Setup files in the pass for this tool
        virtual void setup(Time& runTime);

        //- Return functionObjects for this passItem
        virtual List<dictionary> functionObjects()
        {
            // By default we will return no function objects; this can be
            // overwritten by the specific pass items
            return List<dictionary>(0);
        };

        //- set the wire length
        void setWireLength(scalar wireLength)
        {
            wireLength_ = wireLength;
        }

        //- get the wire length
        scalar wireLength()
        {
            return wireLength_;
        }

        //- set the wire axial resolution
        void setWireAxialResolution(scalar wireAxialResolution)
        {
            wireAxialResolution_ = wireAxialResolution;
        }

        //- get the wire axial resolution
        scalar wireAxialResolution()
        {
            return wireAxialResolution_;
        }

        //- get pass nr
        const int& passNr()
        {
            return passNr_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
