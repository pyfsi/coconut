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
    processPass

Description
    Class holding data and functions for a process pass

SourceFiles
    processPass.C
    processPassTemplates.C

Author
    Philip Cardiff, UCD.
    Peter De Jaeger, Bekaert.

\*---------------------------------------------------------------------------*/

#ifndef processPass_H
#define processPass_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "tmp.H"
#include "autoPtr.H"
#include "passItem.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "wirePassItem.H"
#include "IOobjectList.H"
#include "fvMesh.H"
#include "dataHandler.H"
#include "boundaryConditionsContainer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Class forward declaration
    class processLine;

/*---------------------------------------------------------------------------* \
                         Class processPass Declaration
\*---------------------------------------------------------------------------*/

class processPass
{
    // Private data

        //- Reference to the parent class (process line)
        const processLine* parentPtr_;

        // Name of the process pass
        const word name_;

        // Reference to the pass dictionary
        const dictionary dict_;

        // controlDict to allow creation of the time object
        dictionary controlDict_;

        // Time for the pass case
        autoPtr<Time> runTimePtr_;

        // Mesh for the pass case
        mutable autoPtr<fvMesh> meshPtr_;

        //- List of pass items i.e. tools and wire
        mutable PtrList<passItem> passItems_;

        //- Pass number
        const int passNr_;

        //- Take reference to the process line data container
        mutable dataHandler* dataContainerPtr_;

    // Private Member Functions

        //- Read the mesh from the case and set the mesh pointer
        void readMesh() const;

        //- Return reference to passItems
        PtrList<passItem>& passItems();

        //- Return const reference to passItems
        const PtrList<passItem>& passItems() const;

        //- Make the passItems
        void makePassItems() const;

        //- Merge pass item meshes
        void mergePassItemMeshes();

        //- Merge pass item material dictionaries
        void mergePassItemMaterialProperties();

        //- Merge pass item boundary conditions
        void mergePassItemBoundaryConditions();

        //- Check if the pass items have fields; if found, then merge them with
        //  the base mesh fields (or create them if not found in the base mesh)
        void mergePassItemFields(const label passItemIndex);

        //- Create the material indicator field
        void writeMaterialIndicatorField();

        //- Create the zones for the topological mesh modifiers
        void createMeshModifierZones();

        //- Write createPatchDict in case of an axi-symmetric simulation
        void writeCreateWedgePatchDict();

        //- Add functions objects to the controlDict
        void addFunctionsObjectsToControlDict();

        //- Write gnuplot script
        void writeGnuplotScript
        (
            const word patchName,
            const word dataType
        );

        //- Read pass custom settings and write them to the pass files
        void overwritePassDefaults();

        //- Generate PDF graphs using gnuplot
        void generateGnuplotPDFs();

        //- Generate images using ParaView
        void generateParaViewImages();

        //- Merge passItem fields with the base mesh fields
        template<class Type>
        void MergeVolFields
        (
            const IOobjectList& baseObjects,
            const IOobjectList& passItemObjects,
            const fvMesh& baseMesh,
            const fvMesh& passItemMesh
        );

        //- For when the rollSetup is specified, create the rollSetup meshDict
        //  in the system folder of the pass
        void generateRollSetupMeshDict(const word& rollSetup);

        //- Return the passItem dict
        const dictionary& passItemsDict() const
        {
            return dict().subDict("passItems");
        }

        //- Check if a roller is valid for a given rollSetup and symmetry
        Switch rollValid
        (
            const word& rollSetup,
            const word& processLineSymm,
            const word& passItemRollType
        ) const;

        //- Disallow copy construct
        processPass(const processPass&);

        //- Disallow default bitwise assignment
        void operator=(const processPass&);


public:

    //- Runtime type information
    TypeName("processPass");


    // Constructors

        //- Construct from dictionary
        processPass
        (
            const processLine& parent,
            const word& name,
            const dictionary& dict,
            const int& passNr,
            dataHandler* dataContainerPtr
        );


    // Destructor

        virtual ~processPass();


    // Member Functions

        //- Const reference to the parent class (process pass)
        const processLine& parent() const;

        //- Return const reference to pass number
        const int& passNr() const
        {
            return passNr_;
        }

        //- Return reference to the dataContainer
        dataHandler& dataContainer() const
        {
            return *dataContainerPtr_;
        }

        //- Return const reference to name of the pass
        const word& name() const
        {
            return name_;
        }

        //- Return const reference to process pass dictionary
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Return reference to runtime
        Time& runTime();

        //- Return const reference to runtime
        const Time& runTime() const;

        //- Return reference to mesh
        fvMesh& mesh();

        //- Return const reference to mesh
        const fvMesh& mesh() const;

        //- Return const reference to the pass item with the given name
        const passItem& passItemObject(const word& passItemName) const;

        //- Return reference to the wire pass item in the list of items
        wirePassItem& wire();

        //- Setup the pass
        void setup();

        //- Run the pass
        void run();

        //- Perform post processing
        void postProcess();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "processPassTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
