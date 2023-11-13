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
    wirePassItem

Description
    wirePassItem class

SourceFiles
    wirePassItem.C

Author
    Philip Cardiff, UCD.
    Peter De Jaeger, Bekaert.

\*---------------------------------------------------------------------------*/

#ifndef wirePassItem_H
#define wirePassItem_H

#include "passItem.H"
#include "boundaryConditionsContainer.H"
#include "dataHandler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class wirePassItem Declaration
\*---------------------------------------------------------------------------*/

class wirePassItem
:
    public passItem
{
private:

    // Private Data

        //- Switch if mesh needs to be generated
        Switch createMesh_;

        //- Switch if DU boundary conditions dictonary is generated
        Switch boundaryConditionsDUGenerated_;

        //- Switch if T boundary conditions dictonary is generated
        Switch boundaryConditionsTGenerated_;

        //- wire downstream patch name
        word wireDownstreamPatchName_;

        //- wire upstream patch name
        word wireUpstreamPatchName_;

        //- wire contact patch name
        word wireContactPatchName_;

        //- wire shadow patch name
        word wireContactShadowPatchName_;

        //- wire shadow patch friction coefficients list
        scalarList frictionCoefficientsAxial_;

        //- wire shadow patch names list
        wordList shadowPatchNames_;

        //- wire XZ symmetry plane patch name
        word wireSymmetryYPatchName_;

        //- wire XY symmetry plane patch name
        word wireSymmetryZPatchName_;

        //- wire front patch name
        word wireFrontPatchName_;

        //- wire back patch name
        word wireBackPatchName_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        wirePassItem(const wirePassItem&);

        //- Disallow default bitwise assignment
        void operator=(const wirePassItem&);


public:

    //- Runtime type information
    TypeName("wire");


    // Constructors

        //- Construct from name of dictionary
        wirePassItem
        (
            const processPass& parent,
            const word& name,
            const dictionary& dict,
            const int& passNr,
            dataHandler* dataContainerPtr
        );


    // Destructor

        virtual ~wirePassItem();


    // Member Functions

        //- Write mesh dict
        void writeMeshDict(Time& runTime, fileName geometryFile);

        //- Return const reference to createMesh switch
        const Switch& createMesh() const
        {
            return createMesh_;
        }

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

        //- Populate shadow contact patch list
        void addShadowPatchName(word shadowPatchName);

        //- Populate friction coefficients of the shadow patches
        void addFrictionCoefficientAxial(scalar frictionCoefficientAxial);

        //- Return reference to wire shadow contact patch list
        wordList& shadowPatchNames();

        //- Translate the passItem mesh by the given vector
        virtual void positionMesh(const Time& runTime, const fvMesh& wireMesh);

        //- Return functionObjects for this passItem
        virtual List<dictionary> functionObjects();

        //- Write out the current deformed wire profile i.e. the coordinates of
        // the boundray points on the wire downstream patch
        virtual void writeCurrentProfile(const Time& runTime);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
