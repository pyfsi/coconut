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
    boundaryConditionsContainer

Description
    Process line class in charge of passes.

SourceFiles
    boundaryConditionsContainer.C

Author
    Philip Cardiff, UCD.
    Peter De Jaeger, Bekaert.

\*---------------------------------------------------------------------------*/

#ifndef boundaryConditionsContainer_H
#define boundaryConditionsContainer_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "scalar.H"
#include "vector.H"
#include "scalarList.H"
#include "dataHandler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class boundaryConditionsContainer Declaration
\*---------------------------------------------------------------------------*/

class boundaryConditionsContainer
{
    // Private data

        //- Take reference to the passItem boundary conditions input dictionary
        dictionary& dict_;

        //- Take reference to the dataContainer
        dataHandler* dataContainerPtr_;

private:

    // Private Member Functions

        //- Print default boundary condition setup info
        void printDefaultBCInfo(word patchName, word BCtype);

        //- Print copy boundary condition setup info
        void printCopyBCInfo(word patchName);

        //- Check if input boundary dictionary will be copied
        bool copyDict(word patchName);

        //- Disallow default bitwise copy construct
        boundaryConditionsContainer(const boundaryConditionsContainer&);

        //- Disallow default bitwise assignment
        void operator=(const boundaryConditionsContainer&);


public:

    //- Runtime type information
    TypeName("boundaryConditionsContainer");


    // Constructors

        //- Construct from dictionary
        boundaryConditionsContainer
        (
            dictionary& dict,
            dataHandler* dataContainerPtr
        );


    // Destructor

        virtual ~boundaryConditionsContainer();


    // Member Functions

        //- Return reference to the dataContainer
        dataHandler& dataContainer()
        {
            return *dataContainerPtr_;
        }

        // Momentum boundary conditions
        //- Add solidSymmetry boundary condition
        void solidSymmetry(word patchName);

        //- Add fixedTangentialNormalPressure boundary condition
        //  Optional sub-directory for timeVsPressure files
        void fixedTangentialNormalPressure
        (
            const word& patchName,
            const word& subDirectory = ""
        );

        //- Add fixedDisplacementZeroShear boundary condition
        //  Optional sub-directory for timeVsDisp files
        void fixedDisplacementZeroShear
        (
            const word& patchName,
            const word& subDirectory = ""
        );

        //- Add solidContactSlave boundary condition
        void solidContactSlave
        (
            word patchName, word shadowPatchName
        );

        //- Add solidContactMaster boundary condition
        void solidContactMaster
        (
            const word& patchName,
            const word& shadowPatchName,
            const scalar frictionCoefficientAxial,
            const scalar targetAveragePenetration,
            const scalar minPenaltyScale,
            const scalar maxPenaltyScale
        );

        //- Add solidContactMaster boundary condition for 1-many contact
        void solidContactMaster
        (
            word patchName,
            wordList shadowPatchNames,
//            word frictionLaw,
            scalarList frictionCoefficients
        );

        //- Add solidRigidContact boundary condition for 1-many contact
        void solidRigidContact
        (
            word patchName,
            dictionary triSurfacesDict,
//            word frictionLaw,
            scalarList frictionCoefficients
        );

        //- Add solidTraction boundary condition
        void solidTraction(word patchName);

        //- Add fixedDisplacement boundary condition
        void fixedDisplacement(word patchName);

        //- Add solidWedge boundary condition
        void solidWedge(word patchName);

        //- Add roller boundary condition
        void roller
        (
            word rollerType,
            word rollerPatchName,
            vector rollerPosition,
            scalar rollGap
        );

        //- Add roller boundary condition
        void freeRoller
        (
            word rollerType,
            word rollerPatchName,
            vector rollerPosition,
            scalar rollGap
        );

        // Thermal boundary conditions
        //- Add fixedTemperatureGradient boundary condition
        void fixedTemperatureGradient(word patchName);

        //- Add fixedTemperature boundary condition
        void fixedTemperature(word patchName);

        //- Add thermalContactSlave boundary condition
        void thermalContactSlave
        (
            word patchName, word shadowPatchName
        );

        //- Add thermalContactMaster boundary condition
        void thermalContactMaster
        (
            word patchName, word shadowPatchName
        );

        //- Add thermalContactMaster boundary condition
        void thermalContactMaster
        (
            word patchName, wordList shadowPatchNames
        );

        //- Add thermalSymmetry boundary condition
        void thermalSymmetry(word patchName);

        //- Add thermalConvection boundary condition
        void thermalConvection(word patchName);

        //- Add thermalWedge boundary condition
        void thermalWedge(word patchName);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
