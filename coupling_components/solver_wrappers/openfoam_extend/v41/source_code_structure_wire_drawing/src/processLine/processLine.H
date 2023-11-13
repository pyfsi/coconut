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
    processLine

Description
    Process line class in charge of passes.

SourceFiles
    processLine.C

Author
    Philip Cardiff, UCD.
    Peter De Jaeger, Bekaert.

\*---------------------------------------------------------------------------*/

#ifndef processLine_H
#define processLine_H

#include "foamTime.H"
#include "IOdictionary.H"
#include "typeInfo.H"
#include "processPass.H"
#include "PtrList.H"
#include "dataHandler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class processLine Declaration
\*---------------------------------------------------------------------------*/

class processLine
{
    // Private data

        //- Name of the process line
        const word processLineName_;

        //- Pseudo controlDict for the processLine
        mutable dictionary controlDict_;

        //- Pseudo runTime for the processLine
        mutable Time* runTimePtr_;

        //- Process line dictionary
        mutable IOdictionary* dictPtr_;

        //- List of passes
        PtrList<processPass> passes_;

        //- List of inter-passes
        //  These control passing information between passes
        // Idea for the future
        //PtrList<interPass> interPasses_;

        //- dataHandler object for inter-pass information exchange
        dataHandler dataContainer_;

private:

    // Private Member Functions

        //- Make the runTime
        void makeRunTime() const;

        //- Return a reference to the runTime
        const Time& runTime() const;

        //- Make the dictionary
        void makeDict() const;

        //- Return a reference to the process line dictionary
        dictionary& dict() const;

        //- Make the passes
        void makePasses();

        //- Return references to passes
        PtrList<processPass>& passes();

        //- Copy the wire region mesh and fields from the latest time in the
        //  fromPass to the constant and 0 directories of the toPass
        void copyWireMeshAndFields
        (
            const label fromPassIndex, const label toPassIndex
        );

        //- Remesh the wire at the end of a given pass
        //  where we give the time object for the pass of interest as an
        //  argument
        void remeshWire(const label passID);

        //- Disallow default bitwise copy construct
        processLine(const processLine&);

        //- Disallow default bitwise assignment
        void operator=(const processLine&);


public:

    //- Runtime type information
    TypeName("processLine");


    // Constructors

        //- Construct from name of dictionary
        processLine(const fileName& processLineDir, const bool restart);

    // Destructor

        virtual ~processLine();


    // Member Functions

        //- Return reference to the dataContainer
        dataHandler& dataContainer()
        {
            return dataContainer_;
        }

        //- Return a reference to the general settings sub-dictionary
        dictionary& generalSettings();

        //- Setup the process line but do not run the cases
        void setup();

        //- Setup and run the process line
        void setupAndRun();

        //- Perform required post-processing steps
        void postProcess();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
