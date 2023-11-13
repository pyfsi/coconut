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
    dataHandler

Description
    Process line class in charge of general data.

SourceFiles
    dataHandler.C

Author
    Philip Cardiff, UCD.
    Peter De Jaeger, Bekaert.

\*---------------------------------------------------------------------------*/

#ifndef dataHandler_H
#define dataHandler_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "scalar.H"
#include "vector.H"
#include "scalarList.H"

#include <stdio.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class dataHandler Declaration
\*---------------------------------------------------------------------------*/

class dataHandler
{
    // Private data

        //- dataContainer dictionary
        dictionary dataContainerDict_;

        //- Process line progress dict to allow process line restarts
        IOdictionary processLineProgressDict_;

private:

    // Private Member Functions


        //- Disallow default bitwise copy construct
//        dataHandler(const dataHandler&);

        //- Disallow default bitwise assignment
        void operator=(const dataHandler&);


public:

    //- Runtime type information
    TypeName("dataHandler");


    // Constructors

        //- Constructor using Time from the processLine object
        dataHandler(const Time& runTime);


    // Destructor

        virtual ~dataHandler();


    // Member Functions

        //- Return reference to dataContainerDict_
        dictionary& dataContainerDict()
        {
            return dataContainerDict_;
        }

        //- Add entry to data container dictionary
        template<class T>
        void addKey(const word keyName, const T keyValue);

        //- Add pass specific entry to data container dictionary
        template<class T>
        void addKey(const word keyName, const T keyValue, const int passNr);

        //- Interface function to the returnKey templates in returnKey
        // namespace
        template<class T>
        T returnKey(const word keyName);

        //- Interface function to the returnKey templates in returnKey
        // namespace
        template<class T>
        T returnKey(const word keyName, const int passNr);

        //- Print data entry from the dataHandler dictionary
        void printData();

        //- Run system command
        void runSystemCommand
        (
            const word& systemCommand,
            const word& logName,
            const word& functionCalling,
            const bool checkForSystemError = true
        );

        //- Return reference to the progress dict
        IOdictionary& processLineProgressDict()
        {
            return processLineProgressDict_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        namespace returnKey Functions
\*---------------------------------------------------------------------------*/

namespace returnKey
{

    //- Return entry from data container dictionary
    template<class T>
    T returnKey(const word keyName, Foam::dataHandler* dataContainerP);

    //- Specialisation for scalar data
    template<> inline
    scalar returnKey(const word keyName, Foam::dataHandler* dataContainerP);

    //- Return pass specific entry from data container dictionary
    template<class T>
    T returnKey
    (
        const word keyName,
        Foam::dataHandler* dataContainerP,
        const int passNr
    );

    //- Specialisation for scalar data
    template<> inline
    scalar returnKey
    (
        const word keyName,
        Foam::dataHandler* dataContainerP,
        const int passNr
    );

} // End namespace returnKey


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#ifdef NoRepository
#   include "dataHandlerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
