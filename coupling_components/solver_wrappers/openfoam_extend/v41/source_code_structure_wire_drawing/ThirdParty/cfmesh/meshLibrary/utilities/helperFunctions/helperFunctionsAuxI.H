/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "helperFunctionsAux.H"

#include <unistd.h>
#include <stdio.h>

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace help
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

inline int runProcess(const char* progName, const std::string logName)
{
    FILE *handle = popen(progName, "r");

    if( handle == 0 )
    {
        return 1;
    }

    FILE* filePtr = NULL;
    if( logName != "" )
        filePtr = fopen(logName.c_str(), "wa");

    char buf[64];
    size_t readn;
    while( (readn = fread(buf, 1, sizeof(buf), handle)) > 0 )
    {
        fwrite(buf, 1, readn, filePtr?filePtr:stdout);

        if( filePtr )
        {
            //- flush data to a file
            fflush(filePtr);
        }
    }

    return pclose(handle);
}

inline void checkFile(const fileName& file, const word& output)
{
    if (fileStat(file).isValid() == 0)
    {
        FatalError
            << output << nl
            << file
            << abort(FatalError);
    }
}


inline double readScalarFromDict(const word& keyword, const word& subDictionary, bool abortOnError, const fileName& file, const fileName& location)
{
    double value=-123456789;

    fileName currDir = cwd();

    if (fileStat(currDir/location/file).isValid() == 0)
    {
       if (abortOnError)
       {
           FatalError
               << "Error in readScalarFromDict --> tried to read from non-existent file/location:" << nl
               << currDir/location/file
               << abort(FatalError);
       }
       else
           return value;
    }

    Time runTime
    (
        Time::controlDictName,
        currDir.path(),
        currDir.name(),
        "system",
        "constant",
        false
    );

    IOdictionary fileDict
    (
        IOobject
        (
            file,
            location,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if( fileDict.found(subDictionary) && fileDict.isDict(subDictionary) )
    {
        const dictionary& dict = fileDict.subDict(subDictionary);
        if(dict.found(keyword))
        {
            value=readScalar(dict.lookup(keyword));
        }
        else
        {
            if (abortOnError)
            {
                FatalError
                    << "Keyword "<< keyword << " not found in sub-dictionary " << subDictionary << " in:" << nl
                    << currDir/location/file
                    << abort(FatalError);
            }
            else
                return value;
        }
    }
    else
    {
        if (abortOnError)
        {
            FatalError
                << "Sub-dictionary "<< subDictionary << " not found in:" << nl
                << currDir/location/file
                << abort(FatalError);
        }
        else
            return value;
    }

    return value;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
