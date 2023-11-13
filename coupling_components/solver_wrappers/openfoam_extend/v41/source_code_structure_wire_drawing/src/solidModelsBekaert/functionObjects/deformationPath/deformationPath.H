/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    deformationPath

Description
    FunctionObject which extracts the deformation gradient history for
    every cell at every timestep. This will create a separate file for
    each cell so caution should be taken with large meshes.
    
Author
    Michael Clancy UCD
    Philip Cardiff UCD

SourceFiles
    deformationPath.C

\*---------------------------------------------------------------------------*/

#ifndef deformationPath_H
#define deformationPath_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

#include <sstream>
#include <string>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class deformationPath Declaration
\*---------------------------------------------------------------------------*/

class deformationPath
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;


        //- Reference to the mesh and the cells we
        // are interested in tabulating the deformation
        // gradient for
        const fvMesh& mesh_;

        //- the number of cells in the entire mesh
        const label nCells_;

        labelList cellsOfInterest_;

    // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        deformationPath
        (
            const deformationPath&
        );

        //- Disallow default bitwise assignment
        void operator=(const deformationPath&);


public:

    //- Runtime type information
    TypeName("deformationPath");


    // Constructors

        //- Construct from components
        deformationPath
        (
            const word& name,
            const Time&,
            const dictionary&
        );

    // Destructors
    //virtual ~deformationPath();


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
#if FOAMEXTEND > 40
        virtual bool execute(const bool forceWrite);
#else
        virtual bool execute();
#endif

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
