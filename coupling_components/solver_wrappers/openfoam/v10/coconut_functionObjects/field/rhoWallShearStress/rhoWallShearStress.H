/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::functionObjects::rhoWallShearStress

Description
    Calculates and write the shear-stress at wall patches as
    the volVectorField field 'rhoWallShearStress'.

        \f[
            Stress = R \dot n
        \f]

    where
    \vartable
        R       | stress tensor
        n       | patch normal vector (into the domain)
    \endvartable

    The shear-stress symmetric tensor field is retrieved from the turbulence
    model.  All wall patches are included by default; to restrict the
    calculation to certain patches, use the optional 'patches' entry.

    Example of function object specification:
    \verbatim
    rhoWallShearStress1
    {
        type        rhoWallShearStress;
        libs        ("libfieldFunctionObjects.so");
        ...
        patches     (".*Wall");
    }
    \endverbatim

Usage
    \table
        Property | Description               | Required    | Default value
        type     | type name: rhoWallShearStress | yes        |
        patches  | list of patches to process | no         | all wall patches
    \endtable

Note
    Writing field 'rhoWallShearStress' is done by default, but it can be overridden
    by defining an empty \c objects list. For details see writeLocalObjects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    rhoWallShearStress.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_rhoWallShearStress_H
#define functionObjects_rhoWallShearStress_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class rhoWallShearStress Declaration
\*---------------------------------------------------------------------------*/

class rhoWallShearStress
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{
    // Private Data

        //- The name of the phase
        word phaseName_;

        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the shear-stress
        tmp<volVectorField> calcShearStress
        (
            const volSymmTensorField& tau
        );


public:

    //- Runtime type information
    TypeName("rhoWallShearStress");


    // Constructors

        //- Construct from Time and dictionary
        rhoWallShearStress
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        rhoWallShearStress(const rhoWallShearStress&) = delete;


    //- Destructor
    virtual ~rhoWallShearStress();


    // Member Functions

        //- Read the rhoWallShearStress data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the wall shear-stress
        virtual bool execute();

        //- Write the wall shear-stress
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rhoWallShearStress&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
