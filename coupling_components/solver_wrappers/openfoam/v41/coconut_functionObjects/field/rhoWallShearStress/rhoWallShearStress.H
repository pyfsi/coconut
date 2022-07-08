/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

Group
    grpForcesFunctionObjects

Description
    This function object evaluates and outputs the shear stress at wall
    patches.  The result is written as a volVectorField to time directories as
    field 'rhoWallShearStress'

        \f[
            Stress = R \dot n
        \f]

    where
    \vartable
        R       | stress tensor
        n       | patch normal vector (into the domain)
    \endvartable

    The shear stress (symmetrical) tensor field is retrieved from the
    turbulence model.  All wall patches are included by default; to restrict
    the calculation to certain patches, use the optional 'patches' entry.

    Example of function object specification:
    \verbatim
    rhoWallShearStress1
    {
        type        rhoWallShearStress;
        libs ("libfieldFunctionObjects.so");
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

See also
    Foam::functionObject
    Foam::functionObjects::writeFiles
    Foam::functionObjects::pressureTools
    Foam::functionObjects::timeControl

SourceFiles
    rhoWallShearStress.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_rhoWallShearStress_H
#define functionObjects_rhoWallShearStress_H

#include "writeFiles.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class rhoWallShearStress Declaration
\*---------------------------------------------------------------------------*/

class rhoWallShearStress
:
    public writeFiles
{
protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the shear stress
        void calcShearStress
        (
            const fvMesh& mesh,
            const volSymmTensorField& Reff,
            volVectorField& shearStress
        );


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        rhoWallShearStress(const rhoWallShearStress&);

        //- Disallow default bitwise assignment
        void operator=(const rhoWallShearStress&);


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


    //- Destructor
    virtual ~rhoWallShearStress();


    // Member Functions

        //- Read the rhoWallShearStress data
        virtual bool read(const dictionary&);

        //- Calculate the wall shear-stress
        virtual bool execute();

        //- Write the wall shear-stress
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
