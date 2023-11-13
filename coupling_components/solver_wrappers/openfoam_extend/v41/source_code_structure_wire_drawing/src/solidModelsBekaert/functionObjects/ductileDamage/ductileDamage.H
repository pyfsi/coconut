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
    ductileDamage

Description
    FunctionObject that outputs damage fields based on a number of
    different ductile damage/failure models, namely:
        - Lemaitre
        - Cockcroft-Lathem

    It assumed that the stress tensor is called "sigma" or "sigmaCauchy".

    @verbatim
    ductileDamage1
    {
        // Type of functionObject
        type ductileDamage;

        // Where to load it from (if not already in solver)
        // Note: this is not required if the solver already loads this library.
        functionObjectLibs ("libsolidModelsBekaert.so");

        // Optional: specifiy mesh region for multi-region solvers, defaults to
        // region0
        region    region0;

        // Lemaitre damage parameters
        Lemaitre
        {
            // Optional: if not specified, all parameters will default to values
            // taken from Masse thesis
            h                 h [ 1 -1 -2 0 0 0 0 ] 0.2;
            S0                S0 [ 1 -1 -2 0 0 0 0 ] 5.23;
            b                 b [ 0 0 0 0 0 0 0 ] 3.14;
            epsilonD          epsilonD [ 0 0 0 0 0 0 0 ] 0.1;
        }

        // Cockcroft-Latham parameter
        Cockcroft-Latham
        {
            // Cockcroft-Latham parameter that indicates that maximum plastic
            // energy that can be abosrbed before damage occurs
            // More details can be found in the Masse thesis
            Vd                Vd [ 1 -1 -2 0 0 0 0 ] 1.0;
        }
    );
    @endverbatim

Author
    Philip Cardiff UCD

SourceFiles
    ductileDamage.C

\*---------------------------------------------------------------------------*/

#ifndef ductileDamage_H
#define ductileDamage_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ductileDamage Declaration
\*---------------------------------------------------------------------------*/

class ductileDamage
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& runTime_;

        //- Reference to mesh region
        const fvMesh& mesh_;

        //- Accumulated Lemaitre damage internal field
        scalarField oldDamageLemaitreI_;

        //- Accumulated Cockcroft-Latham damage internal field
        scalarField oldDamageCockcroftLathamI_;

        // Lemaitre damage parameters
        // For definition of these parameters, see:
        // Masse, Study and optimization of high carbon steel flat wires, PhD
        // thesis, 2010.
        const dimensionedScalar h_;
        const dimensionedScalar S0_;
        const dimensionedScalar b_;
        const dimensionedScalar epsilonD_;

        // Cockcroft-Latham damage parameter
        // Indicates that amount of plastic energy that can be abosrbed before
        // damage
        const dimensionedScalar cockcroftLathamVd_;

    // Private Member Functions

        //- Write data
        bool writeData();

        //- Write Lemaitre damage field
        bool writeDamageLemaitre();

        //- Write Cockcroft-Latham damage field
        bool writeDamageCockcroftLatham();

        //- Return the equivalent plastic strain field
        //  This is to correctly account for the possibility of multiple
        //  materials
        tmp<volScalarField> epsilonPEq() const;

        //- Return the increment of equivalent plastic strain field
        //  This is to correctly account for the possibility of multiple
        //  materials
        tmp<volScalarField> DEpsilonPEq() const;

        //- Disallow default bitwise copy construct
        ductileDamage
        (
            const ductileDamage&
        );

        //- Disallow default bitwise assignment
        void operator=(const ductileDamage&);


public:

    //- Runtime type information
    TypeName("ductileDamage");


    // Constructors

        //- Construct from components
        ductileDamage
        (
            const word& name,
            const Time&,
            const dictionary&
        );


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
