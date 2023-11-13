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
    dataFriction

Description
    Friction law which interpolates friction coefficient from experimental data
    based on entrainment speed, slide-roll-ration, pressure and temperature.

SourceFiles
    dataFriction.C

Author
    Vanja Skuric, UNIZAG FSB

\*---------------------------------------------------------------------------*/

#ifndef dataFriction_H
#define dataFriction_H

#include "frictionLaw.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class dataFriction Declaration
\*---------------------------------------------------------------------------*/

class dataFriction
:
    public frictionLaw
{
    // Private data

        //- Reference to fvMesh
        const fvMesh& mesh_;

        //- Friction law dict
        dictionary frictionLawDict_;

        //- Coulomb friction coefficient
        scalar frictionCoeff_;

        //- List of parameters for the specific measured data
        //- Order of the parameters is based on the parameter priority
        //- specified in the interpolation dictionary
        scalarListList parameterList_;

        //- Primary parameter list
        DynamicList<scalar> primaryParameters_;

        //- Secondary parameter list divided in primary parameter groups
        DynamicList<DynamicList<scalar> > secondaryParameters_;

        //- Tertiary parameter list divided in primary and secondary parameter
        //- groups
        DynamicList<DynamicList<DynamicList<scalar> > > tertiaryParameters_;

        //- Data list indices divided the same way as tertiary parameters
        DynamicList<DynamicList<DynamicList<scalar> > > dataIndices_;

        //- Experimental data (friction coefficients) data list
        List<List<Tuple2<scalar, scalar> > > frictionCoeffList_;

        //- Use nearest property value bool list
        boolList useNearest_;

        //- Data properties prioritiy
        label SRRPriority_;
        label pressurePriority_;
        label temperaturePriority_;

        //- Property values if properties are fixed
        Tuple2<bool, scalar> fixedSRR_;
        Tuple2<bool, scalar> fixedPressure_;
        Tuple2<bool, scalar> fixedTemperature_;

        //- Temperature field pointer
        const volScalarField* TPtr_;

    // Private Member Functions

        //- Extrapolate friction coefficient from measured data
        void calcFrictionCoeff
        (
            const scalar entrSpeed,
            const scalar SRR,
            const scalar pressure,
            const scalar temperature
        );

        //- Get high and low values and their indices from list around the give
        //- value
        //- [0][0] = low value index
        //- [0][1] = low value
        //- [1][0] = high value index
        //- [1][1] = high value
        List<Tuple2<label, scalar> > highLowValues
        (
            const scalar value,
            const scalarList& valueList,
            const bool useNearest
        );

        //- Lookup temperature field
        void lookupTemperatureField();

        //- Return temperature field
        const volScalarField& T()
        {
            if (!TPtr_)
            {
                lookupTemperatureField();
            }

            return *TPtr_;
        }

        //- Disallow default bitwise copy construct
        //dataFriction(const dataFriction&);

        //- Disallow default bitwise assignment
        void operator=(const dataFriction&);

public:

    //- Runtime type information
    TypeName("data");

    // Static data members


    // Constructors

        //- Construct from dictionary
        dataFriction
        (
            const word& name,
            const frictionContactModel& fricModel,
            const dictionary& dict
        );

        //- Construct as a copy
        dataFriction(const dataFriction& fricLaw);

        //- Construct and return a clone
        virtual autoPtr<frictionLaw> clone() const
        {
            return
                autoPtr<frictionLaw>
                (
                    new dataFriction(*this)
                );
        }



    // Destructor

        virtual ~dataFriction();


    // Member Functions

        //- Return traction that would cause slipping
        virtual scalar slipTraction
        (
            const scalar contactPressure,         // Contact pressure
            const vector& faceSlip,               // Slip vector
            const vector& slaveFaceVelocity,      // Velocity of slave face
            const vector& masterFaceVelocity,     // Velocity of master face
            const label slavePatchIndex,          // Slave patch index
            const label faceIndex                 // Local slave face ID
        );

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
