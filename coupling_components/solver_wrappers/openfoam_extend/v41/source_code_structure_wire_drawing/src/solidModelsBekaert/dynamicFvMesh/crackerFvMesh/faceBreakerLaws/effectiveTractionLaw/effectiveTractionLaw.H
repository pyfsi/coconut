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
    effectiveTractionLaw

Description
    Failure class that breaks internal mesh faces and coupled faces when the
    effective traction fraction exceeds unity.

    The effectve traction fraction is the effective traction normalised with
    respect to the material strength, defined as:

        tracFracEff = (tN/tNC)^2 + (tS/tSC)^2

    where
        tNC is the strength is opening
        tSC is the strength is shear
        tN is the normal traction for tN > 0 and is zero for tN < 0
        tS is the shear traction

    Only one face is broken per call where the face with the maximum traction
    fraction is chosen if multiple faces exceed unity.

    After a break, all volFields values on the new faces are updated by
    averaging from the pervious owner and neighbour cells.

    surfaceFields will also be mapped if there is a corresponding volField.

Author
    Philip Cardif UCD

SourceFiles
    effectiveTractionLaw.C

\*---------------------------------------------------------------------------*/

#ifndef effectiveTractionLaw_H
#define effectiveTractionLaw_H

#include "faceBreakerLaw.H"
#include "crackPathLimiter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class effectiveTractionLaw Declaration
\*---------------------------------------------------------------------------*/

class effectiveTractionLaw
:
    public faceBreakerLaw
{
    // Private data

        // Patch ID of the cohesive patch
        mutable label* cohesivePatchIDPtr_;

        // Allow cracks along coupled faces
        const Switch allowCoupledFaces_;

        // List of internal faces to break
        mutable labelList* facesToBreakPtr_;

        // List of internal faces to normal flips
        mutable boolList* facesToBreakFlipPtr_;

        // List of coupled faces to break
        mutable labelList* coupledFacesToBreakPtr_;

        //- Crack path limiter
        mutable autoPtr<crackPathLimiter> pathLimiterPtr_;


    // Private Member Functions

        //- Calculate the cohesive patch ID
        void calcCohesivePatchID() const;

        //- Return cohesive patch ID
        label cohesivePatchID() const;

        //- Calculate all faces to break
        void calcAllFacesToBreak() const;

        //- Disallow default bitwise copy construct
        effectiveTractionLaw(const effectiveTractionLaw&);

        //- Disallow default bitwise assignment
        void operator=(const effectiveTractionLaw&);


public:

    //- Runtime type information
    TypeName("effectiveTraction");

    // Static data members


    // Constructors

        //- Construct from dictionary
        effectiveTractionLaw
        (
            const word& name,
            fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~effectiveTractionLaw();


    // Member Functions

        // Return list of internal faces to break
        virtual const labelList& facesToBreak() const;

        // Return list of internal faces to break normal flips
        virtual const boolList& facesToBreakFlip() const;

        // Return list of coupled faces to break
        virtual const labelList& coupledFacesToBreak() const;

        //- Delete demand driven data
        virtual void clearOut();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
