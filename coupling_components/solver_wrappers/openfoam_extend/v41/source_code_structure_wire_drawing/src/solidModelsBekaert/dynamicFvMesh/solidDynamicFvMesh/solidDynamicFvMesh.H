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
    Foam::solidDynamicFvMesh

Description
    Moving mesh for updated Lagrangian solid solver.

SourceFiles
    solidDynamicFvMesh.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef solidDynamicFvMesh_H
#define solidDynamicFvMesh_H

#include "dynamicFvMesh.H"
#include "pointMesh.H"
#include "newLeastSquaresVolPointInterpolation.H"
#include "meshSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidDynamicFvMesh Declaration
\*---------------------------------------------------------------------------*/

class solidDynamicFvMesh
:
    public dynamicFvMesh
{
    // Private Data Members

        //- Dynamic mesh dictionary
        dictionary dict_;

        //- Point mesh
        pointMesh pMesh_;

        //- Point displacement field on pMesh
        mutable pointVectorField* pointDUPtr_;

        //- Interpolor for vol to point
        newLeastSquaresVolPointInterpolation volToPointInterp_;

        //- Global face zones
        mutable labelList* globalFaceZonesPtr_;

        //- Global to local face zone point mapping
        mutable labelListList* globalToLocalFaceZonePointMapPtr_;

        //- Optional mesh smoother
        autoPtr<meshSmoother> meshSmootherPtr_;

        //- Mesh smoothing frequency i.e. 1 means every time-step, 2 means
        //  every 2 time-stpes
        int meshSmoothFreq_;


    // Private Member Functions

        // Make the point displacement increment field
        void makePointDU() const;

        // Return the point displacement increment field
        pointVectorField& pointDU();

        //- Find global face zones
        void calcGlobalFaceZones() const;

        //- Find global to local point addressing
        void calcGlobalToLocalFaceZonePointMap() const;

        //- Disallow default bitwise copy construct
        solidDynamicFvMesh(const solidDynamicFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const solidDynamicFvMesh&);


public:

    //- Runtime type information
    TypeName("solidDynamicFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit solidDynamicFvMesh(const IOobject& io);


    // Destructor

        virtual ~solidDynamicFvMesh();


    // Member Functions

        //- Get global face zones list
        const labelList& globalFaceZones() const;

        //- Get global to local face zone point addressing
        const labelListList& globalToLocalFaceZonePointMap() const;
  
        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
