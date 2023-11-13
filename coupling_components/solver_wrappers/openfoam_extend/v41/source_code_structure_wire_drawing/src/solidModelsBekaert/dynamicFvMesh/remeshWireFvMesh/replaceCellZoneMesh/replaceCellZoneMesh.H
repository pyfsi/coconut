/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    Foam::replaceCellZoneMesh

Description
    Remesh the given cellZone and create a new mesh containing the newly meshed
    cellZone and the remainder of the original mesh.

Author
    Philip Cardiff, UCD

SourceFiles
    replaceCellZoneMesh.C

\*---------------------------------------------------------------------------*/

#ifndef replaceCellZoneMesh_H
#define replaceCellZoneMesh_H

#include "fvMesh.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declarations
    class mergePolyMesh;

/*---------------------------------------------------------------------------* \
                      Class replaceCellZoneMesh Declaration
\*---------------------------------------------------------------------------*/

class replaceCellZoneMesh
{
    // Private Data

        //- Store a reference to the oldMesh (may be decomposed)
        fvMesh& oldMesh_;

        //- Old global mesh
        mutable autoPtr<fvMesh> oldGlobalMeshPtr_;

        //- Cell zone name
        const word cellZoneName_;

        //- Upstream patch name
        const word& upstreamPatchName_;

        //- Upstream faceZone name
        const word& upstreamFaceZoneName_;

        //- Downstream patch name
        const word& downstreamPatchName_;

        //- Downstream faceZone name
        const word& downstreamFaceZoneName_;

        //- Store a reference to the remesh dictionary
        const dictionary& dict_;

        //- New cellZone mesh
        mutable autoPtr<mergePolyMesh> newCellZoneMeshPtr_;

        //- New global mesh (including new cellZone mesh)
        mutable autoPtr<mergePolyMesh> newGlobalMeshPtr_;

        //- Store the maxCellSize for remeshing
        mutable scalar maxCellSize_;

        //- Are we running in parallel
        mutable bool parRun_;

        //- Global time (point to serial case duirng parallel runs)
        mutable autoPtr<Time> globalTimePtr_;

    // Private Member Functions

        //- Create the old global mesh
        virtual void createOldGlobalMesh() const;

        //- Create the new cellZone mesh
        virtual void createNewCellZoneMesh() const;

        //- Create the new global mesh (including new cellZone mesh)
        virtual void createNewGlobalMesh() const;

        //- Calculate the maxCellSize
        virtual void calcMaxCellSize() const;

        //- Write an STL of the boundary to disk
        const word writeBoundarySTL(const fvMesh& mesh) const;

        //- Suspend OpenFOAM MPI
        virtual void suspendMPI() const;

        //- Resume OpenFOAM MPI
        virtual void resumeMPI() const;

        //- Force all procs to synchronise
        virtual void syncMPI() const;

        //- Make the global time
        virtual void makeGlobalTime() const;

        //- Return a const-reference to the global time
        virtual const Time& globalTime() const;

        //- Disallow default bitwise copy construct
        replaceCellZoneMesh(const replaceCellZoneMesh&);

        //- Disallow default bitwise assignment
        void operator=(const replaceCellZoneMesh&);


public:

    //- Runtime type information
    TypeName("replaceCellZoneMesh");

    // Constructor

        //- Construct from database
        explicit replaceCellZoneMesh
        (
            fvMesh& mesh,
            const word& cellZoneName,
            const word& upstreamPatchName,
            const word& upstreamFaceZoneName,
            const word& downstreamPatchName,
            const word& downstreamFaceZoneName,
            const dictionary& dict
        );


    // Destructor

        virtual ~replaceCellZoneMesh();


    // Member Functions

        // Access functions

            //- Return const reference access to old mesh
            virtual const fvMesh& oldMesh() const
            {
                return oldMesh_;
            }

            //- Return const reference access to old global mesh
            virtual const fvMesh& oldGlobalMesh() const;

            //- Return non-const reference access to old global mesh
            virtual fvMesh& oldGlobalMesh();

            //- Return non-const reference access to new cellZone mesh
            virtual const polyMesh& newCellZoneMesh() const;

            //- Return non-const reference access to new global mesh
            virtual const polyMesh& newGlobalMesh() const;

            //- Return const reference access to old global pointMesh
            virtual const pointMesh& oldGlobalPointMesh() const
            {
                return pointMesh::New(oldGlobalMesh());
            }

            //- Return const reference access to new global pointMesh
            virtual const pointMesh& newGlobalPointMesh() const
            {
                return pointMesh::New(newGlobalMesh());
            }

            //- Return the maxCellSize
            virtual scalar maxCellSize() const;

        // Evaluation functions

            //- Clear out all demand driven data
            virtual void clearOut();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #   include "replaceCellZoneMeshTemplates.C"
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
