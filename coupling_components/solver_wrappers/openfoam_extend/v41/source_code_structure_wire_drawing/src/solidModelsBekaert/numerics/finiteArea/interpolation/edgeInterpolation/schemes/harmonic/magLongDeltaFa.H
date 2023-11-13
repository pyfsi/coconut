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
    Foam::magLongDeltaFa

Description
    Long distance to the neighbour via edge centre

SourceFiles
    magLongDeltaFa.C

\*---------------------------------------------------------------------------*/

#ifndef magLongDeltaFa_H
#define magLongDeltaFa_H

#include "MeshObject.H"
#include "faMesh.H"
#include "edgeFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                          Class magLongDeltaFa Declaration
\*---------------------------------------------------------------------------*/

class magLongDeltaFa
:
    public MeshObject<faMesh, magLongDeltaFa>
{
    // Private data

        //- Magnitude of the normal distances to the edge
        mutable edgeScalarField* magLongDeltaFaPtr_;


    // Private member functions

        //- Construct long distance
        void makeMagLongDistance() const;

        //- Construct long distance for patch
        tmp<scalarField> calcMagLongDistance
        (
            const label patchi
        ) const;

        //- Clear data
        void clearOut() const;


public:

    // Declare name of the class and its debug switch
    TypeName("magLongDeltaFa");


    // Constructors

        //- Construct given an faMesh
        explicit magLongDeltaFa(const faMesh&);


    // Destructor

        virtual ~magLongDeltaFa();


    // Member functions

        //- Return reference to the magnitude of the long distance
        const edgeScalarField& magDelta() const;

        //- Return reference to the magnitude of the long distance
        //- for patch i
        const scalarField& magDelta(const label patchi) const;

        //- Update after mesh motion:
        //  Delete data when the mesh moves
        virtual bool movePoints() const;

        //- Update after topo change:
        //  Delete data when mesh changes
        virtual bool updateMesh(const mapPolyMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
