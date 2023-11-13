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
    wireStreamlines

Description
    Class in charging of finding the wire streamlines, which flow from the
    wire upstream patch to the wire downstream patch.

    We assume that the wire always flows in the positive x direction, then
    the upstream patch will be the patch with the most negative x coordinates
    on average; alternatively, the upstreamPatchName can be provided.

SourceFiles
    wireStreamlines.C

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef wireStreamlines_H
#define wireStreamlines_H

#include "polyMesh.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class wireStreamlines Declaration
\*---------------------------------------------------------------------------*/

class wireStreamlines
{
    // Private data

        //- Const reference to the wire mesh
        const polyMesh& mesh_;

        //- Upstream patch index
        label upstreamPatchID_;

        //- Point streamlines
        //  A list of point IDs for the streamlines flowing from the upstream
        //  patch to the downstream patch
        mutable List< SLList<label> >* pointStreamlinesPtr_;

        //- Cell streamlines
        //  A list of cells IDs for the streamlines flowing from the upstream
        //  patch to the downstream patch
        mutable List< SLList<label> >* cellStreamlinesPtr_;


private:

    // Private Member Functions

        //- Make the point streamlines
        void makePointStreamlines() const;

        //- Make the cell streamlines
        void makeCellStreamlines() const;

        //- Disallow default bitwise copy construct
        wireStreamlines(const wireStreamlines&);

        //- Disallow default bitwise assignment
        void operator=(const wireStreamlines&);


public:

    //- Runtime type information
    TypeName("wireStreamlines");


    // Constructors

        //- Construct from the wire mesh
        wireStreamlines(const polyMesh& wireMesh);

        //- Construct from the wire mesh and name of wire upstream patch
        wireStreamlines
        (
            const polyMesh& wireMesh,
            const word& wireUpstreamPatchName
        );

    // Destructor

        virtual ~wireStreamlines();


    // Member Functions

        //- Return a reference to the point streamlines
        List< SLList<label> >& pointStreamlines();

        //- Return a const reference to the point streamlines
        const List< SLList<label> >& pointStreamlines() const;

        //- Return a reference to the cell streamlines
        List< SLList<label> >& cellStreamlines();

        //- Return a const reference to the cell streamlines
        const List< SLList<label> >& cellStreamlines() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
