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
    Foam::cartesian2DFvMesh

Description
    Dynamic mesh class that remeshes the domain using cfMesh cartesian2DMesh
    and globalFaceZones are recreated.

SourceFiles
    cartesian2DFvMesh.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef cartesian2DFvMesh_H
#define cartesian2DFvMesh_H

#include "dynamicFvMesh.H"
#include "boundBox.H"
#include "globalFaceZone.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                           Class cartesian2DFvMesh Declaration
\*---------------------------------------------------------------------------*/

class cartesian2DFvMesh
:
    public dynamicFvMesh
{
    // Private Data Members

        //- Settings dictionary
        dictionary dict_;

        //- cartesianMesh settings
        const scalar maxCellSize_;

        //- Optional bounding box to restrict size of face zones
        const boundBox bb_;

        //- globalFaceZone class to create face zones in parallel
        globalFaceZone gfz_;

        //- Field indicating boundary faces currently in a global face zone
        //  This field is for visualisation e.g. when checking the bounding box
        //  size
        volScalarField indicator_;

        // Optional: enable/disable update of indicator field
        Switch updateIndicator_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        cartesian2DFvMesh(const cartesian2DFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const cartesian2DFvMesh&);


public:

    //- Runtime type information
    TypeName("cartesian2DFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit cartesian2DFvMesh(const IOobject& io);


    // Destructor

        virtual ~cartesian2DFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
