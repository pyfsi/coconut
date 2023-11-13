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
    Foam::resizeGlobalFaceZonesFvMesh

Description
    Dynamic mesh class that creates global face zones in parallel for each
    contact patch.
    The globalFaceZone class does all the work.

SourceFiles
    resizeGlobalFaceZonesFvMesh.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef resizeGlobalFaceZonesFvMesh_H
#define resizeGlobalFaceZonesFvMesh_H

#include "dynamicFvMesh.H"
#include "boundBox.H"
#include "globalFaceZone.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                           Class resizeGlobalFaceZonesFvMesh Declaration
\*---------------------------------------------------------------------------*/

class resizeGlobalFaceZonesFvMesh
:
    public dynamicFvMesh
{
    // Private Data Members

        //- Settings dictionary
        dictionary dict_;

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
        resizeGlobalFaceZonesFvMesh(const resizeGlobalFaceZonesFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const resizeGlobalFaceZonesFvMesh&);


public:

    //- Runtime type information
    TypeName("resizeGlobalFaceZonesFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit resizeGlobalFaceZonesFvMesh(const IOobject& io);


    // Destructor

        virtual ~resizeGlobalFaceZonesFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
