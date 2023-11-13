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
    Foam::dynamicRefineGlobalFaceZoneFvMesh

Description
    Dynamic mesh class that derives from dynamicRefineFvMesh with the added
    capability to regenerating the global face zones on the fly in parallel
    from patches of type solidContact.
    Also optional bounding box to limit the face zones.

SourceFiles
    dynamicRefineGlobalFaceZoneFvMesh.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef dynamicRefineGlobalFaceZoneFvMesh_H
#define dynamicRefineGlobalFaceZoneFvMesh_H

#include "dynamicRefineFvMesh.H"
#include "volFields.H"
#include "boundBox.H"
#include "globalFaceZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicRefineGlobalFaceZoneFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicRefineGlobalFaceZoneFvMesh
:
    public dynamicRefineFvMesh
{
    // Private Data Members

        //- Settings dictionary
        dictionary dict_;

        //- Discretisation residual error field for slecting cells to refine
        volScalarField residualError_;

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
        dynamicRefineGlobalFaceZoneFvMesh
        (
            const dynamicRefineGlobalFaceZoneFvMesh&
        );

        //- Disallow default bitwise assignment
        void operator=(const dynamicRefineGlobalFaceZoneFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicRefineGlobalFaceZoneFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicRefineGlobalFaceZoneFvMesh(const IOobject& io);


    // Destructor

        virtual ~dynamicRefineGlobalFaceZoneFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
