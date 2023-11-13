/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Class
    rollerSurfaceCutter

Description


SourceFiles
    rollerSurfaceCutter.C

\*---------------------------------------------------------------------------*/

#ifndef rollerSurfaceCutter_H
#define rollerSurfaceCutter_H

#include "plane.H"
#include "PtrList.H"
#include "List.H"
#include "DynList.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class triSurf;

/*---------------------------------------------------------------------------*\
                    Class rollerSurfaceCutter Declaration
\*---------------------------------------------------------------------------*/

class rollerSurfaceCutter
{
    // Private data
        //- reference to the surface mesh
        triSurf& rollSurf_;

        //- const reference to the cutting planes
        const plane& plane_;

        //- patch name of the triangles in the plane
        const word patchName_;

        //- classify vertices according to their position
        List<direction> pointType_;

    // Private member functions
        //- classify surface vertices
        bool classifyVertices();

        //- create new triangles from cut edges
        void newTrianglesFromCutEdges(const std::map<label, label>& cutEdges);

        //- sort plane points in a cyclic order
        void sortPlanePoints(DynList<label>& planePoints) const;

        //- patch index for a created plane
        label patchForPlane() const;

        //- create intersection points and triangles
        void createIntersections();

        //- remove triangles and vertices that are not needed any more
        void removeOutsideTriangles();

        //- Disallow default bitwise copy construct
        rollerSurfaceCutter(const rollerSurfaceCutter&);

        //- Disallow default bitwise assignment
        void operator=(const rollerSurfaceCutter&);

    // Enumerators

        enum pointTypes_
        {
            NONE = 0,
            OUTSIDE = 1,
            INSIDE = 4,
            ATPLANE = 16,
        };

public:

    // Constructors

        //- Construct from triSurf and planes
        rollerSurfaceCutter
        (
            triSurf& rollSurf,
            const plane& cutPlane,
            const word patchName = ""
        );


    // Destructor

        ~rollerSurfaceCutter();


    // Member Functions
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
