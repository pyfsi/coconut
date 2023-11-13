/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
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

Description
    Calulate the face centres and areas.

    Calculate the centre by breaking the face into triangles using the face
    centre and area-weighted averaging their centres.  This method copes with
    small face-concavity.

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "helperFunctions.H"

# ifdef USE_OMP
# include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcFaceCentresAndAreas() const
{
    if( faceCentresPtr_ || faceAreasPtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcFaceCentresAndAreas() const")
            << "Face centres or face areas already calculated"
            << abort(FatalError);
    }

    const pointFieldPMG& points = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    faceCentresPtr_ = new vectorLongList(faces.size());
    vectorLongList& fCtrs = *faceCentresPtr_;

    faceAreasPtr_ = new vectorLongList(faces.size());
    vectorLongList& fAreas = *faceAreasPtr_;

    makeFaceCentresAndAreas(points, fCtrs, fAreas);
}

void polyMeshGenAddressing::makeFaceCentresAndAreas
(
    const pointFieldPMG& points,
    vectorLongList& fCtrs,
    vectorLongList& fAreas
) const
{
    const faceListPMG& faces = mesh_.faces();

    # ifdef USE_OMP
    # pragma omp parallel for if( faces.size() > 1000 ) schedule(guided, 100)
    # endif
    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        fCtrs[faceI] = help::faceCentre(points, f);
        fAreas[faceI] = help::faceAreaVector(points, f);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorLongList& polyMeshGenAddressing::faceCentres() const
{
    if( !faceCentresPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const vectorField& polyMeshGenAddressing::faceCentres() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcFaceCentresAndAreas();
    }

    return *faceCentresPtr_;
}

const vectorLongList& polyMeshGenAddressing::faceAreas() const
{
    if( !faceAreasPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const vectorField& polyMeshGenAddressing::faceAreas() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcFaceCentresAndAreas();
    }

    return *faceAreasPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
