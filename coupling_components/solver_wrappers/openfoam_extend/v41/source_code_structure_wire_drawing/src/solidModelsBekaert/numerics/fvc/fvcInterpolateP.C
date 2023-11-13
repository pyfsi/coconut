/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fvcInterpolateP.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField<Type, pointPatchField, pointMesh>& pf
)
{
    const fvMesh& mesh =
        pf.db().parent().objectRegistry::lookupObject<fvMesh>
        (
            pf.mesh().mesh().name()
        );

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                pf.name() + "f",
                pf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("0", pf.dimensions(), pTraits<Type>::zero)
        )
    );
    Field<Type>& sfI = tsf().internalField();

    const vectorField& points = mesh.points();
    const faceList& faces = mesh.faces();

    const Field<Type>& pfI = pf.internalField();

    forAll(sfI, faceI)
    {
        const face& curFace = faces[faceI];

        // Calculate face centre value as weighted average of facePoint values
        sfI[faceI] = curFace.average(points, pfI);
    }

    forAll(tsf().boundaryField(), patchI)
    {
        // Directly use the vol field face values
        tsf().boundaryField()[patchI] = vf.boundaryField()[patchI];
    }

    tsf().correctBoundaryConditions();

    return tsf;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >interpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>& pf
)
{
    // PC Fix 06-April-2016: lookup region name from pf mesh
    const fvMesh& mesh =
        pf.db().parent().objectRegistry::lookupObject<fvMesh>
        (
            pf.mesh().mesh().name()
        );

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                pf.name() + "f",
                pf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("0", pf.dimensions(), pTraits<Type>::zero)
        )
    );
    Field<Type>& sfI = tsf().internalField();

    const vectorField& points = mesh.points();
    const faceList& faces = mesh.faces();

    const Field<Type>& pfI = pf.internalField();

    forAll(sfI, faceI)
    {
        const face& curFace = faces[faceI];

        // Calculate face centre value as weighted average of facePoint values
        sfI[faceI] = curFace.average(points, pfI);
    }

    forAll(tsf().boundaryField(), patchI)
    {
        Field<Type>& patchS = tsf().boundaryField()[patchI];

        forAll(patchS, faceI)
        {
            label globalFaceID = mesh.boundaryMesh()[patchI].start() + faceI;

            const face& curFace = mesh.faces()[globalFaceID];

            patchS[faceI] = curFace.average(points, pfI);
        }
    }

    tsf().correctBoundaryConditions();

    return tsf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
