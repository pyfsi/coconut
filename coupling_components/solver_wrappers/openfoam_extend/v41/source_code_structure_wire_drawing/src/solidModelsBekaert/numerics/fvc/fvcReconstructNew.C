/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvcReconstructNew.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcSurfaceIntegrate.H"
//#include "extrapolatedCalculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// template<class Type>
// tmp
// <
//     GeometricField
//     <
//         typename outerProduct<vector,Type>::type, fvPatchField, volMesh
//     >
// >
// reconstructNew
// (
//     const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
// )
// {
//     typedef typename outerProduct<vector, Type>::type GradType;

//     const fvMesh& mesh = ssf.mesh();

//     surfaceVectorField SfHat(mesh.Sf()/mesh.magSf());

//     tmp< GeometricField<GradType, fvPatchField, volMesh> > treconField
//     (
//         new GeometricField<GradType, fvPatchField, volMesh>
//         (
//             IOobject
//             (
//                 "volIntegrate("+ssf.name()+')',
//                 ssf.instance(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::NO_WRITE
//             ),
//             inv(surfaceSum(SfHat*mesh.Sf()))&surfaceSum(SfHat*ssf),
//             //extrapolatedCalculatedFvPatchField<GradType>::typeName
//             calculatedFvPatchField<GradType>::typeName
//         )
//     );

//     treconField.ref().correctBoundaryConditions();

//     return treconField;
// }


// Taken from http://bugs.openfoam.org/file_download.php?file_id=684&type=bug
template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
reconstructNew
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    surfaceVectorField faceVols
    (
        //mesh.Sf()/(mesh.magSf()*mesh.nonOrthDeltaCoeffs())
        mesh.Sf()/(mesh.magSf()*mesh.deltaCoeffs())
    );

    const surfaceScalarField& weights = mesh.weights();

    volTensorField T
    (
        IOobject
        (
            "T",
            ssf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("0", dimVolume, pTraits<tensor>::zero),
        zeroGradientFvPatchField<tensor>::typeName
        //"zeroGradient"
    );

    GeometricField <GradType, fvPatchField, volMesh> A
    (
        IOobject
        (
            "A",
            ssf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<GradType>
        (
            "0",
            dimLength * ssf.dimensions(),
            pTraits<GradType>::zero
        ),
        zeroGradientFvPatchField<GradType>::typeName
        //"zeroGradient"
    );

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    forAll(faceVols, facei)
    {
        T[owner[facei]] +=
            mesh.Sf()[facei]*faceVols[facei]*(1.0 - weights[facei]);
        T[neighbour[facei]] +=
            mesh.Sf()[facei]*faceVols[facei]*weights[facei];

        A[owner[facei]] += faceVols[facei]*ssf[facei]*(1.0 - weights[facei]);
        A[neighbour[facei]] += faceVols[facei]*ssf[facei]*weights[facei];
    }

    forAll(faceVols.boundaryField(), patchi)
    {
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
        const vectorField& pFaceVols = faceVols.boundaryField()[patchi];
        const scalarField& pWeights = weights.boundaryField()[patchi];
        const Field<Type>& pssf = ssf.boundaryField()[patchi];
        const unallocLabelList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        forAll(faceVols.boundaryField()[patchi], facei)
        {
            //PDJ 01-02-2018: removed const for addT and addA
            tensor addT = pSf[facei]*pFaceVols[facei];
            GradType addA = pFaceVols[facei]*pssf[facei];

            if (faceVols.boundaryField()[patchi].coupled())
            {
                addT *= (1.0 - pWeights[facei]);
                addA *= (1.0 - pWeights[facei]);
            }

            T[pFaceCells[facei]] += addT;
            A[pFaceCells[facei]] += addA;
        }
    }

    T.correctBoundaryConditions();
    A.correctBoundaryConditions();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > treconField
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                "volIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            inv(T) & A,
            zeroGradientFvPatchField<GradType>::typeName
        )
    );

    treconField().correctBoundaryConditions();

    return treconField;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
reconstructNew
(
    const tmp< GeometricField<Type, fvsPatchField, surfaceMesh> >& tssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp< GeometricField<GradType, fvPatchField, volMesh> > tvf
    (
        fvc::reconstruct(tssf())
    );
    tssf.clear();
    return tvf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
