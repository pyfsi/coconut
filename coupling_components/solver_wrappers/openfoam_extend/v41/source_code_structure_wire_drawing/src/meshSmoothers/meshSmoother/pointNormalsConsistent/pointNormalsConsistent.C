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

#include "pointNormalsConsistent.H"
#include "pointMesh.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::pointNormalsConsistent
(
    PtrList<pointField>& patchPointNormals,
    const polyMesh& mesh
)
{
    // Calculate weight for each point on each patch
    // Each point will have a weight for each of its neighbour faces (i.e.
    // pointFaces)

    // Initialise point weights
    List<scalarListList> weights(mesh.boundaryMesh().size());
    forAll(weights, patchI)
    {
        weights[patchI].setSize(mesh.boundaryMesh()[patchI].nPoints());
    }


    // Calculate normals for each patch
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (mesh.boundaryMesh()[patchI].coupled())
        {
            continue;
        }

        // Initialise sum of weights
        // We need a different field for each patch as each patch calculates
        // point normals indepedently
        pointScalarField sumWeights
        (
            IOobject
            (
                "sumWeights_pointNormalsConsistent",
                mesh.polyMesh::instance(),
                mesh
            ),
            pointMesh::New(mesh),
            dimensionedScalar("zero", dimless, 0.0)
        );

        const polyPatch& ppatch = mesh.boundaryMesh()[patchI];
        const labelListList& pointFaces = ppatch.pointFaces();
        const labelList& meshPoints = ppatch.meshPoints();

        // Initialise weights for current patch
        scalarListList& curWeights = weights[patchI];
        forAll(curWeights, pointI)
        {
            curWeights[pointI].setSize(pointFaces[pointI].size());
        }

        scalarField& curSumWeightsI = sumWeights.internalField();

        // Set un-normalised weights to 1.0 i.e. arithmetic average
        forAll(curWeights, pI)
        {
            scalarList& pw = curWeights[pI];
            const labelList& curPointFaces = pointFaces[pI];

            forAll(curPointFaces, pfI)
            {
                // Set the un-normalised point weight
                pw[pfI] = 1.0;

                // Increment the weights sum
                curSumWeightsI[meshPoints[pI]] += 1.0;
            }
        }

        // Sync weights in parallel
        forAll(sumWeights.boundaryField(), pI)
        {
            if (sumWeights.boundaryField()[pI].coupled())
            {
                sumWeights.boundaryField()[pI].initAddField();
            }
        }

        forAll(sumWeights.boundaryField(), pI)
        {
            if (sumWeights.boundaryField()[pI].coupled())
            {
                sumWeights.boundaryField()[pI].addField
                (
                    sumWeights.internalField()
                );
            }
        }

        // Normalise the weights
        forAll(curWeights, pI)
        {
            scalarList& pw = curWeights[pI];

            forAll(pw, pfI)
            {
                pw[pfI] /= curSumWeightsI[meshPoints[pI]];
            }
        }

        // Calculate consistent point normals

        pointVectorField pointNormalField
        (
            IOobject
            (
                "pointNormalField_" + ppatch.name(),
                mesh.polyMesh::instance(),
                mesh
            ),
            pointMesh::New(mesh),
            dimensionedVector("zero", dimless, vector::zero)
        );

        const pointField& fn = ppatch.faceNormals();
        forAll(curWeights, pointI)
        {
            point& curNormal = pointNormalField[meshPoints[pointI]];
            const labelList& curFaces = pointFaces[pointI];

            forAll(curFaces, faceI)
            {
                curNormal += curWeights[pointI][faceI]*fn[curFaces[faceI]];
            }
        }

        // Add contributions from other processors
        forAll(pointNormalField.boundaryField(), pI)
        {
            if (pointNormalField.boundaryField()[pI].coupled())
            {
                pointNormalField.boundaryField()[pI].initAddField();
            }
        }

        forAll(pointNormalField.boundaryField(), pI)
        {
            if (pointNormalField.boundaryField()[pI].coupled())
            {
                pointNormalField.boundaryField()[pI].addField
                (
                    pointNormalField.internalField()
                );
            }
        }

        // Update coupled and constrained boundaries
        pointNormalField.correctBoundaryConditions();

        // Copy normals field into patchPointNormals list
        patchPointNormals[patchI] =
            pointNormalField.boundaryField()[patchI].patchInternalField();
    }
}


// ************************************************************************* //
