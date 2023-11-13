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

Global
    axisymmetricSolid

Description
    By default OpenFOAM will solve in all 3 directions for axisymmetric cases,
    to allow for swirl; however, axisymmetric solid models are truly 2-D so this
    file will disable the wedge direction making the model 2-D.

    Note: const_cast is used to change polyMesh.solutionD()

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/


label nWedgePatches = 0;
vector wedgeDirVec = vector::zero;

forAll(mesh.boundaryMesh(), patchI)
{
    if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
    {
        const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>
            (
                mesh.boundaryMesh()[patchI]
            );

        nWedgePatches++;
        wedgeDirVec += cmptMag(wpp.centreNormal());

        // Make sure that solidWedge is used NOT wedge
        if (DU.boundaryField()[patchI].type() != "solidWedge")
        {
            FatalError
                << "solidWedge should be used on DU boundary wedge patches"
                << " as non-orthogonal is important!"
                << abort(FatalError);
        }
    }
}

reduce(nWedgePatches, maxOp<label>());

if (nWedgePatches)
{
    if (nWedgePatches != 2)
    {
        FatalError
            << "For axisymmetric cases, there should be exactly two wedges "
            << "patches!" << abort(FatalError);
    }

    // PC, 22-Nov-18
    // Disabling the 3rd direction slows convergence so we won't do it
    // Info<< nl << "Axisymmetric case: disabling the 3rd direction" << nl
    //     << endl;

    // // We will const_cast as it is currently the tidiest way, until polyMesh
    // // is modified or gives write access to solutionD
    // Vector<label>& solD = const_cast<Vector<label>&>(mesh.solutionD());

    // reduce(wedgeDirVec, sumOp<vector>());

    // wedgeDirVec /= mag(wedgeDirVec);

    // for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    // {
    //     if (wedgeDirVec[cmpt] > 1e-6)
    //     {
    //         solD[cmpt] = -1;
    //     }
    //     else
    //     {
    //         solD[cmpt] = 1;
    //     }
    // }
}


// Check all the face normals are in the same direction
forAll(mesh.boundaryMesh(), patchI)
{
    if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
    {
        // Unit face normals on processor
        const vectorField nf = mesh.boundaryMesh()[patchI].faceNormals();

        if (nf.size() == 0)
        {
            FatalError
                << "There are no faces on the wedge patch "
                << mesh.boundaryMesh()[patchI].name() << " on this processor:"
                << nl << "Every processor should have at least one face on "
                << "each wedge patch"
                << abort(FatalError);
        }

        // Check that all the wedge face normal point in the same direction

        vector firstFaceNOnMasterProc = vector::zero;

        if (Pstream::master())
        {
            firstFaceNOnMasterProc = nf[0];
        }

        // Sync in parallel so that all processors have the master vector
        reduce(firstFaceNOnMasterProc, sumOp<vector>());

        forAll(nf, faceI)
        {
            if ((nf[faceI] & firstFaceNOnMasterProc) < 0)
            {
                FatalError
                    << "On wedge patch "
                    << mesh.boundaryMesh()[patchI].name() << " there are at "
                    << "least two faces with unit normals in the opposite "
                    << "directions" << nl
                    << "Please check that the wedge patches are correctly "
                    << "defined"
                    << abort(FatalError);
            }
        }
    }
}
