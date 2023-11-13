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

\*---------------------------------------------------------------------------*/

#include "cartesian2DFvMesh.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "mySolidContactFvPatchVectorField.H"
#include "meshTriangulation.H"
#include "cartesian2DMeshGenerator.H"
#include "triSurf.H"
#include "triSurfaceDetectFeatureEdges.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cartesian2DFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh, cartesian2DFvMesh, IOobject
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::cartesian2DFvMesh::cartesian2DFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    maxCellSize_(readScalar(dict_.lookup("maxCellSize"))),
    bb_
    (
        dict_.lookupOrDefault<boundBox>
        (
            "boundBox",
            boundBox(vector::min, vector::max)
        )
    ),
    gfz_(*this, bb_),
    indicator_
    (
        IOobject
        (
            "globalFaceZones",
            time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    updateIndicator_
    (
        dict_.lookupOrDefault<Switch>("updateIndicatorField", true)
    )
{
    Info<< "Creating cartesian2DFvMesh" << nl
        << "    boundBox: " << bb_ << nl
        << "    updateIndicatorField: " << updateIndicator_ << endl;

    // Create directory for triSurfaces
    Info<< "    Creating triSurfaces to place the STLs in" << endl;
    rmDir(time().path()/"triSurfaces");
    mkDir(time().path()/"triSurfaces");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cartesian2DFvMesh::~cartesian2DFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::cartesian2DFvMesh::update()
{
    // 2-D cartesian mesher cannot be run in parallel
    if (Pstream::parRun())
    {
        FatalErrorIn("bool Foam::cartesian2DFvMesh::update()")
            << "cartesian2DMeshGenerator cannot be run in parallel"
            << abort(FatalError);
    }

    if (time().timeIndex() != 0)
    {
        // Create a STL surface from the boundary of the current mesh
        // meshTriangulation triangulates all the faces including internal faces

        Info<< "Triangulating entire mesh" << endl;

        meshTriangulation triangulatedFaces
        (
            *this,                             // polyMesh
            -1,                                // internal faces patch
            boolList(this->nCells(), true),    // included cell
            true                               // face centre decomposition
        );

        // Write all triangles to STL (including internal faces)
        // const word surfaceName =
        //     "meshTotal_" + Foam::name(time().timeIndex()) + ".stl";
        // Info<< "Writing surface mesh to " << surfaceName << endl;
        // triangulatedFaces.write(time().path()/"triSurfaces"/surfaceName);


        // Create triangulated surface of just the boundary

        // Mark all boundary faces
        const label nFaces = triangulatedFaces.size();
        const label nIntFaces = triangulatedFaces.nInternalFaces();

        boolList includedFaces(nFaces, false);
        labelList boundaryPointMap(triangulatedFaces.nPoints(), -1);
        labelList boundaryFaceMap(nFaces, -1);

        const labelList& faceMap = triangulatedFaces.faceMap();

        for (int triFaceI = nIntFaces; triFaceI < nFaces; triFaceI++)
        {
            const label faceID = faceMap[triFaceI];

            // Include the face if it is not on an empty patch
            const label patchID = boundaryMesh().whichPatch(faceID);

            if (boundaryMesh()[patchID].type() != emptyPolyPatch::typeName)
            {
                includedFaces[triFaceI] = true;
            }
        }

        triSurface triangulatedBoundaryFaces =
            triangulatedFaces.subsetMesh
            (
                includedFaces,
                boundaryPointMap,
                boundaryFaceMap
            );

        const word boundarySurfaceName =
            time().path()/"triSurfaces"/"meshBoundary_"
          + Foam::name(time().timeIndex()) + ".stl";

        Info<< "Writing boundary surface mesh to " << boundarySurfaceName
            << endl;

        triangulatedBoundaryFaces.write(boundarySurfaceName);


        // Detect edges on the STL file and convert to FMS file
        triSurf stlSurface(boundarySurfaceName);

        const scalar tol = 45.0;
        triSurfaceDetectFeatureEdges edgeDetector(stlSurface, tol);
        edgeDetector.detectFeatureEdges();

        const word boundarySurfaceNameFms =
            time().path()/"triSurfaces"/"meshBoundary_"
          + Foam::name(time().timeIndex()) + ".fms";

        Info << "Writing : " << boundarySurfaceNameFms << endl;

        stlSurface.writeSurface(boundarySurfaceNameFms);


        // Write out mesh dict for cartesianMesh to use the new
        // surface mesh

        WarningIn("bool Foam::cartesian2DFvMesh::update()")
            << "The meshDict will be overwritten" << endl;

        IOdictionary meshDict
        (
            IOobject
            (
                "meshDict",
                time().system(),
                time(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        );

        // Add settings to mesh dict
        meshDict.add("maxCellSize", maxCellSize_);
        meshDict.add
        (
            "surfaceFile", "triSurfaces/meshBoundary_"
          + Foam::name(time().timeIndex()) + ".fms"
        );

        Info<< "Writing system/meshDict" << endl;
        meshDict.regIOobject::write();


        // Create 2-D cartesian mesh
        // cfMesh writes the mesh to constant so we will make a copy of the
        // orginal mesh and move it back after moving the new mesh

        // Make a copy of the original mesh
        rmDir(time().path()/time().constant()/"polyMesh.org");
        cp
        (
            time().path()/time().constant()/"polyMesh",
            time().path()/time().constant()/"polyMesh.org"
        );

        Info<< "Running cartesian2DMesh" << endl;
        cartesian2DMeshGenerator cmg(this->time());

        // Write mesh
        Info<< "Writing new cartesian mesh" << endl;
        cmg.writeMesh();

        // Create time directory
        mkDir(time().path()/time().timeName());

        // Move polyMesh directory from constant to the time directory
        mv
        (
            time().path()/time().constant()/"polyMesh",
            time().path()/time().timeName()/"polyMesh"
        );

        // Move original mesh back
        mv
        (
            time().path()/time().constant()/"polyMesh.org",
            time().path()/time().constant()/"polyMesh"
        );


        Info<< nl << "Number of cells in old mesh: " << nCells() << endl;

        // Re-read the mesh from disc
        readUpdate();
        Info<< nl << "Number of cells in new mesh: " << nCells() << endl;



        // Stop here while testing
        FatalErrorIn("bool Foam::cartesian2DFvMesh::update()")
            << "testing: stop" << abort(FatalError);

        return false;

        // // We will remove all current face zones
        // // This is to avoid problems with faces being present in multiple
        // // face zones
        // gfz_.removeFaceZones();

        // // Lookup DU field
        // const volVectorField& DU = this->lookupObject<volVectorField>("DU");

        // forAll(DU.boundaryField(), patchI)
        // {
        //     const word pType = DU.boundaryField()[patchI].type();

        //     if (pType == mySolidContactFvPatchVectorField::typeName)
        //     {
        //         gfz_.createFaceZone(patchI);

        //         mySolidContactFvPatchVectorField& contactPatch =
        //             const_cast<mySolidContactFvPatchVectorField&>
        //             (
        //                 refCast<const mySolidContactFvPatchVectorField>
        //                 (
        //                     DU.boundaryField()[patchI]
        //                 )
        //             );

        //         // Force contact patch addressing to be re-calculated
        //         contactPatch.clearOut();
        //     }
        // }

        // // Reorder each face zone so that the faces are in the same order on all
        // // procs; this allows parallel syncing of zone data
        // forAll(this->faceZones(), zoneI)
        // {
        //     gfz_.reorderFaceZone(zoneI);
        // }

        // // Set indicator field for visualisation of current zones
        // if (updateIndicator_)
        // {
        //     // Reset field to zero
        //     indicator_ = dimensionedScalar("zero", dimless, 0.0);

        //     const label nActiveFaces = this->nFaces();
        //     forAll(this->faceZones(), zoneI)
        //     {
        //         const faceZone& curZone = this->faceZones()[zoneI];

        //         forAll(curZone, fI)
        //         {
        //             const label faceID = curZone[fI];

        //             if (faceID < nActiveFaces)
        //             {
        //                 const label patchID =
        //                     this->boundaryMesh().whichPatch(faceID);

        //                 if (patchID != -1)
        //                 {
        //                     const label start =
        //                         this->boundaryMesh()[patchID].start();
        //                     indicator_.boundaryField()[patchID][faceID - start] =
        //                         1.0;
        //                 }
        //             }
        //         }
        //     }

        //     // Set internal cell values
        //     scalarField& indicatorI = indicator_.internalField();
        //     forAll(indicator_.boundaryField(), patchI)
        //     {
        //         const word pType = indicator_.boundaryField()[patchI].type();

        //         if (!this->boundaryMesh()[patchI].coupled() && pType != "empty")
        //         {
        //             scalarField& pIndicator = indicator_.boundaryField()[patchI];
        //             const unallocLabelList& faceCells =
        //                 this->boundaryMesh()[patchI].faceCells();

        //             forAll(faceCells, faceI)
        //             {
        //                 indicatorI[faceCells[faceI]] =
        //                     max
        //                     (
        //                         indicatorI[faceCells[faceI]],
        //                         pIndicator[faceI]
        //                     );
        //             }
        //         }
        //     }
        // }
    }

    return true;
}


// ************************************************************************* //
