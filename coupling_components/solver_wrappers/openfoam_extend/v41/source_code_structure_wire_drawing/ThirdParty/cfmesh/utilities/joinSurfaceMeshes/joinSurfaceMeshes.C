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

Description
    cfMesh utility to merge the supplied list of patches onto a single
    patch.

Author
    Ivor Clifford <ivor.clifford@psi.ch>

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "autoPtr.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "demandDrivenData.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "triSurfaceCleanupDuplicates.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mergeSurface(triSurf& mergedSurf, const triSurf& surf)
{
    const label nOrigPoints = mergedSurf.nPoints();
    const label nOrigTriangles = mergedSurf.size();
    const label nOrigFeatureEdges = mergedSurf.nFeatureEdges();

    triSurfModifier sMod(mergedSurf);
    pointField& points = sMod.pointsAccess();
    LongList<labelledTri>& trias = sMod.facetsAccess();
    edgeLongList& featureEdges = sMod.featureEdgesAccess();

    //- copy surface points
    points.setSize(nOrigPoints+surf.nPoints());

    const pointField& sPoints = surf.points();
    forAll(sPoints, pI)
        points[nOrigPoints+pI] = sPoints[pI];

    //- check patch addressing
    std::map<word, label> patchLabel;
    std::map<word, word> patchType;
    forAll(mergedSurf.patches(), patchI)
    {
        const word pName = mergedSurf.patches()[patchI].name();
        patchLabel[pName] = patchI;
        patchType[pName] = mergedSurf.patches()[patchI].geometricType();
    }

    //- copy triangles
    std::set<word> newPatchNames;
    trias.setSize(nOrigTriangles+surf.size());
    forAll(surf, tI)
    {
        labelledTri& newTri = trias[nOrigTriangles+tI];

        const labelledTri& tri = surf[tI];
        newTri[0] = nOrigPoints + tri[0];
        newTri[1] = nOrigPoints + tri[1];
        newTri[2] = nOrigPoints + tri[2];

        const word pName = surf.patches()[tri.region()].name();
        if( patchLabel.find(pName) == patchLabel.end() )
        {
            const label patchID = patchLabel.size();
            patchLabel[pName] = patchID;
            patchType[pName] = surf.patches()[tri.region()].geometricType();
            newPatchNames.insert(pName);
        }

        newTri.region() = patchLabel[pName];
    }

    //- insert new patches
    if( newPatchNames.size() )
    {
        geometricSurfacePatchList& patches = sMod.patchesAccess();

        const label nOrigPatches = patches.size();
        patches.setSize(nOrigPatches+newPatchNames.size());

        forAllConstIter(std::set<word>, newPatchNames, it)
        {
            const label patchI = patchLabel[*it];
            const word pType = patchType[*it];

            patches[patchI].name() = *it;
            patches[patchI].geometricType() = pType;
        }
    }

    //- copy feature edges
    const edgeLongList& sFeatureEdges = surf.featureEdges();
    featureEdges.setSize(nOrigFeatureEdges+sFeatureEdges.size());
    forAll(featureEdges, feI)
    {
        edge& mfe = featureEdges[nOrigFeatureEdges+feI];
        const edge& fe = sFeatureEdges[feI];

        mfe.start() = nOrigPoints + fe.start();
        mfe.end() = nOrigPoints + fe.end();
    }

    //- transfer point subsets
    DynList<label> subsetIds;
    surf.pointSubsetIndices(subsetIds);
    forAll(subsetIds, i)
    {
        const label sId = subsetIds[i];
        const word sName = surf.pointSubsetName(sId);

        label newId = mergedSurf.pointSubsetIndex(sName);
        if( newId < 0 )
            newId = mergedSurf.addPointSubset(sName);

        labelLongList indices;
        surf.pointsInSubset(sId, indices);

        forAll(indices, elI)
            mergedSurf.addPointToSubset(newId, nOrigPoints+indices[elI]);
    }

    //- transfer facet subsets
    subsetIds.clear();
    surf.facetSubsetIndices(subsetIds);
    forAll(subsetIds, i)
    {
        const label sId = subsetIds[i];
        const word sName = surf.facetSubsetName(sId);

        label newId = mergedSurf.facetSubsetIndex(sName);
        if( newId < 0 )
            newId = mergedSurf.addFacetSubset(sName);

        labelLongList indices;
        surf.facetsInSubset(sId, indices);

        forAll(indices, elI)
            mergedSurf.addFacetToSubset(newId, nOrigTriangles+indices[elI]);
    }

    //- transfer feature edge subsets
    subsetIds.clear();
    surf.edgeSubsetIndices(subsetIds);
    forAll(subsetIds, i)
    {
        const label sId = subsetIds[i];
        const word sName = surf.edgeSubsetName(sId);

        label newId = mergedSurf.edgeSubsetIndex(sName);
        if( newId < 0 )
            newId = mergedSurf.addEdgeSubset(sName);

        labelLongList indices;
        surf.edgesInSubset(sId, indices);

        forAll(indices, elI)
            mergedSurf.addEdgeToSubset(newId, nOrigFeatureEdges+indices[elI]);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("output surface file");
    argList::validOptions.insert("surfaceNames", "list of surface files");

    argList args(argc, argv);

    fileName outFileName = args.args()[1];

    //- check if the surface has the extension
    if( outFileName == outFileName.lessExt() )
        outFileName += ".fms";

    DynList<fileName> surfacesToMerge;
    if( args.options().found("surfaceNames") )
    {
        std::string fileNames = args.options()["surfaceNames"];

        //- remove the ( and ) characters
        if( fileNames[0] == token::BEGIN_LIST )
            fileNames.replace(0, 1, "");
        if( fileNames[fileNames.size()-1] == token::END_LIST )
            fileNames.replace(fileNames.size()-1, 1, "");

        //- read file names and append them to the list
        fileName fName;
        fName.reserve(100);
        for(unsigned charI=0;charI<fileNames.size();++charI)
        {
            if( fileNames[charI] == token::SPACE )
            {
                //- add the file name to the list and clear it
                if( fName.size() )
                    surfacesToMerge.append(fName);

                fName.clear();
            }
            else
            {
                //- add the current character to the current file name
                fName += fileNames[charI];
            }
        }

        //- append the last filename to the list
        if( fName.size() )
            surfacesToMerge.append(fName);
    }

    //- create the resulting surface
    triSurf mergedSurfaces;

    // Read original surfaces
    forAll(surfacesToMerge, surfI)
    {
        const fileName surfName = surfacesToMerge[surfI];

        Info << "Adding surface " << surfName << endl;

        const triSurf surf(surfName);

        mergeSurface(mergedSurfaces, surf);
    }

    //- merge indentical entities
    Info << "Creating octree" << endl;
    meshOctree octree(mergedSurfaces);

    Info << "Refining octree" << endl;
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(20, 20);

    Info << "Merging identities" << endl;
    triSurfaceCleanupDuplicates(octree).mergeIdentities();

    // Write new surface mesh
    Info << "Writting surface" << endl;
    mergedSurfaces.writeSurface(outFileName);

    Info << "Surface written to " << outFileName <<  endl;

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
