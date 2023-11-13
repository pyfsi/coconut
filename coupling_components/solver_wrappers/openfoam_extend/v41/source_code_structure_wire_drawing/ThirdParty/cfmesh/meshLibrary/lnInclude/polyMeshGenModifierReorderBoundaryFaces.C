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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenModifier.H"
#include "helperFunctions.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::reorderBoundaryFaces()
{
    Info << "Reordering boundary faces " << endl;

    if( Pstream::parRun() )
        reorderProcBoundaryFaces();

    faceListPMG& faces = mesh_.faces_;
    cellListPMG& cells = mesh_.cells_;

    labelLongList newFaceLabel(faces.size());

    label nIntFaces(0), nBndFaces(0);

    label nActiveFaces = faces.size();
    if( Pstream::parRun() )
        nActiveFaces = mesh_.procBoundaries()[0].patchStart();

    const labelLongList& neighbour = mesh_.neighbour();

    labelList nIntFacesAtTask, nBndFacesAtTask;

    faceList newFaces(faces.size());

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1) nowait
        # endif
        forAll(newFaceLabel, faceI)
            newFaceLabel[faceI] = -1;

        //- calculate the number of internal and boundary faces
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 100) \
        reduction(+:nIntFaces,nBndFaces)
        # endif
        for(label faceI=0;faceI<nActiveFaces;++faceI)
        {
            if( neighbour[faceI] >= 0 )
            {
                ++nIntFaces;
            }
            else
            {
                ++nBndFaces;
            }
        }

        //- find internal and boundary faces that are not
        # ifdef USE_OMP
        const label nTasks = 5 * omp_get_num_threads();
        # else
        const label nTasks = 1;
        # endif

        //- count the number of remaining internal faces and the number
        //- of remaining boundary faces
        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            nIntFacesAtTask.setSize(nTasks);
            nBndFacesAtTask.setSize(nTasks);

            for(label taskI=0;taskI<nTasks;++taskI)
            {
                # ifdef USE_OMP
                # pragma omp task default(shared) firstprivate(taskI)
                # endif
                {
                    const label cSize = max(nActiveFaces / nTasks, 1);

                    //- count the number of internal faces
                    label& nif = nIntFacesAtTask[taskI];
                    nif = 0;

                    //- count the number of boundary faces
                    label& nbf = nBndFacesAtTask[taskI];
                    nbf = 0;

                    const label sf = taskI * cSize;
                    label ef = min(sf + cSize, nActiveFaces);
                    if( taskI == (nTasks-1) )
                        ef = nActiveFaces;

                    for(label faceI=sf;faceI<ef;++faceI)
                    {
                        if( neighbour[faceI] >= 0 )
                        {
                            ++nif;
                        }
                        else
                        {
                            ++nbf;
                        }
                    }
                }
            }
        }

        //- set new face labels
        # ifdef USE_OMP
        # pragma omp single
        # endif
        {
            for(label taskI=0;taskI<nTasks;++taskI)
            {
                # ifdef USE_OMP
                # pragma omp task default(shared) firstprivate(taskI)
                # endif
                {
                    const label cSize = max(nActiveFaces / nTasks, 1);
                    const label sf = taskI * cSize;
                    label ef = min(sf + cSize, nActiveFaces);
                    if( taskI == (nTasks-1) )
                        ef = nActiveFaces;

                    label startInt(0), startBnd(nIntFaces);

                    for(label i=0;i<taskI;++i)
                    {
                        startInt += nIntFacesAtTask[i];
                        startBnd += nBndFacesAtTask[i];
                    }

                    for(label faceI=sf;faceI<ef;++faceI)
                    {
                        if( neighbour[faceI] >= 0 )
                        {
                            newFaceLabel[faceI] = startInt++;
                        }
                        else
                        {
                            newFaceLabel[faceI] = startBnd++;
                        }
                    }
                }
            }
        }

        //- set new labels of processor faces
        if( Pstream::parRun() )
        {
            # ifdef USE_OMP
            # pragma omp for schedule(static, 1)
            # endif
            for(label faceI=nActiveFaces;faceI<faces.size();++faceI)
                newFaceLabel[faceI] = faceI;
        }

        //- copy faces to a new list
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(faces, faceI)
            newFaces[newFaceLabel[faceI]].transfer(faces[faceI]);

        //- copy faces back to the original list
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1) nowait
        # endif
        forAll(faces, faceI)
            faces[faceI].transfer(newFaces[faceI]);

        //- renumber cells
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 40)
        # endif
        forAll(cells, cellI)
        {
            cell& c = cells[cellI];

            forAll(c, fI)
                c[fI] = newFaceLabel[c[fI]];
        }
    }

    //- re-create boundary data
    PtrList<boundaryPatch>& boundaries = mesh_.boundaries_;
    if( (boundaries.size() > 1) && (newFaceLabel.size() != 0) )
    {
        if( newFaceLabel[boundaries[0].patchStart()] > nIntFaces )
        {
            //- create an additional patch for internal faces
            PtrList<boundaryPatch> newPatches(boundaries.size()+1);

            HashSet<word> extraPatches;
            if( mesh_.metaData().found("extraPatches") )
                extraPatches =
                    HashSet<word>(mesh_.metaData().lookup("extraPatches"));

            //- check for the existence of the default name
            label counter(0);
            word pName;
            bool found;
            do
            {
                found = false;

                pName = "defaultFaces";
                if( counter )
                    pName += "_"+help::labelToText(counter);

                ++counter;

                forAll(boundaries, patchI)
                    if( boundaries[patchI].patchName() == pName )
                        found = true;

                if( !found )
                    extraPatches.insert(pName);

            } while( found );

            //- write extra patches
            if( extraPatches.size() )
            {
                mesh_.metaData().add("extraPatches", extraPatches, true);
            }

            //- create a new patch
            newPatches.set
            (
                0,
                new boundaryPatch
                (
                    pName,
                    "patch",
                    boundaries[0].patchStart()-nIntFaces,
                    nIntFaces
                )
            );

            //- update existing patches
            forAll(boundaries, patchI)
            {
                boundaryPatch& patch = boundaries[patchI];
                const label start = patch.patchStart();
                const label size = patch.patchSize();

                label nFacesInPatch(0);
                for(label fI=0;fI<size;++fI)
                {
                    if( newFaceLabel[start+fI] >= nIntFaces )
                        ++nFacesInPatch;
                }

                newPatches.set
                (
                    patchI+1,
                    new boundaryPatch
                    (
                        boundaries[patchI].patchName(),
                        boundaries[patchI].patchType(),
                        nFacesInPatch,
                        newFaceLabel[start]
                    )
                );
            }

            //- copy patches
            boundaries.setSize(newPatches.size());

            forAll(boundaries, patchI)
                boundaries.set(patchI, new boundaryPatch(newPatches[patchI]));
        }
        else
        {
            //- additional patch is not needed
            //- change the starting label of the existing patches
            label startBndFace(nIntFaces);
            forAll(boundaries, patchI)
            {
                boundaryPatch& patch = boundaries[patchI];
                const label start = patch.patchStart();
                const label size = patch.patchSize();

                label nFacesInPatch(0);
                for(label fI=0;fI<size;++fI)
                {
                    if( newFaceLabel[start+fI] >= nIntFaces )
                        ++nFacesInPatch;
                }

                patch.patchStart() = startBndFace;
                patch.patchSize() = nFacesInPatch;

                startBndFace += patch.patchSize();
            }
        }
    }
    else if( boundaries.size() == 0 )
    {
        boundaries.setSize(1);

        boundaries.set
        (
            0,
            new boundaryPatch
            (
                "defaultFaces",
                "patch",
                nBndFaces,
                nIntFaces
            )
        );
    }
    else
    {
        boundaries[0].patchStart() = nIntFaces;
        boundaries[0].patchSize() = nBndFaces;
    }

    //- remove boundary faces
    newFaces.clear();

    //- update face subsets
    mesh_.updateFaceSubsets(newFaceLabel);

    //- delete invalid data
    mesh_.clearOut();
    this->clearOut();

    Info << "Finished reordering boundary faces" << endl;
}

void polyMeshGenModifier::reorderProcBoundaryFaces()
{
    PtrList<processorBoundaryPatch>& procBoundaries = mesh_.procBoundaries_;
    if( procBoundaries.size() == 0 )
    {
        Warning << "Processor " << Pstream::myProcNo() << " has no "
            << "processor boundaries!" << endl;
        return;
    }

    //- check if there exist any internal or ordinary bnd faces
    //- which appear after processor bnd faces. Move those faces before
    //- the processor boundary
    const label origProcStart = procBoundaries[0].patchStart();
    label nProcFaces(0);
    forAll(procBoundaries, patchI)
        nProcFaces += procBoundaries[patchI].patchSize();

    faceListPMG& faces = mesh_.faces_;
    cellListPMG& cells = mesh_.cells_;

    const label shift = faces.size() - (origProcStart + nProcFaces);
    if( shift == 0 )
        return;
    if( shift < 0 )
        FatalErrorIn
        (
            "void polyMeshGenModifier::reorderProcBoundaryFaces()"
        ) << "Missing some faces!" << abort(FatalError);

    labelLongList newFaceLabel(faces.size(), -1);

    //- faces added after processor boundaries should be moved up front
    faceList facesAtEnd(shift);
    label counter(0);
    for(label faceI=(origProcStart + nProcFaces);faceI<faces.size();++faceI)
    {
        facesAtEnd[counter].transfer(faces[faceI]);
        newFaceLabel[faceI] = origProcStart + counter;
        ++counter;
    }

    //- shift proc faces
    forAllReverse(procBoundaries, patchI)
    {
        const label start = procBoundaries[patchI].patchStart();
        const label end = start + procBoundaries[patchI].patchSize();

        //- set patch start to the new value
        procBoundaries[patchI].patchStart() += shift;

        for(label faceI=end-1;faceI>=start;--faceI)
        {
            faces[faceI+shift].transfer(faces[faceI]);
            newFaceLabel[faceI] = faceI + shift;
        }
    }

    //- store faces taken from the end
    forAll(facesAtEnd, fI)
    {
        faces[origProcStart+fI].transfer(facesAtEnd[fI]);
    }

    //- set correct patch size
    PtrList<boundaryPatch>& boundaries = mesh_.boundaries_;
    if( boundaries.size() == 1 )
    {
        boundaries[0].patchSize() =
            procBoundaries[0].patchStart() - boundaries[0].patchStart();
    }
    else
    {
        const label start = boundaries[0].patchStart();

        boundaries.clear();
        boundaries.setSize(1);
        boundaries.set
        (
            0,
            new boundaryPatch
            (
                "defaultFaces",
                "patch",
                procBoundaries[0].patchStart() - start,
                start
            )
        );
    }

    //- renumber cells
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(cells, cellI)
    {
        cell& c = cells[cellI];

        forAll(c, fI)
            if( newFaceLabel[c[fI]] != -1 )
                c[fI] = newFaceLabel[c[fI]];
    }

    //- update face subsets
    mesh_.updateFaceSubsets(newFaceLabel);

    //- delete invalid data
    mesh_.clearOut();
    this->clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
