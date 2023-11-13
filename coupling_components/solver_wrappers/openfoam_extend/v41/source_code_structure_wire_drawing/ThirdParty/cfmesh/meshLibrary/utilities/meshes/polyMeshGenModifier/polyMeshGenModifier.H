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
    polyMeshGenModifier

Description
    Modifier for polyMeshGen

SourceFiles
    polyMeshGenModifier.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGenModifier_H
#define polyMeshGenModifier_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGen.H"
#include "boolList.H"
#include "VRWGraph.H"
#include "demandDrivenData.H"

namespace Foam
{

// Forward declarations
class VRWGraphList;
class polyMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class polyMeshGenModifier
{
    // Private data
        //- reference to the mesh
        polyMeshGen& mesh_;

        //- helper data for adding cells
        VRWGraph* pointFacesPtr_;

    // Private member functions
        //- calculate and return point faces
        inline void calculatePointFaces()
        {
            const faceListPMG& faces = mesh_.faces();

            pointFacesPtr_ = new VRWGraph();
            VRWGraph& pointFaces = *pointFacesPtr_;

            pointFaces.reverseAddressing(mesh_.points().size(), faces);
        }

        //- re-order positions of processor boundary faces
        //- they should comea immediately after the internal faces
        void reorderProcBoundaryFaces();

protected:

        VRWGraph& pointFaces()
        {
            if( !pointFacesPtr_ )
                calculatePointFaces();

            return *pointFacesPtr_;
        }

public:

    // Constructors
        //- Construct from the reference to the mesh
        polyMeshGenModifier(polyMeshGen& mesh)
        :
            mesh_(mesh),
            pointFacesPtr_(NULL)
        {}

    // Destructor
        ~polyMeshGenModifier()
        {}

    // Member functions
        //- access to mesh points
        inline pointFieldPMG& pointsAccess();

        //- access to mesh faces
        inline faceListPMG& facesAccess();

        //- access to cells
        inline cellListPMG& cellsAccess();

        //- access to processor boundary data
        inline PtrList<processorBoundaryPatch>& procBoundariesAccess();

        //- access to boundary data
        inline PtrList<boundaryPatch>& boundariesAccess();

        //- functions which change the mesh
        //- change position of a point
        inline const point& movePoint(const label pointI, const point& p);

        //- backup all points
        inline void backupPoints();

        //- backup points in a list
        template<class labelListType>
        inline void backupPoints(const labelListType&);

        //- backup a point
        inline void backupPoint(const label pointI);

        //- set the coordinate of a backup
        inline void setBackupPoint(const label pointI, const point& pOrig);

        //- revert all points
        inline void revertPoints();

        //- revert a point
        inline void revertPoint(const label pointI);

        //- remove the point backup
        inline void clearPointBackup(const label pointI);

        //- clear backup of points
        inline void clearPointBackup();

        //- lock all points
        inline void lockPoints();

        //- lock a point
        inline void lockPoint(const label pointI);

        //- unlock all points
        inline void unlockPoints();

        //- unlock a point
        inline void unlockPoint(const label pointI);

    // Complex mesh modifiers
        //- reorder boundary faces
        void reorderBoundaryFaces();

        //- remove unused vertices
        void removeUnusedVertices();

        //- remove faces
        void removeFaces(const boolList& removeFace);

        //- remove duplicate faces from the mesh
        void removeDuplicateFaces();

        //- remove cells
        void removeCells
        (
            const boolList& removeCell,
            const bool removeProcFaces = true
        );

        //- add cells (vertices must be added)
        void addCells(const LongList<faceList>& cellFaces);
        void addCells(const VRWGraphList& cellFaces);
        void addCell(const faceList& cellFaces);

        //- merge another mesh into the current one
        void addMesh(const polyMeshGen&);

        //- merge a polyMesh into the current one
        void addMesh(const polyMesh&);

        //- replace the boundary with new boundary faces
        void replaceBoundary
        (
            const wordList& patchNames,
            const VRWGraph& boundaryFaces,
            const labelLongList& faceOwners,
            const labelLongList& facePatches
        );

        //- add additional faces into processor patches
        void addProcessorFaces
        (
            const VRWGraph& procFaces,
            const labelLongList& facePatches
        );

        //- add new processor patch and return its label
        label addProcessorPatch(const label otherProcLabel);

        //- remove empty processor patch
        bool removeEmptyProcessorPatches();

        //- add buffer cells needed for exporting the mesh in the format
        //- required by some solvers
        void addBufferCells();

        //- zip up topologically open cells
        void zipUpCells();

        //- reorder the cells and faces to reduce the matrix bandwidth
        void renumberMesh();

        //- generate zones from subsets
        void generateCellZonesFromSubsets(const DynList<label>& indices);
        void generateCellZonesFromSubsets();
        void generateFaceZonesFromSubsets(const DynList<label>& indices);
        void generateFaceZonesFromSubsets();
        void generatePointZonesFromSubsets(const DynList<label>& indices);
        void generatePointZonesFromSubsets();

        //- clear out unnecessary data (pointFacesPtr_);
        inline void clearOut()
        {
            deleteDemandDrivenData(pointFacesPtr_);
        }

        //- clear out all allocated data
        inline void clearAll()
        {
            clearOut();
            mesh_.clearOut();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyMeshGenModifierI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
