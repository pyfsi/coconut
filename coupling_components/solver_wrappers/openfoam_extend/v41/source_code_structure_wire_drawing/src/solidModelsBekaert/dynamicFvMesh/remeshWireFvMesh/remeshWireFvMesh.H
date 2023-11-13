/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

Class
    Foam::remeshWireFvMesh

Description
    Work in progress.

    Re-mesh the wire and map the fields.

Author
    Philip Cardiff, UCD
    Peter De Jaeger, Bekaert

SourceFiles
    remeshWireFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef remeshWireFvMesh_H
#define remeshWireFvMesh_H

#include "layerAdditionFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Forward declaration
    class IOobjectList;
    class meshToMesh;
    class geometricFieldContainer;

/*---------------------------------------------------------------------------* \
                      Class remeshWireFvMesh Declaration
\*---------------------------------------------------------------------------*/

class remeshWireFvMesh
:
    public layerAdditionFvMesh
{
    // Private Data

        //- remeshWire switch
        const Switch remeshWire_;

        //- remeshWire frequency
        const int remeshFrequency_;

        //- Keep intermediate meshing and mapping fields
        const Switch keepIntermediateMeshesAndFields_;

        //- Store lists of geometric fields for global mesh during remeshing in
        //  parallel
        autoPtr<geometricFieldContainer> geomFieldsPtr_;


    // Nested helper classes

        //- Store mesh information: used to store info of old mesh
        class meshInfo
        {
            public:

                // Public Data
                const label nPoints_;
                const label nFaces_;
                const label nCells_;
                labelList patchStarts_;
                labelList patchNMeshPoints_;

                // Constructor
                meshInfo(const polyMesh& mesh)
                :
                nPoints_(mesh.nPoints()),
                nFaces_(mesh.nFaces()),
                nCells_(mesh.nCells()),
                patchStarts_(mesh.boundaryMesh().size(), -1),
                patchNMeshPoints_(patchStarts_.size(), -1)
                {
                    forAll(patchStarts_, patchI)
                    {
                        patchStarts_[patchI] =
                            mesh.boundaryMesh()[patchI].start();
                        patchNMeshPoints_[patchI] =
                            mesh.boundaryMesh()[patchI].nPoints();
                    }
                }
        };

    // Private Member Functions

        //- Return const-access to the remesh dict
        const dictionary& remeshDict() const
        {
            return dict().subDict(type() + "Coeffs");
        }

        //- Return non-const access to the remesh dict
        dictionary& remeshDict()
        {
            return dict().subDict(type() + "Coeffs");
        }

        //- Should the wire be re-meshed this time-step?
        bool checkRemeshWire() const;

        //- Write all fields of the given mesh to disk
        void writeMeshFields(const fvMesh& mesh) const;

        //- Map fields from old to new mesh
        //  This function calls mapConsistentMesh
        void mapFieldsAndOverwriteTargetMesh
        (
            const polyMesh& meshSource,
            fvMesh& meshTarget
        );

        //- Map fields from old to new mesh (taken from mapFields.C)
        void mapConsistentMesh
        (
            const fvMesh& meshSource,
            const fvMesh& meshTarget
        );

        //- Correct the material index fields
        void correctMaterialIndexFields(const fvMesh& mesh) const;

        //- Update boundary conditions
        void updateBoundaryConditionsAndInterp
        (
            const fvMesh& mesh,
            const remeshWireFvMesh::meshInfo& oldMeshInfo,
            const bool onlyMasterWritesFields = true
        ) const;

        //- Copy the zones from meshSource to meshTarget
        void copyZones(const fvMesh& meshSource, fvMesh& meshTarget) const;

        //- Decompopse the global mesh and place into the local proc meshes
        void decomposeGlobalMesh
        (
            const fvMesh& globalMesh,
            fvMesh& localMesh,
            const meshInfo& oldLocalMeshInfo
        );

        //- Replace the meshTarget fields with the meshSource fields
        virtual void replaceMeshFields
        (
            const fvMesh& meshTarget,
            const fvMesh& meshSource,
            const bool writeFields
        );

        //- Map vol fields between meshes, where fields are read from disk
        //  This is directly based on MapConsistentVolFields in the mapFields
        //  utility
        template<class Type>
        void newMapConsistentVolFields
        (
            const IOobjectList& objects,
            const meshToMesh& meshToMeshInterp
        );

        //- Map vol fields between meshes, where source fields are read from
        //  memory
        template<class Type>
        void newMapConsistentVolFields
        (
            const meshToMesh& meshToMeshInterp
        );

        //- Replace oldMesh fields with newMesh fields
        template
        <
            class Type,
            template<class> class PatchFieldType,
            class MeshType,
            class FieldMeshType
        >
        void ReplaceOldMeshFieldsWithNewMeshFields
        (
            const fvMesh& oldMesh,
            const fvMesh& newMesh,
            const FieldMeshType& oldFieldMesh,
            const FieldMeshType& newFieldMesh,
            const bool writeFields
        );

        //- Replace oldMesh point fields with newMesh point fields
        template<class Type>
        void ReplaceOldMeshPointFieldsWithNewMeshPointFields
        (
            const fvMesh& oldMesh,
            const fvMesh& newMesh,
            const pointMesh& oldPMesh,
            const pointMesh& newPMesh,
            const bool writeFields
        );

        //- Write all objects in the mesh to disk
        template
        <
            class Type, template<class> class PatchFieldType, class MeshType
        >
        void WriteAllObjects(const fvMesh& mesh) const;

        //- Disallow default bitwise copy construct
        remeshWireFvMesh(const remeshWireFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const remeshWireFvMesh&);


public:

    //- Runtime type information
    TypeName("remeshWireFvMesh");


    // Constructors

        //- Construct from database
        explicit remeshWireFvMesh(const IOobject& io);

    // Destructor

        virtual ~remeshWireFvMesh();

    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "remeshWireFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
