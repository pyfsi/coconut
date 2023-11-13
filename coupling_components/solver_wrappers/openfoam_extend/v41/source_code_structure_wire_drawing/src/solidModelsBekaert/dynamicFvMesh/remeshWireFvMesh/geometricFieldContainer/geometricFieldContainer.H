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
    Foam::geometricFieldContainer

Description
    This container will reads all geometric fields (volField, pointField,
    surfaceField) from disk for the current mesh at the current time; it will
    then store these fields in pointer lists and register them to the mesh
    objectRegistry.

    This class is useful for constructing a global mesh and fields from disk
    during parallel running.

Author
    Philip Cardiff, UCD

SourceFiles
    geometricFieldContainer.C

\*---------------------------------------------------------------------------*/

#ifndef geometricFieldContainer_H
#define geometricFieldContainer_H

#include "fvMesh.H"
#include "GeometricFields.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------* \
                      Class geometricFieldContainer Declaration
\*---------------------------------------------------------------------------*/

class geometricFieldContainer
{
    // Private Data

        //- Store a reference to the mesh
        const fvMesh& mesh_;

        //- Store a reference to the pointMesh
        const pointMesh& pMesh_;

        //- Lists of volFields
        PtrList<volScalarField> volScalarFields_;
        PtrList<volVectorField> volVectorFields_;
        PtrList<volTensorField> volTensorFields_;
        PtrList<volSymmTensorField> volSymmTensorFields_;
        PtrList<volDiagTensorField> volDiagTensorFields_;
        PtrList<volSphericalTensorField> volSphericalTensorFields_;

        //- Lists of surfaceFields
        PtrList<surfaceScalarField> surfaceScalarFields_;
        PtrList<surfaceVectorField> surfaceVectorFields_;
        PtrList<surfaceTensorField> surfaceTensorFields_;
        PtrList<surfaceSymmTensorField> surfaceSymmTensorFields_;
        PtrList<surfaceDiagTensorField> surfaceDiagTensorFields_;
        PtrList<surfaceSphericalTensorField> surfaceSphericalTensorFields_;

        //- Lists of pointFields
        PtrList<pointScalarField> pointScalarFields_;
        PtrList<pointVectorField> pointVectorFields_;
        PtrList<pointTensorField> pointTensorFields_;
        PtrList<pointSymmTensorField> pointSymmTensorFields_;
        PtrList<pointDiagTensorField> pointDiagTensorFields_;
        PtrList<pointSphericalTensorField> pointSphericalTensorFields_;


    // Private Member Functions

        //- Read fields and set store in list
        template
        <
            class Type,
            template<class> class PatchFieldType,
            class MeshType,
            class FieldMeshType
        >
        void readFields
        (
            PtrList< GeometricField<Type, PatchFieldType, MeshType> >& ptrList,
            const IOobjectList& objects,
            const FieldMeshType& mesh
        );

        //- Disallow default bitwise copy construct
        geometricFieldContainer(const geometricFieldContainer&);

        //- Disallow default bitwise assignment
        void operator=(const geometricFieldContainer&);


public:

    //- Runtime type information
    TypeName("geometricFieldContainer");

    // Constructor

        //- Construct from database
        explicit geometricFieldContainer
        (
            const fvMesh& mesh
        );


    // Destructor

        virtual ~geometricFieldContainer();


    // Member Functions

        // Access functions

            //- Return const reference access to mesh
            virtual const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return const reference access to volScalarFields
            virtual const PtrList<volScalarField>& volScalarFields() const
            {
                return volScalarFields_;
            }

            //- Return const reference access to volVectorFields
            virtual const PtrList<volVectorField>& volVectorFields() const
            {
                return volVectorFields_;
            }

            //- Return const reference access to volTensorFields
            virtual const PtrList<volTensorField>& volTensorFields() const
            {
                return volTensorFields_;
            }

            //- Return const reference access to volSymmTensorFields
            virtual const PtrList<volSymmTensorField>&
            volSymmTensorFields() const
            {
                return volSymmTensorFields_;
            }

            //- Return const reference access to volDiagTensorFields
            virtual const PtrList<volDiagTensorField>&
            volDiagTensorFields() const
            {
                return volDiagTensorFields_;
            }

            //- Return const reference access to volSphericalTensorFields
            virtual const PtrList<volSphericalTensorField>&
            volSphericalTensorFields() const
            {
                return volSphericalTensorFields_;
            }

            //- Return const reference access to surfaceScalarFields
            virtual const PtrList<surfaceScalarField>&
            surfaceScalarFields() const
            {
                return surfaceScalarFields_;
            }

            //- Return const reference access to surfaceVectorFields
            virtual const PtrList<surfaceVectorField>&
            surfaceVectorFields() const
            {
                return surfaceVectorFields_;
            }

            //- Return const reference access to surfaceTensorFields
            virtual const PtrList<surfaceTensorField>&
            surfaceTensorFields() const
            {
                return surfaceTensorFields_;
            }

            //- Return const reference access to surfaceSymmTensorFields
            virtual const PtrList<surfaceSymmTensorField>&
            surfaceSymmTensorFields() const
            {
                return surfaceSymmTensorFields_;
            }

            //- Return const reference access to surfaceDiagTensorFields
            virtual const PtrList<surfaceDiagTensorField>&
            surfaceDiagTensorFields() const
            {
                return surfaceDiagTensorFields_;
            }

            //- Return const reference access to surfaceSphericalTensorFields
            virtual const PtrList<surfaceSphericalTensorField>&
            surfaceSphericalTensorFields() const
            {
                return surfaceSphericalTensorFields_;
            }

            //- Return const reference access to pointScalarFields
            virtual const PtrList<pointScalarField>& pointScalarFields() const
            {
                return pointScalarFields_;
            }

            //- Return const reference access to pointVectorFields
            virtual const PtrList<pointVectorField>& pointVectorFields() const
            {
                return pointVectorFields_;
            }

            //- Return const reference access to pointTensorFields
            virtual const PtrList<pointTensorField>& pointTensorFields() const
            {
                return pointTensorFields_;
            }

            //- Return const reference access to pointSymmTensorFields
            virtual const PtrList<pointSymmTensorField>&
            pointSymmTensorFields() const
            {
                return pointSymmTensorFields_;
            }

            //- Return const reference access to pointDiagTensorFields
            virtual const PtrList<pointDiagTensorField>&
            pointDiagTensorFields() const
            {
                return pointDiagTensorFields_;
            }

            //- Return const reference access to pointSphericalTensorFields
            virtual const PtrList<pointSphericalTensorField>&
            pointSphericalTensorFields() const
            {
                return pointSphericalTensorFields_;
            }

        // Evaluation functions

            //- Read fields
            virtual void readFields();

            //- Clear out all demand driven data
            virtual void clearOut();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "geometricFieldContainerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
