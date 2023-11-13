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

\*---------------------------------------------------------------------------*/

#include "geometricFieldContainer.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(geometricFieldContainer, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::geometricFieldContainer::geometricFieldContainer
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    pMesh_(pointMesh::New(mesh_)),
    volScalarFields_(),
    volVectorFields_(),
    volTensorFields_(),
    volSymmTensorFields_(),
    volDiagTensorFields_(),
    volSphericalTensorFields_(),
    surfaceScalarFields_(),
    surfaceVectorFields_(),
    surfaceTensorFields_(),
    surfaceSymmTensorFields_(),
    surfaceDiagTensorFields_(),
    surfaceSphericalTensorFields_(),
    pointScalarFields_(),
    pointVectorFields_(),
    pointTensorFields_(),
    pointSymmTensorFields_(),
    pointDiagTensorFields_(),
    pointSphericalTensorFields_()
{
    readFields();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::geometricFieldContainer::~geometricFieldContainer()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::geometricFieldContainer::readFields()
{
    if (volScalarFields_.size())
    {
        FatalErrorIn(type() + "::readFields()")
            << "Fields have already been read!" << abort(FatalError);
    }

    Info<< "Reading fields for mesh = " << mesh_.name()
        << " at time = " << mesh_.time().timeName() << endl;

    // Read the names of objects written to the current time step
    const IOobjectList objects(mesh_, mesh_.time().timeName());

    // Read fields and store them in lists

    readFields(volScalarFields_, objects, mesh_);
    readFields(volVectorFields_, objects, mesh_);
    readFields(volTensorFields_, objects, mesh_);
    readFields(volSymmTensorFields_, objects, mesh_);
    readFields(volDiagTensorFields_, objects, mesh_);
    readFields(volSphericalTensorFields_, objects, mesh_);

    readFields(surfaceScalarFields_, objects, mesh_);
    readFields(surfaceVectorFields_, objects, mesh_);
    readFields(surfaceTensorFields_, objects, mesh_);
    readFields(surfaceSymmTensorFields_, objects, mesh_);
    readFields(surfaceDiagTensorFields_, objects, mesh_);
    readFields(surfaceSphericalTensorFields_, objects, mesh_);

    readFields(pointScalarFields_, objects, pMesh_);
    readFields(pointVectorFields_, objects, pMesh_);
    readFields(pointTensorFields_, objects, pMesh_);
    readFields(pointSymmTensorFields_, objects, pMesh_);
    readFields(pointDiagTensorFields_, objects, pMesh_);
    readFields(pointSphericalTensorFields_, objects, pMesh_);
}

void Foam::geometricFieldContainer::clearOut()
{
    volScalarFields_.clear();
    volVectorFields_.clear();
    volTensorFields_.clear();
    volSymmTensorFields_.clear();
    volDiagTensorFields_.clear();
    volSphericalTensorFields_.clear();

    surfaceScalarFields_.clear();
    surfaceVectorFields_.clear();
    surfaceTensorFields_.clear();
    surfaceSymmTensorFields_.clear();
    surfaceDiagTensorFields_.clear();
    surfaceSphericalTensorFields_.clear();

    pointScalarFields_.clear();
    pointVectorFields_.clear();
    pointTensorFields_.clear();
    pointSymmTensorFields_.clear();
    pointDiagTensorFields_.clear();
    pointSphericalTensorFields_.clear();
}

// ************************************************************************* //
