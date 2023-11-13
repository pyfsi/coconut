/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "passItem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(passItem, 0);
defineRunTimeSelectionTable(passItem, dictionary);


// * * * * * * * * * *  Private Member Functions  * * * * * * * * * * * * * * //

void passItem::readMesh(const Time& runTime) const
{
    if (!meshPtr_.empty())
    {
        FatalErrorIn("void Foam::passItem::readMesh()")
            << "pointer already set" << abort(FatalError);
    }

    if (debug)
    {
        Info<< "        Reading the mesh for region " << name() << endl;
    }

    meshPtr_.set
    (
        new fvMesh
        (
            IOobject
            (
                name(),
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

passItem::passItem
(
    const processPass& parent,
    const word& name,
    const dictionary& dict,
    const int& passNr,
    dataHandler* dataContainerPtr
)
:
    parentPtr_(&parent),
    name_(name),
    dict_(dict),
    passNr_(passNr),
    wireLength_(),
    wireAxialResolution_(),
    dataContainerPtr_(dataContainerPtr)
{}


// * * * * * * * * * * * * * * * * Public Members  * * * * * * * * * * * * * //


const processPass& passItem::parent() const
{
    return *parentPtr_;
};


fvMesh& passItem::mesh(const Time& runTime)
{
    // Force mesh to be re-read
    meshPtr_.clear();

    if (meshPtr_.empty())
    {
        readMesh(runTime);
    }

    return meshPtr_();
}


const fvMesh& passItem::mesh(const Time& runTime) const
{
    // Force mesh to be re-read
    meshPtr_.clear();

    if (meshPtr_.empty())
    {
        readMesh(runTime);
    }

    // Re-read mesh from disk in case it has changed
    meshPtr_->readUpdate();

    return meshPtr_();
}


void passItem::positionMesh(const Time& runTime, const fvMesh& wireMesh)
{
    // Lookup translation vector from the dict
    const vector translation =
        meshInputDict().lookupOrDefault<vector>
        (
            "positioningVector",
            vector::zero
        );

    positionMesh(runTime, translation);
}


void passItem::positionMesh(const Time& runTime, const vector& translation)
{
    if (mag(translation) > SMALL)
    {
        // Change to pass directory
        chDir(runTime.path());


        // We will use the transformPoints utility to translate the sub-meshes

        dataContainer().runSystemCommand
        (
            word
            (
                "transformPoints -region " + name()
                + " -translate "
                + "\"("
                + Foam::name(translation.x()) + " "
                + Foam::name(translation.y()) + " "
                + Foam::name(translation.z())
                + ")\""
            ),
            "log.transformPoints_" + name(),
            "void passItem::positionMesh"
        );
    }
}


void passItem::setup(Time& runTime)
{
    // Create passItem subMesh directories
    mkDir(runTime.path()/"constant"/name());
    mkDir(runTime.path()/"system"/name());
    mkDir(runTime.path()/"0"/name());

    // In case of a roller passItem, there can be multiple rollers in one item.
    // How to deal with this, thinking of meshing, setting BCs, material etc...

    // Setup passItem mesh
    setupMesh(runTime);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// ************************************************************************* //
