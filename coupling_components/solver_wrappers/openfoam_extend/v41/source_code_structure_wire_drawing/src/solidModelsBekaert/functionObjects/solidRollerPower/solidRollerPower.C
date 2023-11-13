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

\*----------------------------------------------------------------------------*/

#include "solidRollerPower.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "rollerFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidRollerPower, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidRollerPower,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidRollerPower::writeData()
{
    if (patchFound_)
    {
        const fvMesh& mesh = time_.lookupObject<fvMesh>("region0");

        if (mesh.foundObject<volSymmTensorField>(stressName_))
        {
            // Calculate the net external power at the specified patch
            // externalPower = surface integral of force dotted with velocity

            // Lookup the Cauchy stress tensor field
            const symmTensorField& sigma =
                mesh.lookupObject<volSymmTensorField>
                (
                    stressName_
                ).boundaryField()[historyPatchID_];

            // Lookup the displacement increment patch field and cast to a
            // roller boundary condition
            const rollerFvPatchVectorField& roller =
                refCast<const rollerFvPatchVectorField>
                (
                    mesh.lookupObject<volVectorField>
                    (
                        "DU"
                    ).boundaryField()[historyPatchID_]
                );


            // Calculate moment arms

            vectorField patchC = mesh.C().boundaryField()[historyPatchID_];

            const vector axis = roller.rotationAxis();

            const vector pointOnAxis =
                roller.initialRotationOrigin()
              + roller.currentAxisDisplacement();

            vectorField r = patchC - pointOnAxis;
            r -= axis*(axis & r);


            // Calculate torque

            const vectorField force =
                mesh.Sf().boundaryField()[historyPatchID_] & sigma;

            const vector torque = gSum(r ^ force);


            // Calculate rotational power
            // P = torque & (rotationAxis*rotationalVelocity)
            // where rotationalVelocity is in rad/s

            const vector rotationalVelocity =
                axis*roller.currentRpm()*2.0*mathematicalConstant::pi/60.0;

            const scalar power = torque & rotationalVelocity;

            if (Pstream::master())
            {
                historyFilePtr_()
                    << time_.time().value() << " " << power << endl;
            }
        }
        else
        {
            InfoIn(this->name() + " function object constructor")
                << stressName_ << " field not found!" << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidRollerPower::solidRollerPower
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    historyPatchID_(-1),
    patchFound_(false),
    stressName_
    (
        dict.found("stressName")
      ? word(dict.lookup("stressName"))
      : word("sigmaCauchy")
    ),
    historyFilePtr_(NULL)
{
    Info<< "Creating " << this->name() << " function object" << endl;

    word historyPatchName("notSpecified");
    if (dict.found("historyPatch"))
    {
        dict.lookup("historyPatch")
            >> historyPatchName;
    }
    else
    {
        WarningIn(this->name() + " function object constructor")
            << "solidRollerPower: historyPatch not specified" << endl;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>("region0");

    historyPatchID_ = mesh.boundaryMesh().findPatchID(historyPatchName);
    if (historyPatchID_ == -1)
    {
        WarningIn(this->name() + " function object constructor")
            << "history patch " << historyPatchName << " not found"
            << endl;
    }
    else
    {
        patchFound_ = true;
    }

    // Create history file if not already created
    if (historyFilePtr_.empty() && patchFound_)
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                time_.timeName(mesh.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"history"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"history"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            // Open new file at start up
            historyFilePtr_.reset
            (
                new OFstream
                (
                    historyDir/"solidRollerPower" + historyPatchName + ".dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time" << " "
                    << "rotationalPower" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidRollerPower::start()
{
    return writeData();
}


#if FOAMEXTEND > 40
bool Foam::solidRollerPower::execute(const bool forceWrite)
#else
bool Foam::solidRollerPower::execute()
#endif
{
    return writeData();
}


bool Foam::solidRollerPower::read(const dictionary& dict)
{
    return true;
}

// ************************************************************************* //
