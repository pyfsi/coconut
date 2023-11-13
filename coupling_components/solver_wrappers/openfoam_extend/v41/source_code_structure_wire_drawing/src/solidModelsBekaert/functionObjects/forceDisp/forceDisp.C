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

#include "forceDisp.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(forceDisp, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        forceDisp,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::forceDisp::writeData()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>("region0");

    if (mesh.foundObject<volSymmTensorField>("sigma"))
    {
        const symmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>
            (
                "sigma"
            ).boundaryField()[historyPatchID_];

        const vectorField& U =
            mesh.lookupObject<volVectorField>
            (
                "U"
            ).boundaryField()[historyPatchID_];

        // Patch force
        vector force = gSum(mesh.Sf().boundaryField()[historyPatchID_] & sigma);

        // Arithmetic average disp and force on patch
        vector avDisp = average(U);

        if (Pstream::master())
        {
            historyFilePtr_()
                << time_.time().value() << " "
                    << avDisp.x() << " "
                    << avDisp.y() << " "
                    << avDisp.z() << " "
                    << force.x() << " "
                    << force.y() << " "
                    << force.z();

            historyFilePtr_() << endl;
        }
    }
    else
    {
        Info<< "sigma not found" << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::forceDisp::forceDisp
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
    historyFilePtr_(NULL)
{
    Info << "Creating " << this->name() << " function object." << endl;

    word historyPatchName("notSpecified");
    if (dict.found("historyPatch"))
    {
        dict.lookup("historyPatch") >> historyPatchName;
    }
    else
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "historyPatch not specified."
            << abort(FatalError);
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>("region0");

   historyPatchID_ = mesh.boundaryMesh().findPatchID(historyPatchName);
    if (historyPatchID_ == -1)
    {
        FatalErrorIn(this->name() + " function object constructor")
            << "history patch " << historyPatchName << " not found."
            << abort(FatalError);
    }

    // Create history file if not already created
    if (historyFilePtr_.empty())
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
            historyFilePtr_.reset(new OFstream(historyDir/"history.dat"));

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time" << " "
                    << "dispX" << " " << "dispY" << " "
                    << "dispZ" << " "
                    << "forceX" << " " << "forceY" << " "
                    << "forceZ";

                historyFilePtr_() << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::forceDisp::start()
{
    return writeData();
}


#if FOAMEXTEND > 40
bool Foam::forceDisp::execute(const bool forceWrite)
#else
bool Foam::forceDisp::execute()
#endif
{
    return writeData();
}


bool Foam::forceDisp::read(const dictionary& dict)
{
    return true;
}

// ************************************************************************* //
