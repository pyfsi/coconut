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

#include "solidDisp.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidDisp, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidDisp,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidDisp::writeData()
{
    if (patchFound_)
    {
        const fvMesh& mesh =
            time_.lookupObject<fvMesh>("region0");

        if (mesh.foundObject<volVectorField>("U"))
        {
            const vectorField& U =
                mesh.lookupObject<volVectorField>
                (
                    "U"
                ).boundaryField()[historyPatchID_];

            const scalar minX = gMin(U.component(vector::X));
            const scalar minY = gMin(U.component(vector::Y));
            const scalar minZ = gMin(U.component(vector::Z));

            const scalar maxX = gMax(U.component(vector::X));
            const scalar maxY = gMax(U.component(vector::Y));
            const scalar maxZ = gMax(U.component(vector::Z));

            const vector avDisp = gAverage(U);

            if (Pstream::master())
            {
                historyFilePtr_()
                    << time_.time().value()
                    << " " << minX
                    << " " << minY
                    << " " << minZ
                    << " " << maxX
                    << " " << maxY
                    << " " << maxZ
                    << " " << avDisp.x()
                    << " " << avDisp.y()
                    << " " << avDisp.z()
                    << endl;
            }
        }
        else
        {
            InfoIn(this->name() + " function object constructor")
                << "U not found" << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidDisp::solidDisp
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
    historyFilePtr_(NULL)
{
    Info<< "Creating " << this->name() << " function object" << endl;

    word historyPatchName("notSpecified");

    if (dict.found("historyPatch"))
    {
        dict.lookup("historyPatch") >> historyPatchName;
    }
    else
    {
        WarningIn(this->name() + " function object constructor")
            << "solidDisp: historyPatch not specified" << endl;
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
                    historyDir/"solidDisp"+historyPatchName+".dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time"
                    << " " << "minX"
                    << " " << "minY"
                    << " " << "minZ"
                    << " " << "maxX"
                    << " " << "maxY"
                    << " " << "maxZ"
                    << " " << "avX"
                    << " " << "avY"
                    << " " << "avZ"
                    << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidDisp::start()
{
    return writeData();
}


#if FOAMEXTEND > 40
bool Foam::solidDisp::execute(const bool forceWrite)
#else
bool Foam::solidDisp::execute()
#endif
{
    return writeData();
}


bool Foam::solidDisp::read(const dictionary& dict)
{
    return true;
}

// ************************************************************************* //
