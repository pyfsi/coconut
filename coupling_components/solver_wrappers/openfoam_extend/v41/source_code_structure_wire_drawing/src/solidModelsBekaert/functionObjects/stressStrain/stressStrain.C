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

#include "stressStrain.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stressStrain, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        stressStrain,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::stressStrain::writeData()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>("region0");

    if (mesh.foundObject<volSymmTensorField>(sigmaName_))
    {
        const symmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>
            (
                sigmaName_
            ).boundaryField()[historyPatchID_];

        const symmTensorField& epsilon =
            mesh.lookupObject<volSymmTensorField>
            (
                epsilonName_
            ).boundaryField()[historyPatchID_];

        // Arithmetic average stress and strain on patch
        symmTensor avStress = average(sigma);
        symmTensor avStrain = average(epsilon);

        if (Pstream::master())
        {
            historyFilePtr_()
                << time_.time().value() << " "
                    << avStrain.xx() << " "
                    << avStrain.xy() << " "
                    << avStrain.xz() << " "
                    << avStrain.yy() << " "
                    << avStrain.yz() << " "
                    << avStrain.zz() << " "
                    << avStress.xx() << " "
                    << avStress.xy() << " "
                    << avStress.xz() << " "
                    << avStress.yy() << " "
                    << avStress.yz() << " "
                    << avStress.zz();

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

Foam::stressStrain::stressStrain
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
    sigmaName_
    (
        dict.found("stressName")
        ? word(dict.lookup("stressName"))
        : word("sigma")
    ),
    epsilonName_
    (
        dict.found("strainName")
        ? word(dict.lookup("strainName"))
        : word("epsilon")
    ),
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
            historyFilePtr_.reset
                (
                    new OFstream(historyDir/"stressStrain"+historyPatchName+".dat")
                );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time" << " "
                    << "strainXX" << " " << "strainXY" << " "
                    << "strainXZ" << " " << "strainYY" << " "
                    << "strainYZ" << " " << "strainZZ" << " "
                    << "stressXX" << " " << "stressXY" << " "
                    << "stressXZ" << " " << "stressYY" << " "
                    << "stressYZ" << " " << "stressZZ";

                historyFilePtr_() << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::stressStrain::start()
{
    return writeData();
}


#if FOAMEXTEND > 40
bool Foam::stressStrain::execute(const bool forceWrite)
#else
bool Foam::stressStrain::execute()
#endif
{
    return writeData();
}


bool Foam::stressStrain::read(const dictionary& dict)
{
    return true;
}

// ************************************************************************* //
