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

#include "volAvgStress.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(volAvgStress, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        volAvgStress,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::volAvgStress::writeData()
{
  const fvMesh& mesh =
    time_.lookupObject<fvMesh>("region0");

  if (mesh.foundObject<volSymmTensorField>(stressName_))
    {
      const symmTensorField& sigma =
	mesh.lookupObject<volSymmTensorField>
	(
	 stressName_
	 ).internalField();

      symmTensor avStress = gSum(mesh.V()*sigma)/gSum(mesh.V());

      if (Pstream::master())
	{
	  historyFilePtr_()
	    << time_.time().value()
	    << " " << avStress.xx() << " " << avStress.xy()
	    << " " << avStress.xz() << " " << avStress.yy()
	    << " " << avStress.yz() << " " << avStress.zz();

	  historyFilePtr_() << endl;
	}
    }
  else
    {
      InfoIn(this->name() + " function object constructor")
	<< stressName_ << " not found" << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::volAvgStress::volAvgStress
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    stressName_
    (
        dict.found("stressName")
        ? word(dict.lookup("stressName"))
        : word("sigma")
    ),
    historyFilePtr_(NULL)
{
    Info<< "Creating " << this->name() << " function object" << endl;


    const fvMesh& mesh =
        time_.lookupObject<fvMesh>("region0");


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
                    new OFstream
                    (
                        historyDir/"volAvgStress.dat"
                    )
                );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time" << " "
                    << "stressXX" << " " << "stressXY" << " "
                    << "stressXZ" << " " << "stressYY" << " "
                    << "stressYZ" << " " << "stressZZ" << " " << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::volAvgStress::start()
{
    return writeData();
}


#if FOAMEXTEND > 40
bool Foam::volAvgStress::execute(const bool forceWrite)
#else
bool Foam::volAvgStress::execute()
#endif
{
    return writeData();
}


bool Foam::volAvgStress::read(const dictionary& dict)
{
    return true;
}

// ************************************************************************* //
