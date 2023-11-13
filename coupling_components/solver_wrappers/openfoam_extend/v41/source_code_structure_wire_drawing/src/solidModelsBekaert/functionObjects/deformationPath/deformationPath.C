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

#include "deformationPath.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(deformationPath, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        deformationPath,
        dictionary
     );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::deformationPath::writeData()
{
    const tensorField& defGradI =
        mesh_.lookupObject<volTensorField>
        (
         "F"
         ).internalField();

    // Loop over each cell, creating a new OFstream object
    // each time.

    // The OFstream for the last cell is automatically closed upon object deletion
    // the old method which had all the cells open at once and pointed to by
    // historyFilePtrList_ does not behave as expected with increased mesh sizes
    // for unknown reasons.
      
    forAll(cellsOfInterest_, cellI)
    {

      fileName historyDir;
      word startTimeName =
	time_.timeName(mesh_.time().startTime().value());

      historyDir = time_.path()/"history"/startTimeName/"deformationPaths";


      word deformationPathString = "deformationPath_" + Foam::name(cellI) + ".dat";
      fileName deformationPathFile(deformationPathString);

      // The ios_base::app argument indicates we wish to append to an existing file
      OFstream historyFile(historyDir/deformationPathFile, ios_base::app);

      historyFile.setEof(); // jump to the end

      // note: only the master cells will be recorded with this,
      // a more general approach is needed which also avoids separate
      // processors editing the same file
      if (Pstream::master())
        {
	  historyFile()
	    << "    ( "
	    << time_.time().value()
	    << "    (" << defGradI[cellI].xx() << " " << defGradI[cellI].xy()
	    << " " << defGradI[cellI].xz() << " " << defGradI[cellI].yx()
	    << " " << defGradI[cellI].yy() << " " << defGradI[cellI].yz()
	    << " " << defGradI[cellI].zx() << " " << defGradI[cellI].zy()
	    << " " << defGradI[cellI].zz() << ") )";

	  historyFile() << endl;

	  if ( time_.time().value() == time_.endTime().value() )
            {
	      historyFile() << ")";
            }
        }

    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::deformationPath::deformationPath
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    mesh_(time_.lookupObject<fvMesh>("region0")),
    nCells_(mesh_.nCells()),
    cellsOfInterest_
    (
     dict.found("cellZone")
     ? (labelList)mesh_.cellZones()[mesh_.cellZones().findZoneID(word(dict.lookup("cellZone")))]
     : labelList()
    )

{
    Info<< "Creating " << this->name() << " function object" << endl;

    // If a specific region isnt chosen set the cellsOfInterest_ to all
    // cells
    if (!dict.found("cellZone"))
      {
	Info<< "No cellZone region specified, defaulting to all " << nCells_
	    << " cells" << endl;
	cellsOfInterest_.setSize(nCells_);
	for (int cellI=0; cellI<nCells_; cellI++)
	  {
	    cellsOfInterest_[cellI] = cellI;
	  }
      }
    else
      {
	Info<< "Saving deformation gradient history for all " << cellsOfInterest_.size()
	    << " cells in mesh zone " << word(dict.lookup("cellZone")) << endl;
      }
    
    forAll(cellsOfInterest_, cellI)
    {
        fileName historyDir;

        word startTimeName =
            time_.timeName(mesh_.time().startTime().value());

        if (Pstream::parRun())
        {
	  FatalError << "This function object currently only works in serial"
		     << abort(FatalError);
        }
        else
        {
            historyDir = time_.path()/"history"/startTimeName/"deformationPaths";
        }

        word deformationPathString = "deformationPath_" + Foam::name(cellI) + ".dat";
        fileName deformationPathFile(deformationPathString);

	OFstream historyFile(historyDir/deformationPathFile);

        if (Pstream::master())
        {
            // Create directory if does not exist.
            mkDir(historyDir);

            historyFile
                << "(" << endl;
        }
    }
}

// * * * * * * * * * * * * * * *   Destructors     * * * * * * * * * * * * * //

// Foam::deformationPath::~deformationPath()
// {
// }

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::deformationPath::start()
{
    return writeData();
}


#if FOAMEXTEND > 40
bool Foam::deformationPath::execute(const bool forceWrite)
#else
bool Foam::deformationPath::execute()
#endif
{
    return writeData();
}


bool Foam::deformationPath::read(const dictionary& dict)
{
    return true;
}

// ************************************************************************* //
