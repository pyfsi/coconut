/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
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

Application
    extractDeformationGradientHistory

Description
    This post-processing utility is limited to running on process line cases where there is a wire
    mesh zone, and an inlet patch named "wireUpstream"

    The controlDict startTime should be set to "latestTime"

    A tensor interpolation table is created along the wire axis (assumed to be the X direction)
    for each boundary cell on wireUpstream. This may assumption may not be valid if sharp direction changes
    in the X direction are experienced.

    The tensor interpolation series are intended for use in RVE simulations.

    NOTE: This utility requires the mesh to be created by the extrusion of a slice (ordered and neat)

Author
    Michael Clancy UCD
    Philip Cardiff UCD
    Peter De Jaeger Bekaert


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "newLeastSquaresVolPointInterpolation.H"
#include "twoDPointCorrector.H"
#include "symmetryPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Total deformation gradient
    Info << "runTime.timeName(): " << runTime.timeName() << endl;
    
    volTensorField F
    (
         IOobject
         (
             "F",
             runTime.timeName(),
             mesh,
             IOobject::MUST_READ,
             IOobject::AUTO_WRITE
          ),
         mesh
     );

    const tensorField& defGradI = F.internalField();
    
    label wireUpstreamID = mesh.boundaryMesh().findPatchID("wireUpstream");
    label wireDownstreamID = mesh.boundaryMesh().findPatchID("wireDownStream1");

    Info<< "patchID for wireUpstream: " << wireUpstreamID << endl;

    labelList downstreamCells = mesh.boundaryMesh()[wireDownstreamID].faceCells();

    
    forAll(mesh.boundaryMesh()[wireUpstreamID].faceCells(), faceCellI)
    {

        label cellID = mesh.boundaryMesh()[wireUpstreamID].faceCells()[faceCellI];
        Info<< "Creating deformation gradient history file for faceCell " << cellID << endl;


        fileName historyDir = "./deformationPaths";
        mkDir(historyDir);

        word deformationPathString = "deformationPath_" + Foam::name(cellID) + ".dat";
        fileName deformationPathFile(deformationPathString);

        // The ios_base::app argument indicates we wish to append to an existing file
        // not necessary here
        OFstream historyFile(historyDir/deformationPathFile);//, ios_base::app);
        //historyFile.setEof(); // jump to the end
        historyFile()
            << "defGrad.xx()" << "," << "defGradI.xy()"
            << "," << "defGradI.xz()" << "," << "defGradI.yx()"
            << "," << "defGradI.yy()" << "," << "defGradI.yz()"
            << "," << "defGradI.zx()" << "," << "defGradI.zy()"
            << "," << "defGradI.zz()" << ","
            << "xLoc" << ","
            << "yLoc" << ","
            << "zLoc\n";

        
        
        bool hasXNeighbour = true;
        int idx = 0;
        do {
            if (Pstream::master())
            {
                historyFile()
                    << defGradI[cellID].xx() << "," << defGradI[cellID].xy()
                    << "," << defGradI[cellID].xz() << "," << defGradI[cellID].yx()
                    << "," << defGradI[cellID].yy() << "," << defGradI[cellID].yz()
                    << "," << defGradI[cellID].zx() << "," << defGradI[cellID].zy()
                    << "," << defGradI[cellID].zz() << ","
                    << mesh.C()[cellID].x() << ","
                    << mesh.C()[cellID].y() << ","
                    << mesh.C()[cellID].z();


                historyFile() << endl;
            }


            // check if the cell is on the downstream patch and exit if so
            forAll(downstreamCells, downstreamCellI)
            {
                if (cellID == downstreamCells[downstreamCellI])
                    hasXNeighbour = false;
            }

            // get the neighbouring cell with the greatest x distance
            scalar maxXDist = 0;
            label nextCellID = -1;
            if (hasXNeighbour)
            {
                forAll(mesh.cellCells()[cellID], i)
                {
                    label neighbourCellID = mesh.cellCells()[cellID][i];
                    scalar cellDist = mesh.C()[neighbourCellID].x() - mesh.C()[cellID].x();
                    if (cellDist > maxXDist)
                    {
                        maxXDist = cellDist;
                        nextCellID = neighbourCellID;
                    }
                }
                cellID = nextCellID;
            }

            idx++;

        } while (hasXNeighbour);

    }

    
    Info<< "\nEnd\n" << endl;

    return(0);
}


// ************************************************************************* //
