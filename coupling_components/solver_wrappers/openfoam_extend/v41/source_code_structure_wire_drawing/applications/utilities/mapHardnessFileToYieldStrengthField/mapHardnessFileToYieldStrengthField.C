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

Application
    mapHardnessFileToYieldStrengthField

Description
    Reads in the a text file with a hardness map, converts this to yield
    strength, and interpolates this as a field to the mesh.

    It is assumed that the hardness map file takes the following format:
    
        # horizontalCoordinate verticalCoordinate VickersHardness
        0.1    0.3    301
        -0.2   0.0    350
        ...

    where the coordinates are in mm.

    Interpolation is performed using radial basis functions

Author
   Philip Cardiff, UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include <memory>
#include "RBFInterpolation.H"
#include "TPSFunction.H"
#include <cstdlib>

using namespace rbf;
//using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("hardnessFileName");
    argList::validArgs.append("patchName");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    // Read arguments
    const word hardnessFileName(args.additionalArgs()[0]);
    const word patchName(args.additionalArgs()[1]);

    // Read the hardness map file
    Info<< "Reading hardness map file: " << hardnessFileName << endl;
    IFstream hardnessFile(hardnessFileName);

    // Check the patch name
    const label patchID = mesh.boundaryMesh().findPatchID(patchName);
    const fvPatch& fpatch = mesh.boundary()[patchID];
    const vectorField& fpatchCf = fpatch.Cf();
    const scalar avgFpatchX = average(fpatchCf.component(Foam::vector::X));

    // Skip first line (assumed to be a header line)
    string tempString;
    hardnessFile.getLine(tempString);

    // Read the data lines
    int nData = 0;
    DynamicList<Foam::vector> hardnessCoords(10000);
    DynamicList<scalar> yieldStrength(10000);
    while (hardnessFile.getLine(tempString))
    {
        // Convert to a list of words
        wordList values(fileName(tempString).components(' '));
        //Info<< "values " << values << endl;

        if (values.size() < 3)
        {
            FatalError
                << "Line " << int(nData + 1) << " has less than three values"
                << abort(FatalError);
        }

        // We assume the out of plane direction is x
        hardnessCoords.append
        (
            Foam::vector
            (
                avgFpatchX,
                0.001*strtod(values[1].c_str(), NULL),
                0.001*strtod(values[0].c_str(), NULL)
            )
        );

        // Multiply Vickers hardness by 3 to approximately convert to yield
        // strength in MPa and by 1e6 to convert to Pa
        yieldStrength.append(1e6*3*strtod(values[2].c_str(), NULL));
        nData++;
    }
    Info<< "Number of data points: " << nData << endl;
    hardnessCoords.shrink();
    yieldStrength.shrink();

    const vectorField controlPoints(hardnessCoords); 
    const scalarField controlValues(yieldStrength); 

    //Info<< "Control points = " << controlPoints << endl;

    // Create an RBF interpolator fram hardness points to patch face centres
    Info<< "Creating the RBF function" << endl;
    std::shared_ptr<RBFFunctionInterface> rbfFunction =
        std::shared_ptr<RBFFunctionInterface>(new TPSFunction());

    Info<< "Creating the RBF interpolator" << endl;
    RBFInterpolation rbfInterp(rbfFunction);

    matrix controlX(controlPoints.size(), 3);
    matrix dataX(fpatchCf.size(), 3);

    forAll(controlPoints, pointI)
    {
        controlX(pointI, 0) = controlPoints[pointI].x();
        controlX(pointI, 1) = controlPoints[pointI].y();
        controlX(pointI, 2) = controlPoints[pointI].z();
    }

    forAll(fpatchCf, faceI)
    {
        dataX(faceI, 0) = fpatchCf[faceI].x();
        dataX(faceI, 1) = fpatchCf[faceI].y();
        dataX(faceI, 2) = fpatchCf[faceI].z();
    }

    // Compute the RBFs
    Info<< "Computing the RBF" << endl;
    rbfInterp.compute(controlX, dataX);

    // Create yield strength file
    Info<< "Creating the sigmaY field" << endl;
    volScalarField sigmaY
    (
        IOobject
        (
            "sigmaY",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    );

    // Interpolate the yield strength to the patch
    Info<< "Interpolating the field" << endl;

    matrix controlY(controlValues.size(), 1);
    matrix dataY(fpatchCf.size(), 1);

    forAll(controlValues, pointI)
    {
        controlY(pointI, 0) = controlValues[pointI];
    }

    rbfInterp.interpolate(controlY, dataY);

    const labelList& faceCells = fpatch.faceCells();
    forAll(fpatchCf, faceI)
    {
        sigmaY.boundaryField()[patchID][faceI] = dataY(faceI, 0);
        sigmaY.internalField()[faceCells[faceI]] = dataY(faceI, 0);
    }

    // Write the field
    Info<< "Writing sigmaY field" << endl;
    sigmaY.write();

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
