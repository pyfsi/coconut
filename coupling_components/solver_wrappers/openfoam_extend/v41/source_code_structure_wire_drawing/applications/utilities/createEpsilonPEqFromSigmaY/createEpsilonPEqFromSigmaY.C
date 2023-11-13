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
    createEpsilonPEqFromSigmaY

Description
    Reads in the sigmaY (yield stress) field and plasticStrainVsYieldStress file.
    An equivalent epsilonPEq (equivalent plastic strain) field is then created
    and written.

Author
   Philip Cardiff, UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("plasticStrainVsYieldStressFileName");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    // Read arguments
    const word epsilonPEqVsSigmaYFileName(args.additionalArgs()[0]);

    // Read the sigmaY field
    Info<< "Reading sigmaY field "<< nl << endl;
    volScalarField sigmaY
    (
        IOobject
        (
            "sigmaY",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Create the epsilonPEq field
    Info<< "Creating epsilonPEq field "<< nl << endl;
    volScalarField epsilonPEq
    (
        IOobject
        (
            "epsilonPEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    // Read the plasticStrainVsYieldStress file
    IFstream inFile(epsilonPEqVsSigmaYFileName);
    List<Tuple2<scalar, scalar> > inData;
    inFile
        >> inData;

    // We will now reverse the order of the XY data ie. switch X with Y and
    // re-write the new table to disk
    forAll(inData, i)
    {
        Swap(inData[i].first(), inData[i].second());
    }
    fileName sigmaYToEpsilonPEqFileName("sigmaYToEpsilonPEq");
    {
        OFstream outFile(sigmaYToEpsilonPEqFileName);
        outFile
            << inData;
    }

    // Read in the new XY data relating sigmaY to epsilonPEq
    interpolationTable<scalar> sigmaYToEpsilonPEq(sigmaYToEpsilonPEqFileName);

    // Set epsilonPEq internal values
    Info<< "Setting epsilonPEq based on sigmaY "<< nl << endl;
    scalarField& epsilonPEqI = epsilonPEq.internalField();
    scalarField& sigmaYI = sigmaY.internalField();
    forAll(epsilonPEqI, cellI)
    {
        if (sigmaYI[cellI] > SMALL)
        {
            epsilonPEqI[cellI] = sigmaYToEpsilonPEq(sigmaYI[cellI]);
        }
    }

    // Set epsilonPEq boundary values
    forAll(epsilonPEq.boundaryField(), patchI)
    {
        scalarField& epsilonPEqP = epsilonPEq.boundaryField()[patchI];
        scalarField& sigmaYP = sigmaY.boundaryField()[patchI];
        forAll(epsilonPEqP, faceI)
        {
            if (sigmaYP[faceI] > SMALL)
            {
                epsilonPEqP[faceI] = sigmaYToEpsilonPEq(sigmaYP[faceI]);
            }
        }
    }

    // Write the field
    Info<< "Writing epsilonPEq field" << endl;
    epsilonPEq.write();

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
