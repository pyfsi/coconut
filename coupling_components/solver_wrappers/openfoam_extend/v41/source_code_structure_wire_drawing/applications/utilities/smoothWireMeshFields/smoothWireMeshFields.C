/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    smoothWireMeshFields

Description
    A utility to smooth the wire mesh and fields in the flow direction.

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wireStreamlines.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "IOobjectList.H"
#include "readFields.H"
#include "SetAverage.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    // The new wire length must be specified
    argList::validArgs.append("newWireLength");

#   include "setRootCase.H"
#   include "createTime.H"

    // Read the new wire length
    const scalar newWireLength
    (
        readScalar(IStringStream(args.additionalArgs()[0])())
    );


    // Smooth the wire fields first, then the mesh

    Info<< nl << "Smoothing the wire fields" << endl;

    // Set runTime to the latest time-step so that the new fields will be
    // read from and written to the latest time-step
    const instantList& times = runTime.times();

    // Index of the latest time-step
    const label latestTimeID = times.size() - 1;

    // Set the previous time to be the latest time-step
    runTime.setTime(times[latestTimeID], latestTimeID);

    Info<< "Reading the wire mesh from time: " << runTime.value() << endl;

    // Read the wire mesh
    fvMesh mesh
    (
        IOobject
        (
            "wire",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // Create a wire streamline object
    wireStreamlines streamlines(mesh);

    // Get the cell streamlines
    const List< SLList<label> >& cellStreamlines =
        streamlines.cellStreamlines();

    // Find the end of the wire (most positive x)
    const scalar wireMaxX = max(mesh.points().component(vector::X));

    if (wireMaxX < SMALL)
    {
        FatalError
            << "Wire smoothing cannot be performed because the wire does not"
            << " have any region with an x cordinate greater than x = 0"
            << abort(FatalError);
    }

    // Calculate x coordinate limits for the averaging range
    // as middle of the wire that has passed by x = 0
    const scalar averagingRangeMinX = 0.5*wireMaxX - 0.1*wireMaxX;
    const scalar averagingRangeMaxX = 0.5*wireMaxX + 0.1*wireMaxX;
    
    Info << "    averagingRangeMinX, averagingRangeMaxX: "
         << averagingRangeMinX << ", " << averagingRangeMaxX << endl;

    // Create a list of streamlines including only the cells in the averaging
    // range
    List< SLList<label> > cellStreamlinesInRange(cellStreamlines.size());

    const vectorField& C = mesh.C().internalField();

    forAll(cellStreamlines, streamlineI)
    {
        const labelList curStreamline = cellStreamlines[streamlineI];

        forAll(curStreamline, cI)
        {
            // Check if cell is in the averaging range

            const label cellID = curStreamline[cI];

            const scalar nextCellX = C[cellID].x();

            if
            (
                nextCellX > averagingRangeMinX && nextCellX < averagingRangeMaxX
            )
            {
                cellStreamlinesInRange[streamlineI].append(cellID);
            }
        }
    }

    // Check if cells were found
    forAll(cellStreamlinesInRange, streamlineI)
    {
        if (cellStreamlinesInRange[streamlineI].size() == 0)
        {
            FatalError
                << "There were no cells found within the smooth range:"
                << " try making the 'range' larger" << exit(FatalError);
        }
    }

    Info<< "There are " << cellStreamlinesInRange[0].size()
        << " cells in the field-smoothing range" << endl;


    // Read the volFields into the objectRegistry
    // These fields will then be averaged along the streamlines within the
    // averaging range

    // Search for list of objects for this time
    Info<< "Reading wire fields from time: " << runTime.timeName() << endl;
    IOobjectList objects(mesh, runTime.timeName());

    PtrList<volScalarField> volScalarFields;
    readFields(mesh, objects, volScalarFields);

    PtrList<volVectorField> volVectorFields;
    readFields(mesh, objects, volVectorFields);

    PtrList<volSphericalTensorField> volSphericalTensorFields;
    readFields(mesh, objects, volSphericalTensorFields);

    PtrList<volSymmTensorField> volSymmTensorFields;
    readFields(mesh, objects, volSymmTensorFields);

    PtrList<volTensorField> volTensorFields;
    readFields(mesh, objects, volTensorFields);


    // Set fields to be average value along streamline within range
    forAll(cellStreamlines, cI)
    {
        const labelList streamCells = cellStreamlines[cI];
        const labelList streamCellsInRange = cellStreamlinesInRange[cI];

        SetAverage<scalar>(mesh, streamCells, streamCellsInRange);
        SetAverage<vector>(mesh, streamCells, streamCellsInRange);
        SetAverage<tensor>(mesh, streamCells, streamCellsInRange);
        SetAverage<symmTensor>(mesh, streamCells, streamCellsInRange);
        SetAverage<diagTensor>(mesh, streamCells, streamCellsInRange);
        SetAverage<sphericalTensor>(mesh, streamCells, streamCellsInRange);
    }


    // Smooth the wire mesh


    Info<< nl << "Smoothing the wire mesh" << endl;

    // Take a reference to the mesh points
    const pointField& points = mesh.points();

    // Get the point streamlines
    const List< SLList<label> >& pointStreamlines =
        streamlines.pointStreamlines();

    // We will calculate the new smooth wire points
    pointField newPoints = mesh.allPoints();

    // We consistently assume the wire flows in the positive x direction on
    // average
    const vector flowDir = vector(1, 0, 0);

    // Start x coordinate of the streamline for the new mesh
    const scalar newStartPoint = -newWireLength/2.0;

    // We will calculate the average point coordinate in the averaging range
    vectorField streamlineAveragePoint =
        vectorField(pointStreamlines.size(), vector::zero);

    // We will calculate the number of points in the averaging range
    labelList nStreamlineAveragePoint =
        labelList(pointStreamlines.size(), 0);

    // Loop through all streamlines
    forAll(pointStreamlines, streamlineI)
    {
        const labelList streamPoints = pointStreamlines[streamlineI];

        // Move along the points on this streamline
        forAll(streamPoints, pI)
        {
            // Check if the point is in the averaging range
            const label pointID = streamPoints[pI];
            const vector& curPoint = points[pointID];
            const scalar curPointX = flowDir & curPoint;

            if
            (
                curPointX > averagingRangeMinX && curPointX < averagingRangeMaxX
            )
            {
                // Take the radial component
                streamlineAveragePoint[streamlineI] +=
                    (I - sqr(flowDir)) & curPoint;
                nStreamlineAveragePoint[streamlineI]++;
            }
        }
    }

    if (min(nStreamlineAveragePoint) == 0)
    {
        FatalError
            << "No points found in the specified range:"
            << " try making the 'range' larger!" << abort(FatalError);
    }

    // Calculate the average
    forAll(nStreamlineAveragePoint, streamlineI)
    {
        streamlineAveragePoint[streamlineI]
            /= nStreamlineAveragePoint[streamlineI];
    }


    // Calculate the new mesh spacing
    const scalar meshSpacing = newWireLength/(pointStreamlines[0].size() - 1);


    // Set the new points on each streamlines
    forAll(pointStreamlines, streamlineI)
    {
        const labelList streamPoints = pointStreamlines[streamlineI];

        // Move along the points on this streamline
        forAll(streamPoints, pI)
        {
            const label pointID = streamPoints[pI];

            // Update the point coordinates, consisting of an axial and a radial
            // component, as well as a translation offset
            newPoints[pointID] =
                (flowDir*(newStartPoint + pI*meshSpacing))
              + ((I - sqr(flowDir)) & streamlineAveragePoint[streamlineI]);
        }
    }

    // Move the mesh points
    mesh.movePoints(newPoints);
    mesh.moving(false);
    mesh.changing(false);
    mesh.setPhi().writeOpt() = IOobject::NO_WRITE;

    // Write the smoothed wire mesh
    Info<< "Writing the smooth wire mesh and fields" << endl;
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
