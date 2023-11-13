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

\*---------------------------------------------------------------------------*/

#include "wireStreamlines.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wireStreamlines, 0);
}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void Foam::wireStreamlines::makePointStreamlines() const
{
    if (pointStreamlinesPtr_)
    {
        FatalErrorIn("void Foam::wireStreamlines::makePointStreamlines() const")
            << "pointer already set" << abort(FatalError);
    }

    // We assume the wire is flowing in the positive x direction
    const vector flowDir = vector(1, 0, 0);

    // Take a reference to the wire mesh
    const polyMesh& mesh = mesh_;
    const pointField& points = mesh.points();
    const labelListList& pointPoints = mesh.pointPoints();

    // Take a reference to the wire upstream polyPatch
    const polyPatch& ppatch = mesh.boundaryMesh()[upstreamPatchID_];
    const labelList& meshPoints = ppatch.meshPoints();

    // Set the pointer list size to the number of points on the upstream patch
    pointStreamlinesPtr_ = new List< SLList<label> >(ppatch.nPoints());
    List< SLList<label> >& pointStreamlines = *pointStreamlinesPtr_;

    // We will set each point streamline
    // Each point streamline consists of a list of pointIDs along the wire
    // starting at a point on the the upstream patch

    forAll(pointStreamlines, streamlineI)
    {
        const label pointID = meshPoints[streamlineI];

        // Insert first point
        pointStreamlines[streamlineI].insert(pointID);

        SLList<label> pointsToCheck;
        pointsToCheck.insert(pointID);

        do
        {
            const label pointToCheckID = pointsToCheck.removeHead();
            const vector& curPointToCheck = points[pointToCheckID];

            // Find pointPoint in most aligned with flow direction

            const labelList& curPointPoints = pointPoints[pointToCheckID];

            label nextPointID = -1;
            scalar bestDir = 0;

            forAll(curPointPoints, ppI)
            {
                const label pointPointID = curPointPoints[ppI];

                // Vector from the current point to the current pointPoint
                // neighbour
                vector localDir = points[pointPointID] - curPointToCheck;
                localDir /= mag(localDir);

                // Calculate the dot product with the wire flow direction
                scalar b = localDir & flowDir;

                // Find the vectors that is most aligned with thw wire flow
                // direction
                if (b > bestDir && b > 0.7071)
                {
                    bestDir = b;
                    nextPointID = pointPointID;
                }
            }

            // If we found a viable point, then we store it and add it to the
            // points to be checked
            if (nextPointID != -1)
            {
                // Add the next point to the streamline
                pointStreamlines[streamlineI].append(nextPointID);

                // Add the next point to be checked
                pointsToCheck.insert(nextPointID);
            }
        } while (pointsToCheck.size());
    }
}


void Foam::wireStreamlines::makeCellStreamlines() const
{
    if (cellStreamlinesPtr_)
    {
        FatalErrorIn("void Foam::wireStreamlines::makeCellStreamlines() const")
            << "pointer already set" << abort(FatalError);
    }

    // We assume the wire is flowing in the positive x direction
    const vector flowDir = vector(1, 0, 0);

    // Take a reference to the wire mesh
    const polyMesh& mesh = mesh_;
    const vectorField& C = mesh.cellCentres();
    const labelListList& cellCells = mesh.cellCells();

    // Take a reference to the wire upstream polyPatch
    const polyPatch& ppatch = mesh.boundaryMesh()[upstreamPatchID_];
    const labelList& faceCells = ppatch.faceCells();

    cellStreamlinesPtr_ = new List< SLList<label> >(ppatch.size());

    List< SLList<label> >& cellStreamlines = *cellStreamlinesPtr_;

    forAll(cellStreamlines, streamlineI)
    {
        const label cellID = faceCells[streamlineI];

        // Insert first point
        cellStreamlines[streamlineI].insert(cellID);

        SLList<label> cellsToCheck;
        cellsToCheck.insert(cellID);

        do
        {
            const label cellToCheckID = cellsToCheck.removeHead();
            const vector& curCellToCheck = C[cellToCheckID];

            // Find cellCell in most aligned with flow direction

            const labelList& curCellCells = cellCells[cellToCheckID];

            label nextCellID = -1;
            scalar bestDir = 0;

            forAll(curCellCells, ppI)
            {
                const label cellCellID = curCellCells[ppI];

                vector localDir = C[cellCellID] - curCellToCheck;
                localDir /= mag(localDir);

                const scalar b = flowDir & localDir;

                if (b > bestDir && b > 0.7071)
                {
                    bestDir = b;
                    nextCellID = cellCellID;
                }
            }

            if (nextCellID != -1)
            {
                // Add the next cell
                cellStreamlines[streamlineI].append(nextCellID);
                cellsToCheck.insert(nextCellID);
            }
        } while (cellsToCheck.size());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wireStreamlines::wireStreamlines(const polyMesh& wireMesh)
    :
    mesh_(wireMesh),
    upstreamPatchID_(-1),
    pointStreamlinesPtr_(NULL),
    cellStreamlinesPtr_(NULL)
{
    // Find the upstream patch

    // We assume that the wire always flows in the positive x direction, then
    // the upstream patch will be the patch with the most negative x coordinates
    // on average
    // CAREFUL: this method can fail if the die upstream patch has a more
    // positive X coordinate

    scalar minX = GREAT;

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const scalar patchAverageX =
            average(mesh_.boundaryMesh()[patchI].faceCentres()).x();

        if (patchAverageX < minX)
        {
            minX = patchAverageX;
            upstreamPatchID_ = patchI;
        }
    }

    if (upstreamPatchID_ == -1)
    {
        FatalErrorIn
        (
            "Foam::wireStreamlines::wireStreamlines"
            "(const polyMesh& wireMesh)"
        )   << "Something went wrong finding the wire upstream patch!"
            << abort(FatalError);
    }

    Info<< "    " << type() << ": wire upstream patch: "
        << mesh_.boundaryMesh()[upstreamPatchID_].name()
        << endl;
}


Foam::wireStreamlines::wireStreamlines
(
    const polyMesh& wireMesh,
    const word& wireUpstreamPatchName
)
    :
    mesh_(wireMesh),
    upstreamPatchID_(-1),
    pointStreamlinesPtr_(NULL),
    cellStreamlinesPtr_(NULL)
{
    // Find the upstream patch

    upstreamPatchID_ = mesh_.boundaryMesh().findPatchID(wireUpstreamPatchName);

    if (upstreamPatchID_ == -1)
    {
        FatalErrorIn
        (
            "Foam::wireStreamlines::wireStreamlines"
            "(const polyMesh& wireMesh, const word& wireUpstreamPatchName)"
        )   << "Wire upstream patch not found: " << wireUpstreamPatchName
            << abort(FatalError);
    }

    Info<< "    " << type() << ": wire upstream patch: "
        << mesh_.boundaryMesh()[upstreamPatchID_].name()
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::wireStreamlines::~wireStreamlines()
{}


// * * * * * * * * * * * * * * * Public Member Functions  * * *  * * * * * * //


Foam::List< Foam::SLList<Foam::label> >&
Foam::wireStreamlines::pointStreamlines()
{
    if (!pointStreamlinesPtr_)
    {
        makePointStreamlines();
    }

    return *pointStreamlinesPtr_;
}


const Foam::List< Foam::SLList<Foam::label> >&
Foam::wireStreamlines::pointStreamlines() const
{
    if (!pointStreamlinesPtr_)
    {
        makePointStreamlines();
    }

    return *pointStreamlinesPtr_;
}


Foam::List< Foam::SLList<Foam::label> >&
Foam::wireStreamlines::cellStreamlines()
{
    if (!cellStreamlinesPtr_)
    {
        makeCellStreamlines();
    }

    return *cellStreamlinesPtr_;
}


const Foam::List< Foam::SLList<Foam::label> >&
Foam::wireStreamlines::cellStreamlines() const
{
    if (!cellStreamlinesPtr_)
    {
        makeCellStreamlines();
    }

    return *cellStreamlinesPtr_;
}


// ************************************************************************* //
