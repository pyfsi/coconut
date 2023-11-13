/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "nagataPatch.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nagataPatch, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nagataPatch::calcFaceInterpCoeffs() const
{
    if (faceInterpCoeffsPtr_.valid())
    {
        FatalErrorIn(type() + "::calcFaceInterpCoeffs()")
            << "pointer already set" << abort(FatalError);
    }

    // Required mesh data
    const List<labelledTri>& faces = triSurf_.localFaces();
    const pointField& points = triSurf_.localPoints();
    const List< List<vector> >& faceCurvatureCoeffs =
        this->faceCurvatureCoeffs();

    if (debug > 1)
    {
        Info<< "localPoints: " << points << nl << nl
            << "faces: " << faces << nl << endl;
    }

    // Initialise the field of coefficients: there are six vectors for each face
    faceInterpCoeffsPtr_.set(new List< List<vector> >(triSurf_.size()));
    List< List<vector> >& faceInterpCoeffs = faceInterpCoeffsPtr_();

    forAll(faceInterpCoeffs, faceI)
    {
        List<vector>& curFaceInterpCoeffs = faceInterpCoeffs[faceI];
        const List<vector>& curFaceCurvatureCoeffs = faceCurvatureCoeffs[faceI];

        curFaceInterpCoeffs.setSize(6, vector::zero);

        const face& curFace = faces[faceI];

        const point& point0 = points[curFace[0]];
        const point& point1 = points[curFace[1]];
        const point& point2 = points[curFace[2]];

        // c1
        const vector& c1 = curFaceCurvatureCoeffs[0];

        // c2
        const vector& c2 = curFaceCurvatureCoeffs[1];

        // c3
        const vector& c3 = curFaceCurvatureCoeffs[2];

        // c00
        curFaceInterpCoeffs[0] = point0;

        // c10
        curFaceInterpCoeffs[1] = point1 - point0 - c1;

        // c01
        curFaceInterpCoeffs[2] = point2 - point1 + c1 - c3;

        // c11
        curFaceInterpCoeffs[3] = c3 - c1 - c2;

        // c20
        curFaceInterpCoeffs[4] = c1;

        // c02
        curFaceInterpCoeffs[5] = c2;
    }
}


const Foam::List< Foam::List<Foam::vector> >&
Foam::nagataPatch::faceInterpCoeffs() const
{
    if (faceInterpCoeffsPtr_.empty())
    {
        calcFaceInterpCoeffs();
    }

    return faceInterpCoeffsPtr_();
}


void Foam::nagataPatch::calcFaceCurvatureCoeffs() const
{
    if (faceCurvatureCoeffsPtr_.valid())
    {
        FatalErrorIn(type() + "::calcFaceCurvatureCoeffs()")
            << "pointer already set" << abort(FatalError);
    }

    // Required mesh data
    const List<labelledTri>& faces = triSurf_.localFaces();
    const pointField& points = triSurf_.localPoints();
    const pointField& pointNormals = triSurf_.pointNormals();
    const boolList& faceOnEdgeOfPatch = this->faceOnEdgeOfPatch();

    if (debug)
    {
        Info<< nl << "points and normals" << endl;
        forAll(points, pI)
        {
            Info<< "p: " << points[pI]
                << ", n: " << pointNormals[pI] << endl;
        }

        Info<< "faceNormals: " << triSurf_.faceNormals() << endl;
        Info<< "pointNormals: " << triSurf_.pointNormals() << endl;
    }

    // Initialise the field of coefficients: there are three vectors for each
    // face
    faceCurvatureCoeffsPtr_.set(new List< List<vector> >(triSurf_.size()));
    List< List<vector> >& faceCurvatureCoeffs = faceCurvatureCoeffsPtr_();

    forAll(faceCurvatureCoeffs, faceI)
    {
        List<vector>& curFaceCurvatureCoeffs = faceCurvatureCoeffs[faceI];

        curFaceCurvatureCoeffs.setSize(3, vector::zero);

        const face& curFace = faces[faceI];

        const label point0ID = curFace[0];
        const label point1ID = curFace[1];
        const label point2ID = curFace[2];

        const point& x0 = points[point0ID];
        const point& x1 = points[point1ID];
        const point& x2 = points[point2ID];

        // Take a copy of the normals as we may update them
        vector n0 = pointNormals[point0ID];
        vector n1 = pointNormals[point1ID];
        vector n2 = pointNormals[point2ID];

        if (correctPointNormalsOnPatchEdges_)
        {
            // Update point normal if they are on the edge of a patch/region
            // so that they are calculate using face normal in the same
            // patch/region. This allows us to keep sharp feature lines at the
            // edge of patches
            if (faceOnEdgeOfPatch[faceI])
            {
                const label curRegionID = triSurf_[faceI].region();
                updatePointNormalUsingFacesNormalInRegion
                (
                    n0, point0ID, curRegionID, triSurf_
                );
                updatePointNormalUsingFacesNormalInRegion
                (
                    n1, point1ID, curRegionID, triSurf_
                );
                updatePointNormalUsingFacesNormalInRegion
                (
                    n2, point2ID, curRegionID, triSurf_
                );
            }
        }

        // c1
        curvatureCoeff(curFaceCurvatureCoeffs[0], x0, x1, n0, n1);

        // c2
        curvatureCoeff(curFaceCurvatureCoeffs[1], x1, x2, n1, n2);

        // c3
        curvatureCoeff(curFaceCurvatureCoeffs[2], x0, x2, n0, n2);

        if (debug > 1)
        {
            Info<< "x0: " << x0 << nl
                << "x1: " << x1 << nl
                << "x2: " << x2 << nl
                << "n0: " << n0 << nl
                << "n1: " << n1 << nl
                << "n2: " << n2 << nl
                << "c1: " << curFaceCurvatureCoeffs[0] << nl
                << "c2: " << curFaceCurvatureCoeffs[1] << nl
                << "c3: " << curFaceCurvatureCoeffs[2] << nl
                << endl;
        }
    }
}


const Foam::List< Foam::List<Foam::vector> >&
Foam::nagataPatch::faceCurvatureCoeffs() const
{
    if (faceCurvatureCoeffsPtr_.empty())
    {
        calcFaceCurvatureCoeffs();
    }

    return faceCurvatureCoeffsPtr_();
}


void Foam::nagataPatch::localCoordinates
(
    scalar& eta, scalar& xi, const point& pt, const label triFaceID
)
{
    const labelledTri& curFace = triSurf_.localFaces()[triFaceID];
    const pointField& localPoints = triSurf_.localPoints();

    // Physical points
    const vector& x0 = localPoints[curFace[0]];
    const vector& x1 = localPoints[curFace[1]];
    const vector& x2 = localPoints[curFace[2]];

    // Calculate transformation matrix from physical space to parent space,
    // where x0 is at (0, 0), x1 is at (1, 0) and x2 is at (1, 1). The two local
    // coordinates are eta and xi
    const vector col3 = (x1 - x0)^(x2 - x1);
    const tensor A
    (
        x1.x() - x0.x(), x2.x() - x1.x(), col3.x(),
        x1.y() - x0.y(), x2.y() - x1.y(), col3.y(),
        x1.z() - x0.z(), x2.z() - x1.z(), col3.z()
    );

    // We could probably do this more efficiently but this is fine for now
    // Better: express in terms of analytical inverse
    const vector localCoords = inv(A) & (pt - x0);

    eta = localCoords.x();
    xi = localCoords.y();

    if (eta < -0.001 || eta > 1.001 || xi < -0.001 || xi > 1.001)
    {
        WarningIn(type() + "::localCoordinates(...)")
            << "eta and/or xi are less than 0.0 or greater than 1.0" << nl
            << "limiting eta and xi!" << nl
            << "eta = " << eta << ", xi = " << xi << ", z = " << localCoords.z()
            << nl << "A & localCoords = " << (A & localCoords)
            << nl << "inv(A) & pt = " << (inv(A) & pt) << endl;

        eta = max(min(eta, 1.0), 0.0);
        xi = max(min(xi, 1.0), 0.0);
    }
    else if
    (
        eta < -SMALL || eta > (1.0 + SMALL) || xi < -SMALL || xi > (1.0 + SMALL)
    )
    {
        // No need for a warning: just limit
        eta = max(min(eta, 1.0), 0.0);
        xi = max(min(xi, 1.0), 0.0);
    }
}


void Foam::nagataPatch::curvatureCoeff
(
    vector& coeff,
    const vector& x0,
    const vector& x1,
    const vector& n0,
    const vector& n1
) const
{
    const scalar a = n0 & n1;
    const scalar magA = mag(a);

    // Unit vector along the edge
    const vector b = (x1 - x0)/mag(x1 - x0);

    // In the singular case, we will default to a linear interpolation

    bool singularCase = false;

    if (magA > (1.0 - SMALL) || magA < (-1 + SMALL))
    {
        singularCase = true;
    }
    else if ((n0 & b)*(n1 & b) > 0.0)
    {
        singularCase = true;
    }
    else if
    (
        (mag(n0 & b) < epsilon1_ || mag(n1 & b) < epsilon1_)
     && mag((n0 & b) + (n1 & b)) > epsilon2_
    )
    {
        singularCase = true;
    }

    // Set curvature coefficient
    if (singularCase)
    {
        coeff = vector::zero;
    }
    else
    {
        RectangularMatrix<scalar> term1(3, 2);
        term1[0][0] = n0[0];
        term1[1][0] = n0[1];
        term1[2][0] = n0[2];
        term1[0][1] = n1[0];
        term1[1][1] = n1[1];
        term1[2][1] = n1[2];

        if (debug > 1)
        {
            Info<< "        term1: " << term1 << endl;
        }

        RectangularMatrix<scalar> term2(2, 2);
        term2[0][0] = 1;
        term2[1][0] = -a;
        term2[0][1] = -a;
        term2[1][1] = 1;

        if (debug > 1)
        {
            Info<< "        term2: " << term2 << endl
                << "        term1*term2: " << multiply(term1, term2) << endl;
        }

        RectangularMatrix<scalar> term3(2, 1);
        term3[0][0] = n0 & (x1 - x0);
        term3[1][0] = -n1 & (x1 - x0);

        if (debug > 1)
        {
            Info<< "        term3: " << term3 << endl
                << "        multiply(multiply(term1, term2), term3): "
                << multiply(multiply(term1, term2), term3) << endl;
        }

        RectangularMatrix<scalar> res =
            (1.0/(1.0 - a*a))*multiply(multiply(term1, term2), term3);

        if (debug > 1)
        {
            Info<< "        res: " << res << endl;
        }

        coeff.x() = res[0][0];
        coeff.y() = res[1][0];
        coeff.z() = res[2][0];
    }
}


Foam::RectangularMatrix<Foam::scalar> Foam::nagataPatch::multiply
(
    const Foam::RectangularMatrix<Foam::scalar>& a,
    const Foam::RectangularMatrix<Foam::scalar>& b
) const
{
    if (a.m() != b.n())
    {
        FatalErrorIn(type() + "::multiply(...")
            << "Matrices are incomplatible for multiplication"
            << abort(FatalError);
    }

    RectangularMatrix<scalar> result(a.n(), b.m(), 0.0);

    for (int i = 0; i < a.n(); i++)
    {
        for (int j = 0; j < b.m(); j++)
        {
            for (int k = 0; k < a.m(); k++)
            {
                result[i][j] += a[i][k]*b[k][j];
            }
        }
    }

    return result;
}


void Foam::nagataPatch::updatePointNormalUsingFacesNormalInRegion
(
    vector& n,
    const label pointID,
    const label curRegionID,
    const triSurface& triSurf
) const
{
    if (debug > 1)
    {
        Info<< "normal before: " << n << endl;
    }

    // Reset normal to zero
    n = vector::zero;

    // Take a reference to the pointFaces
    const labelList& curPointFaces = triSurf.pointFaces()[pointID];

    // Calculate the new normal using neighbour faces in the given region
    label nFacesOnPatch = 0;
    forAll(curPointFaces, pfI)
    {
        const label faceID = curPointFaces[pfI];

        if (triSurf[faceID].region() == curRegionID)
        {
            n += triSurf.faceNormals()[faceID];
            nFacesOnPatch++;
        }
    }

    // Normalise the normal
    n /= nFacesOnPatch;

    if (debug > 1)
    {
        Info<< "normal after: " << n << endl;
    }
}


const Foam::boolList& Foam::nagataPatch::faceOnEdgeOfPatch() const
{
    if (faceOnEdgeOfPatchPtr_.empty())
    {
        calcFaceOnEdgeOfPatch();
    }

    return faceOnEdgeOfPatchPtr_();
}


void Foam::nagataPatch::calcFaceOnEdgeOfPatch() const
{
    if (faceOnEdgeOfPatchPtr_.valid())
    {
        FatalErrorIn(type() + "::calcFaceOnEdgeOfPatch()")
            << "pointer already set" << abort(FatalError);
    }

    faceOnEdgeOfPatchPtr_.set(new boolList(triSurf_.size(), false));

    boolList& faceOnEdgeOfPatch = faceOnEdgeOfPatchPtr_();

    const labelListList& faceFaces = triSurf_.faceFaces();

    forAll(triSurf_, faceI)
    {
        if (faceOnEdgeOfPatch[faceI])
        {
            continue;
        }

        const label curFaceRegionID = triSurf_[faceI].region();

        const labelList& curFaceFaces = faceFaces[faceI];

        forAll(curFaceFaces, ffI)
        {
            const label neiFaceID = curFaceFaces[ffI];
            const label neiFaceRegionID = triSurf_[neiFaceID].region();

            // If two neighbouring faces share different region IDs then they
            // are both on the edge of a region/patch
            if (curFaceRegionID != neiFaceRegionID)
            {
                faceOnEdgeOfPatch[faceI] = true;
                faceOnEdgeOfPatch[neiFaceID] = true;
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::nagataPatch::nagataPatch
(
    const triSurface& triSurf,
    const scalar epsilon1,
    const scalar epsilon2,
    const bool correctPointNormalsOnPatchEdges
)
:
    triSurf_(triSurf),
    epsilon1_(epsilon1),
    epsilon2_(epsilon2),
    correctPointNormalsOnPatchEdges_(correctPointNormalsOnPatchEdges),
    faceInterpCoeffsPtr_(),
    faceCurvatureCoeffsPtr_(),
    faceOnEdgeOfPatchPtr_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nagataPatch::~nagataPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::nagataPatch::projectPoint
(
    point& projectedPoint,
    const List<vector>& faceCoeffs,
    const scalar eta,
    const scalar xi
)
{
    // Set the coefficients for this face
    const vector& c00 = faceCoeffs[0];
    const vector& c10 = faceCoeffs[1];
    const vector& c01 = faceCoeffs[2];
    const vector& c11 = faceCoeffs[3];
    const vector& c20 = faceCoeffs[4];
    const vector& c02 = faceCoeffs[5];

    // Calculate the projected point
    projectedPoint =
        c00 + c10*eta + c01*xi + c11*eta*xi + c20*eta*eta + c02*xi*xi;

    if (debug > 1)
    {
        InfoIn(type() + "::projectPoint(...)")
            << "projectedPoint: " << projectedPoint << nl
            << "eta: " << eta << nl
            << "xi: " << xi << nl
            << "c00: " << c00 << nl
            << "c10: " << c10 << nl
            << "c01: " << c01 << nl
            << "c11: " << c11 << nl
            << "c20: " << c20 << nl
            << "c02: " << c02 << nl
            << "c00: " << c00 << nl
            << "c10*eta: " << c10*eta << nl
            << "c01*xi: " << c01*xi << nl
            << "c11*eta*xi: " << c11*eta*xi << nl
            << "c20*eta*eta: " << c20*eta*eta << nl
            << "c02*xi*xi: " << c02*xi*xi << nl
            << endl;
    }
}


Foam::tmp<Foam::pointField> Foam::nagataPatch::projectPoints
(
    const pointField& origPoints,
    const labelList& triFaceIDs
)
{
    // Prepare the return field
    tmp<pointField> tprojectedPoints
    (
        new pointField(origPoints.size(), vector::zero)
    );
    pointField& projectedPoints = tprojectedPoints();

    // Lookup the face coefficients
    const List< List<vector> >& faceCoeffs = faceInterpCoeffs();

    // Local coordinates
    scalar eta = 0;
    scalar xi = 0;

    forAll(projectedPoints, pointI)
    {
        const point& origPoint = origPoints[pointI];
        const label faceID = triFaceIDs[pointI];

        if (debug > 1)
        {
            InfoIn(type() + "::projectPoints(...)")
                << "originalPoint: " << origPoint << endl;
        }

        // Calculate local coordinates within the face for the orignal point
        localCoordinates(eta, xi, origPoint, faceID);

        // Project point
        projectPoint(projectedPoints[pointI], faceCoeffs[faceID], eta, xi);
    }

    return tprojectedPoints;
}


Foam::triSurface Foam::nagataPatch::refineAndProjectPatches
(
    const wordList& patchesToProject,
    const int nRefinementLevels
)
{
    // Find the region indices of the patches to project
    labelList patchesToProjectIDs(patchesToProject.size(), -1);
    forAll(patchesToProjectIDs, pI)
    {
        forAll(triSurf_.patches(), patchI)
        {
            if (patchesToProject[pI] == triSurf_.patches()[patchI].name())
            {
                patchesToProjectIDs[pI] = patchI;
                break;
            }
        }

        if (patchesToProjectIDs[pI] == -1)
        {
            FatalError
                << "Patch " << patchesToProject[pI] << " not found!"
                << abort(FatalError);
        }
    }

    // Make faces for refinement and projection
    DynamicList<label> includedFaces(triSurf_.size());
    label nFaces = 0;
    forAll(patchesToProjectIDs, pI)
    {
        forAll(triSurf_, faceI)
        {
            if (triSurf_[faceI].region() == patchesToProjectIDs[pI])
            {
                // A face can only be in one region so it will not be added
                // twice
                includedFaces.append(faceI);
                nFaces++;
            }
        }
    }
    includedFaces.shrink();

    Info<< "Refining " << nFaces << " faces out of " << triSurf_.size()
        << endl;

    // Create refined surface
    // This surface contains new points added at the end of the list
    Info<< "Refining the surface: .";
    triSurface* refinedSurfPtr =
        new triSurface
        (
            triSurfaceTools::redGreenRefine
            (
                triSurf_,
                includedFaces
            )
        );

    for (int i = 1; i < nRefinementLevels; i++)
    {
        Info<< ".";

        // Update faces for refinement and projection list
        triSurface& refinedSurf = *refinedSurfPtr;
        includedFaces = DynamicList<label>(refinedSurf.size());
        forAll(patchesToProjectIDs, pI)
        {
            forAll(refinedSurf, faceI)
            {
                if (refinedSurf[faceI].region() == patchesToProjectIDs[pI])
                {
                    // A face can only be in one region so it will not be added
                    // twice
                    includedFaces.append(faceI);
                }
            }
        }
        includedFaces.shrink();

        triSurface* tmpRefinedSurfPtr = refinedSurfPtr;
        refinedSurfPtr =
            new triSurface
            (
                triSurfaceTools::redGreenRefine
                (
                    *tmpRefinedSurfPtr,
                    includedFaces
                )
            );

        delete tmpRefinedSurfPtr;
    }
    triSurface& refinedSurf = *refinedSurfPtr;
    Info<< endl;

    if (debug)
    {
        // Write refined surface
        fileName refinedSurfName("refined.stl");
        Info<< nl << "Writing " << refinedSurfName << endl;
        refinedSurf.write(refinedSurfName);
    }

    // Create a map for the new points

    if (debug)
    {
        Info<< nl << "Performing search to create point-to-face map" << endl;
    }

    triSurfaceSearch surfSearch(triSurf_);
    const vector span = vector(0.01*boundBox(triSurf_.points()).span());
    //const vector span = vector(boundBox(triSurf_.points()).span());
    pointField addedPoints
    (
        refinedSurf.nPoints() - triSurf_.nPoints(), vector::zero
    );
    forAll(addedPoints, pI)
    {
        addedPoints[pI] = refinedSurf.points()[triSurf_.nPoints() + pI];
    }

    const labelList addedPointsFaceMap =
        surfSearch.calcNearestTri(addedPoints, span);


    // Calculate new positions for the added points

    if (debug)
    {
        Info<< nl << "Projecting points" << endl;
    }

    const pointField newPoints = projectPoints(addedPoints, addedPointsFaceMap);

    // Fix edges shared with patches that were not projected. To do this, we
    // must first find points on patches that were not projected
    labelList patchesNotProjectedIDs
    (
        triSurf_.patches().size() - patchesToProjectIDs.size(), -1
    );

    label pnfID = 0;
    forAll(refinedSurf.patches(), patchI)
    {
        bool patchFound = false;
        forAll(patchesToProjectIDs, pI)
        {
            if (patchesToProjectIDs[pI] == patchI)
            {
                patchFound = true;
                break;
            }
        }

        if (!patchFound)
        {
            patchesNotProjectedIDs[pnfID++] = patchI;
        }
    }

    if (debug)
    {
        forAll(patchesNotProjectedIDs, pI)
        {
            Info<< "Not projecting patch: "
                << refinedSurf.patches()[patchesNotProjectedIDs[pI]].name()
                << endl;
        }
    }

    boolList pointsProjected(refinedSurf.nPoints(), true);
    forAll(patchesNotProjectedIDs, pI)
    {
        forAll(refinedSurf, faceI)
        {
            if (refinedSurf[faceI].region() == patchesNotProjectedIDs[pI])
            {
                // The points on this face should not be projected
                forAll(refinedSurf[faceI], pointI)
                {
                    pointsProjected[refinedSurf[faceI][pointI]] = false;
                }
            }
        }
    }

    pointField updatedPoints = refinedSurf.points();
    forAll(newPoints, pI)
    {
        const label pointID = triSurf_.nPoints() + pI;
        if (pointsProjected[pointID])
        {
            updatedPoints[pointID] = newPoints[pI];
        }
    }
    refinedSurf.movePoints(updatedPoints);

    if (debug)
    {
        forAll(refinedSurf.patches(), pI)
        {
            Info<< "refinedSurf.patches()[" << pI << "]: "
                << "name = " << refinedSurf.patches()[pI].name()
                << ", index = " << refinedSurf.patches()[pI].index() << endl;
        }
    }

    return refinedSurf;
}


void Foam::nagataPatch::removeEmptyPacthes(triSurface& triSurf) const
{
    // Make a copy of the patches
    geometricSurfacePatchList patchesCopy = triSurf.patches();

    // Mark all patches to keep
    boolList keepPatch(patchesCopy.size(), false);
    forAll(triSurf, faceI)
    {
        keepPatch[triSurf[faceI].region()] = true;
    }

    if (debug)
    {
        Info<< "Patches to keep:" << endl;
    }

    int nPatchesToKeep = 0;
    labelList patchMap(patchesCopy.size(), -1);
    forAll(keepPatch, patchI)
    {
        if (keepPatch[patchI])
        {
            patchMap[patchI] = nPatchesToKeep;

            nPatchesToKeep++;

            if (debug)
            {
                Info<< "    patch = " << patchI
                    << ", name = " << patchesCopy[patchI].name() << endl;
            }
        }
    }


    // Resize and set the patches
    triSurf.patches().clear();
    triSurf.patches().setSize(nPatchesToKeep);
    label newPatchID = 0;
    forAll(keepPatch, pI)
    {
        if (keepPatch[pI])
        {
            triSurf.patches()[newPatchID] =
                geometricSurfacePatch
                (
                    patchesCopy[pI].geometricType(),
                    patchesCopy[pI].name(),
                    newPatchID
                );

            newPatchID++;
        }
    }

    // Update face regions
    forAll(triSurf, faceI)
    {
        triSurf[faceI].region() = patchMap[triSurf[faceI].region()];
    }

    if (debug)
    {
        Info<< "New surface patches are: " << triSurf.patches() << endl;
    }
}


void Foam::nagataPatch::clearOut()
{
    faceInterpCoeffsPtr_.clear();
    faceCurvatureCoeffsPtr_.clear();
    faceOnEdgeOfPatchPtr_.clear();
}


// ************************************************************************* //
