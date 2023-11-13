#include "DxfFileToTriSurfConverter.H"

#include "cubicBSpline.H"
#include "DxfFile.H"
#include "triSurfModifier.H"

#include <cmath>
#include <iterator>

namespace Foam
{

// ============================================================================
// class DxfFileToTriSurfConverter
// ============================================================================

// ----------------------------------------------------------------------------
// Public member functions
// ----------------------------------------------------------------------------

DxfFileToTriSurfConverter::DxfFileToTriSurfConverter
(
    const DxfFile& file,
    const scalar&  angleTol,
    const scalar&  discretisationLength
)
    : file_(file)
    , angleTol_(angleTol)
    , discretisationLength_(discretisationLength)
{ }


void DxfFileToTriSurfConverter::Convert
(
    triSurf& surf
)
{
    auto& entities = file_.GetEntities();

    // Iterating through all entities and adding them to the surface
    for (auto it = entities.begin(); it != entities.end(); ++it)
    {
        AddEntityToTriSurf(it, surf);
    }

    ScaleToMeters(surf);
}


// ----------------------------------------------------------------------------
// Private member functions
// ----------------------------------------------------------------------------


void DxfFileToTriSurfConverter::AddArcToTriSurf
(
    const DxfEntityArc& arc,
    triSurf& surf
)
{
    triSurfModifier sMod(surf);

    auto& centre     = arc.GetCentre();
    auto& radius     = arc.GetRadius();
    auto  angleUnits = file_.GetAngleUnits();
    auto  startAngle = arc.GetStartAngle();
    auto  endAngle   = arc.GetEndAngle();
    auto  normal     = arc.GetNormal();

    // transformation vectors to the world coordinte system
    const vector Ax = arc.Ax(normal);
    const vector Ay = arc.Ay(normal);

    // centre in the world coordinate system
    const point Wc = centre.x() * Ax + centre.y() * Ay + centre.z() * normal;

    if (angleUnits == DxfFile::AngleUnits::Gradians)
    {
        startAngle *= 0.9;
        endAngle   *= 0.9;
    }
    else if (angleUnits == DxfFile::AngleUnits::Radians)
    {
        startAngle *= 180.0 / M_PI;
        endAngle   *= 180.0 / M_PI;
    }

    pointField& points = sMod.pointsAccess();

    if (startAngle > endAngle)
    {
        endAngle += 360.0;
    }

    label nDivisions = std::ceil((endAngle - startAngle) / angleTol_);
    const label nPoints = nDivisions + 1;

    label pointI = points.size();

    points.setSize(pointI + nPoints);

    //- calculate the step in radians
    const scalar startRad = startAngle * M_PI / 180.0;
    const scalar angleStep =
        ((endAngle - startAngle) / nDivisions) * M_PI / 180.0;

    //- add points into the surface
    label pId = sMod.surface().pointSubsetIndex("_corners_");
    if( pId < 0 )
        pId = sMod.surface().addPointSubset("_corners_");

    sMod.surface().addPointToSubset(pId, pointI);
    sMod.surface().addPointToSubset(pId, pointI+nPoints-1);

    for(label pI=0;pI<nPoints;++pI)
    {
        const scalar phi = startRad+pI*angleStep;
        const point p = Wc + radius * (Ax * std::cos(phi) + Ay * std::sin(phi));

        points[pointI+pI] = p;
    }

    //- add edges into the surface
    for(label eI=1;eI<nPoints;++eI)
    {
        sMod.featureEdgesAccess().append(edge(pointI+eI-1, pointI+eI));
    }
}


void DxfFileToTriSurfConverter::AddCircleToTriSurf
(
    const DxfEntityCircle& circle,
    triSurf& surf
)
{
    triSurfModifier sMod(surf);

    auto& centre = circle.GetCentre();
    auto& radius = circle.GetRadius();
    auto  normal = circle.GetNormal();

    // transformation vectors to the world coordinate system
    const vector Ax = circle.Ax(normal);
    const vector Ay = circle.Ay(normal);

    // centre in the world coordinate system
    const point Wc = centre.x() * Ax + centre.y() * Ay + centre.z() * normal;

    pointField& points = sMod.pointsAccess();

    const label nPoints = std::ceil(360.0 / angleTol_);

    const label pointI = points.size();

    points.setSize(pointI + nPoints);

    //- calculate the step in radians
    const scalar angleStep = (360.0 / nPoints) * M_PI / 180.0;

    //- add points into the surface
    for(label pI=0;pI<nPoints;++pI)
    {
        const scalar phi = angleStep * pI;
        const point p = Wc + radius * (Ax * std::cos(phi) + Ay * std::sin(phi));

        points[pointI+pI] = p;
    }

    //- add edges into the surface
    sMod.featureEdgesAccess().append(edge(points.size()-1, pointI));
    for(label eI=1;eI<nPoints;++eI)
    {
        sMod.featureEdgesAccess().append(edge(pointI+eI-1, pointI+eI));
    }
}


void DxfFileToTriSurfConverter::AddEllipseToTriSurf
(
    const DxfEntityEllipse& ellipse,
    triSurf& surf
)
{
    triSurfModifier sMod(surf);

    auto& centre     = ellipse.GetCentre();
    auto& endMajor   = ellipse.GetEndMajor();
    auto& normal     = ellipse.GetNormal();
    auto& axisRatio  = ellipse.GetAxisRatio();
    auto  angleUnits = file_.GetAngleUnits();
    auto  startAngle = ellipse.GetStartAngle();
    auto  endAngle   = ellipse.GetEndAngle();


    if (angleUnits == DxfFile::AngleUnits::Gradians)
    {
        startAngle *= 0.9;
        endAngle   *= 0.9;
    }
    else if (angleUnits == DxfFile::AngleUnits::Radians)
    {
        startAngle *= 180.0 / M_PI;
        endAngle   *= 180.0 / M_PI;
    }

    vector b = axisRatio * (normal ^ endMajor);

    pointField& points = sMod.pointsAccess();

    const scalar radTol = angleTol_ * M_PI / 180.0;
    label nDivisions = std::ceil((endAngle - startAngle) / radTol);
    const label nPoints = nDivisions + 1;

    label pointI = points.size();

    points.setSize(pointI + nPoints);

    //- calculate the step in radians
    if( startAngle > endAngle )
    {
        endAngle += 360.0;
    }

    const scalar angleStep = ((endAngle - startAngle) / nDivisions);

    //- add points into the surface
    label pId = sMod.surface().pointSubsetIndex("_corners_");
    if( pId < 0 )
        pId = sMod.surface().addPointSubset("_corners_");

    sMod.surface().addPointToSubset(pId, pointI);
    sMod.surface().addPointToSubset(pId, pointI+nPoints-1);

    for(label pI=0;pI<nPoints;++pI)
    {
        const scalar phi = startAngle + pI * angleStep;
        const point p = centre + endMajor * std::cos(phi) + b * std::sin(phi);

        points[pointI+pI] = p;
    }

    //- add edges into the surface
    for(label eI=1;eI<nPoints;++eI)
    {
        sMod.featureEdgesAccess().append(edge(pointI+eI-1, pointI+eI));
    }
}


void DxfFileToTriSurfConverter::AddInsertToTriSurf
(
    const DxfEntityInsert& insert,
    triSurf& surf
)
{
    auto& blocks = file_.GetBlocks();

    // Iteratating through all blocks
    for (auto it = blocks.begin(); it != blocks.end(); ++it)
    {
        const DxfBlock& block = **it;

        // If current block name is same as block name in INSERT entity
        if (block.GetName() == insert.GetBlockName())
        {
            // Adding block to edge subset
            auto blockIndex = surf.addEdgeSubset(block.GetName());

            // Number of edges in surface before adding entities to surface
            auto numOfEdgesBefore = surf.nFeatureEdges();

            auto& entities = block.GetEntities();

            // Iterating through all entities in the current block
            // and adding them to the surface
            for (auto it = entities.begin(); it != entities.end(); ++it)
            {
                AddEntityToTriSurf(it, surf);
            }

            // Adding new edges to edge subset
            for (auto i = numOfEdgesBefore; i != surf.nFeatureEdges(); ++i)
            {
                surf.addEdgeToSubset(blockIndex, i);
            }
        }
    }
}


void DxfFileToTriSurfConverter::AddLineToTriSurf
(
    const DxfEntityLine& line,
    triSurf& surf
)
{
    triSurfModifier sMod(surf);

    pointField&   points    = sMod.pointsAccess();
    edgeLongList& featEdges = sMod.featureEdgesAccess();

    const point& lineStartPoint = line.GetStart();
    const point& lineEndPoint   = line.GetEnd();

    const vector lineVector = lineEndPoint - lineStartPoint;
    const scalar lineLength = mag(lineVector);

    const label numOfMidpoints = lineLength / (discretisationLength_ + VSMALL);

    const label numOfSublines = numOfMidpoints + 1;

    const vector sublineVector = lineVector / numOfSublines;

    // ------------------------------------------------------------------------
    // Refining line
    // ------------------------------------------------------------------------

    std::vector<point> refinedLinePoints;

    // Inserting starting point
    refinedLinePoints.push_back(lineStartPoint);

    for (auto i = 0; i != numOfMidpoints; ++i)
    {
        // Inserting midpoints
        const point midpoint = refinedLinePoints.back() + sublineVector;
        refinedLinePoints.push_back(midpoint);
    }

    // Inserting ending point
    refinedLinePoints.push_back(lineEndPoint);

    // ------------------------------------------------------------------------
    // Inserting refined line into mesh
    // ------------------------------------------------------------------------

    const label oldNumOfPoints = points.size();
    const label newNumOfPoints = oldNumOfPoints + refinedLinePoints.size();

    points.setSize(newNumOfPoints);

    label pId = sMod.surface().pointSubsetIndex("_corners_");

    if (pId < 0)
    {
        pId = sMod.surface().addPointSubset("_corners_");
    }

    label pointI = oldNumOfPoints;

    // Inserting first point
    points[pointI] = refinedLinePoints.front();
    surf.addPointToSubset(pId, pointI);

    ++pointI;

    for
    (
        auto i = decltype(refinedLinePoints.size())(1);
        i != refinedLinePoints.size();
        ++i
    )
    {
        // Inserting other points inside the loop
        points[pointI] = refinedLinePoints[i];
        surf.addPointToSubset(pId, pointI);

        // Inserting edges
        featEdges.append(edge(pointI, pointI - 1));

        ++pointI;
    }
}


void DxfFileToTriSurfConverter::AddPolyLineToTriSurf
(
    const DxfEntityPolyLine& polyLine,
    triSurf& surf
)
{
    triSurfModifier sMod(surf);

    pointField&   points    = sMod.pointsAccess();
    edgeLongList& featEdges = sMod.featureEdgesAccess();

    const vector normal    = polyLine.GetNormal();
    const scalar elevation = polyLine.GetElevation();

    // Transformation vectors to the world coordinate system
    const vector Ax = polyLine.Ax(normal);
    const vector Ay = polyLine.Ay(normal);

    const point Wc = elevation * normal;

    const auto& xCoords = polyLine.GetVerticesX();
    const auto& yCoords = polyLine.GetVerticesY();

    // ------------------------------------------------------------------------
    // Calculating poly line points
    // ------------------------------------------------------------------------

    std::vector<point> polyLinePoints;

    for
    (
        auto i = decltype(polyLine.GetNumOfVertices())(0);
        i != polyLine.GetNumOfVertices();
        ++i
    )
    {
        polyLinePoints.push_back(Wc + xCoords[i] * Ax + yCoords[i] * Ay);
    }

    if (polyLinePoints.size() == 0)
    {
        return;
    }

    // ------------------------------------------------------------------------
    // Refining poly line
    // ------------------------------------------------------------------------

    std::vector<point> refinedPolyLinePoints;

    // Inserting starting point
    refinedPolyLinePoints.push_back(polyLinePoints.front());

    // Inserting other points inside this loop
    for
    (
        auto i = decltype(polyLinePoints.size())(1);
        i < polyLinePoints.size();
        ++i
    )
    {
        const point& segmentStartPoint = polyLinePoints[i - 1];
        const point& segmentEndPoint   = polyLinePoints[i];

        const vector segmentVector = segmentEndPoint - segmentStartPoint;
        const scalar segmentLength = mag(segmentVector);

        const label numOfMidpoints =
            segmentLength / (discretisationLength_ + VSMALL);

        const label numOfSubsegments = numOfMidpoints + 1;

        const vector subsegmentVector = segmentVector / numOfSubsegments;

        for (auto i = 0; i != numOfMidpoints; ++i)
        {
            // Inserting midpoints
            const point midpoint =
                refinedPolyLinePoints.back() + subsegmentVector;

            refinedPolyLinePoints.push_back(midpoint);
        }

        // Inserting ending point. This point is also starting point of the
        // following poly line segment.
        refinedPolyLinePoints.push_back(segmentEndPoint);
    }

    // ------------------------------------------------------------------------
    // Inserting refined poly line into mesh
    // ------------------------------------------------------------------------

    const label oldNumOfPoints = points.size();
    const label newNumOfPoints = oldNumOfPoints + refinedPolyLinePoints.size();

    points.setSize(newNumOfPoints);

    label pId = surf.pointSubsetIndex("_corners_");

    if (pId < 0)
    {
        pId = surf.addPointSubset("_corners_");
    }

    label pointI = oldNumOfPoints;

    // Inserting first point
    points[pointI] = refinedPolyLinePoints.front();
    surf.addPointToSubset(pId, pointI);

    ++pointI;

    for
    (
        auto i = decltype(refinedPolyLinePoints.size())(1);
        i != refinedPolyLinePoints.size();
        ++i
    )
    {
        // Inserting other points
        points[pointI] = refinedPolyLinePoints[i];
        surf.addPointToSubset(pId, pointI);

        // Inserting edges
        featEdges.append(edge(pointI, pointI - 1));

        ++pointI;
    }

    // If poly line is closed
    if (polyLine.GetFlag() == 1)
    {
        featEdges.append(edge(oldNumOfPoints, newNumOfPoints - 1));
    }
}


void DxfFileToTriSurfConverter::AddSplineToTriSurf
(
    const DxfEntitySpline& dxfEntitySpline,
    triSurf& surf
)
{
    const scalar tol = cos(angleTol_ * M_PI / 180.0);


    LongList<point> splinePoints(dxfEntitySpline.GetNumOfControlPoints());

    auto itX = dxfEntitySpline.GetControlPointsX().begin();
    auto itY = dxfEntitySpline.GetControlPointsY().begin();
    auto itZ = dxfEntitySpline.GetControlPointsZ().begin();

    for (auto i = 0; i != dxfEntitySpline.GetNumOfControlPoints(); ++i)
    {
        splinePoints[i] = point(*itX, *itY, *itZ);

        itX = std::next(itX);
        itY = std::next(itY);
        itZ = std::next(itZ);
    }



    cubicBSpline spline(splinePoints, "cubicBSpline");

    spline.setTangentAtStartingPoint(dxfEntitySpline.GetStartTangent());
    spline.setTangentAtEndPoint(dxfEntitySpline.GetEndTangent());

    label nDivisions = 100;
    scalar dt;
    bool finished;
    do
    {
        finished = true;
        dt = 1.0 / nDivisions;
        for (scalar t = 0.0; t < 1.0; )
        {
            if (t + 2.0 * dt > 1.0)
            {
                break;
            }

            const point p0 = spline.evaluate(t);
            const point p1 = spline.evaluate(t+dt);
            const point p2 = spline.evaluate(t+2.0*dt);

            vector v1 = p1 - p0;
            v1 /= (mag(v1) + VSMALL);
            vector v2 = p2 - p1;
            v2 /= (mag(v2) + VSMALL);

            if (mag(v2 & v1) < tol)
            {
                //- increase the number of subdivisions
                nDivisions *= 2;
                finished = false;
                break;
            }

            t += dt;
        }
    }
    while (!finished);



    triSurfModifier sMod(surf);

    pointField& points = sMod.pointsAccess();

    label pointI = points.size();

    points.setSize(pointI + nDivisions + 1);

    label pId = sMod.surface().pointSubsetIndex("_corners_");
    if (pId < 0)
    {
        pId = sMod.surface().addPointSubset("_corners_");
    }

    sMod.surface().addPointToSubset(pId, pointI);
    sMod.surface().addPointToSubset(pId, points.size() - 1);

    scalar t = 0.0;
    for (auto i = 0; i < nDivisions; ++i)
    {
        points[pointI + i] = spline.evaluate(t);
        t += dt;
    }
    points[pointI + nDivisions] = spline.evaluate(1.0);

    edgeLongList& edges = sMod.featureEdgesAccess();

    if (dxfEntitySpline.GetFlag() == 1)
    {
        edges.append(edge(pointI + nDivisions, pointI));
    }

    for (auto i = 1; i < nDivisions; ++i)
    {
        edges.append(edge(pointI + i - 1, pointI + i));
    }
}


void DxfFileToTriSurfConverter::ScaleToMeters(triSurf& surf)
{
    auto scalingFactor = 1.0;

    switch (file_.GetUnits())
    {
        case DxfFile::Units::Unitless:     scalingFactor = 1.0;      break;
        case DxfFile::Units::Inches:       scalingFactor = 0.0254;   break;
        case DxfFile::Units::Feet:         scalingFactor = 0.3048;   break;
        case DxfFile::Units::Miles:        scalingFactor = 1609.34;  break;
        case DxfFile::Units::Millimeters:  scalingFactor = 0.001;    break;
        case DxfFile::Units::Centimeters:  scalingFactor = 0.01;     break;
        case DxfFile::Units::Meters:       scalingFactor = 1.0;      break;
        case DxfFile::Units::Kilometers:   scalingFactor = 1000.0;   break;
        case DxfFile::Units::Microinches:  scalingFactor = 2.54e-8;  break;
        case DxfFile::Units::Mils:         scalingFactor = 2.54e-5;  break;
        case DxfFile::Units::Yards:        scalingFactor = 0.9144;   break;
        case DxfFile::Units::Angstroms:    scalingFactor = 1e-10;    break;
        case DxfFile::Units::Nanometers:   scalingFactor = 1e-9;     break;
        case DxfFile::Units::Microns:      scalingFactor = 1e-6;     break;
        case DxfFile::Units::Decimeters:   scalingFactor = 0.1;      break;
        case DxfFile::Units::Decameters:   scalingFactor = 10;       break;
        case DxfFile::Units::Hectometers:  scalingFactor = 100;      break;
        case DxfFile::Units::Gigameters:   scalingFactor = 1e9;      break;

        case DxfFile::Units::AstronomicalUnits:
            scalingFactor = 149597870700;
            break;

        case DxfFile::Units::LightYears:
            scalingFactor = 9.4607e15;
            break;

        case DxfFile::Units::Parsecs:
            scalingFactor = 3.0857e16;
            break;

        default:
            break;
    }

    triSurfModifier sMod(surf);
    pointField& points = sMod.pointsAccess();
    forAll(points, pI)
        points[pI] *= scalingFactor;
}

} // End namespace Foam
