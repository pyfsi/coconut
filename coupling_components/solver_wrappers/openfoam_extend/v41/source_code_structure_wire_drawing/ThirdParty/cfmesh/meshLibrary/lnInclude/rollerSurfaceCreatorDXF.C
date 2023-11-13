/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "rollerSurfaceCreatorDXF.H"
#include "demandDrivenData.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "triSurfaceCopyParts.H"
#include "DxfFile.H"
#include "DxfFileParser.H"
#include "DxfFileToTriSurfConverter.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "triSurfaceCleanupDuplicates.H"

#include "addToRunTimeSelectionTable.H"

#include <map>

//#define DEBUGSurfaceCreator

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rollerSurfaceCreatorDXF, 0);
addToRunTimeSelectionTable
(
    rollerSurfaceCreator,
    rollerSurfaceCreatorDXF,
    dictionary
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rollerSurfaceCreatorDXF::parseGeometryFromDXF()
{
    //- generate triangulated surface
    surfPtr_ = new triSurf();

    //
    // Opening input file
    //

    std::ifstream ifs(fName_, std::ios::binary);
    if (!ifs.is_open())
    {
        Info << "ERROR: Cannot open file \"" << fName_ << "\"" << endl;
        std::abort();
    }

    //
    // Parsing input file
    //

    DxfFile file;

    try
    {
        DxfFileParser parser;
        parser.Parse(ifs, file);
    }
    catch (const std::exception& e)
    {
        Info << "ERROR: Error while parsing DXF file \"" << fName_ << "\": "
             << e.what() << endl;
        std::abort();
    }

    if (file.GetUnits() != DxfFile::Units::Millimeters)
    {
        Warning << "THE UNITS IN THE FILE \"" << fName_ << "\" ARE NOT mm."
                << " CONTINUING WITH THE ASSUMPTION THEY ARE PROVIDED IN mm!!"
                << endl;

        file.SetUnits(DxfFile::Units::Millimeters);
    }

    //
    // Converting DXF file to triSurf
    //

    try
    {
        DxfFileToTriSurfConverter converter(file, 0.1);
        converter.Convert(*surfPtr_);
    }
    catch (const std::exception& e)
    {
        Info << "ERROR: Error while converting DXF file "
             << "\"" << fName_ << "\" to triSurf: " << e.what() << endl;
        std::abort();
    }

    Switch useContact(false);
    if( dict_.parent().found("contact") )
    {
        dict_.parent().readIfPresent("contact", useContact);

        if( useContact )
        {
            const label sId = surfPtr_->edgeSubsetIndex(rollPosition_);

            if( sId < 0 )
            {
                FatalError << "dxfFile " << fName_ << " does not contain "
                    << "the block " << rollPosition_
                    << ". Please create a block for each roll and start again."
                    << exit(FatalError);
            }

            //- the dxf contains blocks. Use the parts belonging to this roll
            triSurfaceCopyParts copyParts(*surfPtr_);
            wordList parts(1);
            parts[0] = rollPosition_;

            triSurf* copyPtr = copyParts.copySurface(parts);
            delete surfPtr_;
            surfPtr_ = copyPtr;
        }
    }

    if( !useContact )
    {
        if( file.GetBlocks().size() != 0 )
        {
            if( surfPtr_->edgeSubsetIndex(rollPosition_) >= 0 )
            {
                //- the dxf contains blocks.
                //- Use the parts belonging to this roll
                Warning << "dxf file contains block " << rollPosition_
                     << ". Assuming this block contains the geometry"
                     << " of the roll" << endl;

                triSurfaceCopyParts copyParts(*surfPtr_);
                wordList parts(1);
                parts[0] = rollPosition_;

                triSurf* copyPtr = copyParts.copySurface(parts);
                delete surfPtr_;
                surfPtr_ = copyPtr;
            }
            else
            {
                Warning << "dxf file " << fName_
                     << " does not contain the block "
                     << rollPosition_ << ". Assuming all entities are part"
                     << " of the roll at position " << rollPosition_ << endl;
            }
        }
    }

    //- merge duplicate vertices in the dxf file
    mergeDuplicates();

    # ifdef DEBUGSurfaceCreator
    triSurfaceExtrude2DEdges ext(*surfPtr_);
    const triSurf* extrudedPtr = ext.extrudeSurface();
    extrudedPtr->writeSurface("extrudedSurf.stl");
    deleteDemandDrivenData(extrudedPtr);
    # endif

    triSurfModifier sMod(*surfPtr_);

    pointField& points = sMod.pointsAccess();

    if( points.size() == 0 )
    {
        FatalErrorIn
        (
            "void rollerSurfaceCreatorDXF::parseGeometryFromDXF()"
        ) << "Surface " << fName_ << " does not contain any points!!!"
          << exit(FatalError);
    }

    //- check if the surface is in x-y space
    forAll(points, pointI)
    {
        if( mag(points[pointI].z()) > SMALL )
        {
            FatalErrorIn
            (
                "void rollerSurfaceCreatorDXF::parseGeometryFromDXF()"
            ) << "Surface " << fName_ << " is not in x-y space"
              << exit(FatalError);
        }
    }

    //- transform the surface from the x-y coordinates, as defined by Bekaert,
    //- into the lower half of the x-y space used for defining
    //- the profile from points
    tensor transformation(tensor::zero);

    if( rollPosition_ == "topRoll" )
    {
        transformation.xx() = 1.0;
        transformation.yy() = -1.0;
    }
    else if( rollPosition_ == "bottomRoll" )
    {
        transformation.xx() = 1.0;
        transformation.yy() = 1.0;
    }
    else if( rollPosition_ == "leftRoll" )
    {
        //- flip the y coordinate
        forAll(points, pI)
            points[pI].y() *= -1.0;

        transformation.xy() = -1.0;
        transformation.yx() = 1.0;
    }
    else if( rollPosition_ == "rightRoll" )
    {
        transformation.xy() = 1.0;
        transformation.yx() = -1.0;
    }
    else if( rollPosition_ == "bottomLeftRoll" )
    {
        transformation.xx() = transformation.yy() = 0.5;
        transformation.yx() = sqrt(3.0) / 2.0;
        transformation.xy() = -1.0 * transformation.yx();
    }
    else if( rollPosition_ == "bottomRightRoll" )
    {
        transformation.xx() = transformation.yy() = 0.5;
        transformation.yx() = -sqrt(3.0) / 2.0;
        transformation.xy() = -1.0 * transformation.yx();
    }
    else
    {
        FatalErrorIn
        (
            "void rollerSurfaceCreatorDXF::parseGeometryFromDXF()"
        ) << "Roll position " << rollPosition_ << " is not a valid one"
          << exit(FatalError);
    }

    //- transform the points into default position for meshing
    //- lower half of the x-y plane
    forAll(points, pointI)
    {
        points[pointI] = transformation & points[pointI];
    }

    # ifdef DEBUGSurfaceCreator
    extrudedPtr = ext.extrudeSurface();
    extrudedPtr->writeSurface("extrudedAfterTransformation.stl");
    deleteDemandDrivenData(extrudedPtr);
    Info << "Number of edges after transformation "
         << surfPtr_->featureEdges().size() << endl;
    # endif

    //- find min and max x coordinates, and max y coordinate
    scalar maxY(-VGREAT);
    std::map<scalar, std::map<scalar, std::set<label> > > xyToIndex;

    forAll(points, pointI)
    {
        const scalar currX = points[pointI].x();
        const scalar currY = points[pointI].y();

        xyToIndex[currX][currY].insert(pointI);

        maxY = max(maxY, currY);
    }

    //- maxY respresent the radial shift of the roll from the origin
    rollRadialShift_ = -maxY;

    //- transform the profile such that maxY becomes 0.0
    forAll(points, pointI)
        points[pointI].y() -= maxY;
    maxY = 0.0;

    scalar minX = xyToIndex.begin()->first;
    scalar maxX = xyToIndex.rbegin()->first;

    //- find edges at min and max x coordinates
    DynList<label> edgesAtMinX, edgesAtMaxX;

    const edgeLongList& surfEdges = surfPtr_->featureEdges();

    forAll(surfEdges, seI)
    {
        const edge& se = surfEdges[seI];
        const point& s = points[se.start()];
        const point& e = points[se.end()];

        if( mag(s.x() - minX) < SMALL && mag(e.x() - minX) < SMALL )
        {
            edgesAtMinX.append(seI);
        }

        if( mag(s.x() - maxX) < SMALL && mag(e.x() - maxX) < SMALL )
        {
            edgesAtMaxX.append(seI);
        }
    }

    //- set the roll width to the requested value
    rollAxialShift_ = 0.5 * (minX + maxX);

    if
    (
        (minX < (rollAxialShift_ - 0.5 * rollWidth_ - SMALL)) ||
        (maxX > (rollAxialShift_ + 0.5 * rollWidth_ + SMALL))
    )
    {
        FatalError << "The specified roll profile for " << rollPosition_
            << " is smaller than the requested roll width!"
            << "Width of the specified profile " << (maxX-minX) << nl
            << "Requested width " << rollWidth_ << nl
            << "Axial shift from zero " << rollAxialShift_ << nl
            << "Please check your settings!" << exit(FatalError);
    }
    else
    {
        if( useContact && (rollWidth_ > (maxX - minX + SMALL)) )
        {
            FatalError << "It is not allowed to set the rollWidth"
                << " larger than the width of the drawing when meshing"
                << " a case with contact. This often generated overlapping"
                << " meshes!"
                << " Desired roll width for the roll " << rollPosition_
                << " is " << (maxX - minX)
                << exit(FatalError);
        }

        //- setting correct roll width
        minX = rollAxialShift_ - 0.5 * rollWidth_;
        maxX = rollAxialShift_ + 0.5 * rollWidth_;
    }

    //- move the points with min y coordinate onto the correct position
    const scalar minY = maxY - 0.5 * (outerDiameter_ - innerDiameter_);

    label start(-1), end(-1);

    scalar maxYCoord(-VGREAT);
    if( edgesAtMinX.size() == 0 )
    {
        //- the edge does not exist, create additional vertex and the edge
        label minXPointI = *xyToIndex.begin()->second.rbegin()->second.begin();

        points[minXPointI].x() = minX;

        if( xyToIndex.begin()->second.size() > 1 )
        {
            Warning << "Found more than one point with min x coordinate"
                    << endl;
        }

        //- find the maximum y coordinate on the line at minX
        maxYCoord = points[minXPointI].y();

        //- create a missing edge
        edgesAtMinX.append(surfPtr_->nFeatureEdges());

        start = points.size();
        surfPtr_->appendFeatureEdge(edge(minXPointI, start));

        surfPtr_->appendVertex(point(minX, minY, 0));
    }
    else
    {
        //- find the node with min y coordinate and set it to minY
        scalar minYCoord(VGREAT);
        label minYPointI(-1);
        forAll(edgesAtMinX, i)
        {
            const edge& e = surfEdges[edgesAtMinX[i]];

            //- find the maximum y coordinate on the line at minX
            maxYCoord = max(maxYCoord, points[e.start()].y());
            maxYCoord = max(maxYCoord, points[e.end()].y());

            if( points[e.start()].y() < minYCoord )
            {
                points[e.start()].x() = minX;
                minYCoord = points[e.start()].y();
                minYPointI = e.start();
            }

            if( points[e.end()].y() < minYCoord )
            {
                points[e.end()].x() = minX;
                minYCoord = points[e.end()].y();
                minYPointI = e.end();
            }
        }

        points[minYPointI].y() = minY;
        start = minYPointI;
    }

    //- check the properties of the roll: the inner diameter shall not be
    //- greater than the min diameter of the contact part of the roll
    if( maxYCoord < minY )
    {
        FatalError
            << "The settings for the roll at position " << rollPosition_
            << " are not valid!" << nl
            << "Requested inner diameter " << innerDiameter_ << nl
            << "Max allowed inner diameter "
            << (outerDiameter_-2.0 * mag(maxYCoord)) << nl
            << "Please check your settings!" << exit(FatalError);
    }

    maxYCoord = -VGREAT;
    if( edgesAtMaxX.size() == 0 )
    {
        //- the edge does not exist, create vertex and edge
        label maxXPointI = *xyToIndex.rbegin()->second.rbegin()->second.begin();
        points[maxXPointI].x() = maxX;

        if( xyToIndex.rbegin()->second.size() != 1 )
        {
            Warning << "Found more than one point with max x coordinate"
                    << endl;
        }

        //- find the maximum y coordinate on the line at maxX
        maxYCoord = points[maxXPointI].y();

        //- create a missing edge
        edgesAtMaxX.append(surfPtr_->nFeatureEdges());

        end = points.size();
        surfPtr_->appendFeatureEdge(edge(maxXPointI, end));

        //- create new vertex
        surfPtr_->appendVertex(point(maxX, minY, 0));
    }
    else
    {
        //- find the node with min y coordinate and set it to minY
        scalar minYCoord(VGREAT);
        label minYPointI(-1);
        forAll(edgesAtMaxX, i)
        {
            const edge& e = surfEdges[edgesAtMaxX[i]];

            //- find the maximum y coordinate on the line at maxX
            maxYCoord = max(maxYCoord, points[e.start()].y());
            maxYCoord = max(maxYCoord, points[e.end()].y());


            if( points[e.start()].y() < minYCoord )
            {
                points[e.start()].x() = maxX;
                minYCoord = points[e.start()].y();
                minYPointI = e.start();
            }

            if( points[e.end()].y() < minYCoord )
            {
                points[e.start()].x() = maxX;
                minYCoord = points[e.end()].y();
                minYPointI = e.end();
            }
        }

        points[minYPointI].y() = minY;
        end = minYPointI;
    }

    //- check the properties of the roll: the inner diameter shall not be
    //- greater than the min diameter of the contact part of the roll
    if( maxYCoord < minY )
    {
        FatalError
            << "The settings for the roll at position " << rollPosition_
            << " are not valid!" << nl
            << "Requested inner diameter " << innerDiameter_ << nl
            << "Max allowed inner diameter "
            << (outerDiameter_-2.0 * mag(maxYCoord)) << nl
            << "Please check your settings!" << exit(FatalError);
    }

    //- refine edges at min and max x coordinates
    if( useContact )
    {
        edgeLongList& fEdges = sMod.featureEdgesAccess();

        while( edgesAtMinX.size() < 10 )
        {
            //- recursively subdivide edges
            forAllReverse(edgesAtMinX, i)
            {
                edge& fe = fEdges[edgesAtMinX[i]];

                const point c = fe.centre(points);

                //- append a new edge
                edgesAtMinX.append(fEdges.size());
                fEdges.append(edge(points.size(), fe.end()));

                //- modify the existing edge
                fe.end() = points.size();

                //- append new vertex
                surfPtr_->appendVertex(c);
            }
        }

        while( edgesAtMaxX.size() < 10 )
        {
            //- recursively subdivide edges
            forAllReverse(edgesAtMaxX, i)
            {
                edge& fe = fEdges[edgesAtMaxX[i]];

                const point c = fe.centre(points);

                //- append a new edge
                edgesAtMaxX.append(fEdges.size());
                fEdges.append(edge(points.size(), fe.end()));

                //- modify the existing edge
                fe.end() = points.size();

                //- append new vertex
                surfPtr_->appendVertex(c);
            }
        }

        Info << "Subdivision completed" << endl;
    }

    # ifdef DEBUGSurfaceCreator
    {
        triSurfaceExtrude2DEdges ext(*surfPtr_);
        const triSurf* extrudedPtr = ext.extrudeSurface();
        extrudedPtr->writeSurface("extrudedAfterMinMax.stl");
        deleteDemandDrivenData(extrudedPtr);
        Info << "Number of edges " << surfEdges.size() << endl;
        std::exit(0);
    }
    # endif

    //- create subsets that shall become patches
    const label contactId =
        surfPtr_->addEdgeSubset
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERCONTACT
            )
        );

    label rollerToAirId(-1);
    scalar contactMinX = -VGREAT;
    scalar contactMaxX = VGREAT;

    if( contactWidth_ > SMALL )
    {
        rollerToAirId =
            surfPtr_->addEdgeSubset
            (
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERTOAIR
                )
            );

        contactMinX = rollAxialShift_ - 0.5 * contactWidth_;
        contactMaxX = rollAxialShift_ + 0.5 * contactWidth_;

        //- cut the geometry to have get the exact contact area length
        cutContactArea(*surfPtr_);
    }

    const label backId =
        surfPtr_->addEdgeSubset
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERBACK
            )
        );
    const label frontId =
        surfPtr_->addEdgeSubset
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERFRONT
            )
        );
    const label axisId =
        surfPtr_->addEdgeSubset
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERAXIS
            )
        );

    scalar minXContact(VGREAT), maxXContact(-VGREAT);
    scalar xMinPointI(-1), xMaxPointI(-1);

    forAll(surfEdges, eI)
    {
        const edge& e = surfEdges[eI];

        bool foundEdge(false);

        if( edgesAtMinX.contains(eI) )
        {
            //- min and max contact points
            if( points[e.start()].x() < minXContact )
            {
                xMinPointI = e.start();
                minXContact = points[e.start()].x();
            }

            if( points[e.start()].x() > maxXContact )
            {
                xMaxPointI = e.start();
                maxXContact = points[e.start()].x();
            }

            if( points[e.end()].x() < minXContact )
            {
                xMinPointI = e.end();
                minXContact = points[e.end()].x();
            }

            if( points[e.end()].x() > maxXContact )
            {
                xMaxPointI = e.end();
                maxXContact = points[e.end()].x();
            }

            surfPtr_->addEdgeToSubset(backId, eI);
            foundEdge = true;
        }

        if( edgesAtMaxX.contains(eI) )
        {
            //- min and max contact points
            if( points[e.start()].x() < minXContact )
            {
                xMinPointI = e.start();
                minXContact = points[e.start()].x();
            }

            if( points[e.start()].x() > maxXContact )
            {
                xMaxPointI = e.start();
                maxXContact = points[e.start()].x();
            }

            if( points[e.end()].x() < minXContact )
            {
                xMinPointI = e.end();
                minXContact = points[e.end()].x();
            }

            if( points[e.end()].x() > maxXContact )
            {
                xMaxPointI = e.end();
                maxXContact = points[e.end()].x();
            }

            surfPtr_->addEdgeToSubset(frontId, eI);
            foundEdge = true;
        }

        if( !foundEdge )
        {
            //- min and max contact points
            if( points[e.start()].x() < minXContact )
            {
                xMinPointI = e.start();
                minXContact = points[e.start()].x();
            }

            if( points[e.start()].x() > maxXContact )
            {
                xMaxPointI = e.start();
                maxXContact = points[e.start()].x();
            }

            if( points[e.end()].x() < minXContact )
            {
                xMinPointI = e.end();
                minXContact = points[e.end()].x();
            }

            if( points[e.end()].x() > maxXContact )
            {
                xMaxPointI = e.end();
                maxXContact = points[e.end()].x();
            }

            if
            (
                (rollerToAirId >= 0) &&
                (
                    (
                        (points[e.end()].x() < contactMinX) ||
                        (points[e.start()].x() < contactMinX)
                    ) ||
                    (
                        (points[e.end()].x() > contactMaxX) ||
                        (points[e.start()].x() > contactMaxX)
                    )
                )
            )
            {
                surfPtr_->addEdgeToSubset(rollerToAirId, eI);
            }
            else
            {
                surfPtr_->addEdgeToSubset(contactId, eI);
            }
        }
    }

    //- move the points in the contact patch to match the rest of the mesh
    points[xMinPointI].x() = minX;
    points[xMaxPointI].x() = maxX;
    forAll(edgesAtMinX, i)
    {
        const edge& e = surfEdges[edgesAtMinX[i]];
        points[e.start()].x() = minX;
        points[e.end()].x() = minX;
    }
    forAll(edgesAtMaxX, i)
    {
        const edge& e = surfEdges[edgesAtMaxX[i]];
        points[e.start()].x() = maxX;
        points[e.end()].x() = maxX;
    }

    //- create an edge between the last two points
    surfPtr_->addEdgeToSubset(axisId, surfPtr_->featureEdges().size());
    surfPtr_->appendFeatureEdge(edge(start, end));

    //- mark corner nodes that shall be extruded to feature edges
    std::map<label, std::set<label> > pointSubsets;
    DynList<label> eSubsets;
    surfPtr_->edgeSubsetIndices(eSubsets);
    forAll(eSubsets, i)
    {
        const label sId = eSubsets[i];

        labelLongList sEdges;
        surfPtr_->edgesInSubset(sId, sEdges);

        forAll(sEdges, j)
        {
            const edge& fe = surfPtr_->featureEdges()[sEdges[j]];

            pointSubsets[fe.start()].insert(sId);
            pointSubsets[fe.end()].insert(sId);
        }
    }

    const label cornerId = surfPtr_->addPointSubset("_corners_");
    for
    (
        std::map<label, std::set<label> >::const_iterator it =
            pointSubsets.begin();
        it!=pointSubsets.end();
        ++it
    )
    {
        if( it->second.size() > 1 )
        {
            surfPtr_->addPointToSubset(cornerId, it->first);
        }
    }

    # ifdef DEBUGSurfaceCreator
    //- is every edge in a subset
    labelList edgeInSubset(surfEdges.size(), -1);

    DynList<label> esIds;
    surfPtr_->edgeSubsetIndices(esIds);
    forAll(esIds, i)
    {
        labelLongList edgesInSubset;
        surfPtr_->edgesInSubset(esIds[i], edgesInSubset);

        forAll(edgesInSubset, j)
            edgeInSubset[edgesInSubset[j]] = esIds[i];
    }

    forAll(edgeInSubset, i)
    {
        if( edgeInSubset[i] < 0 )
        {
            FatalError << "Edge " << i << " is not assigned to any subset"
                       << abort(FatalError);
        }
    }

    extrudedPtr = ext.extrudeSurface();
    extrudedPtr->writeSurface("surfAfterSubsets.stl");
    deleteDemandDrivenData(extrudedPtr);

    Info << "Number of edges after subsets " << surfEdges.size() << endl;
    # endif

    //- extrude lines into a 2D ribbon
    triSurfaceExtrude2DEdges extruder(*surfPtr_);
    triSurf* extrudedSurfPtr = new triSurf();
    extruder.extrudeSurface(*extrudedSurfPtr);

    deleteDemandDrivenData(surfPtr_);
    surfPtr_ = extrudedSurfPtr;

    //- the first and the last point differ in the z coordinate
    transformationMatrix_.yy() = 1.0;
    transformationMatrix_.xz() = 1.0;
    transformationMatrix_.zx() = -1.0;

    const tensor invTransform = inv(transformationMatrix_);
    pointField& pts = triSurfModifier(*surfPtr_).pointsAccess();

    forAll(pts, pI)
    {
        point& p = pts[pI];

        p = invTransform & p;

        p.y() += outerDiameter_ / 2.0;
    }

    # ifdef DEBUGSurfaceCreator
    surfPtr_->writeSurface("surf.fms");
    surfPtr_->writeSurface("surf.stl");
    Info << "Surface generation finished " << endl;
    ::exit(0);
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollerSurfaceCreatorDXF::rollerSurfaceCreatorDXF
(
    const word rollPosition,
    const direction symm,
    const rollingMillPatchNamesHandler& patchHandler,
    const dictionary& dict,
    const scalar tol
)
:
    rollerSurfaceCreator
    (
        rollPosition,
        symm,
        patchHandler,
        dict,
        tol
    )
{}

rollerSurfaceCreatorDXF::rollerSurfaceCreatorDXF
(
    const rollerSurfaceCreatorDXF& creator
)
:
    rollerSurfaceCreator(creator)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rollerSurfaceCreatorDXF::generateGeometry()
{
    if( !dict_.readIfPresent("dxfFile", fName_) )
    {
        FatalErrorIn
        (
            "void rollerSurfaceCreatorDXF::generateGeometry()"
        ) << "dxfFile is not present in dictionary "
          << dict_ << exit(FatalError);
    }

    if( fName_.ext() == "dxf" || fName_.ext() == "DXF" )
    {
        //- geometry is given as a dxf file
        if( rollPosition_.empty() )
        {
            FatalErrorIn
            (
                "void rollerSurfaceCreatorDXF::generateGeometry()"
            ) << "Invalid roll position " << rollPosition_ << exit(FatalError);
        }

        parseGeometryFromDXF();

        detectIndependentRegions();

        positionProfileInAxialDirection();
    }
    else
    {
        FatalErrorIn
        (
            "void rollerSurfaceCreatorDXF::generateGeometry()"
        ) << "Geometry is not given in a dxf file " << fName_
          << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
