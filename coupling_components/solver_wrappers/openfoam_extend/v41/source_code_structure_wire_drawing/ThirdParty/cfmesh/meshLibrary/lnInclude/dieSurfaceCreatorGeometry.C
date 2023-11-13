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

#include "dieSurfaceCreatorGeometry.H"
#include "rollingMillMesh.H"
#include "demandDrivenData.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "polyMeshGenModifier.H"
#include "splineBase.H"
#include "cubicBSpline.H"
#include "IFstream.H"
#include "boundBox.H"
#include "addToRunTimeSelectionTable.H"

//#define DEBUGDie

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(dieSurfaceCreatorGeometry, 0);
addToRunTimeSelectionTable
(
    dieSurfaceCreator,
    dieSurfaceCreatorGeometry,
    dictionary
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dieSurfaceCreatorGeometry::parseGeometryFile()
{

    //- exit in case the geometry is an empty string
    if( fName_ == "" )
        FatalErrorIn
        (
            "void dieSurfaceCreatorGeometry::"
            "parseGeometryFile()"
        ) << "Input geometry is not given!" << exit(FatalError);

    IFstream file(fName_);

    while( !file.eof() )
    {
        //- read the next token
        token t;
        file.read(t);

        if( t == token::BEGIN_LIST )
        {
            file.putBack(t);

            //- read a characteristic point
            point p;
            file >> p;

            // PDJ20/02/2019: Points in list are assumed in m instead of mm
            //characteristicPoints_.append(0.001 * p);
            characteristicPoints_.append(p);
        }
        else if( t.isWord() )
        {
            const word w = t.wordToken();

            if( w == "mode" )
            {
                word mode;
                file >> mode;
                
                interpolationType_ = mode;
            }
            else if( w == "dirFront" )
            {
                file >> dirFront_;
            }
            else if( w == "dirBack" )
            {
                file >> dirBack_;
            }
            else
            {
                Warning << "Unknown keyword " << w << " found!" << endl;
            }
        }
        else if( file.good() )
        {
            Warning << "Unknown token " << t.info() << " found!" << endl;
        }
    }

    // PDJ20/02/2019: interpolationType_ defaulted to b-spline.
//    if( !dict_.readIfPresent("mode", interpolationType_) )
//    {
//        FatalError << "mode is not defined in dictionary" << dict_
//             << exit(FatalError);
//    }

    if( dict_.readIfPresent("dirFront", dirFront_) )
    {
        Info << "Gradient on the front side " << dirFront_ << endl;
    }

    if( dict_.readIfPresent("dirBack", dirBack_) )
    {
        Info << "Gradient on the back side " << dirBack_ << endl;
    }

    Info << "Profile calculation mode " << interpolationType_ << endl;
    Info << "Characteristic points " << characteristicPoints_ << endl;
}

autoPtr<splineBase>
dieSurfaceCreatorGeometry::createProfile() const
{
    forAll(characteristicPoints_, pI)
    {
        if( mag(characteristicPoints_[pI].z()) > SMALL )
        {
            FatalErrorIn
            (
                "autoPtr<splineBase> dieSurfaceCreatorGeometry"
                "::createProfile() const"
            ) << "Profile is not in the x-y plane" << exit(FatalError);
        }

        if( characteristicPoints_[pI].y() < -SMALL )
        {
            FatalErrorIn
            (
                "autoPtr<splineBase> dieSurfaceCreatorGeometry"
                "::createProfile() const"
            ) << "Point " << characteristicPoints_[pI]
              << " has negative y value" << exit(FatalError);
        }
    }

    //- create the profile
    autoPtr<splineBase> splinePtr;

    if( interpolationType_ == "Bezier" )
    {
        splinePtr = splineBase::New(characteristicPoints_, "bezierSpline");
    }
    else if( interpolationType_ == "b-spline" )
    {
        splinePtr = splineBase::New(characteristicPoints_, "cubicBSpline");
    }
    else if( interpolationType_ == "b-splineWithTangents")
    {
        cubicBSpline* bSplinePtr =
            new cubicBSpline(characteristicPoints_, "cubicBSpline");

        bSplinePtr->setTangentAtStartingPoint(dirFront_);
        bSplinePtr->setTangentAtEndPoint(dirBack_);

        splinePtr = bSplinePtr->clone();
        deleteDemandDrivenData(bSplinePtr);
    }
    else
    {
        splinePtr = splineBase::New(characteristicPoints_, interpolationType_);
    }

    return splinePtr;
}

void dieSurfaceCreatorGeometry::createTriangulatedSurface()
{
    //- tranform the points to the y-z plane
    autoPtr<splineBase> splinePtr = createProfile();

    //- create the spline interpolating a given set of points
    LongList<point> curvePoints;
    splinePtr->createPolyLine(geometryTol_, curvePoints);

    if( curvePoints.size() < 2 )
    {
        FatalErrorIn
        (
            "void dieSurfaceCreatorGeometry"
            "::createTriangulatedSurface()"
        ) << "The number of generates vertices is smaller than 2"
          << exit(FatalError);
    }

    if( curvePoints[curvePoints.size()-1].x() > curvePoints[0].x() )
    {
        profilePoints_.setSize(curvePoints.size());

        forAll(curvePoints, pI)
            profilePoints_[pI] = curvePoints[pI];
    }
    else
    {
        profilePoints_.setSize(curvePoints.size());

        label counter(0);
        forAllReverse(curvePoints, pI)
            profilePoints_[counter++] = curvePoints[pI];
    }

    //- create surface mesh
    triSurf surf;
    triSurfModifier sMod(surf);

    //- create surface points
    pointField& pts = sMod.pointsAccess();
    pts.transfer(profilePoints_);

    //- create feature edges
    for(label i=1;i<pts.size();++i)
        surf.appendFeatureEdge(edge(i-1, i));

    scalar minX(VGREAT), maxX(-VGREAT);
    label minId(-1), maxId(-1);

    forAll(surf.points(), pI)
    {
        const point& p = surf.points()[pI];

        if( mag(p.z()) > SMALL )
        {
            FatalErrorIn
            (
                "void dieSurfaceCreatorDXF::parseDXFFile()"
            ) << "Geometry is not in the x-y plane!" << exit(FatalError);
        }

        if( p.x() < minX )
        {
            minX = p.x();
            minId = pI;
            inletDiameter_ = 2.0 * p.y();
        }

        if( p.x() > maxX )
        {
            maxX = p.x();
            maxId = pI;

            outletDiameter_ = 2.0 * p.y();
        }
    }

//    if( !dict_.readIfPresent("outerDiameter", outerDiameter_) )
//    {
//        FatalError << "outerDiameter is not present in dictionary " << dict_
//            << exit(FatalError);
//    }

    //- the number of points before extrusion
//    const label nPts = pts.size();

        DynList<label> edgeSubsetIDs;
    surf.edgeSubsetIndices(edgeSubsetIDs);
    forAll(edgeSubsetIDs, i)
        surf.removeEdgeSubset(edgeSubsetIDs[i]);
    edgeSubsetIDs.clear();

    edgeSubsetIDs.append
    (
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEWIRE
            )
        )
    );

    forAll(surf.featureEdges(), eI)
        surf.addEdgeToSubset(edgeSubsetIDs.lastElement(), eI);

    //- create downstream patch
    const label nPointsBefore = surf.nPoints();
    const label nEdgesBefore = surf.nFeatureEdges();

    //- add additional vertices
    surf.appendVertex(point(maxX, 0.5 * outerDiameter_, 0.0));
    surf.appendVertex(point(minX, 0.5 * outerDiameter_, 0.0));

    //- create missing edges
    surf.appendFeatureEdge(edge(maxId, nPointsBefore));
    surf.appendFeatureEdge(edge(nPointsBefore, nPointsBefore+1));
    surf.appendFeatureEdge(edge(nPointsBefore+1, minId));

    //- add new edges into subsets
    //- downstream
    edgeSubsetIDs.append
    (
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEDOWNSTREAM
            )
        )
    );

    surf.addEdgeToSubset(edgeSubsetIDs.lastElement(), nEdgesBefore);

    //- die to housing
    edgeSubsetIDs.append
    (
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEHOUSING
            )
        )
    );

    surf.addEdgeToSubset(edgeSubsetIDs.lastElement(), nEdgesBefore+1);

    //- die upstream
    edgeSubsetIDs.append
    (
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEUPSTREAM
            )
        )
    );

    surf.addEdgeToSubset(edgeSubsetIDs.lastElement(), nEdgesBefore+2);

    //- extrude feature edges to a 2D surface
    triSurfaceExtrude2DEdges extruder(surf);

    surfPtr_ = new triSurf();
    extruder.extrudeSurface(*surfPtr_);

    geometricSurfacePatchList& sPatches =
        triSurfModifier(*surfPtr_).patchesAccess();

    forAll(sPatches, patchI)
        sPatches[patchI].geometricType() = "patch";

    Info << "Number of surface triangles " << surfPtr_->size() << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dieSurfaceCreatorGeometry::
dieSurfaceCreatorGeometry
(
    const rollingMillPatchNamesHandler& patchHandler,
    const dictionary& dict,
    const scalar tol
)
:
    dieSurfaceCreator(patchHandler, dict, tol),
    fName_(),
    characteristicPoints_(),
    interpolationType_(),
    profilePoints_(),
    dirFront_(),
    dirBack_(),
    inletDiameter_(0.0),
    outletDiameter_(0.0),
    outerDiameter_(0.0)
{
    //- convert the extension to lowercase
    if( dict.found("geometryFile") )
    {
        //PDJ12/02/2019: read interpolationType from wireMeshDict.
        // Default to b-spline.
        if( dict.found("interpolationType") )
        {
            dict.readIfPresent("interpolationType", interpolationType_);
            Info<< "interpolationType is: " << interpolationType_ << endl;
        }
        else
        {
            interpolationType_ = "b-spline";
        }

        if( !dict_.readIfPresent("outerDiameter", outerDiameter_) )
        {
            FatalError << "outerDiameter is not present in dictionary " << dict_
                << exit(FatalError);
        }

        dict.readIfPresent("geometryFile", fName_);

        parseGeometryFile();

        createTriangulatedSurface();
    }
    else
    {
        FatalError << "geometryFile is not present in dictionary " << dict
            << exit(FatalError);
    }
}

dieSurfaceCreatorGeometry::
dieSurfaceCreatorGeometry
(
    const dieSurfaceCreatorGeometry& die
)
:
    dieSurfaceCreator
    (
        die.patchHandler_,
        die.dict_,
        die.geometryTol_
    ),
    fName_(die.fName_),
    characteristicPoints_(die.characteristicPoints_),
    interpolationType_(die.interpolationType_),
    profilePoints_(die.profilePoints_),
    dirFront_(die.dirFront_),
    dirBack_(die.dirBack_),
    inletDiameter_(die.inletDiameter_),
    outletDiameter_(die.outletDiameter_),
    outerDiameter_(die.outerDiameter_)
{}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar dieSurfaceCreatorGeometry::inletDiameter() const
{
    const pointField& pts = surfPtr_->points();

    scalar minX(VGREAT);
    label minPoint(-1);

    forAll(pts, i)
    {
        if( pts[i].x() < minX )
        {
            minX = pts[i].x();
            minPoint = i;
        }
    }

    if( minPoint != -1 )
        return 2. * pts[minPoint].y();

    return 0.0;
}

scalar dieSurfaceCreatorGeometry::outletDiameter() const
{
    const pointField& pts = surfPtr_->points();

    scalar maxX(VGREAT);
    label maxPoint(-1);

    forAll(pts, i)
    {
        if( pts[i].x() > maxX )
        {
            maxX = pts[i].x();
            maxPoint = i;
        }
    }

    if( maxPoint != -1 )
        return 2. * pts[maxPoint].y();

    return 0.0;
}

scalar dieSurfaceCreatorGeometry::outerDiameter() const
{
    scalar d(0.0);
    dict_.readIfPresent("outerDiameter", d);

    const pointField& pts = surfPtr_->points();
    forAll(pts, i)
        d = max(d, 2.0 * pts[i].y());

    return d;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
