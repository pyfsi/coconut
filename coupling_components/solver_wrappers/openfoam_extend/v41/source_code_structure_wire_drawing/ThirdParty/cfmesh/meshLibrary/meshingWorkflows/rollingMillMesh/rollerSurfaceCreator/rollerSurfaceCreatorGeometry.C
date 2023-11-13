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

#include "rollerSurfaceCreatorGeometry.H"
#include "demandDrivenData.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "splineBase.H"
#include "cubicBSpline.H"
#include "IFstream.H"
#include "boundBox.H"

#include "addToRunTimeSelectionTable.H"

//#define DEBUGSurfaceCreator

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rollerSurfaceCreatorGeometry, 0);
addToRunTimeSelectionTable
(
    rollerSurfaceCreator,
    rollerSurfaceCreatorGeometry,
    dictionary
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rollerSurfaceCreatorGeometry::parseGeometryFile()
{
    //- exit in case the geometry is an empty string
    if( fName_ == "" )
        FatalErrorIn
        (
            "void rollerSurfaceCreatorGeometry::parseGeometryFile()"
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

            characteristicPoints_.append(0.001 * p);
        }
        else if( file.good() )
        {
            Warning << "Unknown token " << t.info() << " found!" << endl;
        }
    }

    Info << "Profile calculation mode " << interpolationType_ << endl;
    Info << "Characteristic points " << characteristicPoints_ << endl;

        //- find the coordinate with all zeros
    bool allZeroX(true);

    forAll(characteristicPoints_, pI)
    {
        const point& p = characteristicPoints_[pI];

        if( mag(p.x()) > SMALL )
            allZeroX = false;
    }

    if( !allZeroX )
    {
        FatalErrorIn
        (
            "void rollerSurfaceCreatorGeometry::checkPointLocations()"
        ) << "The profile is not in the y-z plane" << exit(FatalError);
    }
}

autoPtr<splineBase> rollerSurfaceCreatorGeometry::createRollerProfile() const
{
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

void rollerSurfaceCreatorGeometry::cutContactArea(LongList<point>& points)
{
    if( contactWidth_ < 0.0 || points.size() == 0 )
        return;

    const scalar r = 0.5 * contactWidth_;

    LongList<point> pCopy;

    point p = points[0];
    pCopy.append(points[0]);

    for(label i=1;i<points.size();++i)
    {
        const point& p1 = points[i];

        //- positive side
        if( (p.z() < r) && (p1.z() > r) )
        {
            pCopy.append
            (
                point
                (
                    0.0,
                    (r-p.z()) * (p1.y()-p.y()) / (p1.z()-p.z()) + p.y(),
                    r
                )
            );
        }
        else if( (p1.z() < r) && (p.z() > r) )
        {
            pCopy.append
            (
                point
                (
                    0.0,
                    (r-p1.z()) * (p.y()-p1.y()) / (p.z()-p1.z()) + p1.y(),
                    r
                )
            );
        }

        //- negative side
        if( (p.z() < -r) && (p1.z() > -r) )
        {
            pCopy.append
            (
                point
                (
                    0.0,
                    (-r-p.z()) * (p1.y()-p.y()) / (p1.z()-p.z()) + p.y(),
                    -r
                )
            );
        }
        else if( (p1.z() < -r) && (p.z() > -r) )
        {
            pCopy.append
            (
                point
                (
                    0.0,
                    (-r-p1.z()) * (p.y()-p1.y()) / (p.z()-p1.z()) + p1.y(),
                    -r
                )
            );
        }

        //- update p
        p = p1;
        pCopy.append(p1);
    }

    points.transfer(pCopy);
}

void rollerSurfaceCreatorGeometry::createTriangulatedSurface()
{
    //- tranform the points to the y-z plane
    autoPtr<splineBase> splinePtr = createRollerProfile();

    //- create the spline interpolating a given set of points
    LongList<point> curvePoints;
    splinePtr->createPolyLine(geometryTol_, curvePoints);

    //- create additional points to obey range constraints
    cutContactArea(curvePoints);

    //- translate the points such that topmost point is at y=0.0
    scalar maxY(-VGREAT);
    forAll(curvePoints, i)
        maxY = max(curvePoints[i].y(), maxY);

    forAll(curvePoints, i)
        curvePoints[i].y() -= maxY;

    //- create surface mesh
    triSurf surf;
    triSurfModifier sMod(surf);

    //- patch names
    const label contactId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERCONTACT
            )
        );

    label rollerToAirId(-1);
    if( contactWidth_ > SMALL )
    {
        rollerToAirId =
            surf.addEdgeSubset
            (
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERTOAIR
                )
            );
    }

    const label frontId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERFRONT
            )
        );

    const label backId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERBACK
            )
        );

    const label axisId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERAXIS
            )
        );

    //- create surface points
    pointField& pts = sMod.pointsAccess();

    pts.setSize(2 + curvePoints.size());
    pts[0] = vector(0.0, 0.5 * innerDiameter_, characteristicPoints_[0].z());
    pts[1] =
        vector
        (
            0.0,
            0.5 * innerDiameter_,
            characteristicPoints_[characteristicPoints_.size()-1].z()
        );

    label counter(2);
    forAll(curvePoints, pI)
    {
        pts[counter] = curvePoints[pI];
        pts[counter].y() += outerDiameter_ / 2.0;

        ++counter;
    }

    //- axis region
    surf.appendFeatureEdge(edge(0, 1));
    surf.addEdgeToSubset(axisId, 0);

    //- front region
    surf.appendFeatureEdge(edge(1, pts.size()-1));
    surf.addEdgeToSubset(frontId, 1);

    //- back region
    surf.appendFeatureEdge(edge(0, 2));
    surf.addEdgeToSubset(backId, 2);

    //- contact region
    for(label pI=3;pI<pts.size();++pI)
    {
        const edge e(pI-1, pI);
        surf.appendFeatureEdge(e);

        if
        (
            (rollerToAirId >= 0) &&
            (
                (
                    (pts[e.start()].z() < -0.5 * contactWidth_) ||
                    (pts[e.end()].z() < -0.5 * contactWidth_)
                ) ||
                (
                    (pts[e.start()].z() > 0.5 * contactWidth_) ||
                    (pts[e.end()].z() > 0.5 * contactWidth_)
                )
            )
        )
        {
            //- add to rollerToAir patch
            surf.addEdgeToSubset(rollerToAirId, pI);
        }
        else
        {
            //- add to the contact patch
            surf.addEdgeToSubset(contactId, pI);
        }
    }

    //- transform points into the x-y plane
    forAll(pts, pI)
        pts[pI] = (transformationMatrix_ & pts[pI]);

    //- extrude feature edges to a 2D surface
    triSurfaceExtrude2DEdges extruder(surf);

    surfPtr_ = new triSurf();
    extruder.extrudeSurface(*surfPtr_);

    //- create patches
    triSurfModifier surfMod(*surfPtr_);

    const tensor inverseTransformation = inv(transformationMatrix_);
    pointField& surfPoints = surfMod.pointsAccess();
    forAll(surfPoints, pI)
    {
        point& p = surfPoints[pI];

        p = (inverseTransformation & p);
    }

    //- create feature edges
    const label nPts = surf.nPoints();
    surfPtr_->appendFeatureEdge(edge(0, nPts));
    surfPtr_->appendFeatureEdge(edge(1, 1+nPts));
    surfPtr_->appendFeatureEdge(edge(2, 2+nPts));
    surfPtr_->appendFeatureEdge(edge(nPts-1, 2*nPts-1));

    //- add feature edges needed in case of a polyLine interpolation
    //- other interpolation modes are differentiable
    if( interpolationType_ == "polyLine" )
    {
        for(label i=1;i<characteristicPoints_.size()-1;++i)
            surfPtr_->appendFeatureEdge(edge(2+i, 2+i+nPts));
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollerSurfaceCreatorGeometry::rollerSurfaceCreatorGeometry
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
    ),
    characteristicPoints_(),
    interpolationType_("polyLine"),
    dirBack_(vector::zero),
    dirFront_(vector::zero)
{}

rollerSurfaceCreatorGeometry::rollerSurfaceCreatorGeometry
(
    const rollerSurfaceCreatorGeometry& creator
)
:
    rollerSurfaceCreator(creator),
    characteristicPoints_(creator.characteristicPoints_),
    interpolationType_(creator.interpolationType_),
    dirBack_(creator.dirBack_),
    dirFront_(creator.dirFront_)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rollerSurfaceCreatorGeometry::generateGeometry()
{
    //- geometry file
    if( !dict_.readIfPresent("geometryFile", fName_) )
    {
        FatalErrorIn
        (
            "void rollerSurfaceCreatorGeometry::generateGeometry()"
        ) << "geometryFile is not specified in dictionary "
          << dict_ << exit(FatalError);
    }

    //- interpolation mode
    if( !dict_.readIfPresent("mode", interpolationType_) )
    {
        FatalErrorIn
        (
            "void rollerSurfaceCreatorGeometry::generateGeometry()"
        ) << "mode is not specified in dictionary "
          << dict_ << exit(FatalError);
    }

    //- tangent vectors at the front and back sides
    dict_.readIfPresent("dirBack", dirBack_);
    dict_.readIfPresent("dirFront", dirFront_);

    //- geometry given from a set of points
    parseGeometryFile();

    createTriangulatedSurface();

    detectIndependentRegions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
