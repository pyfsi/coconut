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

#include "rollerSurfaceCreator.H"
#include "demandDrivenData.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "polyMeshGenModifier.H"
#include "splineBase.H"
#include "cubicBSpline.H"
#include "IFstream.H"
#include "boundBox.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "triSurfaceCleanupDuplicates.H"
#include "triSurfaceRemoveFacets.H"
#include "triSurfaceChecks.H"
#include "rollingMillMesh.H"

#include <map>

//#define DEBUGSurfaceCreator

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(rollerSurfaceCreator, 0);
defineRunTimeSelectionTable(rollerSurfaceCreator, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rollerSurfaceCreator::createTransformationMatrix()
{
    Info << "Calculating transformation matrix" << endl;

    //- the first and the last point differ in the z coordinate
    transformationMatrix_.yy() = 1.0;
    transformationMatrix_.xz() = 1.0;
    transformationMatrix_.zx() = -1.0;

    Info << "Finished calculating transformation matrix" << endl;
}

void rollerSurfaceCreator::mergeDuplicates()
{
    //- extruded extruded surface first
    triSurfaceExtrude2DEdges extruder(*surfPtr_);

    const triSurf* extrudedPtr = extruder.extrudeSurface();

    //- cleanup duplicates from the extruded surface
    meshOctree octree(*extrudedPtr);
    meshOctreeCreator(octree).createOctreeWithRefinedBoundary(15, 50);

    triSurfaceCleanupDuplicates cleaner(octree);
    cleaner.mergeIdentities();

    //- find the average z to take the vertices at z=0, only
    const pointField& points = extrudedPtr->points();

    scalar avgZ = (min(points).z() + max(points).z()) / 2;

    deleteDemandDrivenData(surfPtr_);

    labelList newLabel(points.size(), -1);
    label nNewPoints(0);

    const triSurf& es = *extrudedPtr;

    deleteDemandDrivenData(surfPtr_);
    surfPtr_ = new triSurf();
    forAll(es, tI)
    {
        const labelledTri& lt = es[tI];

        DynList<label> e;
        forAll(lt, i)
        {
            if( points[lt[i]].z() < avgZ )
            {
                if( newLabel[lt[i]] == -1 )
                    newLabel[lt[i]] = nNewPoints++;

                e.append(newLabel[lt[i]]);
            }
        }

        if( e.size() == 2 )
        {
            surfPtr_->appendFeatureEdge(edge(e[0], e[1]));
        }
    }

    triSurfModifier sMod(*surfPtr_);
    pointField& pts = sMod.pointsAccess();
    pts.setSize(nNewPoints);

    forAll(points, pI)
    {
        if( newLabel[pI] >= 0 )
            pts[newLabel[pI]] = points[pI];
    }

    //- delete extruded surface
    deleteDemandDrivenData(extrudedPtr);

    # ifdef DEBUGSurfaceCreator
    forAll(pts, pI)
    {
        for(label pJ=pI+1;pJ<pts.size();++pJ)
        {
            if( mag(pts[pI] - pts[pJ]) < 1e-6 )
                Info << "Duplicate vertices in surface "
                     << pI << " coord " << pts[pI]
                     << " and " << pJ << " coord " << pts[pJ] << endl;
        }
    }
    # endif
}

void rollerSurfaceCreator::detectIndependentRegions() const
{
    const label nIndependentRegions =
        triSurfaceChecks::checkDisconnectedParts(*surfPtr_);

    if( nIndependentRegions != 1 )
    {
        surfPtr_->writeSurface("surfWithInvalidRegions.fms");

        Warning << "The geometry of a roll in file " << fName_
            << " does not consist of a single closed region."
            << " It is written into surfWithInvalidRegions.fms." << endl;

        labelLongList facetsInRegion;
        triSurfaceChecks::checkDisconnectedParts(*surfPtr_, facetsInRegion);

        const LongList<labelledTri>& trias = surfPtr_->facets();
        const geometricSurfacePatchList& patches = surfPtr_->patches();
        label mainRegion(-1);

        const word backName =
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERBACK
            );
        const word frontName =
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERFRONT
            );
        const word axisName =
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERAXIS
            );

        forAll(trias, tI)
        {
            const label patchI = trias[tI].region();
            const word pName = patches[patchI].name();

            if
            (
                (pName.find(backName) != word::npos) ||
                (pName.find(frontName) != word::npos) ||
                (pName.find(axisName) != word::npos)
            )
            {
                mainRegion = facetsInRegion[tI];
            }
        }

        if( mainRegion < 0 )
        {
            FatalErrorIn
            (
                "void rollerSurfaceCreator::detectIndependentRegions() const"
            ) << "Cannot detect region of interest" << exit(FatalError);
        }

        const label rId = surfPtr_->addFacetSubset("removeFacets");
        forAll(facetsInRegion, i)
        {
            if( facetsInRegion[i] != mainRegion )
                surfPtr_->addFacetToSubset(rId, i);
        }

        triSurfaceRemoveFacets remover(*surfPtr_);
        remover.selectFacetsInSubset("removeFacets");
        remover.removeFacets();

        surfPtr_->removeFacetSubset(rId);
        surfPtr_->writeSurface("surfAfterCleanup.fms");
    }
}

void rollerSurfaceCreator::cutContactArea(triSurf& surf) const
{
    if( contactWidth_ < 0.0 )
        return;

    const pointField& points = surf.points();

    triSurfModifier sMod(surf);
    edgeLongList& featureEdges = sMod.featureEdgesAccess();

    const scalar pxmax = rollAxialShift_ + 0.5 * contactWidth_;
    const scalar pxmin = rollAxialShift_ - 0.5 * contactWidth_;

    forAll(featureEdges, eI)
    {
        edge& e = featureEdges[eI];

        scalar p0x = points[e.start()].x();
        scalar p1x = points[e.end()].x();
        scalar p0y = points[e.start()].y();
        scalar p1y = points[e.end()].y();

        //- cut on the positive side of the x axis
        if( (p1x > pxmax) && (p0x < pxmax) )
        {
            const scalar py = (pxmax - p0x) * ((p1y - p0y) / (p1x - p0x)) + p0y;

            const label s = e.start();
            e.start() = points.size();
            featureEdges.append(edge(s, points.size()));

            surf.appendVertex(point(pxmax, py, 0.0));
        }
        else if( (p0x > pxmax) && (p1x < pxmax) )
        {
            const scalar py = (pxmax - p1x) * ((p0y - p1y) / (p0x - p1x)) + p1y;

            const label s = e.start();
            e.start() = points.size();
            featureEdges.append(edge(s, points.size()));

            surf.appendVertex(point(pxmax, py, 0.0));
        }

        //- cut on the negative side
        p0x = points[e.start()].x();
        p1x = points[e.end()].x();
        p0y = points[e.start()].y();
        p1y = points[e.end()].y();

        if( (p1x > pxmin) && (p0x < pxmin) )
        {
            const scalar py = (pxmin - p0x) * (p1y - p0y) / (p1x - p0x) + p0y;

            const label s = e.start();
            e.start() = points.size();
            featureEdges.append(edge(s, points.size()));

            surf.appendVertex(point(pxmin, py, 0.0));
        }
        else if( (p0x > pxmin) && (p1x < pxmin) )
        {
            const scalar py = (pxmin - p1x) * (p0y - p1y) / (p0x - p1x) + p1y;

            const label s = e.start();
            e.start() = points.size();
            featureEdges.append(edge(s, points.size()));

            surf.appendVertex(point(pxmin, py, 0.0));
        }
    }
}

void rollerSurfaceCreator::positionProfileInAxialDirection()
{
    bool centreProfile(false);
    if( !dict_.readIfPresent("centreInAxialDirection", centreProfile) )
        return;
    if( !centreProfile )
        return;

    Info << "Translating roll at position " << rollPosition_
         << " to the centre of the coordinate system" << endl;

    const boundBox bb(surfPtr_->points());

    if( rollPosition_ == "topRoll" || rollPosition_ == "bottomRoll" )
    {
        const vector disp(0., 0., -0.5 * (bb.max().z() + bb.min().z()));

        Info << "Displacement " << disp << endl;

        triSurfModifier(*surfPtr_).pointsAccess() += disp;
    }
    else if( rollPosition_ == "leftRoll" || rollPosition_ == "rightRoll" )
    {
        const vector disp(0., -0.5 * (bb.max().y() + bb.min().y()), 0.);

        Info << "Displacement " << disp << endl;

        triSurfModifier(*surfPtr_).pointsAccess() += disp;
    }

    Info << "Finished positioning roll at position " << rollPosition_ << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollerSurfaceCreator::rollerSurfaceCreator
(
    const word rollPosition,
    const direction symm,
    const rollingMillPatchNamesHandler& patchHandler,
    const dictionary& dict,
    const scalar tol
)
:
    patchHandler_(patchHandler),
    fName_(),
    dict_(dict),
    typeOfSymmetry_(symm),
    geometryTol_(tol),
    transformationMatrix_(symmTensor::zero),
    innerDiameter_(-SMALL),
    outerDiameter_(-SMALL),
    rollWidth_(-SMALL),
    rollAxialShift_(0.0),
    rollRadialShift_(0.0),
    surfPtr_(NULL),
    rollPosition_(rollPosition),
    rollerAxis_(),
    contactWidth_(-1.0)
{
    if( !dict_.readIfPresent("rollWidth", rollWidth_) )
    {
        FatalError
          << "rollWidth is not present in dictionary " << dict_
          << exit(FatalError);
    }

    if( !dict_.readIfPresent("innerDiameter", innerDiameter_) )
    {
        FatalError
          << "innerDiameter is not present in dictionary " << dict_
          << exit(FatalError);
    }

    if( !dict_.readIfPresent("outerDiameter", outerDiameter_) )
    {
        FatalError
          << "outerDiameter is not present in dictionary " << dict_
          << exit(FatalError);
    }

    dict_.readIfPresent("contactWidth", contactWidth_);

    createTransformationMatrix();
}

rollerSurfaceCreator::rollerSurfaceCreator(const rollerSurfaceCreator& creator)
:
    patchHandler_(creator.patchHandler_),
    fName_(creator.fName_),
    dict_(creator.dict_),
    typeOfSymmetry_(creator.typeOfSymmetry_),
    geometryTol_(creator.geometryTol_),
    transformationMatrix_(creator.transformationMatrix_),
    innerDiameter_(creator.innerDiameter_),
    outerDiameter_(creator.outerDiameter_),
    rollWidth_(creator.rollWidth_),
    rollAxialShift_(creator.rollAxialShift_),
    rollRadialShift_(creator.rollRadialShift_),
    surfPtr_(NULL),
    rollPosition_(creator.rollPosition_),
    rollerAxis_(creator.rollerAxis_),
    contactWidth_(creator.contactWidth_)
{
    if( creator.surfPtr_ )
        surfPtr_ = new triSurf(*creator.surfPtr_);
}

rollerSurfaceCreator::~rollerSurfaceCreator()
{
    deleteDemandDrivenData(surfPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const scalar& rollerSurfaceCreator::geometricTolerance() const
{
    return geometryTol_;
}

label rollerSurfaceCreator::symmetryType() const
{
    return typeOfSymmetry_;
}

const word& rollerSurfaceCreator::rollPosition() const
{
    return rollPosition_;
}

const vector& rollerSurfaceCreator::rotationAxis() const
{
    return rollerAxis_;
}

point rollerSurfaceCreator::rotationOrigin() const
{
    return vector::zero;
}

const scalar& rollerSurfaceCreator::innerRollDiameter() const
{
    return innerDiameter_;
}

const scalar& rollerSurfaceCreator::outerRollDiameter() const
{
    return outerDiameter_;
}

const scalar& rollerSurfaceCreator::rollWidth() const
{
    return rollWidth_;
}

const scalar& rollerSurfaceCreator::axialShift() const
{
    return rollAxialShift_;
}

const scalar& rollerSurfaceCreator::radialShift() const
{
    return rollRadialShift_;
}

const scalar& rollerSurfaceCreator::contactWidth() const
{
    return contactWidth_;
}

const triSurf& rollerSurfaceCreator::surface() const
{
    return *surfPtr_;
}

const triSurf* rollerSurfaceCreator::transformedSurface() const
{
    triSurf* surfPtr = new triSurf(*surfPtr_);

    triSurfModifier sMod(*surfPtr);

    pointField& pts = sMod.pointsAccess();

    forAll(pts, pI)
        pts[pI] = (transformationMatrix_ & pts[pI]);

    return surfPtr;
}

const triSurf* rollerSurfaceCreator::transformedSurfaceForContact() const
{
    triSurf* surfPtr = new triSurf(*surfPtr_);

    triSurfModifier sMod(*surfPtr);

    pointField& pts = sMod.pointsAccess();

    forAll(pts, pI)
    {
        pts[pI].y() -= rollRadialShift_;
        pts[pI] = (transformationMatrix_ & pts[pI]);
    }

    return surfPtr;
}

void rollerSurfaceCreator::transformToPosition(polyMeshGen& mesh) const
{
    polyMeshGenModifier meshModifier(mesh);
    pointFieldPMG& points = meshModifier.pointsAccess();

    //- calculate the inverse tranformation matrix
    const tensor backwardTransformation = inv(transformationMatrix_);

    //- transform the points back to the original coordinate system
    forAll(points, pI)
        points[pI] = (backwardTransformation & points[pI]);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
