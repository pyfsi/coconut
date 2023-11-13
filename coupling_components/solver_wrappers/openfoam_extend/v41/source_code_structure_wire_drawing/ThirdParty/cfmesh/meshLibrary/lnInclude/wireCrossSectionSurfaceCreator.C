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

#include "wireCrossSectionSurfaceCreator.H"
#include "rollingMillMesh.H"
#include "demandDrivenData.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "polyMeshGenModifier.H"
#include "splineBase.H"
#include "cubicBSpline.H"
#include "DxfFile.H"
#include "DxfFileParser.H"
#include "DxfFileToTriSurfConverter.H"
#include "IFstream.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void wireCrossSectionSurfaceCreator::createTransformationMatrix()
{
    Info << "Calculating transformation matrix" << endl;

    //- find the coordinate with all zeros
    bool allZeroX(true);

    bool allPosY(true), allNegY(true);
    bool allPosZ(true), allNegZ(true);

    forAll(characteristicPoints_, pI)
    {
        const point& p = characteristicPoints_[pI];

        if( mag(p.x()) > SMALL )
            allZeroX = false;

        if( p.y() > SMALL )
            allNegY = false;
        if( p.y() < -SMALL )
            allPosY = false;
        if( p.z() > SMALL )
            allNegZ = false;
        if( p.z() < -SMALL )
            allPosZ = false;
    }

    if( !allZeroX )
    {
        FatalErrorIn
        (
            "void wireCrossSectionSurfaceCreator::checkPointLocations()"
        ) << "The profile is not in the y-z plane" << exit(FatalError);
    }

    if( allPosZ )
        typeOfSymmetry_ |= POSITIVE_Z;

    if( allNegZ )
        typeOfSymmetry_ |= NEGATIVE_Z;

    if( allPosY )
        typeOfSymmetry_ |= POSITIVE_Y;

    if( allNegY )
        typeOfSymmetry_ |= NEGATIVE_Y;

    //- cleanup vertices to avoid problems with symmetric setup
    if( typeOfSymmetry_ )
    {
        if( typeOfSymmetry_ & (POSITIVE_Y|NEGATIVE_Y) )
        {
            LongList<point> pCopy;
            forAll(characteristicPoints_, i)
            {
                bool hasNext(false), hasPrev(false);
                if( (i+1) < characteristicPoints_.size() )
                {
                    if( mag(characteristicPoints_[i+1].y()) > SMALL )
                        hasNext = true;
                }
                if( (i-1) >= 0 )
                {
                    if( mag(characteristicPoints_[i-1].y()) > SMALL )
                        hasPrev = true;
                }

                if( hasPrev || hasNext )
                    pCopy.append(characteristicPoints_[i]);
            }

            if( pCopy.size() < characteristicPoints_.size() )
            {
                Info << "Filtered "
                     << (characteristicPoints_.size()-pCopy.size())
                     << " points" << endl;

                characteristicPoints_ = pCopy;

                Info << "First vertex " << characteristicPoints_[0] << endl;
                Info << "Last vertex "
                     << characteristicPoints_[characteristicPoints_.size()-1]
                     << endl;

                //- set the y coordinate fo the first and the last point to 0
                characteristicPoints_[0].y() = 0.0;
                characteristicPoints_[characteristicPoints_.size()-1].y() = 0.0;
            }
        }

        if( typeOfSymmetry_ & (POSITIVE_Z|NEGATIVE_Z) )
        {
            LongList<point> pCopy;
            forAll(characteristicPoints_, i)
            {
                bool hasNext(false), hasPrev(false);
                if( (i+1) < characteristicPoints_.size() )
                {
                    if( mag(characteristicPoints_[i+1].z()) > SMALL )
                        hasNext = true;
                }
                if( (i-1) >= 0 )
                {
                    if( mag(characteristicPoints_[i-1].z()) > SMALL )
                        hasPrev = true;
                }

                if( hasPrev || hasNext )
                    pCopy.append(characteristicPoints_[i]);
            }

            if( pCopy.size() < characteristicPoints_.size() )
            {
                Info << "Filtered "
                     << (characteristicPoints_.size()-pCopy.size())
                     << " points" << endl;

                characteristicPoints_ = pCopy;

                Info << "First vertex " << characteristicPoints_[0] << endl;
                Info << "Last vertex "
                     << characteristicPoints_[characteristicPoints_.size()-1]
                     << endl;

                //- set the z coordinate fo the first and the last point to 0
                characteristicPoints_[0].z() = 0.0;
                characteristicPoints_[characteristicPoints_.size()-1].z() = 0.0;
            }
        }
    }

    //- transformation matrix from the y-z to the x-y plane
    transformationMatrix_ = tensor::zero;
    transformationMatrix_.yy() = 1.0;
    transformationMatrix_.xz() = 1.0;
    transformationMatrix_.zx() = -1.0;

    Info << "Finished calculating transformation matrix" << endl;
}

void wireCrossSectionSurfaceCreator::parseGeometryFile()
{
    //- exit in case the geometry is an empty string
    if( fName_ == "" )
        FatalErrorIn
        (
            "void wireCrossSectionSurfaceCreator::parseGeometryFile()"
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

            // PDJ12/02/2019: Points in list are assumed in m instead of mm
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

    Info << "Profile calculation mode " << interpolationType_ << endl;
    Info << "Characteristic points " << characteristicPoints_ << endl;
}

autoPtr<splineBase> wireCrossSectionSurfaceCreator::createProfile() const
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

void wireCrossSectionSurfaceCreator::createTriangulatedSurface()
{
    //- tranform the points to the y-z plane
    autoPtr<splineBase> splinePtr = createProfile();

    //- create the spline interpolating a given set of points
    LongList<point> curvePoints;
    splinePtr->createPolyLine(geometryTol_, curvePoints);

    //- create surface mesh
    triSurf surf;
    triSurfModifier sMod(surf);

    //- check the symmetry of the setup
    label counter(0);

    if( typeOfSymmetry_ & POSITIVE_Y )
        ++counter;
    if( typeOfSymmetry_ & NEGATIVE_Y )
        ++counter;
    if( typeOfSymmetry_ & POSITIVE_Z )
        ++counter;
    if( typeOfSymmetry_ & NEGATIVE_Z )
        ++counter;

    //- create feature edges from points in this subset
    label pId = sMod.surface().pointSubsetIndex("_corners_");
    if( pId < 0 )
        pId = sMod.surface().addPointSubset("_corners_");

    //- create surface points
    pointField& pts = sMod.pointsAccess();
    pts.setSize(curvePoints.size());

    forAllReverse(curvePoints, pI)
    {
        pts[pI] = curvePoints[pI];

        //- create feature edges
        if( interpolationType_ == "polyLine" )
            sMod.surface().addPointToSubset(pId, pI);
    }

    //- create feature edges that are extruded into triangles
    const label contactId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForWire
            (
                rollingMillPatchNamesHandler::WIRECONTACT
            )
        );

    DynList<word> symmPatchNames;

    if( counter )
    {
        const label nPts = pts.size();

        label originI(-1);
        label yPlaneId(-1), zPlaneId(-1);
        if( counter > 1 )
        {
            //- there exists a corner at the origin of the system
            originI = pts.size();
            surf.appendVertex(vector::zero);

            yPlaneId =
                surf.addEdgeSubset
                (
                    patchHandler_.patchNameForWire
                    (
                        rollingMillPatchNamesHandler::WIRESYMMY
                    )
                );
            symmPatchNames.append(surf.edgeSubsetName(yPlaneId));

            zPlaneId =
                surf.addEdgeSubset
                (
                    patchHandler_.patchNameForWire
                    (
                        rollingMillPatchNamesHandler::WIRESYMMZ
                    )
                );
            symmPatchNames.append(surf.edgeSubsetName(zPlaneId));
        }

        label nEdges(0);
        for(label i=1;i<nPts;++i)
        {
            surf.appendFeatureEdge(edge(i-1, i));
            surf.addEdgeToSubset(contactId, nEdges++);
        }

        if( originI >= 0 )
        {
            //- two symmetry planes detected
            surf.appendFeatureEdge(edge(originI-1, originI));
            if( mag(pts[originI-1].y()) < mag(pts[0].y()) )
            {
                surf.addEdgeToSubset(yPlaneId, nEdges++);
            }
            else
            {
                surf.addEdgeToSubset(zPlaneId, nEdges++);
            }

            surf.appendFeatureEdge((edge(0, originI)));
            if( mag(pts[0].z()) < mag(pts[originI-1].z()) )
            {
                surf.addEdgeToSubset(zPlaneId, nEdges++);
            }
            else
            {
                surf.addEdgeToSubset(yPlaneId, nEdges++);
            }
        }
        else
        {
            //- only one symmetry plane detected
            if( mag(pts[0].y()) < mag(pts[0].z()) )
            {
                //- point is in the y=0 plane
                yPlaneId =
                    surf.addEdgeSubset
                    (
                        patchHandler_.patchNameForWire
                        (
                            rollingMillPatchNamesHandler::WIRESYMMY
                        )
                    );
                symmPatchNames.append(surf.edgeSubsetName(yPlaneId));

                surf.appendFeatureEdge(edge(0, nPts-1));
                surf.addEdgeToSubset(yPlaneId, nEdges++);
            }
            else
            {
                //- point is in the z=0 plane
                zPlaneId =
                    surf.addEdgeSubset
                    (
                        patchHandler_.patchNameForWire
                        (
                            rollingMillPatchNamesHandler::WIRESYMMZ
                        )
                    );
                symmPatchNames.append(surf.edgeSubsetName(zPlaneId));

                surf.appendFeatureEdge(edge(0, nPts-1));
                surf.addEdgeToSubset(zPlaneId, nEdges++);
            }
        }
    }
    else
    {
        label nEdges(0);
        forAll(pts, eI)
        {
            surf.appendFeatureEdge(edge(eI, pts.fcIndex(eI)));
            surf.addEdgeToSubset(contactId, nEdges++);
        }
    }

    //- transform points into the x-y plane
    forAll(pts, pI)
    {
        pts[pI] = (transformationMatrix_ & pts[pI]);
    }

    //- extrude feature edges to a 2D surface
    triSurfaceExtrude2DEdges extruder(surf);

    surfPtr_ = new triSurf();
    extruder.extrudeSurface(*surfPtr_);

    //- set correct patch types at symmetry planes
    triSurfModifier surfMod(*surfPtr_);
    geometricSurfacePatchList& patches = surfMod.patchesAccess();

    forAll(patches, patchI)
    {
        if( symmPatchNames.contains(patches[patchI].name()) )
            patches[patchI].geometricType() = "symmetryPlane";
    }

    //- back transform the surface to the y-z space
    const tensor inverseTransformation = inv(transformationMatrix_);
    pointField& surfPoints = surfMod.pointsAccess();
    forAll(surfPoints, pI)
    {
        point& p = surfPoints[pI];

        p = (inverseTransformation & p);
    }
}

void wireCrossSectionSurfaceCreator::parseDXFFile()
{
    //- transformation matrix from the y-z to the x-y plane
    transformationMatrix_ = tensor::zero;
    transformationMatrix_.yy() = 1.0;
    transformationMatrix_.xz() = 1.0;
    transformationMatrix_.zx() = -1.0;

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

    triSurf surf;

    try
    {
        DxfFileToTriSurfConverter converter(file, 0.1);
        converter.Convert(surf);
    }
    catch (const std::exception& e)
    {
        Info << "ERROR: Error while converting DXF file "
             << "\"" << fName_ << "\" to triSurf: " << e.what() << endl;
        std::abort();
    }

    //- extrude feature edges to a 2D surface
    triSurfaceExtrude2DEdges extruder(surf);

    surfPtr_ = new triSurf();
    extruder.extrudeSurface(*surfPtr_);

    //- create patches
    triSurfModifier surfMod(*surfPtr_);
    geometricSurfacePatchList& patches = surfMod.patchesAccess();
    LongList<labelledTri>& triangles = surfMod.facetsAccess();

    //- transform the surface into the y-z plane
    const tensor inverseTransformation = inv(transformationMatrix_);
    pointField& surfPoints = surfMod.pointsAccess();
    forAll(surfPoints, pI)
    {
        point& p = surfPoints[pI];

        p = (inverseTransformation & p);
    }

    if( wireProfileCentering_ )
    {
        //- calculate bounding box and position its centre at origin
        boundBox bb(surfPoints);
        vector disp(0.0, bb.midpoint().y(), bb.midpoint().z());
        Info << "Bounding box centre " << disp << endl;
        surfPoints -= disp;
    }

    //- set correct patches
    forAll(triangles, triI)
    {
        triangles[triI].region() = 0;
    }

    //- no symmetry
    patches.setSize(1);
    patches[0].name() =
        patchHandler_.patchNameForWire
        (
            rollingMillPatchNamesHandler::WIRECONTACT
        );
    patches[0].geometricType() = "patch";
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

wireCrossSectionSurfaceCreator::wireCrossSectionSurfaceCreator
(
    const dictionary& wireDict,
    const rollingMillPatchNamesHandler& patchHandler,
    const scalar tol
)
:
    fName_(),
    patchHandler_(patchHandler),
    geometryTol_(tol),
    transformationMatrix_(symmTensor::zero),
    characteristicPoints_(),
    interpolationType_(),
    surfPtr_(NULL),
    typeOfSymmetry_(NONE),
    dirFront_(),
    dirBack_(),
    wireProfileCentering_(false)
{
    if( wireDict.found("dxfFile") )
    {
        wireDict.readIfPresent("dxfFile", fName_);
    }
    else if( wireDict.found("geometryFile") )
    {
        wireDict.readIfPresent("geometryFile", fName_);
    }
    else if( wireDict.found("surfaceFile") )
    {
        wireDict.readIfPresent("surfaceFile", fName_);

        surfPtr_ = new triSurf(fName_);
        return;
    }

    //PDJ12/02/2019: read interpolationType from wireMeshDict.
    // Default to b-spline.
    if( wireDict.found("interpolationType") )
    {
        wireDict.readIfPresent("interpolationType", interpolationType_);
    }
    else
    {
        interpolationType_ = "b-spline";
    }

    //- check if the profile shall be centred
    wireDict.readIfPresent("wireProfileCentering", wireProfileCentering_);

    //- convert the extension to lowercase
    word ext = fName_.ext();
    for(unsigned i=0;i<ext.size();++i)
        ext[i] = tolower(ext[i]);

    if( ext == "dxf" )
    {
        //- generate the surface mesh from a DXF profile
        parseDXFFile();
    }
    else
    {
        parseGeometryFile();

        createTransformationMatrix();

        createTriangulatedSurface();
    }
}

wireCrossSectionSurfaceCreator::~wireCrossSectionSurfaceCreator()
{
    deleteDemandDrivenData(surfPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const triSurf& wireCrossSectionSurfaceCreator::surface() const
{
    return *surfPtr_;
}

const triSurf* wireCrossSectionSurfaceCreator::transformedSurface() const
{
    triSurf* surfPtr = new triSurf(*surfPtr_);

    triSurfModifier sMod(*surfPtr);

    pointField& pts = sMod.pointsAccess();

    forAll(pts, pI)
        pts[pI] = (transformationMatrix_ & pts[pI]);

    return surfPtr;
}

void wireCrossSectionSurfaceCreator::updatePointLocations
(
    polyMeshGen& mesh
) const
{
    Info << "Updating point locations" << endl;

    const faceListPMG& faces = mesh.faces();
    polyMeshGenModifier meshModifier(mesh);
    pointFieldPMG& points = meshModifier.pointsAccess();

    autoPtr<splineBase> splinePtr = createProfile();

    boolList usedPoint(points.size(), false);

    const PtrList<boundaryPatch>& patches = mesh.boundaries();

    forAll(patches, patchI)
    {
        const boundaryPatch& patch = patches[patchI];

        //- check if the patch represents the contact patch
        if( patch.patchName().find("contactPatch") == word::npos )
            continue;

        const label s = patch.patchStart();
        const label e = s + patch.patchSize();

        for(label faceI=s;faceI<e;++faceI)
        {
            const face& f = faces[faceI];

            forAll(f, pI)
            {
                if( usedPoint[f[pI]] )
                    continue;

                usedPoint[f[pI]] = true;

                point& p = points[f[pI]];

                point newP = splinePtr->nearestPointOnSpline(p);

                p.x() = newP.x();
                p.y() = newP.y();
            }
        }
    }

    Info << "Finished updating point locations" << endl;
}

void wireCrossSectionSurfaceCreator::transformToPosition
(
    polyMeshGen& mesh
) const
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
