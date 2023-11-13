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

#include "rollingMillMesh.H"
#include "polyMeshGenModifier.H"
#include "revolve2DMesh.H"
#include "triSurfModifier.H"

#include "boundBox.H"

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const wordHashSet rollGeometryInfo::rollPositions = populateRollPositions();

const tensor rollGeometryInfo::transformationMatrix =
    tensor(0., 0., 1.0, 0., 1.0, 0., -1.0, 0., 0.);

wordHashSet rollGeometryInfo::populateRollPositions()
{
    wordHashSet positions;

    positions.insert("topRoll");
    positions.insert("bottomRoll");
    positions.insert("leftRoll");
    positions.insert("rightRoll");
    positions.insert("bottomLeftRoll");
    positions.insert("bottomRightRoll");

    return positions;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollGeometryInfo::rollGeometryInfo
(
    triSurf* surfacePtr,
    const dictionary& rollDict,
    const word rollPosition
)
:
    surfPtr_(surfacePtr),
    rollDictPtr_(&rollDict),
    outerDiameter_(0.0),
    rollPosition_(rollPosition),
    rollOffset_(""),
    circumGrading_(1.2),
    circumScaling_(1.0),
    circumferentialResolution_(20),
    axialResolution_(10),
    axialTranslation_(0.0),
    radialTranslation_(0.0),
    radialGrading_(1.0),
    axialGrading_(1.0),
    axialScaling_(1.0),
    radialScaling_(1.0),
    offsetX_(0.0),
    radialOffset_(0.0),
    contactWidth_(-1.0),
    extraContactAngle_(10.0),
    symmetryType_(rollingMillGeometryHandler::NONE),
    isClosedSetup_(false),
    centreInAxialDirection_(false)
{
    //- read meshing parameters
    //- position every roll in axial direction
    if( rollDict.found("centreInAxialDirection") )
    {
        centreInAxialDirection_ =
            readBool(rollDict.lookup("centreInAxialDirection"));
    }

    //- extra contact angle
    if( rollDict.found("extraContactAngle") )
    {
        extraContactAngle_ = readScalar(rollDict.lookup("extraContactAngle"));
    }

    if( rollDict.found("circumferentialResolution") )
        circumferentialResolution_ =
            readLabel(rollDict.lookup("circumferentialResolution"));

    if( rollDict.found("circumGrading") )
        circumGrading_ = readScalar(rollDict.lookup("circumGrading"));

    if( rollDict.found("circumScaling") )
        circumScaling_ = readScalar(rollDict.lookup("circumScaling"));

    if( rollDict.found("radialGrading") )
        radialGrading_ = readScalar(rollDict.lookup("radialGrading"));

    if( rollDict.found("radialScaling") )
        radialScaling_ = readScalar(rollDict.lookup("radialScaling"));

    if( rollDict.found("axialGrading") )
        axialGrading_ = readScalar(rollDict.lookup("axialGrading"));

    if( rollDict.found("axialScaling") )
        axialScaling_ = readScalar(rollDict.lookup("axialScaling"));

    //- data for wire straighteners
    //- offset in the direction of the wire
    if( rollDict.found("offsetX") )
    {
        offsetX_ = readScalar(rollDict.lookup("offsetX"));
    }

    //- offset in the radial direction
    if( rollDict.found("radialOffset") )
    {
        radialOffset_ = readScalar(rollDict.lookup("radialOffset"));
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollGeometryInfo::~rollGeometryInfo()
{
    deleteDemandDrivenData(surfPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const triSurf& rollGeometryInfo::surface() const
{
    return *surfPtr_;
}

const dictionary& rollGeometryInfo::rollDict() const
{
    return *rollDictPtr_;
}

word& rollGeometryInfo::rollOffset()
{
    return rollOffset_;
}

const word& rollGeometryInfo::rollOffset() const
{
    return rollOffset_;
}

scalar& rollGeometryInfo::outerDiameter()
{
    return outerDiameter_;
}

const scalar& rollGeometryInfo::outerDiameter() const
{
    return outerDiameter_;
}

scalar& rollGeometryInfo::circumGrading()
{
    return circumGrading_;
}

const scalar& rollGeometryInfo::circumGrading() const
{
    return circumGrading_;
}

scalar& rollGeometryInfo::circumScaling()
{
    return circumScaling_;
}

const scalar& rollGeometryInfo::circumScaling() const
{
    return circumScaling_;
}

label& rollGeometryInfo::circumResolution()
{
    return circumferentialResolution_;
}

label rollGeometryInfo::circumResolution() const
{
    return circumferentialResolution_;
}

label& rollGeometryInfo::axialResolution()
{
    return axialResolution_;
}

label rollGeometryInfo::axialResolution() const
{
    return axialResolution_;
}

scalar& rollGeometryInfo::axialShiftDistance()
{
    return axialTranslation_;
}

const scalar& rollGeometryInfo::axialShiftDistance() const
{
    return axialTranslation_;
}

scalar& rollGeometryInfo::radialShiftDistance()
{
    return radialTranslation_;
}

const scalar& rollGeometryInfo::radialShiftDistance() const
{
    return radialTranslation_;
}

scalar& rollGeometryInfo::radialGrading()
{
    return radialGrading_;
}

const scalar& rollGeometryInfo::radialGrading() const
{
    return radialGrading_;
}

scalar& rollGeometryInfo::axialGrading()
{
    return axialGrading_;
}

const scalar& rollGeometryInfo::axialGrading() const
{
    return axialGrading_;
}

scalar& rollGeometryInfo::radialScaling()
{
    return radialScaling_;
}

const scalar& rollGeometryInfo::radialScaling() const
{
    return radialScaling_;
}

scalar& rollGeometryInfo::axialScaling()
{
    return axialScaling_;
}

const scalar& rollGeometryInfo::axialScaling() const
{
    return axialScaling_;
}

scalar& rollGeometryInfo::offsetX()
{
    return offsetX_;
}

const scalar& rollGeometryInfo::offsetX() const
{
    return offsetX_;
}

scalar& rollGeometryInfo::radialOffset()
{
    return radialOffset_;
}

const scalar& rollGeometryInfo::radialOffset() const
{
    return radialOffset_;
}

scalar& rollGeometryInfo::contactWidth()
{
    return contactWidth_;
}

const scalar& rollGeometryInfo::contactWidth() const
{
    return contactWidth_;
}

scalar& rollGeometryInfo::extraContactAngle()
{
    return extraContactAngle_;
}

scalar rollGeometryInfo::extraContactAngle() const
{
    return extraContactAngle_;
}

direction& rollGeometryInfo::typeOfSymmetry()
{
    return symmetryType_;
}

direction rollGeometryInfo:: typeOfSymmetry() const
{
    return symmetryType_;
}

bool& rollGeometryInfo::isClosedSetup()
{
    return isClosedSetup_;
}

bool rollGeometryInfo::isClosedSetup() const
{
    return isClosedSetup_;
}

bool& rollGeometryInfo::centreInAxialDirection()
{
    return centreInAxialDirection_;
}

bool rollGeometryInfo::centreInAxialDirection() const
{
    return centreInAxialDirection_;
}

const triSurf* rollGeometryInfo::transformedSurface() const
{
    triSurf* surfPtr = new triSurf(*surfPtr_);

    //- tranform the points into the x-y plane
    triSurfModifier sMod(*surfPtr);
    pointField& pts = sMod.pointsAccess();

    forAll(pts, pI)
        pts[pI] = (transformationMatrix & pts[pI]);

    return surfPtr;
}

void rollGeometryInfo::backTransform2DMeshToPosition
(
    polyMeshGen& mesh
) const
{
    polyMeshGenModifier meshModifier(mesh);
    pointFieldPMG& points = meshModifier.pointsAccess();

    //- calculate the inverse tranformation matrix
    const tensor backwardTransformation = inv(transformationMatrix);

    //- transform the points back to the original coordinate system
    forAll(points, pI)
        points[pI] = (backwardTransformation & points[pI]);
}

word rollGeometryInfo::rollPosition() const
{
    return rollPosition_;
}

vector rollGeometryInfo::rotationAxis() const
{
    if( rollPosition_ == "topRoll" )
    {
        //- roll revolves in the positive z direction
        return vector(0., 0., -1.0);
    }
    else if( rollPosition_ == "bottomRoll" )
    {
        //- roll revolves in the negative z direction
        return vector(0., 0., 1.0);
    }
    else if( rollPosition_ == "leftRoll" )
    {
        //- roll revolves in the negative y direction
        return vector(0., -1.0, 0.0);
    }
    else if( rollPosition_ == "rightRoll" )
    {
        //- roll revolves in the positive y direction
        return vector(0., 1.0, 0.);
    }
    else if( rollPosition_ == "bottomLeftRoll" )
    {
        //- roll revolves in the negative y and positive z direction
        return vector(0., -sqrt(3.0)/2.0, 0.5);
    }
    else if( rollPosition_ == "bottomRightRoll" )
    {
        //- roll revolves in the positive y direction and positive z direction
        return vector(0., sqrt(3.0)/2.0, 0.5);
    }

    FatalErrorIn
    (
        "vector rollingMillMesh::rollGeometryInfo::rotationAxis() const"
    ) << "Unknown roll position " << rollPosition_ << exit(FatalError);

    return vector::zero;
}

vector rollGeometryInfo::axialShift() const
{
    if( rollPositions.found(rollPosition_) )
    {
        //- shift the roll in the axial direction
        return rotationAxis() * axialTranslation_;
    }
    else
    {
        //- shift the wire in the axial direction
        return vector(axialTranslation_, 0., 0.);
    }
}

vector rollGeometryInfo::radialShift(const scalar wireDiameter) const
{
    scalar shiftLength = 0.5 * (outerDiameter_ + wireDiameter);
    shiftLength += radialOffset_;
    shiftLength += radialTranslation_;

    if( rollPosition_ == "topRoll" )
    {
        //- roll revolves in the positive z direction
        return vector(0., 1., 0.) * shiftLength;
    }
    else if( rollPosition_ == "bottomRoll" )
    {
        //- roll revolves in the negative z direction
        return vector(0., -1., 0.) * shiftLength;
    }
    else if( rollPosition_ == "leftRoll" )
    {
        //- roll revolves in the negative y direction
        return vector(0., 0., -1.0) * shiftLength;
    }
    else if( rollPosition_ == "rightRoll" )
    {
        //- roll revolves in the positive y direction
        return vector(0., 0., 1.) * shiftLength;
    }
    else if( rollPosition_ == "bottomLeftRoll" )
    {
        //- roll revolves in the negative y and positive z direction
        return vector(0., -0.5, -sqrt(3.0)/2.0) * shiftLength;
    }
    else if( rollPosition_ == "bottomRightRoll" )
    {
        //- roll revolves in the positive y direction and positive z direction
        return vector(0., -0.5 ,sqrt(3.0)/2.0) * shiftLength;
    }

    FatalErrorIn
    (
        "vector rollingMillMesh::rollGeometryInfo::"
        "radialShift(const scalar) const"
    ) << "Unknown roll position " << rollPosition_ << exit(FatalError);

    return vector::zero;
}

point rollGeometryInfo::origin() const
{
    return vector::zero;
}

dictionary rollGeometryInfo::anisotropic2DMeshing() const
{
    //- create a dictionary with anisotropic properties
    dictionary dict("anisotropicSources");

    //- find the axial direction
    const vector axialVec = rotationAxis();

    //- find the radial direction
    vector radialVec = -1.0 * radialShift();
    radialVec /= (mag(radialVec) + VSMALL);

    //- calculate the scaling distance in radial and axial directions
    scalar axialDistance(0.0);
    scalar minRadialDistance(VGREAT);
    scalar maxRadialDistance = 0.5 * outerDiameter_;

    const pointField& pts = surfPtr_->points();
    forAll(pts, pI)
    {
        const point& p = pts[pI];

        axialDistance = max(axialDistance, mag(p & axialVec));

        minRadialDistance = min(minRadialDistance, mag(p & radialVec));
    }

    //- define a dictionary for the radial direction
    dictionary radialPlane("radialGrading");
    radialPlane.add("type", "plane");
    radialPlane.add("normal", -1.0 * (transformationMatrix & radialVec));
    radialPlane.add
    (
        "origin",
        maxRadialDistance * (transformationMatrix & radialVec)
    );
    radialPlane.add("scalingDistance", maxRadialDistance - minRadialDistance);
    radialPlane.add("scalingFactor", radialScaling_);
    radialPlane.add("gradingFactor", radialGrading_);
    radialPlane.add("scaleFromZero", 1);

    dict.add("radialGrading", radialPlane, true);

    //- define a dictionary for the axial direction
    if( contactWidth_ < 0.0 )
    {
        dictionary planeZ("axialGrading");
        planeZ.add("type", "plane");
        planeZ.add("normal", (transformationMatrix & axialVec));
        planeZ.add("origin", vector::zero);
        planeZ.add("scalingDistance", axialDistance);
        planeZ.add("scalingFactor", axialScaling_);
        planeZ.add("gradingFactor", axialGrading_);
        planeZ.add("symmetricScaling", 1);
        planeZ.add("scaleFromZero", 1);

        dict.add("axialGrading", planeZ, true);
    }
    else
    {
        const scalar d = 0.5 * contactWidth_;

        dictionary planeZ("axialGrading1");
        planeZ.add("type", "plane");
        planeZ.add("normal", (transformationMatrix & axialVec));
        planeZ.add("origin", d * (transformationMatrix & axialVec));
        planeZ.add("scalingDistance", mag(axialDistance-d));
        planeZ.add("scalingFactor", axialScaling_);
        planeZ.add("gradingFactor", axialGrading_);
        planeZ.add("symmetricScaling", 0);
        planeZ.add("scaleFromZero", 1);

        dict.add("axialGrading1", planeZ, true);

        dictionary planeZNeg("axialGrading2");
        planeZNeg.add("type", "plane");
        planeZNeg.add("normal", -1.0 * (transformationMatrix & axialVec));
        planeZNeg.add("origin", -d * (transformationMatrix & axialVec));
        planeZNeg.add("scalingDistance", mag(axialDistance-d));
        planeZNeg.add("scalingFactor", axialScaling_);
        planeZNeg.add("gradingFactor", axialGrading_);
        planeZNeg.add("symmetricScaling", 0);
        planeZNeg.add("scaleFromZero", 1);

        dict.add("axialGrading2", planeZNeg, true);
    }

    return dict;
}

void rollGeometryInfo::updateSurfacePatches
(
    const wordList& patchNames,
    const wordList& patchTypes,
    const labelLongList& patchIndices
) const
{
    if( !surfPtr_ )
    {
        FatalError << "Cannot update surface patches because the surface"
            << " does not yet exist" << abort(FatalError);
    }

    if( surfPtr_->size() != patchIndices.size() )
    {
        FatalError << "Inconsistent number of triangles and the surface mesh"
            << " and the specified patch indices" << abort(FatalError);
    }

    if( patchNames.size() != patchTypes.size() )
    {
        FatalError << "Sizes of patchNames and patchTypes are not equal"
            << abort(FatalError);
    }

    triSurfModifier sMod(*surfPtr_);

    geometricSurfacePatchList& patches = sMod.patchesAccess();

    patches.setSize(patchNames.size());

    forAll(patches, patchI)
    {
        patches[patchI].name() = patchNames[patchI];
        patches[patchI].geometricType() = patchTypes[patchI];
    }

    LongList<labelledTri>& triangles = sMod.facetsAccess();

    forAll(triangles, triI)
    {
        triangles[triI].region() = patchIndices[triI];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void rollingMillMesh::generateRollerMeshes()
{
    const PtrList<rollGeometryInfo>& geometries = geomHandler_.rollPositions();

    const bool singleRollerContactPatch =
        geomHandler_.singleRollerContactPatch();

    forAll(geometries, geomI)
    {
        const rollGeometryInfo& geomInfo = geometries[geomI];

        Info << "Generating volume mesh for roll "
             << geomInfo.rollPosition() << endl;

        //- get the surface transformed to the x-y coordinates for 2D meshing
        const triSurf* surfPtr = geomInfo.transformedSurface();

        //- set the anisotropic meshing
        meshDict_.set("anisotropicSources", geomInfo.anisotropic2DMeshing());

        //- generate a 2D mesh of a cross section
        HashSet<const word> contactNames;
        contactNames.insert
        (
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERCONTACT
            )
        );
        crossSectionMeshGenerator sectionMesh
        (
            surfPtr,
            meshDict_,
            geomInfo.rollDict(),
            db_,
            regionName_,
            timeStep_,
            contactNames
        );

        polyMeshGen* meshPtr = sectionMesh.meshPtr();

        //- the surface is not needed any more
        deleteDemandDrivenData(surfPtr);

        //- tranform the mesh back to the original coodinates
        geomInfo.backTransform2DMeshToPosition(*meshPtr);

        Info << "Generating a 3D mesh from the cross-section mesh" << endl;
        revolve2DMesh meshRevolver(*meshPtr);

        //- set the patch at the zero plane
        meshRevolver.setRevolvingPatch("bottomEmptyFaces");

        //- get the circumferential resolution
        meshRevolver.setCircumResolution(geomInfo.circumResolution());

        //- set the grading in the circumferential direction
        meshRevolver.setMaxRatio(geomInfo.circumGrading());

        //- set the origin and the rotation axis
        meshRevolver.setOrigin(geomInfo.origin());
        meshRevolver.setRotationAxis(geomInfo.rotationAxis());

        //- set the local refinement in the contact region
        const scalar contactLength = geomHandler_.contactLength();
        const scalar outerRollDiameter = geomInfo.outerDiameter();

        scalar angle = 2.0 * contactLength / (outerRollDiameter + VSMALL);
        angle = min(angle, 2.0 * M_PI);

        const scalar geometryTolerance = geomHandler_.geometryTolerance();

        const scalar angleStep =
            2.0 * geomInfo.circumScaling() *
            Foam::acos(1.0 - geometryTolerance / (0.5*outerRollDiameter));

        meshRevolver.setIntervalResolution(0.0, angle, angleStep);

        //- rename contact patches
        const scalar angleTol = geomInfo.extraContactAngle() * M_PI / 180.0;
        meshRevolver.setIntervalPatches
        (
            angle+angleTol,
            2.0 * M_PI-angleTol,
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERCONTACT
            ),
            patchHandler_.patchNameForRoll
            (
                rollingMillPatchNamesHandler::ROLLERTOAIR
            )
        );

        if( geomHandler_.areRollsClosed() )
        {
            meshRevolver.setIntervalPatches
            (
                angle+angleTol,
                2.0 * M_PI-angleTol,
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERTOROLLERBACK
                ),
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERTOAIR
                )
            );

            meshRevolver.setIntervalPatches
            (
                angle+angleTol,
                2.0 * M_PI-angleTol,
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERTOROLLERFRONT
                ),
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERTOAIR
                )
            );
        }

        //- create a 3D mesh
        meshRevolver.generateRevolvedMesh();

        //- translate the mesh to the position when the rolls are open
        Info << "Initial clearance "
             << geomHandler_.initialRollGapClearance() << endl;
        const vector translationVector =
            geomInfo.radialShift
            (
                geomHandler_.initialRollGapClearance() > -SMALL ?
                geomHandler_.initialRollGapClearance() :
                0.25 * geomInfo.outerDiameter()
            );
        const scalar offsetX = geomInfo.offsetX();

        polyMeshGenModifier meshModifier(*meshPtr);
        pointFieldPMG& points = meshModifier.pointsAccess();

        scalar minY(VGREAT), maxY(-VGREAT), minZ(VGREAT), maxZ(-VGREAT);
        forAll(points, pI)
        {
            point& p = points[pI];

            p += translationVector;

            p.x() += offsetX;

            minY = min(p.y(), minY);
            maxY = max(p.y(), maxY);
            minZ = min(p.z(), minZ);
            maxZ = max(p.z(), maxZ);
        }

        if
        (
            geomInfo.rollPosition() == "topRoll" ||
            geomInfo.rollPosition() == "bottomRoll"
        )
        {
            Info << geomInfo.rollPosition()
                 << " diameter " << (maxY-minY)
                 << " width " << (maxZ-minZ) << endl;
        }
        else if
        (
            geomInfo.rollPosition() == "leftRoll" ||
            geomInfo.rollPosition() == "rightRoll"
        )
        {
            Info << geomInfo.rollPosition()
                 << " diameter " << (maxZ-minZ)
                 << " width " << (maxY-minY) << endl;
        }

        //- add cells to the zone
        Info << "Generating zones" << endl;
        const word rollerName = geomInfo.rollPosition()+geomInfo.rollOffset();

        label cId;
        if( rollerName != "" )
        {
            cId = meshPtr->addCellZone(rollerName);
        }
        else
        {
            cId = meshPtr->addCellZone("roller");
        }

        forAll(meshPtr->cells(), cellI)
            meshPtr->addCellToZone(cId, cellI);

        //- set the patch names and types
        PtrList<boundaryPatch>& boundaries = meshModifier.boundariesAccess();

        forAll(boundaries, patchI)
        {
            boundaryPatch& patch = boundaries[patchI];

            if
            (
                patch.patchName() ==
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERCONTACT
                )
            )
            {
                label start = patch.patchStart();
                const label size = patch.patchSize();

                //- add faces in the contact patch to a zone
                label fId;
                if( (rollerName != "") && !singleRollerContactPatch )
                {
                    fId = meshPtr->addFaceZone(rollerName+patch.patchName());
                    patch.patchName() = rollerName + patch.patchName();
                }
                else if( singleRollerContactPatch )
                {
                    const word pName =
                        patchHandler_.patchNameForRoll
                        (
                            rollingMillPatchNamesHandler::ROLLERCONTACT
                        );

                    fId = meshPtr->addFaceZone(pName);
                    patch.patchName() = pName;
                }
                else
                {
                    fId = meshPtr->addFaceZone(patch.patchName());
                }

                for(label i=0;i<size;++i)
                    meshPtr->addFaceToZone(fId, start++);
            }
            else if
            (
                patch.patchName() ==
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERTOAIR
                )
            )
            {
                label start = patch.patchStart();
                const label size = patch.patchSize();

                //- add faces in the rollerToAir patch to a zone
                label fId;
                if( (rollerName != "") && !singleRollerContactPatch )
                {
                    fId = meshPtr->addFaceZone(rollerName+patch.patchName());
                    patch.patchName() = rollerName + patch.patchName();
                }
                else if( singleRollerContactPatch )
                {
                    const word pName =
                        patchHandler_.patchNameForRoll
                        (
                            rollingMillPatchNamesHandler::ROLLERTOAIR
                        );

                    fId = meshPtr->addFaceZone(pName);
                    patch.patchName() = pName;
                }
                else
                {
                    fId = meshPtr->addFaceZone(patch.patchName());
                }

                for(label i=0;i<size;++i)
                    meshPtr->addFaceToZone(fId, start++);
            }
            else if
            (
                patch.patchName() !=
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERTOAIR
                )
            )
            {
                patch.patchName() += geomInfo.rollOffset();

                //- add roller name as a prefix
                if( rollerName != "" )
                {
                    patch.patchName() = rollerName + patch.patchName();
                }
            }
        }

        if( globalMeshPtr_ )
        {
            //- add a new roll into the existing mesh
            polyMeshGenModifier(*globalMeshPtr_).addMesh(*meshPtr);

            deleteDemandDrivenData(meshPtr);
        }
        else
        {
            //- insert the mesh as a global mesh
            globalMeshPtr_ = meshPtr;
        }
    }

    Info << "Finished generating a 3D mesh from the cross-section mesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
