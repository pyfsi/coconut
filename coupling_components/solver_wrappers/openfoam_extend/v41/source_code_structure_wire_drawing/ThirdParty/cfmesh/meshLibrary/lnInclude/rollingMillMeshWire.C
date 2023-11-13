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
#include "triSurfModifier.H"
//#include "meshSurfaceDistanceFromGeometry.H"
//#include "meshOctreeModifier.H"
//#include "meshOctreeInsideOutside.H"
#include "extrude2DMesh.H"
#include "wireBlockMeshGenerator.H"

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const tensor wireGeometryInfo::transformationMatrix =
    tensor(0., 0., 1.0, 0., 1.0, 0., -1.0, 0., 0.);

wireGeometryInfo::wireGeometryInfo
(
    const dictionary& wireDict,
    const triSurf* surfacePtr
)
:
    surfPtr_(surfacePtr),
    wireDictPtr_(&wireDict),
    wireDiameter_(-1.0),
    axialGrading_(1.0),
    axialShift_(vector::zero),
    axialResolution_(20),
    isRollingPass_(false),
    isDrawingPass_(false),
    wallPatchName_(),
    symmetryType_(),
    wedgeAngle_(2.5),
    isCoatingPresent_(false)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

wireGeometryInfo::~wireGeometryInfo()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool wireGeometryInfo::isCrossSectionCircular() const
{
    return !surfPtr_;
}

const triSurf& wireGeometryInfo::surface() const
{
    if( !surfPtr_ )
    {
        FatalErrorIn
        (
            "const triSurf& wireGeometryInfo::surface() const"
        ) << "Cross section is circular and surface is not specified"
          << abort(FatalError);
    }

    return *surfPtr_;
}

const dictionary& wireGeometryInfo::wireDict() const
{
    return *wireDictPtr_;
}

scalar& wireGeometryInfo::wireDiameter()
{
    return wireDiameter_;
}

const scalar& wireGeometryInfo::wireDiameter() const
{
    return wireDiameter_;
}

const scalar& wireGeometryInfo::axialGrading() const
{
    return axialGrading_;
}

scalar& wireGeometryInfo::axialGrading()
{
    return axialGrading_;
}

label wireGeometryInfo::axialResolution() const
{
    return axialResolution_;
}

label& wireGeometryInfo::axialResolution()
{
    return axialResolution_;
}

const vector& wireGeometryInfo::axialShift() const
{
    return axialShift_;
}

vector& wireGeometryInfo::axialShift()
{
    return axialShift_;
}

const word& wireGeometryInfo::wallPatchName() const
{
    return wallPatchName_;
}

word& wireGeometryInfo::wallPatchName()
{
    return wallPatchName_;
}

direction wireGeometryInfo::typeOfSymmetry() const
{
    return symmetryType_;
}

direction& wireGeometryInfo::typeOfSymmetry()
{
    return symmetryType_;
}

scalar wireGeometryInfo::wedgeAngle() const
{
    return wedgeAngle_;
}

scalar& wireGeometryInfo::wedgeAngle()
{
    return wedgeAngle_;
}

bool wireGeometryInfo::isCoatingPresent() const
{
    return isCoatingPresent_;
}

bool& wireGeometryInfo::isCoatingPresent()
{
    return isCoatingPresent_;
}

const triSurf* wireGeometryInfo::transformedSurface() const
{
    triSurf* surfPtr = new triSurf(*surfPtr_);

    //- tranform the points into the x-y plane
    triSurfModifier sMod(*surfPtr);
    pointField& pts = sMod.pointsAccess();

    forAll(pts, pI)
        pts[pI] = (transformationMatrix & pts[pI]);

    return surfPtr;
}

void wireGeometryInfo::backTransform2DMeshToPosition(polyMeshGen& mesh) const
{
    polyMeshGenModifier meshModifier(mesh);
    pointFieldPMG& points = meshModifier.pointsAccess();

    //- calculate the inverse tranformation matrix
    const tensor backwardTransformation = inv(transformationMatrix);

    //- transform the points back to the original coordinate system
    forAll(points, pI)
        points[pI] = (backwardTransformation & points[pI]);
}

// * * * * * * * * * * * * Private member functions  * * * * * * * * * * * * //

void rollingMillMesh::generateWireMesh()
{
    const wireGeometryInfo& geomInfo = geomHandler_.wireGeometry();

    //- create a 2D cross section of the wire
    polyMeshGen* meshPtr = NULL;

    if( geomInfo.isCrossSectionCircular() )
    {
        //- block-structured wire topology in case of circular cross-section
        wireBlockMeshGenerator wbm
        (
            db_,
            patchHandler_,
            geomInfo,
            regionName_,
            geomHandler_.geometryTolerance()
        );

        meshPtr = wbm.generateCrossSectionMesh();
    }
    else
    {
        //- get the wire surface transformed into the x-y coordinates
        const triSurf* surfPtr = geomInfo.transformedSurface();

        //- generate a 2D mesh of a cross section
        HashSet<const word> contactNames;
        contactNames.insert
        (
            patchHandler_.patchNameForWire
            (
                rollingMillPatchNamesHandler::WIRECONTACT
            )
        );
        crossSectionMeshGenerator sectionMesh
        (
            surfPtr,
            meshDict_,
            geomInfo.wireDict(),
            db_,
            regionName_,
            timeStep_,
            contactNames,
            true
        );

        meshPtr = sectionMesh.meshPtr();

        //- the surface is not needed any more
        deleteDemandDrivenData(surfPtr);
    }

    //- tranform the mesh back to the original coordinates
    geomInfo.backTransform2DMeshToPosition(*meshPtr);

    //- extrude the mesh in the x direction
    extrude2DMesh extruder(*meshPtr);

    extruder.setExtrusionPatch("bottomEmptyFaces");

    //- set the number of subdivisions
    extruder.setNumberOfSubdivisions(geomInfo.axialResolution());

    //- set the extrusion vector
    extruder.setExtrusionVector
    (
        vector(-geomHandler_.wireLength(), 0., 0.)
    );

    //- set the grading factor
    extruder.setGradingFactor(geomInfo.axialGrading());

    //- generate the extruded mesh
    extruder.generateExtrudedMesh();

    //- translate the mesh a bit
    polyMeshGenModifier meshModifier(*meshPtr);
    pointFieldPMG& points = meshModifier.pointsAccess();

    point minX(VGREAT, VGREAT, VGREAT), maxX(-VGREAT, -VGREAT, -VGREAT);
    forAll(points, pI)
    {
        minX = min(minX, points[pI]);
        maxX = max(maxX, points[pI]);
    }

    //- position the middle of the wire at x = 0
    vector translationVec = -0.5 * vector(maxX.x() + minX.x(), 0., 0.);
    translationVec += geomInfo.axialShift();

    forAll(points, pI)
        points[pI] += translationVec;

    //- rename patches
    forAll(meshModifier.boundariesAccess(), patchI)
    {
        boundaryPatch& patch = meshModifier.boundariesAccess()[patchI];

        if( patch.patchName() == "bottomEmptyFacesCopy" )
        {
            //- correct name for the upstream patch
            patch.patchName() =
                patchHandler_.patchNameForWire
                (
                    rollingMillPatchNamesHandler::WIREUPSTREAM
                );
            patch.patchType() =
                patchHandler_.patchTypeForWire
                (
                    rollingMillPatchNamesHandler::WIREUPSTREAM,
                    rollingMillGeometryHandler::symmetryTypes_
                    (
                        geomInfo.typeOfSymmetry()
                    )
                );

            const label upstreamId = meshPtr->addFaceZone("wireInletFaces");
            const label sId = meshPtr->addFaceSubset("wireInletFaces");
            for(label fI=0;fI<patch.patchSize();++fI)
            {
                meshPtr->addFaceToZone(upstreamId, patch.patchStart()+fI);
                meshPtr->addFaceToSubset(sId, patch.patchStart()+fI);
            }
        }
        else if( patch.patchName() == "bottomEmptyFaces" )
        {
            //- correct patch name for the downstream patch
            patch.patchName() =
                patchHandler_.patchNameForWire
                (
                    rollingMillPatchNamesHandler::WIREDOWNSTREAM
                );
            patch.patchType() =
                patchHandler_.patchTypeForWire
                (
                    rollingMillPatchNamesHandler::WIREDOWNSTREAM,
                    rollingMillGeometryHandler::symmetryTypes_
                    (
                        geomInfo.typeOfSymmetry()
                    )
                );

            const label downId = meshPtr->addFaceZone("wireOutletFaces");
            const label downSetId = meshPtr->addFaceSubset("wireOutletFaces");
            for(label fI=0;fI<patch.patchSize();++fI)
            {
                meshPtr->addFaceToZone(downId, patch.patchStart()+fI);
                meshPtr->addFaceToSubset(downSetId, patch.patchStart()+fI);
            }
        }
        else if
        (
            patch.patchName() ==
            patchHandler_.patchNameForWire
            (
                rollingMillPatchNamesHandler::WIRESYMMY
            )
        )
        {
            patch.patchType() =
                patchHandler_.patchTypeForWire
                (
                    rollingMillPatchNamesHandler::WIRESYMMY,
                    rollingMillGeometryHandler::symmetryTypes_
                    (
                        geomInfo.typeOfSymmetry()
                    )
                );
        }
        else if
        (
            patch.patchName() ==
            patchHandler_.patchNameForWire
            (
                rollingMillPatchNamesHandler::WIRESYMMZ
            )
        )
        {
            patch.patchType() =
                patchHandler_.patchTypeForWire
                (
                    rollingMillPatchNamesHandler::WIRESYMMZ,
                    rollingMillGeometryHandler::symmetryTypes_
                    (
                        geomInfo.typeOfSymmetry()
                    )
                );
        }
        else if
        (
            patch.patchName() ==
            patchHandler_.patchNameForWire
            (
                rollingMillPatchNamesHandler::WIRECONTACT
            )
        )
        {
            //- contact patch name
            label fzId(-1);

            //- there is not coating on the wire
            fzId = meshPtr->addFaceZone(geomInfo.wallPatchName());

            //- create a face zone
            const label start = patch.patchStart();
            const label end = start + patch.patchSize();

            for(label fI=start;fI<end;++fI)
                meshPtr->addFaceToZone(fzId, fI);
        }
    }

    //- add cells to zone wire
    const label cId = meshPtr->addCellZone("wire");
    const label wimId = meshPtr->addCellSubset("wireInletFacesMasterCells");
    const label womId = meshPtr->addCellSubset("wireOutletFacesMasterCells");

    forAll(meshPtr->cells(), cellI)
    {
        meshPtr->addCellToZone(cId, cellI);
        meshPtr->addCellToSubset(wimId, cellI);
        meshPtr->addCellToSubset(womId, cellI);
    }

    globalMeshPtr_ = meshPtr;
}

void rollingMillMesh::generateWireCoatingMesh()
{
    const wireGeometryInfo& geomInfo = geomHandler_.wireGeometry();

    //- check if the coating exists
    if( geomInfo.isCoatingPresent() )
    {
        const dictionary& wireDict = geomInfo.wireDict();

        if( wireDict.found("coatingDict") )
        {
            const dictionary& coatingDict = wireDict.subDict("coatingDict");

            //- thickness of the coating
            scalar coatingThickness(1e-6);
            if
            (
                !coatingDict.readIfPresent("coatingThickness", coatingThickness)
            )
            {
                FatalError << "coatingThickness is not specified in coatingDict"
                    << exit(FatalError);
            }

            //- number of layers in the coating mesh
            label nLayers(1);
            coatingDict.readIfPresent("nSubdivisions", nLayers);

            //- thickness ratio
            scalar thicknessRatio(1.0);
            coatingDict.readIfPresent("thicknessRatio", thicknessRatio);

            //- copy the wire mesh into a new mesh
            polyMeshGen* meshPtr =
                new polyMeshGen(globalMeshPtr_->returnTime());

            polyMeshGenModifier meshModifier(*meshPtr);
            meshModifier.addMesh(*globalMeshPtr_);

            //- extrude coating mesh
            extrude2DMesh coatingExtruder(*meshPtr);
            coatingExtruder.setExtrusionPatch(geomInfo.wallPatchName());
            coatingExtruder.extrudeInPatchNormalDirection(true);
            coatingExtruder.setNumberOfSubdivisions(nLayers);
            coatingExtruder.extrusionLength(coatingThickness);
            coatingExtruder.setGradingFactor(thicknessRatio);
            coatingExtruder.removeOriginalCells(true);
            coatingExtruder.setCellSubsetName("coating");

            coatingExtruder.generateExtrudedMesh();

            //- rename patches in the extruded mesh
            //- patch types remain the same
            const word upstreamName =
                patchHandler_.patchNameForCoating
                (
                    rollingMillPatchNamesHandler::COATINGUPSTREAM
                );
            const word downstreamName =
                patchHandler_.patchNameForCoating
                (
                    rollingMillPatchNamesHandler::COATINGDOWNSTREAM
                );
            forAll(meshModifier.boundariesAccess(), patchI)
            {
                boundaryPatch& patch = meshModifier.boundariesAccess()[patchI];

                //- check all possible patch names and select
                //- the appropriate one. Patch names for wire and coating
                //- MUST have the same ordering in rollingMillPatchNamesHandler
                for(label i=0;i<8;++i)
                {
                    word pName =
                        patchHandler_.patchNameForWire
                        (
                            rollingMillPatchNamesHandler::wirePatchKeys(i)
                        );

                    if( patch.patchName().find(pName) != word::npos )
                    {
                        word newName =
                            patchHandler_.patchNameForCoating
                            (
                                rollingMillPatchNamesHandler::coatingPatchKeys
                                (
                                    i
                                )
                            );

                        if( patch.patchName().find("Copy") != word::npos )
                        {
                            //- this patch is contact between the coating
                            //- and the wire
                            newName =
                                patchHandler_.patchNameForCoating
                                (
                                    rollingMillPatchNamesHandler::COATINGWIRE
                                );
                        }

                        patch.patchName() = newName;

                        if( newName == upstreamName )
                        {
                            //- create coatingInletFaces zone
                            const label fId =
                                meshPtr->addFaceZone("coatingInletFaces");
                            const label sfId =
                                meshPtr->addFaceSubset("coatingInletFaces");

                            const label start = patch.patchStart();
                            for(label fI=0;fI<patch.patchSize();++fI)
                            {
                                meshPtr->addFaceToZone(fId, start+fI);
                                meshPtr->addFaceToSubset(sfId, start+fI);
                            }
                        }
                        else if( newName == downstreamName )
                        {
                            //- create coatingOutletFaces zone
                            const label fId =
                                meshPtr->addFaceZone("coatingOutletFaces");
                            const label sfId =
                                meshPtr->addFaceSubset("coatingOutletFaces");

                            const label start = patch.patchStart();
                            for(label fI=0;fI<patch.patchSize();++fI)
                            {
                                meshPtr->addFaceToZone(fId, start+fI);
                                meshPtr->addFaceToSubset(sfId, start+fI);
                            }
                        }
                    }
                }
            }

            //- add cells to zone coating
            const label cId = meshPtr->addCellZone("coating");
            const label scId =
                meshPtr->addCellSubset("coatingInletFacesMasterCells");
            const label scId1 =
                meshPtr->addCellSubset("coatingOutletFacesMasterCells");
            forAll(meshPtr->cells(), cellI)
            {
                meshPtr->addCellToZone(cId, cellI);
                meshPtr->addCellToSubset(scId, cellI);
                meshPtr->addCellToSubset(scId1, cellI);
            }

            //- rename the contact patch in the wire mesh
            polyMeshGenModifier globalMeshModifier(*globalMeshPtr_);

            forAll(globalMeshModifier.boundariesAccess(), patchI)
            {
                boundaryPatch& patch =
                    globalMeshModifier.boundariesAccess()[patchI];

                if
                (
                    patch.patchName() ==
                    patchHandler_.patchNameForWire
                    (
                        rollingMillPatchNamesHandler::WIRECONTACT
                    )
                )
                {
                    //- change the name associated with WIRECOATING
                    patch.patchName() =
                        patchHandler_.patchNameForWire
                        (
                            rollingMillPatchNamesHandler::WIRECOATING
                        );
                }
            }

            polyMeshGenModifier(*globalMeshPtr_).addMesh(*meshPtr);
            deleteDemandDrivenData(meshPtr);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
