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
#include "meshSurfaceEngine.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "triSurfacePatchManipulator.H"
#include "casingAxialCrossSectionSurfaceCreator.H"
#include "dieCrossSectionSurfaceCreator.H"
#include "dieSurfaceCreator.H"
#include "profiledDieGeometryInterpolator.H"
#include "profiledDieMeshGenerator.H"
#include "revolve2DMesh.H"

#include "helperFunctions.H"

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dieGeometryInfo::dieGeometryInfo
(
    const dictionary& dieDict,
    const rollingMillGeometryHandler& geomHandler,
    const rollingMillPatchNamesHandler& patchHandler
)
:
    dieDictPtr_(&dieDict),
    axialLength_(),
    outerRadius_(),
    inletRadius_(),
    outletRadius_(),
    geometryTolerance_(geomHandler.geometryTolerance()),
    wedgeAngle_(2.5),
    innerSurfPtr_(NULL),
    inletSurfPtr_(NULL),
    outletSurfPtr_(NULL),
    casingSurfPtr_(NULL),
    geomHandler_(geomHandler),
    patchHandler_(patchHandler)
{
    //- check/read geometry tolerance
    if( dieDict.found("geometryTolerance") )
    {
        geometryTolerance_ =
            readScalar(dieDict.lookup("geometryTolerance"));
    }

    if( true )
    {
        autoPtr<dieSurfaceCreator> dieSurfCreator =
            dieSurfaceCreator::New
            (
                patchHandler_,
                dieDict,
                geometryTolerance_
            );

        innerSurfPtr_ = new triSurf(dieSurfCreator->surface());
        inletRadius_ = 0.5 * dieSurfCreator->inletDiameter();
        outletRadius_ = 0.5 * dieSurfCreator->outletDiameter();
        outerRadius_ = 0.5 * dieSurfCreator->outerDiameter();

        //- transfer information about profiled dies
        const DynList<std::pair<scalar, std::shared_ptr<triSurf> > >& sections =
            dieSurfCreator->crossSections();

        crossSections_.setSize(sections.size());
        forAll(sections, i)
        {
            crossSections_.set
            (
                i,
                new std::pair<scalar, std::shared_ptr<triSurf> >
                (
                    sections[i].first,
                    sections[i].second
                )
            );
        }
    }

    //- read the information about the casing
    Switch hasCasing(false);
    if( dieDict.found("casing") )
    {
        hasCasing.readIfPresent("casing", dieDict);

        if( hasCasing )
        {
            const dictionary& casingDict = dieDict.subDict("casingDict");

            //- types of casing (cylindrical)
            autoPtr<casingAxialCrossSectionSurfaceCreator> casingSurfPtr =
                casingAxialCrossSectionSurfaceCreator::New
                (
                    patchHandler_,
                    casingDict,
                    geometryTolerance_
                );

            casingSurfPtr_ = new triSurf(casingSurfPtr->surface());
        }
    }

    //- wedge angle for axi-symmetry
    if( geomHandler_.symmetryType() & rollingMillGeometryHandler::AXISYMMETRIC )
    {
        if( dieDict.found("wedgeAngle") )
        {
            wedgeAngle_ = readScalar(dieDict.lookup("wedgeAngle"));
        }
    }

    if( outletSurfPtr_ )
    {
        //- get the average radius of the outlet profile
        const pointField& pts = outletSurfPtr_->points();
        scalar maxRadius(0.0), minRadius(VGREAT);
        forAll(pts, pI)
        {
            point p = pts[pI];
            p.x() = 0.0;

            scalar r = mag(p);

            maxRadius = max(maxRadius, r);
            minRadius = min(minRadius, r);
        }

        outletRadius_ = 0.5 * (maxRadius + minRadius);
    }

    if( inletSurfPtr_ )
    {
        //- get the average radius of the inlet profile
        const pointField& pts = inletSurfPtr_->points();
        scalar maxRadius(0.0), minRadius(VGREAT);
        forAll(pts, pI)
        {
            point p = pts[pI];
            p.x() = 0.0;

            scalar r = mag(p);

            maxRadius = max(maxRadius, r);
            minRadius = min(minRadius, r);
        }

        inletRadius_ = 0.5 * (maxRadius + minRadius);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dieGeometryInfo::~dieGeometryInfo()
{
    deleteDemandDrivenData(inletSurfPtr_);
    deleteDemandDrivenData(outletSurfPtr_);
    deleteDemandDrivenData(innerSurfPtr_);
    deleteDemandDrivenData(casingSurfPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const dictionary& dieGeometryInfo::dieDict() const
{
    return *dieDictPtr_;
}

scalar dieGeometryInfo::outerDiameter() const
{
    return 2.0 * outerRadius_;
}

scalar dieGeometryInfo::inletDiameter() const
{
    return 2.0 * inletRadius_;
}

bool dieGeometryInfo::isInletCircular() const
{
    return !inletSurfPtr_;
}

bool dieGeometryInfo::isOutletCircular() const
{
    return !outletSurfPtr_;
}

scalar dieGeometryInfo::outletDiameter() const
{
    return 2.0 * outletRadius_;
}

bool dieGeometryInfo::isProfiled() const
{
    return (crossSections_.size() != 0);
}

bool dieGeometryInfo::isSymmetric() const
{
    return geomHandler_.symmetryType();
}

bool dieGeometryInfo::isCasingPresent() const
{
    return (casingSurfPtr_ != NULL);
}

direction dieGeometryInfo::typeOfSymmetry() const
{
    return geomHandler_.symmetryType();
}

bool dieGeometryInfo::isWedge() const
{
    return
    (
        geomHandler_.symmetryType() &
        rollingMillGeometryHandler::AXISYMMETRIC
    );
}

scalar dieGeometryInfo::wedgeAngle() const
{
    return wedgeAngle_;
}

scalar dieGeometryInfo::startCircumAngle() const
{
    const direction sType = geomHandler_.symmetryType();

    if
    (
        sType == rollingMillGeometryHandler::POSY ||
        (
            (sType & rollingMillGeometryHandler::POSY) &&
            (sType & rollingMillGeometryHandler::NEGZ)
        )
    )
    {
        return 1.5 * M_PI;
    }
    else if
    (
        sType == rollingMillGeometryHandler::POSZ ||
        (
            (sType & rollingMillGeometryHandler::POSY) &&
            (sType & rollingMillGeometryHandler::POSZ)
        )
    )
    {
        return 0.0;
    }
    else if
    (
        sType == rollingMillGeometryHandler::NEGY ||
        (
            (sType & rollingMillGeometryHandler::NEGY) &&
            (sType & rollingMillGeometryHandler::POSZ)
        )
    )
    {
        return 0.5 * M_PI;
    }
    else if
    (
        sType == rollingMillGeometryHandler::NEGZ ||
        (
            (sType & rollingMillGeometryHandler::NEGY) &&
            (sType & rollingMillGeometryHandler::NEGZ)
        )
    )
    {
        return M_PI;
    }

    return 0.0;
}

scalar dieGeometryInfo::endCircumAngle() const
{
    const direction sType = geomHandler_.symmetryType();

    if
    (
        sType == rollingMillGeometryHandler::POSY ||
        (
            (sType & rollingMillGeometryHandler::POSY) &&
            (sType & rollingMillGeometryHandler::POSZ)
        )
    )
    {
        return 0.5 * M_PI;
    }
    else if
    (
        sType == rollingMillGeometryHandler::NEGY ||
        (
            (sType & rollingMillGeometryHandler::NEGY) &&
            (sType & rollingMillGeometryHandler::NEGZ)
        )
    )
    {
        return 1.5 * M_PI;
    }
    else if
    (
        sType == rollingMillGeometryHandler::POSZ ||
        (
            (sType & rollingMillGeometryHandler::NEGY) &&
            (sType & rollingMillGeometryHandler::POSZ)
        )
    )
    {
        return M_PI;
    }
    else if
    (
        sType == rollingMillGeometryHandler::NEGZ ||
        (
            (sType & rollingMillGeometryHandler::POSY) &&
            (sType & rollingMillGeometryHandler::NEGZ)
        )
    )
    {
        return 0.0;
    }

    return 2.0 * M_PI;
}

const triSurf& dieGeometryInfo::inletSurface() const
{
    if( !inletSurfPtr_ )
    {
        FatalErrorIn
        (
            "const triSurf& dieGeometryInfo::inletSurface() const"
        ) << "Inlet surface does not exist" << exit(FatalError);
    }

    return *inletSurfPtr_;
}

const triSurf& dieGeometryInfo::outletSurface() const
{
    if( !outletSurfPtr_ )
        FatalErrorIn
        (
            "const triSurf& dieGeometryInfo::outletSurface() const"
        ) << "Outlet surface does not exist" << exit(FatalError);

    return *outletSurfPtr_;
}

const PtrList<std::pair<scalar, std::shared_ptr<triSurf> > >&
dieGeometryInfo::crossSections() const
{
    return crossSections_;
}

const triSurf* dieGeometryInfo::axialCrossSectionSurface() const
{
    if( !innerSurfPtr_ )
    {
        FatalErrorIn
        (
            "const triSurf& dieGeometryInfo::axialCrossSectionSurface() const"
        ) << "die surface does not exist" << exit(FatalError);
    }

    const triSurf* newPtr = new triSurf(*innerSurfPtr_);

    return newPtr;
}

const triSurf* dieGeometryInfo::casingCrossSectionSurface() const
{
    if( !casingSurfPtr_ )
    {
        FatalErrorIn
        (
            "const triSurf* dieGeometryInfo::casingCrossSectionSurface() const"
        ) << "die surface does not exist" << exit(FatalError);
    }

    const triSurf* newPtr = new triSurf(*casingSurfPtr_);

    return newPtr;
}

void rollingMillMesh::generateDieMesh()
{
    const dieGeometryInfo& dieGeom = geomHandler_.dieGeometry();

    const triSurf* surfPtr = dieGeom.axialCrossSectionSurface();

    //- generate a 2D mesh of a cross section
    dictionary dieDict = dieGeom.dieDict();
    if( dieDict.found("casing") )
        dieDict.remove("casing");
    if( dieDict.found("casingDict") )
        dieDict.remove("casingDict");

    //- write axial cross-section to a file
    if( dieDict.found("writeDieProfile") )
    {
        bool writeProfile(false);
        dieDict.readIfPresent("writeDieProfile", writeProfile);

        if( writeProfile )
        {
            surfPtr->writeSurface("dieAxialGeom.stl");
        }
    }

    HashSet<const word> contactNames;
    contactNames.insert
    (
        patchHandler_.patchNameForDie(rollingMillPatchNamesHandler::DIEWIRE)
    );
    crossSectionMeshGenerator sectionMesh
    (
        surfPtr,
        meshDict_,
        dieDict,
        db_,
        regionName_,
        timeStep_,
        contactNames
    );

    polyMeshGen* meshPtr = sectionMesh.meshPtr();

    //- the surface is not needed any more
    deleteDemandDrivenData(surfPtr);

    Info << "Generating a 3D mesh from the cross-section mesh" << endl;
    revolve2DMesh meshRevolver(*meshPtr);

    //- set the patch at the zero plane
    meshRevolver.setRevolvingPatch("bottomEmptyFaces");

    //- get the circumferential resolution
    const scalar inletDiameter = dieGeom.inletDiameter();
    const scalar geometryTolerance = geomHandler_.geometryTolerance();

    scalar angleStep =
        2.0 * Foam::acos(1.0 - geometryTolerance / (0.5*inletDiameter));

    if( dieDict.found("numCellsInCircumferentialDirection") )
    {
        label nCircumferentialDivisions(2);
        dieDict.readIfPresent
        (
            "numCellsInCircumferentialDirection",
            nCircumferentialDivisions
        );

        angleStep = (2.0 * M_PI / nCircumferentialDivisions);
    }

    //- check the revolution angle
    if( dieGeom.isSymmetric() )
    {
        meshRevolver.clearAngleIntervals();

        if( dieGeom.isWedge() )
        {
            //- axi-symmetric setup
            const scalar wedgeAngle = dieGeom.wedgeAngle() * M_PI / 360.0;

            if( angleStep < wedgeAngle )
            {
                Warning << "Wedge angle " << dieGeom.wedgeAngle()
                     << " is greater "
                     << "than " << (angleStep * 180.0 / M_PI)
                     << ", needed to achieve geometry tolerance " << endl;
            }

            angleStep = wedgeAngle;

            meshRevolver.setIntervalResolution(-angleStep, angleStep, 1);
        }
        else
        {
            //- symmetric setup
            meshRevolver.setIntervalResolution
            (
                dieGeom.startCircumAngle(),
                dieGeom.endCircumAngle(),
                angleStep
            );
        }
    }
    else
    {
        meshRevolver.setCircumResolution(angleStep);
    }

    //- set the origin and the rotation axis
    meshRevolver.setOrigin(vector(0, 0., 0.));
    meshRevolver.setRotationAxis(vector(1., 0., 0.));

    //- generate point subsets
    meshRevolver.createPointSubsets();

    //- create a 3D mesh
    meshRevolver.generateRevolvedMesh();

    //- add cells to the zone
    Info << "Generating zones" << endl;
    const label cId = meshPtr->addCellZone("die");
    forAll(meshPtr->cells(), cellI)
        meshPtr->addCellToZone(cId, cellI);

    //- add faces in the contact patch to a zone
    const word contactPatchName =
        patchHandler_.patchNameForDie(rollingMillPatchNamesHandler::DIEWIRE);
    const label fId = meshPtr->addFaceZone(contactPatchName);

    //- set the patch names and types
    PtrList<boundaryPatch>& boundaries =
        polyMeshGenModifier(*meshPtr).boundariesAccess();
    forAll(boundaries, patchI)
    {
        if( boundaries[patchI].patchName() == contactPatchName )
        {
            label start = boundaries[patchI].patchStart();
            const label size = boundaries[patchI].patchSize();

            for(label i=0;i<size;++i)
                meshPtr->addFaceToZone(fId, start++);
        }

        if( dieGeom.isSymmetric() )
        {
            if( boundaries[patchI].patchName() == "defaultFaces" )
            {
                if( dieGeom.isWedge() )
                {
                    boundaries[patchI].patchName() =
                        patchHandler_.patchNameForDie
                        (
                            rollingMillPatchNamesHandler::DIEFRONT
                        );
                    boundaries[patchI].patchType() =
                        patchHandler_.patchTypeForDie
                        (
                            rollingMillPatchNamesHandler::DIEFRONT,
                            rollingMillGeometryHandler::symmetryTypes_
                            (
                                dieGeom.typeOfSymmetry()
                            )
                        );
                }
                else
                {
                    switch( geomHandler_.symmetryType() )
                    {
                        case rollingMillGeometryHandler::NEGY:
                        case rollingMillGeometryHandler::POSY:
                        case
                        (
                            rollingMillGeometryHandler::POSY +
                            rollingMillGeometryHandler::POSZ
                        ):
                        case
                        (
                            rollingMillGeometryHandler::NEGY +
                            rollingMillGeometryHandler::NEGZ
                        ):
                        {
                            boundaries[patchI].patchName() =
                                patchHandler_.patchNameForDie
                                (
                                    rollingMillPatchNamesHandler::DIESYMMY
                                );
                            boundaries[patchI].patchType() =
                                patchHandler_.patchTypeForDie
                                (
                                    rollingMillPatchNamesHandler::DIESYMMY,
                                    rollingMillGeometryHandler::symmetryTypes_
                                    (
                                        dieGeom.typeOfSymmetry()
                                    )
                                );
                        } break;
                        case rollingMillGeometryHandler::NEGZ:
                        case rollingMillGeometryHandler::POSZ:
                        case
                        (
                            rollingMillGeometryHandler::NEGY +
                            rollingMillGeometryHandler::POSZ
                        ):
                        case
                        (
                            rollingMillGeometryHandler::NEGZ +
                            rollingMillGeometryHandler::POSY
                        ):
                        {
                            boundaries[patchI].patchName() =
                                patchHandler_.patchNameForDie
                                (
                                    rollingMillPatchNamesHandler::DIESYMMZ
                                );
                            boundaries[patchI].patchType() =
                                patchHandler_.patchTypeForDie
                                (
                                    rollingMillPatchNamesHandler::DIESYMMZ,
                                    rollingMillGeometryHandler::symmetryTypes_
                                    (
                                        dieGeom.typeOfSymmetry()
                                    )
                                );
                        } break;
                    }
                }
            }
            else if( boundaries[patchI].patchName() == "bottomEmptyFaces" )
            {
                if( dieGeom.isWedge() )
                {
                    boundaries[patchI].patchName() =
                        patchHandler_.patchNameForDie
                        (
                            rollingMillPatchNamesHandler::DIEBACK
                        );
                    boundaries[patchI].patchType() =
                        patchHandler_.patchTypeForDie
                        (
                            rollingMillPatchNamesHandler::DIEBACK,
                            rollingMillGeometryHandler::symmetryTypes_
                            (
                                dieGeom.typeOfSymmetry()
                            )
                        );
                }
                else
                {
                    switch( geomHandler_.symmetryType() )
                    {
                        case rollingMillGeometryHandler::NEGZ:
                        case rollingMillGeometryHandler::POSZ:
                        case
                        (
                            rollingMillGeometryHandler::POSY +
                            rollingMillGeometryHandler::POSZ
                        ):
                        case
                        (
                            rollingMillGeometryHandler::NEGY +
                            rollingMillGeometryHandler::NEGZ
                        ):
                        {
                            boundaries[patchI].patchName() =
                                patchHandler_.patchNameForDie
                                (
                                    rollingMillPatchNamesHandler::DIESYMMZ
                                );
                            boundaries[patchI].patchType() =
                                patchHandler_.patchTypeForDie
                                (
                                    rollingMillPatchNamesHandler::DIESYMMZ,
                                    rollingMillGeometryHandler::symmetryTypes_
                                    (
                                        dieGeom.typeOfSymmetry()
                                    )
                                );
                        } break;
                        case rollingMillGeometryHandler::NEGY:
                        case rollingMillGeometryHandler::POSY:
                        case
                        (
                            rollingMillGeometryHandler::NEGY +
                            rollingMillGeometryHandler::POSZ
                        ):
                        case
                        (
                            rollingMillGeometryHandler::NEGZ +
                            rollingMillGeometryHandler::POSY
                        ):
                        {
                            boundaries[patchI].patchName() =
                                patchHandler_.patchNameForDie
                                (
                                    rollingMillPatchNamesHandler::DIESYMMY
                                );
                            boundaries[patchI].patchType() =
                                patchHandler_.patchTypeForDie
                                (
                                    rollingMillPatchNamesHandler::DIESYMMY,
                                    rollingMillGeometryHandler::symmetryTypes_
                                    (
                                        dieGeom.typeOfSymmetry()
                                    )
                                );
                        } break;
                    }
                }
            }
        }
    }

//    if( !dieGeom.isOutletCircular() )
//    {
//        const triSurf& outletSurf = dieGeom.outletSurface();
//        const VRWGraph& sef = outletSurf.edgeFacets();

//        //- check if there exist any feature edges that shall be preserved
//        std::map<std::pair<label, label>, scalar> featureEdgeAngle;
//        forAll(sef, eI)
//        {
//            if( sef.sizeOfRow(eI) != 2 )
//                continue;

//            const label patch0 = outletSurf[sef(eI, 0)].region();
//            const label patch1 = outletSurf[sef(eI, 1)].region();

//            if( patch0 != patch1 )
//            {
//                const point p =
//                    outletSurf.edges()[eI].centre(outletSurf.points());

//                const vector v = p - point(p.x(), 0., 0.);

//                const scalar r = mag(v) + VSMALL;

//                //- calculate the sine and the cosine
//                const scalar cosAlpha = (v & vector(0., 1., 0.)) / r;
//                const scalar sinAlpha = (v & vector(0., 0., 1.)) / r;

//                //- calculate the angle
//                scalar alpha = acos(cosAlpha);
//                if( sinAlpha < 0.0 )
//                    alpha = 2. * M_PI - alpha;

//                //- save the angle into the map
//                std::pair<label, label> pp
//                (
//                    min(patch0, patch1),
//                    max(patch0, patch1)
//                );
//                featureEdgeAngle[pp] = alpha;
//            }
//        }

//        polyMeshGenModifier meshModifier(*meshPtr);
//        pointFieldPMG& points = meshModifier.pointsAccess();

//        //- the outlet is not circular
//        meshSurfaceEngine mse(*meshPtr);
//        const labelLongList& bp = mse.bp();
//        const labelLongList& facePatch = mse.boundaryFacePatches();
//        const VRWGraph& pFaces = mse.pointFaces();

//        const scalar inletRadius = 0.5 * dieGeom.inletDiameter();
//        const scalar outletRadius = 0.5 * dieGeom.outletDiameter();
//        const scalar outerRadius = 0.5 * dieGeom.outerDiameter();

//        //- calculate displacements for vertices in the radial direction
//        scalar xmin(points[0].x()), xmax(points[0].x());
//        for(label pointI=1;pointI<points.size();++pointI)
//        {
//            xmin = min(xmin, points[pointI].x());
//            xmax = max(xmax, points[pointI].x());
//        }

//        std::map<label, label> planeMovingPoint;

//        DynList<label> subsetIds;
//        meshPtr->pointSubsetIndices(subsetIds);

//        //- find the angles of the planes
//        std::map<label, scalar> planeCircumAngle;

//        if( featureEdgeAngle.size() )
//        {
//            forAll(subsetIds, i)
//            {
//                const word sName = meshPtr->pointSubsetName(subsetIds[i]);

//                if( sName.find("pointsInPlane_") != word::npos )
//                {
//                    const label planeI =
//                        help::textToLabel(sName.substr(14, sName.size()-14));

//                    labelLongList pointsInSubset;
//                    meshPtr->pointsInSubset(subsetIds[i], pointsInSubset);

//                    if( pointsInSubset.size() == 0 )
//                        continue;

//                    const point& p = points[pointsInSubset[0]];
//                    const vector v = p - point(p.x(), 0., 0.);

//                    const scalar r = mag(v) + VSMALL;

//                    //- calculate the sine and the cosine
//                    const scalar cosAlpha = (v & vector(0., 1., 0.)) / r;
//                    const scalar sinAlpha = (v & vector(0., 0., 1.)) / r;

//                    //- calculate the angle
//                    scalar alpha = acos(cosAlpha);
//                    if( sinAlpha < 0.0 )
//                        alpha = 2. * M_PI - alpha;

//                    //- save the angle into the map
//                    planeCircumAngle[planeI] = alpha;
//                }
//            }

//            //- find the best fitting plane for this angle
//            std::set<label> fixedPlanes;
//            for
//            (
//                std::map<std::pair<label, label>, scalar>::const_iterator it=
//                    featureEdgeAngle.begin();
//                it!=featureEdgeAngle.end();
//                ++it
//            )
//            {
//                scalar angleDiff(VGREAT);
//                label bestPlane(-1);

//                for
//                (
//                    std::map<label, scalar>::const_iterator pIt=
//                        planeCircumAngle.begin();
//                    pIt!=planeCircumAngle.end();
//                    ++pIt
//                )
//                {
//                    const scalar diff = mag(it->second - pIt->second);

//                    if( diff < angleDiff )
//                    {
//                        angleDiff = diff;
//                        bestPlane = pIt->first;
//                    }
//                }

//                if( bestPlane != -1 )
//                {
//                    planeCircumAngle[bestPlane] = it->second;
//                    fixedPlanes.insert(bestPlane);
//                }
//            }

//            //- smooth circum angles
//            LongList<scalar> circumAngles;
//            for
//            (
//                std::map<label, scalar>::const_iterator it=
//                    planeCircumAngle.begin();
//                it!=planeCircumAngle.end();
//                ++it
//            )
//            {
//                if( circumAngles.size() <= it->first )
//                    circumAngles.setSize(it->first+1);

//                circumAngles[it->first] = it->second;
//            }

//            const label s = circumAngles.size();

//            for(label iterI=0;iterI<10;++iterI)
//            {
//                forAll(circumAngles, angleI)
//                {
//                    if( fixedPlanes.find(angleI) != fixedPlanes.end() )
//                        continue;

//                    const scalar prevAngle = circumAngles[(angleI-1+s)%s];
//                    scalar nextAngle = circumAngles[(angleI+1)%s];

//                    if
//                    (
//                        (circumAngles[angleI] < prevAngle) ||
//                        (circumAngles[angleI] > nextAngle)
//                    )
//                    {
//                        nextAngle += 2. * M_PI;

//                        circumAngles[angleI] = 0.5 * (prevAngle + nextAngle);
//                        circumAngles[angleI] -= 2. * M_PI;
//                    }
//                    else
//                    {
//                        circumAngles[angleI] = 0.5 * (prevAngle + nextAngle);
//                    }

//                    //- update the angle at this plane
//                    planeCircumAngle[angleI] = circumAngles[angleI];
//                }
//            }
//        }

//        forAll(subsetIds, i)
//        {
//            const word sName = meshPtr->pointSubsetName(subsetIds[i]);

//            if( sName.find("pointsInPlane_") != word::npos )
//            {
//                const label planeI =
//                    help::textToLabel(sName.substr(14, sName.size()-14));

//                labelLongList pointsInSubset;
//                meshPtr->pointsInSubset(subsetIds[i], pointsInSubset);

//                if( pointsInSubset.size() == 0 )
//                    continue;

//                LongList<scalar> axialPos(pointsInSubset.size());
//                LongList<scalar> radialPos(pointsInSubset.size());

//                label movingPoint(-1);

//                forAll(pointsInSubset, pI)
//                {
//                    const label pointI = pointsInSubset[pI];

//                    const label bpI = bp[pointI];

//                    //- check if this point is at the profile
//                    if( bpI >= 0 )
//                    {
//                        DynList<label> pointPatches;
//                        forAllRow(pFaces, bpI, pfI)
//                            pointPatches.appendIfNotIn
//                            (
//                                facePatch[pFaces(bpI, pfI)]
//                            );
//                        if
//                        (
//                            pointPatches.contains(1) &&
//                            pointPatches.contains(2)
//                        )
//                            movingPoint = pointI;
//                    }

//                    const point& p = points[pointI];

//                    axialPos[pI] = (p.x() - xmin) / (xmax - xmin);

//                    const scalar r = sqrt(magSqr(p - point(p.x(), 0., 0.)));

//                    const scalar rmin =
//                        axialPos[pI] * outletRadius +
//                        (1.0 - axialPos[pI]) * inletRadius;

//                    radialPos[pI] = (r - rmin) / (outerRadius - rmin);
//                }

//                //- find the nearest point at the surface mesh
//                if( movingPoint == -1 )
//                    FatalError << "Could not find moving point in subset "
//                        << sName << abort(FatalError);

//                planeMovingPoint[planeI] = movingPoint;

//                //- calculate the direction vector for this plane
//                const point& pMoving = points[movingPoint];
//                const point p0(pMoving.x(), 0., 0.);
//                const vector dirVec =
//                    vector(0., 1., 0.) * cos(planeCircumAngle[planeI]) +
//                    vector(0., 0., 1.) * sin(planeCircumAngle[planeI]);
//                const point p1 = p0 + outerRadius * dirVec;

//                //- calculate the radius of the profile at this angle
//                scalar newOutletRadius(outletRadius);
//                forAll(outletSurf, triI)
//                {
//                    point s;
//                    if( help::triLineIntersection(outletSurf, triI, p0, p1, s) )
//                    {
//                        newOutletRadius = mag(s - point(s.x(), 0., 0.));
//                    }
//                }

//                //- calculate new coordinates and move the point
//                forAll(pointsInSubset, pI)
//                {
//                    const label pointI = pointsInSubset[pI];
//                    point& moveP = points[pointI];

//                    const point axisPoint = point(moveP.x(), 0., 0.);

//                    //- interpolate the min radius
//                    const scalar rmin =
//                        axialPos[pI] * newOutletRadius +
//                        (1.0 - axialPos[pI]) * inletRadius;

//                    //- calculate the new point radius
//                    //- the outer diameter of a die is constant
//                    const scalar newRadius =
//                        radialPos[pI] * outerRadius +
//                        (1.0 - radialPos[pI]) * rmin;

//                    //- move the point to the new location
//                    moveP = axisPoint + newRadius * dirVec;
//                }
//            }
//        }
//    }

    if( dieGeom.isProfiled() )
    {
        //- create a profiled die
        profiledDieMeshGenerator profileMesher
        (
            *meshPtr,
            dieGeom,
            patchHandler_
        );

        profileMesher.generateProfiledDie();
    }
    //- remove all subsets from the mesh
    DynList<label> subsetIds;

    //- remove point subsets
    meshPtr->pointSubsetIndices(subsetIds);

    forAll(subsetIds, i)
        meshPtr->removePointSubset(subsetIds[i]);

    //- remove face subsets
    subsetIds.clear();
    meshPtr->faceSubsetIndices(subsetIds);

    forAll(subsetIds, i)
        meshPtr->removeFaceSubset(subsetIds[i]);

    //- remove cell subsets
    subsetIds.clear();
    meshPtr->cellSubsetIndices(subsetIds);

    forAll(subsetIds, i)
        meshPtr->removeCellSubset(subsetIds[i]);

    //- remove duplicated symmetry patches
    std::map<word, label> patchToId;

    label patchCounter(0);
    forAll(boundaries, patchI)
    {
        const boundaryPatch& patch = boundaries[patchI];

        if( patchToId.find(patch.patchName()) == patchToId.end() )
        {
            patchToId[patch.patchName()] = patchCounter++;
        }
    }

    wordList patchNames(patchCounter), patchTypes(patchCounter);

    labelList newPatchIndex(boundaries.size());
    forAll(boundaries, patchI)
    {
        const boundaryPatch& patch = boundaries[patchI];
        const word& pName = patch.patchName();

        patchNames[patchToId[pName]] = pName;
        patchTypes[patchToId[pName]] = patch.patchType();
        newPatchIndex[patchI] = patchToId[pName];
    }

    VRWGraph newBndFaces;
    labelLongList faceOwner;
    labelLongList facePatch;

    const labelLongList& owner = meshPtr->owner();
    const faceListPMG& faces = meshPtr->faces();
    forAll(boundaries, patchI)
    {
        const boundaryPatch& patch = boundaries[patchI];
        const label start = patch.patchStart();
        const label end = start + patch.patchSize();

        for(label faceI=start;faceI<end;++faceI)
        {
            newBndFaces.appendList(faces[faceI]);
            faceOwner.append(owner[faceI]);
            facePatch.append(newPatchIndex[patchI]);
        }
    }

    polyMeshGenModifier(*meshPtr).replaceBoundary
    (
        patchNames,
        newBndFaces,
        faceOwner,
        facePatch
    );

    forAll(boundaries, patchI)
        boundaries[patchI].patchType() = patchTypes[patchI];

    //- translate die such that the downstream patch is at x = 0
    scalar maxX(-VGREAT);
    pointFieldPMG& points = polyMeshGenModifier(*meshPtr).pointsAccess();
    forAll(points, pointI)
    {
        maxX = max(points[pointI].x(), maxX);
    }

    forAll(points, pointI)
        points[pointI].x() -= maxX;

    if( globalMeshPtr_ )
    {
        //- add a die into the existing mesh
        polyMeshGenModifier(*globalMeshPtr_).addMesh(*meshPtr);

        deleteDemandDrivenData(meshPtr);
    }
    else
    {
        //- insert the mesh as a global mesh
        globalMeshPtr_ = meshPtr;
    }
}

void rollingMillMesh::generateCasingMesh()
{
    const dieGeometryInfo& dieGeom = geomHandler_.dieGeometry();

    //- check if a casing shall be meshed or not
    if( !dieGeom.isCasingPresent() )
        return;

    //- cross section of the casing
    const triSurf* surfPtr = dieGeom.casingCrossSectionSurface();

    const dictionary& dieDict = dieGeom.dieDict();
    if( !dieDict.found("casingDict") || !dieDict.isDict("casingDict") )
    {
        FatalError << "casingDict does nor exist in dieMeshDict "
            << dieDict << exit(FatalError);
    }

    const dictionary& casingDict = dieDict.subDict("casingDict");

    //- generate a 2D mesh of a cross section
    HashSet<const word> contactNames;
    contactNames.insert
    (
        patchHandler_.patchNameForCasing
        (
            rollingMillPatchNamesHandler::CASINGTODIERADIAL
        )
    );
    contactNames.insert
    (
        patchHandler_.patchNameForCasing
        (
            rollingMillPatchNamesHandler::CASINGTODIEAXIAL
        )
    );
    crossSectionMeshGenerator sectionMesh
    (
        surfPtr,
        meshDict_,
        casingDict,
        db_,
        regionName_,
        timeStep_,
        contactNames
    );

    polyMeshGen* meshPtr = sectionMesh.meshPtr();

    //- the surface is not needed any more
    deleteDemandDrivenData(surfPtr);

    Info << "Generating a 3D mesh from the cross-section mesh" << endl;
    revolve2DMesh meshRevolver(*meshPtr);

    //- set the patch at the zero plane
    meshRevolver.setRevolvingPatch("bottomEmptyFaces");

    //- get the circumferential resolution
    const scalar inletDiameter = dieGeom.inletDiameter();
    const scalar geometryTolerance = geomHandler_.geometryTolerance();

    scalar angleStep =
        2.0 * Foam::acos(1.0 - geometryTolerance / (0.5*inletDiameter));

    //- check the revolution angle
    if( dieGeom.isSymmetric() )
    {
        meshRevolver.clearAngleIntervals();

        if( dieGeom.isWedge() )
        {
            //- axi-symmetric setup
            const scalar wedgeAngle = dieGeom.wedgeAngle() * M_PI / 360.0;

            if( angleStep < wedgeAngle )
            {
                Warning << "Wedge angle " << dieGeom.wedgeAngle()
                     << " is greater "
                     << "than " << (angleStep * 180.0 / M_PI)
                     << ", needed to achieve geometry tolerance " << endl;
            }

            angleStep = wedgeAngle;

            meshRevolver.setIntervalResolution(-angleStep, angleStep, 1);
        }
        else
        {
            //- symmetric setup
            meshRevolver.setIntervalResolution
            (
                dieGeom.startCircumAngle(),
                dieGeom.endCircumAngle(),
                angleStep
            );
        }
    }
    else
    {
        meshRevolver.setCircumResolution(angleStep);
    }

    //- set the origin and the rotation axis
    meshRevolver.setOrigin(vector(0, 0., 0.));
    meshRevolver.setRotationAxis(vector(1., 0., 0.));

    //- generate point subsets
    meshRevolver.createPointSubsets();

    //- create a 3D mesh
    meshRevolver.generateRevolvedMesh();

    //- add cells to the zone
    Info << "Generating zones" << endl;
    const label cId = meshPtr->addCellZone("casing");
    forAll(meshPtr->cells(), cellI)
        meshPtr->addCellToZone(cId, cellI);

    //- add faces in the contact patch to a zone
    const word contactPatchName =
        patchHandler_.patchNameForCasing
        (
            rollingMillPatchNamesHandler::CASINGTODIERADIAL
        );
    const label fId = meshPtr->addFaceZone(contactPatchName);

    //- set the patch names and types
    PtrList<boundaryPatch>& boundaries =
        polyMeshGenModifier(*meshPtr).boundariesAccess();
    forAll(boundaries, patchI)
    {
        if( boundaries[patchI].patchName() == contactPatchName )
        {
            label start = boundaries[patchI].patchStart();
            const label size = boundaries[patchI].patchSize();

            for(label i=0;i<size;++i)
                meshPtr->addFaceToZone(fId, start++);
        }

        if( dieGeom.isSymmetric() )
        {
            if( boundaries[patchI].patchName() == "defaultFaces" )
            {
                if( dieGeom.isWedge() )
                {
                    boundaries[patchI].patchName() =
                        patchHandler_.patchNameForCasing
                        (
                            rollingMillPatchNamesHandler::CASINGFRONT
                        );
                    boundaries[patchI].patchType() =
                        patchHandler_.patchTypeForCasing
                        (
                            rollingMillPatchNamesHandler::CASINGFRONT,
                            rollingMillGeometryHandler::symmetryTypes_
                            (
                                dieGeom.typeOfSymmetry()
                            )
                        );
                }
                else
                {
                    switch( geomHandler_.symmetryType() )
                    {
                        case rollingMillGeometryHandler::NEGY:
                        case rollingMillGeometryHandler::POSY:
                        case
                        (
                            rollingMillGeometryHandler::POSY +
                            rollingMillGeometryHandler::POSZ
                        ):
                        case
                        (
                            rollingMillGeometryHandler::NEGY +
                            rollingMillGeometryHandler::NEGZ
                        ):
                        {
                            boundaries[patchI].patchName() =
                                patchHandler_.patchNameForCasing
                                (
                                    rollingMillPatchNamesHandler::CASINGSYMMY
                                );
                            boundaries[patchI].patchType() =
                                patchHandler_.patchTypeForCasing
                                (
                                    rollingMillPatchNamesHandler::CASINGSYMMY,
                                    rollingMillGeometryHandler::symmetryTypes_
                                    (
                                        dieGeom.typeOfSymmetry()
                                    )
                                );
                        } break;
                        case rollingMillGeometryHandler::NEGZ:
                        case rollingMillGeometryHandler::POSZ:
                        case
                        (
                            rollingMillGeometryHandler::NEGY +
                            rollingMillGeometryHandler::POSZ
                        ):
                        case
                        (
                            rollingMillGeometryHandler::NEGZ +
                            rollingMillGeometryHandler::POSY
                        ):
                        {
                            boundaries[patchI].patchName() =
                                patchHandler_.patchNameForCasing
                                (
                                    rollingMillPatchNamesHandler::CASINGSYMMZ
                                );
                            boundaries[patchI].patchType() =
                                patchHandler_.patchTypeForCasing
                                (
                                    rollingMillPatchNamesHandler::CASINGSYMMZ,
                                    rollingMillGeometryHandler::symmetryTypes_
                                    (
                                        dieGeom.typeOfSymmetry()
                                    )
                                );
                        } break;
                    }
                }
            }
            else if( boundaries[patchI].patchName() == "bottomEmptyFaces" )
            {
                if( dieGeom.isWedge() )
                {
                    boundaries[patchI].patchName() =
                        patchHandler_.patchNameForCasing
                        (
                            rollingMillPatchNamesHandler::CASINGBACK
                        );
                    boundaries[patchI].patchType() =
                        patchHandler_.patchTypeForCasing
                        (
                            rollingMillPatchNamesHandler::CASINGBACK,
                            rollingMillGeometryHandler::symmetryTypes_
                            (
                                dieGeom.typeOfSymmetry()
                            )
                        );
                }
                else
                {
                    switch( geomHandler_.symmetryType() )
                    {
                        case rollingMillGeometryHandler::NEGZ:
                        case rollingMillGeometryHandler::POSZ:
                        case
                        (
                            rollingMillGeometryHandler::POSY +
                            rollingMillGeometryHandler::POSZ
                        ):
                        case
                        (
                            rollingMillGeometryHandler::NEGY +
                            rollingMillGeometryHandler::NEGZ
                        ):
                        {
                            boundaries[patchI].patchName() =
                                patchHandler_.patchNameForCasing
                                (
                                    rollingMillPatchNamesHandler::CASINGSYMMZ
                                );
                            boundaries[patchI].patchType() =
                                patchHandler_.patchTypeForCasing
                                (
                                    rollingMillPatchNamesHandler::CASINGSYMMZ,
                                    rollingMillGeometryHandler::symmetryTypes_
                                    (
                                        dieGeom.typeOfSymmetry()
                                    )
                                );
                        } break;
                        case rollingMillGeometryHandler::NEGY:
                        case rollingMillGeometryHandler::POSY:
                        case
                        (
                            rollingMillGeometryHandler::NEGY +
                            rollingMillGeometryHandler::POSZ
                        ):
                        case
                        (
                            rollingMillGeometryHandler::NEGZ +
                            rollingMillGeometryHandler::POSY
                        ):
                        {
                            boundaries[patchI].patchName() =
                                patchHandler_.patchNameForCasing
                                (
                                    rollingMillPatchNamesHandler::CASINGSYMMY
                                );
                            boundaries[patchI].patchType() =
                                patchHandler_.patchTypeForCasing
                                (
                                    rollingMillPatchNamesHandler::CASINGSYMMY,
                                    rollingMillGeometryHandler::symmetryTypes_
                                    (
                                        dieGeom.typeOfSymmetry()
                                    )
                                );
                        } break;
                    }
                }
            }
        }
    }

    //- remove all subsets from the mesh
    DynList<label> subsetIds;

    //- remove point subsets
    meshPtr->pointSubsetIndices(subsetIds);

    forAll(subsetIds, i)
        meshPtr->removePointSubset(subsetIds[i]);

    //- remove face subsets
    subsetIds.clear();
    meshPtr->faceSubsetIndices(subsetIds);

    forAll(subsetIds, i)
        meshPtr->removeFaceSubset(subsetIds[i]);

    //- remove cell subsets
    subsetIds.clear();
    meshPtr->cellSubsetIndices(subsetIds);

    forAll(subsetIds, i)
        meshPtr->removeCellSubset(subsetIds[i]);

    //- remove duplicated symmetry patches
    std::map<word, label> patchToId;

    label patchCounter(0);
    forAll(boundaries, patchI)
    {
        const boundaryPatch& patch = boundaries[patchI];

        if( patchToId.find(patch.patchName()) == patchToId.end() )
        {
            patchToId[patch.patchName()] = patchCounter++;
        }
    }

    wordList patchNames(patchCounter), patchTypes(patchCounter);

    labelList newPatchIndex(boundaries.size());
    forAll(boundaries, patchI)
    {
        const boundaryPatch& patch = boundaries[patchI];
        const word& pName = patch.patchName();

        patchNames[patchToId[pName]] = pName;
        patchTypes[patchToId[pName]] = patch.patchType();
        newPatchIndex[patchI] = patchToId[pName];
    }

    VRWGraph newBndFaces;
    labelLongList faceOwner;
    labelLongList facePatch;

    const labelLongList& owner = meshPtr->owner();
    const faceListPMG& faces = meshPtr->faces();
    forAll(boundaries, patchI)
    {
        const boundaryPatch& patch = boundaries[patchI];
        const label start = patch.patchStart();
        const label end = start + patch.patchSize();

        for(label faceI=start;faceI<end;++faceI)
        {
            newBndFaces.appendList(faces[faceI]);
            faceOwner.append(owner[faceI]);
            facePatch.append(newPatchIndex[patchI]);
        }
    }

    polyMeshGenModifier(*meshPtr).replaceBoundary
    (
        patchNames,
        newBndFaces,
        faceOwner,
        facePatch
    );

    forAll(boundaries, patchI)
        boundaries[patchI].patchType() = patchTypes[patchI];

    if( globalMeshPtr_ )
    {
        //- add a die into the existing mesh
        polyMeshGenModifier(*globalMeshPtr_).addMesh(*meshPtr);

        deleteDemandDrivenData(meshPtr);
    }
    else
    {
        //- insert the mesh as a global mesh
        globalMeshPtr_ = meshPtr;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
