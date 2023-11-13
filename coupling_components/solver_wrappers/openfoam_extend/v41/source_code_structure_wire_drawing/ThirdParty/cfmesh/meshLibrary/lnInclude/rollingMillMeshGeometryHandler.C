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
#include "triSurf.H"
#include "triSurfModifier.H"
#include "triSurfacePatchManipulator.H"
#include "rollerSurfaceCreator.H"
#include "rollerSurfaceCutter.H"
#include "wireCrossSectionSurfaceCreator.H"
#include "boundBox.H"
#include "plane.H"

#include "helperFunctions.H"

//#define DEBUG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const wordHashSet rollingMillGeometryHandler::symmetryTypes =
    rollingMillGeometryHandler::populateSymmetricConfigurations();

wordHashSet rollingMillGeometryHandler::populateSymmetricConfigurations()
{
    wordHashSet validSymmetricConfigurations;

    validSymmetricConfigurations.insert("none");
    validSymmetricConfigurations.insert("posY");
    validSymmetricConfigurations.insert("negY");
    validSymmetricConfigurations.insert("posZ");
    validSymmetricConfigurations.insert("negZ");
    validSymmetricConfigurations.insert("posYposZ");
    validSymmetricConfigurations.insert("negYposZ");
    validSymmetricConfigurations.insert("posYnegZ");
    validSymmetricConfigurations.insert("negYnegZ");
    validSymmetricConfigurations.insert("axiSymmetric");

    return validSymmetricConfigurations;
}

direction rollingMillGeometryHandler::symmetryType() const
{
    return symmetryType_;

    //- old implementation
    bool allPosY(true), allPosZ(true), allNegY(true), allNegZ(true);

    forAll(geometries_, geomI)
    {
        const triSurf& surf = geometries_[geomI].surface();

        const vector translationVec =
            geometries_[geomI].radialShift(wireDiameter_);

        const pointField& points = surf.points();

        forAll(points, pointI)
        {
            const point p = points[pointI] + translationVec;

            if( p.y() > SMALL )
                allNegY = false;
            if( p.y() < -SMALL )
                allPosY = false;
            if( p.z() > SMALL )
                allNegZ = false;
            if( p.z() < -SMALL )
                allPosZ = false;
        }
    }

    direction typeOfSymmetry(NONE);
    if( allNegY )
        typeOfSymmetry |= NEGY;
    if( allPosY )
        typeOfSymmetry |= POSY;
    if( allNegZ )
        typeOfSymmetry |= NEGZ;
    if( allPosZ )
        typeOfSymmetry |= POSZ;

    return typeOfSymmetry;
}

bool rollingMillGeometryHandler::checkHalfSymmetry() const
{
    const direction symmType = symmetryType();

    if( symmType == NEGY )
        return true;
    if( symmType == POSY )
        return true;

    if( symmType == NEGZ )
        return true;
    if( symmType == POSZ )
        return true;

    return false;
}

bool rollingMillGeometryHandler::checkQuarterSymmetry() const
{
    const direction symmType = symmetryType();

    label counter(0);
    if( symmType & NEGY )
        ++counter;
    if( symmType & POSY )
        ++counter;
    if( symmType & NEGZ )
        ++counter;
    if( symmType & POSZ )
        ++counter;

    return (counter == 2);
}

void rollingMillGeometryHandler::writeRollPositions() const
{
    IOdictionary positionsDict
    (
        IOobject
        (
            "rollPositions",
            runTime_.constant(),
            runTime_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    forAll(geometries_, geomI)
    {
        const rollGeometryInfo& geom = geometries_[geomI];

        positionsDict.add(geom.rollPosition(), geom.radialShiftDistance());
    }

    positionsDict.writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED
    );
}

void rollingMillGeometryHandler::transformRollToPosition
(
    triSurf& rollSurf,
    const word& rollPosition,
    const bool flipRoll
) const
{
    Info << "Transforming roll " << rollPosition << " to position" << endl;

    triSurfModifier sMod(rollSurf);
    pointField& pts = sMod.pointsAccess();

    if( flipRoll )
    {
        //- mirror the profile by the x-y plane
        forAll(pts, pI)
            pts[pI].z() *= -1.0;
    }

    if( rollPosition == "topRoll" )
    {
        //- mirror the by the x-z plane
        forAll(pts, pI)
            pts[pI].y() *= -1.0;
    }
    else if( rollPosition == "rightRoll" )
    {
        //- rotate the profile by 90 degrees
        tensor transformationMatrix(tensor::zero);
        transformationMatrix.xx() = 1.0;
        transformationMatrix.yz() = 1.0;
        transformationMatrix.zy() = -1.0;

        forAll(pts, pI)
            pts[pI] = (transformationMatrix & pts[pI]);
    }
    else if( rollPosition == "leftRoll" )
    {
        //- rotate the profile by -90 degrees and flip the y coordinate
        tensor transformationMatrix(tensor::zero);
        transformationMatrix.xx() = 1.0;
        transformationMatrix.zy() = 1.0;
        transformationMatrix.yz() = -1.0;

        forAll(pts, pI)
        {
            pts[pI] = (transformationMatrix & pts[pI]);

            pts[pI].y() *= -1.0;
        }
    }
    else if( rollPosition == "bottomRightRoll" )
    {
        //- rotate the roll by 60 degrees
        tensor transformationMatrix(tensor::zero);
        transformationMatrix.xx() = 1.0;
        transformationMatrix.yy() = cos(-M_PI/3.0);
        transformationMatrix.yz() = -sin(-M_PI/3.0);
        transformationMatrix.zy() = sin(-M_PI/3.0);
        transformationMatrix.zz() = cos(-M_PI/3.0);

        forAll(pts, pI)
            pts[pI] = (transformationMatrix & pts[pI]);
    }
    else if( rollPosition == "bottomLeftRoll" )
    {
        //- rotate the roll by -60 degrees
        tensor transformationMatrix(tensor::zero);
        transformationMatrix.xx() = 1.0;
        transformationMatrix.yy() = cos(M_PI/3.0);
        transformationMatrix.yz() = -sin(M_PI/3.0);
        transformationMatrix.zy() = sin(M_PI/3.0);
        transformationMatrix.zz() = cos(M_PI/3.0);

        forAll(pts, pI)
            pts[pI] = (transformationMatrix & pts[pI]);
    }

    Info << "Finished transforming roll "
         << rollPosition << " to position" << endl;
}

bool rollingMillGeometryHandler::cutGeometry
(
    triSurf& rollSurf,
    const word& rollPosition
) const
{
    //- calculate cutting planes
    PtrList<plane> planes(2);
    DynList<word, 2> patchNames;
    label nPlanes(0);

    if( rollPosition == "topRoll" )
    {
        if( symmetryType_ & NEGY )
            return false;

        if( symmetryType_ & POSZ )
        {
            planes.set(nPlanes, new plane(vector::zero, vector(0., 0., 1.)));
            patchNames.append
            (
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERSYMMZ
                )
            );
            ++nPlanes;
        }
        else if( symmetryType_ & NEGZ )
        {
            planes.set(nPlanes, new plane(vector::zero, vector(0., 0., -1.)));
            patchNames.append
            (
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERSYMMZ
                )
            );

            ++nPlanes;
        }
    }
    else if( rollPosition == "bottomRoll" )
    {
        if( symmetryType_ & POSY )
            return false;

        if( symmetryType_ & POSZ )
        {
            planes.set(nPlanes, new plane(vector::zero, vector(0., 0., 1.)));
            patchNames.append
            (
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERSYMMZ
                )
            );

            ++nPlanes;
        }
        else if( symmetryType_ & NEGZ )
        {
            planes.set(nPlanes, new plane(vector::zero, vector(0., 0., -1.)));
            patchNames.append
            (
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERSYMMZ
                )
            );

            ++nPlanes;
        }
    }
    else if( rollPosition == "leftRoll" )
    {
        if( symmetryType_ & POSZ )
            return false;

        if( symmetryType_ & POSY )
        {
            planes.set(nPlanes, new plane(vector::zero, vector(0., 1., 0.)));
            patchNames.append
            (
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERSYMMY
                )
            );
            ++nPlanes;
        }
        else if( symmetryType_ & NEGY )
        {
            planes.set(nPlanes, new plane(vector::zero, vector(0., -1., 0.)));
            patchNames.append
            (
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERSYMMY
                )
            );

            ++nPlanes;
        }
    }
    else if( rollPosition == "rightRoll" )
    {
        if( symmetryType_ & NEGZ )
            return false;

        if( symmetryType_ & POSY )
        {
            planes.set(nPlanes, new plane(vector::zero, vector(0., 1., 0.)));
            patchNames.append
            (
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERSYMMY
                )
            );

            ++nPlanes;
        }
        else if( symmetryType_ & NEGY )
        {
            planes.set(nPlanes, new plane(vector::zero, vector(0., -1., 0.)));
            patchNames.append
            (
                patchHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERSYMMY
                )
            );

            ++nPlanes;
        }
    }
    else if( rollPosition == "bottomLeftRoll" )
    {
        if( symmetryType_ & (POSY|POSZ) )
        {
            return false;
        }
        else if( symmetryType_ & (NEGY|NEGZ) )
        {
            return true;
        }
    }
    else if( rollPosition == "bottomRightRoll" )
    {
        if( symmetryType_ & (POSY|NEGZ) )
        {
            return false;
        }
        else if( symmetryType_ & (NEGY|POSZ) )
        {
            return true;
        }
    }
    else if( rollPosition == "wire" )
    {
        if( symmetryType_ & POSZ )
        {
            planes.set(nPlanes, new plane(vector::zero, vector(0., 0., 1.)));
            patchNames.append
            (
                patchHandler_.patchNameForWire
                (
                    rollingMillPatchNamesHandler::WIRESYMMZ
                )
            );

            ++nPlanes;
        }
        else if( symmetryType_ & NEGZ )
        {
            planes.set(nPlanes, new plane(vector::zero, vector(0., 0., -1.)));
            patchNames.append
            (
                patchHandler_.patchNameForWire
                (
                    rollingMillPatchNamesHandler::WIRESYMMZ
                )
            );

            ++nPlanes;
        }

        if( symmetryType_ & POSY )
        {
            planes.set(nPlanes, new plane(vector::zero, vector(0., 1., 0.)));
            patchNames.append
            (
                patchHandler_.patchNameForWire
                (
                    rollingMillPatchNamesHandler::WIRESYMMY
                )
            );

            ++nPlanes;
        }
        else if( symmetryType_ & NEGY )
        {
            planes.set(nPlanes, new plane(vector::zero, vector(0., -1., 0.)));
            patchNames.append
            (
                patchHandler_.patchNameForWire
                (
                    rollingMillPatchNamesHandler::WIRESYMMY
                )
            );

            ++nPlanes;
        }
    }

    //- check if any symmetry constraints are imposed
    if( nPlanes == 0 )
        return true;

    planes.setSize(nPlanes);

    for(label i=0;i<nPlanes;++i)
        rollerSurfaceCutter(rollSurf, planes[i], patchNames[i]);

    return rollSurf.size() != 0;
}

void rollingMillGeometryHandler::parseRollDictionary
(
    const dictionary& rollDict,
    const word rollOffsetName
)
{
    //- check the position of the roll
    word position;

    forAllConstIter(wordHashSet, rollGeometryInfo::rollPositions, it)
    {
        if( rollDict.name().find(it.key()) != word::npos )
        {
            position = it.key();
        }
    }

    if( position == "" )
    {
        FatalErrorIn
        (
            "void rollingMillGeometryHandler::"
            "parseRollDictionary(const dictionary&)"
        ) << "Unknown roll position specified in dictionary "
          << rollDict.name() << exit(FatalError);
    }

    Info << "Parsing dictionary for roll " << position << endl;

    //- shall the roll be flipped or not
    bool flip(false);
    if( rollDict.found("flipRoll") )
        flip = readBool(rollDict.lookup("flipRoll"));

    //- create the geometry of a roller
    autoPtr<rollerSurfaceCreator> surfCreatorPtr =
        rollerSurfaceCreator::New
        (
            position,
            rollingMillGeometryHandler::symmetryTypes_(symmetryType_),
            patchHandler_,
            rollDict,
            geometryTolerance_
        );

    surfCreatorPtr->generateGeometry();

    //- set the roll surface
    triSurf* surfPtr = new triSurf(surfCreatorPtr->surface());

    transformRollToPosition(*surfPtr, position, flip);

    bool writeProfile(false);
    rollDict.readIfPresent("writeRollProfile", writeProfile);

    if( writeProfile )
    {
        surfPtr->writeSurface("rollProfile_"+position+".fms");
    }

    if( cutGeometry(*surfPtr, position) )
    {
        const label size = geometries_.size();
        geometries_.setSize(size+1);
        geometries_.set
        (
            size,
            new rollGeometryInfo(surfPtr, rollDict, position)
        );

        geometries_[size].rollOffset() = rollOffsetName;
        geometries_[size].outerDiameter() = surfCreatorPtr->outerRollDiameter();
        geometries_[size].radialShiftDistance() = surfCreatorPtr->radialShift();
        geometries_[size].contactWidth() = surfCreatorPtr->contactWidth();

        if( writeProfile )
        {
            surfPtr->writeSurface("rollProfileAfterSymm_"+position+".fms");
        }
    }

    if( rollDict.found("meshOffset") )
    {
        Info << "Reading offset for rollers" << endl;
        PtrList<entry> offsetRollers(rollDict.lookup("meshOffset"));

        const wordList currentData = rollDict.toc();

        forAll(offsetRollers, offsetI)
        {
            dictionary dict(position);

            dict += offsetRollers[offsetI].dict();

            //- copy setting available in the parent dict
            forAll(currentData, i)
            {
                if( currentData[i] == "meshOffset" )
                    continue;

                if( !dict.found(currentData[i]) )
                {
                    dict.add
                    (
                        rollDict.lookupEntry(currentData[i], false, false)
                    );
                }
            }

            parseRollDictionary(dict, offsetRollers[offsetI].keyword());
        }
    }

    Info << "Finished parsing dictionary for roll " << position << endl;
}

void rollingMillGeometryHandler::parseDictionaries
(
    const dictionary& meshDict
)
{
    Info << "Parsing roll setup from dictionary" << endl;

    //- check if roll setup exists in meshDict
    if( !meshDict.found("rollSetup") )
        FatalErrorIn
        (
            "void rollingMillMesh::rollingMillGeometryHandler::"
            "parseDictionaries(const dictionary& meshDict)"
        ) << "rollSetup keyword does not exist in meshDict!"
          << exit(FatalError);

    const word rollSetup = word(meshDict.lookup("rollSetup"));

    //- read the global geometry tolerance
    if( meshDict.found("geometryTolerance") )
        geometryTolerance_ = readScalar(meshDict.lookup("geometryTolerance"));

    //- read the requested type of symmetry for the given case
    word symmType = "none";

    if( meshDict.found("symmetryType") )
    {
        symmType = word(meshDict.lookup("symmetryType"));

        Info << "Symmetry configuration " << symmType << endl;

        if( !symmetryTypes.found(symmType) )
        {
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "Invalid symmetryType " << symmType
              << " valid symmetric cofigurations are " << symmetryTypes
              << exit(FatalError);
        }

        if( symmType == "none" )
        {
            symmetryType_ = symmetryTypes_::NONE;
        }
        else if( symmType == "axiSymmetric" )
        {
            symmetryType_ = symmetryTypes_::AXISYMMETRIC;

            if( rollSetup != "none" )
            {
                FatalErrorIn
                (
                    "void rollingMillGeometryHandler::parseDictionaries"
                    "(const dictionary&)"
                ) << "axiSymmetric type of symmetry cannot be used"
                  << " for a rolling pass" << exit(FatalError);
            }
        }

        if( symmType.find("posY") != word::npos )
        {
            symmetryType_ |= symmetryTypes_::POSY;
        }

        if( symmType.find("negY") != word::npos )
        {
            symmetryType_ |= symmetryTypes_::NEGY;
        }

        if( symmType.find("posZ") != word::npos )
        {
            symmetryType_ |= symmetryTypes_::POSZ;
        }

        if( symmType.find("negZ") != word::npos )
        {
            symmetryType_ |= symmetryTypes_::NEGZ;
        }
    }

    if( rollSetup == "none" )
    {
        Info << "No rollers considered in this case" << endl;
        isRollerSetup_ = false;

        return;
    }
    else if( rollSetup == "singleRoll" )
    {
        if
        (
            !meshDict.found("singleRollDict") ||
            !meshDict.isDict("singleRollDict")
        )
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "singleRollDict dictionary does not exist in meshDict!"
              << exit(FatalError);

        const dictionary& dict = meshDict.subDict("singleRollDict");

        //- read clearance if it exists
        dict.readIfPresent("initialRollGapClearance", initialRollGapClearance_);

        //- read singleRollerContactPatch
        dict.readIfPresent
        (
            "singleRollerContactPatch",
            singleRollerContactPatch_
        );

        //- read contact length
        if( !dict.readIfPresent("contactLength", contactLength_) )
        {
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "contactLength is not specified in singleRollDict"
              << exit(FatalError);
        }

        //- check the number of roll setup
        const wordList rolls = dict.toc();

        //- check if the roll position exists
        bool found(false);
        forAllConstIter(wordHashSet, rollGeometryInfo::rollPositions, it)
        {
            if( dict.found(it.key()) && dict.isDict(it.key()) )
            {
                if( found )
                {
                    FatalErrorIn
                    (
                        "void rollingMillMesh::rollingMillGeometryHandler::"
                        "parseDictionaries(const dictionary& meshDict)"
                    ) << "More that one roll position"
                      << " specified in singleRollDict" << exit(FatalError);
                }

                parseRollDictionary(dict.subDict(it.key()));
                found = true;
            }
        }

        if( !found )
        {
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "Roll position " << rolls[0] << " is not a valid name."
              << " Valid names are " << rollGeometryInfo::rollPositions
              << exit(FatalError);
        }
    }
    else if( rollSetup == "twoHighRolls" )
    {
        if
        (
            !meshDict.found("twoHighRollsDict") ||
            !meshDict.isDict("twoHighRollsDict")
        )
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "twoHighRollsDict dictionary does not exist in meshDict!"
              << exit(FatalError);

        const dictionary& dict = meshDict.subDict("twoHighRollsDict");

        //- read clearance if it exists
        dict.readIfPresent("initialRollGapClearance", initialRollGapClearance_);

        //- read singleRollerContactPatch
        dict.readIfPresent
        (
            "singleRollerContactPatch",
            singleRollerContactPatch_
        );

        //- read contact length
        if( !dict.readIfPresent("contactLength", contactLength_) )
        {
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "contactLength is not specified in twoHighRollsDict"
              << exit(FatalError);
        }

        //- read the information about the top roll
        if( !(symmetryType_ & NEGY) )
        {
            if( !dict.found("topRoll") || !dict.isDict("topRoll") )
                FatalErrorIn
                (
                    "void rollingMillMesh::rollingMillGeometryHandler::"
                    "parseDictionaries(const dictionary& meshDict)"
                ) << "topRoll dictionary does not exist in twoHighRollsDict!"
                  << exit(FatalError);

            parseRollDictionary(dict.subDict("topRoll"));
        }

        //- read the data about the bottom roll
        if( !(symmetryType_ & POSY) )
        {
            if
            (
                !dict.found("bottomRoll") ||
                !dict.isDict("bottomRoll")
            )
                FatalErrorIn
                (
                    "void rollingMillMesh::rollingMillGeometryHandler::"
                    "parseDictionaries(const dictionary& meshDict)"
                ) << "bottomRoll dictionary does not exist in twoHighRollsDict!"
                  << exit(FatalError);

            parseRollDictionary(dict.subDict("bottomRoll"));
        }

        if( dict.found("contact") )
        {
            Switch useContact(false);
            dict.readIfPresent("contact", useContact);

            if( useContact )
            {
                rollingMillContactHandler(meshDict, geometries_, patchHandler_);

                areRollsClosed_ = true;

                writeRollPositions();
            }
        }
    }
    else if( rollSetup == "edgeRolls" )
    {
        if
        (
            !meshDict.found("edgeRollsDict") ||
            !meshDict.isDict("edgeRollsDict")
        )
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "edgeRollsDict dictionary does not exist in meshDict!"
              << exit(FatalError);

        const dictionary& dict = meshDict.subDict("edgeRollsDict");

        //- read clearance if it exists
        dict.readIfPresent("initialRollGapClearance", initialRollGapClearance_);

        //- read singleRollerContactPatch
        dict.readIfPresent
        (
            "singleRollerContactPatch",
            singleRollerContactPatch_
        );

        //- read contact length
        if( !dict.readIfPresent("contactLength", contactLength_) )
        {
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "contactLength is not specified in edgeRollsDict"
              << exit(FatalError);
        }

        //- read the information about the left roll
        if( !(symmetryType_ & POSZ) )
        {
            if( !dict.found("leftRoll") || !dict.isDict("leftRoll") )
                FatalErrorIn
                (
                    "void rollingMillMesh::rollingMillGeometryHandler::"
                    "parseDictionaries(const dictionary& meshDict)"
                ) << "leftRoll dictionary does not exist in edgeRollsDict!"
                  << exit(FatalError);

            parseRollDictionary(dict.subDict("leftRoll"));
        }

        //- read the data about the right roll
        if( !(symmetryType_ & NEGZ) )
        {
            if
            (
                !dict.found("rightRoll") ||
                !dict.isDict("rightRoll")
            )
                FatalErrorIn
                (
                    "void rollingMillMesh::rollingMillGeometryHandler::"
                    "parseDictionaries(const dictionary& meshDict)"
                ) << "rightRoll dictionary does not exist in edgeRollsDict!"
                  << exit(FatalError);

            parseRollDictionary(dict.subDict("rightRoll"));
        }

        if( dict.found("contact") )
        {
            Switch useContact(false);
            dict.readIfPresent("contact", useContact);

            if( useContact )
            {
                rollingMillContactHandler(meshDict, geometries_, patchHandler_);

                areRollsClosed_ = true;

                writeRollPositions();
            }
        }

    }
    else if( rollSetup == "sideRolls" )
    {
        if
        (
            !meshDict.found("sideRollsDict") ||
            !meshDict.isDict("sideRollsDict")
        )
        {
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "sideRollsDict dictionary does not exist in meshDict!"
              << exit(FatalError);
        }

        const dictionary& dict = meshDict.subDict("sideRollsDict");

        //- read clearance if it exists
        dict.readIfPresent("initialRollGapClearance", initialRollGapClearance_);

        //- read singleRollerContactPatch
        dict.readIfPresent
        (
            "singleRollerContactPatch",
            singleRollerContactPatch_
        );

        //- read contact length
        if( !dict.readIfPresent("contactLength", contactLength_) )
        {
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "contactLength is not specified in sideRollsDict"
              << exit(FatalError);
        }

        //- read the information about the top roll
        if( !(symmetryType_ & NEGY) )
        {
            if( !dict.found("topRoll") || !dict.isDict("topRoll") )
                FatalErrorIn
                (
                    "void rollingMillMesh::rollingMillGeometryHandler::"
                    "parseDictionaries(const dictionary& meshDict)"
                ) << "topRoll dictionary does not exist in sideRollsDict!"
                  << exit(FatalError);

            parseRollDictionary(dict.subDict("topRoll"));
        }

        //- read the data about the bottom roll
        if( !(symmetryType_ & POSY) )
        {
            if
            (
                !dict.found("bottomRoll") ||
                !dict.isDict("bottomRoll")
            )
                FatalErrorIn
                (
                    "void rollingMillMesh::rollingMillGeometryHandler::"
                    "parseDictionaries(const dictionary& meshDict)"
                ) << "bottomRoll dictionary does not exist in sideRollsDict!"
                  << exit(FatalError);

            parseRollDictionary(dict.subDict("bottomRoll"));
        }

        //- read the information about the left roll
        if( !(symmetryType_ & POSZ) )
        {
            if( !dict.found("leftRoll") || !dict.isDict("leftRoll") )
                FatalErrorIn
                (
                    "void rollingMillMesh::rollingMillGeometryHandler::"
                    "parseDictionaries(const dictionary& meshDict)"
                ) << "leftRoll dictionary does not exist in sideRollsDict!"
                  << exit(FatalError);

            parseRollDictionary(dict.subDict("leftRoll"));
        }

        //- read the data about the right roll
        if( !(symmetryType_ & NEGZ) )
        {
            if
            (
                !dict.found("rightRoll") ||
                !dict.isDict("rightRoll")
            )
                FatalErrorIn
                (
                    "void rollingMillMesh::rollingMillGeometryHandler::"
                    "parseDictionaries(const dictionary& meshDict)"
                ) << "bottomRoll dictionary does not exist in sideRollsDict!"
                  << exit(FatalError);

            parseRollDictionary(dict.subDict("rightRoll"));
        }

        if( dict.found("contact") )
        {
            Switch useContact(false);
            dict.readIfPresent("contact", useContact);

            if( useContact )
            {
                rollingMillContactHandler(meshDict, geometries_, patchHandler_);
                areRollsClosed_ = true;

                writeRollPositions();
            }
        }
    }
    else if( rollSetup == "threeRolls" )
    {
        if( (symmetryType_ & POSY) || (symmetryType_ & NEGY) )
        {
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "symmetryType " << symmType << " specified in meshDict"
              << " is not possible for rollSetup " << rollSetup
              << exit(FatalError);
        }

        if
        (
            !meshDict.found("threeRollsDict") ||
            !meshDict.isDict("threeRollsDict")
        )
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "threeRollsDict dictionary does not exist in meshDict!"
              << exit(FatalError);

        const dictionary& dict = meshDict.subDict("threeRollsDict");

        //- read clearance if it exists
        dict.readIfPresent("initialRollGapClearance", initialRollGapClearance_);

        //- read singleRollerContactPatch
        dict.readIfPresent
        (
            "singleRollerContactPatch",
            singleRollerContactPatch_
        );

        //- read contact length
        if( !dict.readIfPresent("contactLength", contactLength_) )
        {
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "contactLength is not specified in threeRollsDict"
              << exit(FatalError);
        }

        if( symmetryType_ & (POSY|NEGY) )
        {
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "Is is not possible to use posY or negY symmetry"
              << " for threeRolls setup!"
              << exit(FatalError);
        }

        //- read the information about the top roll
        if( !dict.found("topRoll") || !dict.isDict("topRoll") )
            FatalErrorIn
            (
                "void rollingMillMesh::rollingMillGeometryHandler::"
                "parseDictionaries(const dictionary& meshDict)"
            ) << "topRoll dictionary does not exist in threeRollsDict!"
              << exit(FatalError);

        parseRollDictionary(dict.subDict("topRoll"));

        //- read the data about the bottom-left roll
        if( !(symmetryType_ & POSZ) )
        {
            if
            (
                !dict.found("bottomLeftRoll") ||
                !dict.isDict("bottomLeftRoll")
            )
                FatalErrorIn
                (
                    "void rollingMillMesh::rollingMillGeometryHandler::"
                    "parseDictionaries(const dictionary& meshDict)"
                ) << "bottomLeftRoll dictionary does not"
                  << " exist in threeRollsDict!"
                  << exit(FatalError);

            parseRollDictionary(dict.subDict("bottomLeftRoll"));
        }

        //- read the data about the bottom-right roll
        if( !(symmetryType_ & NEGZ) )
        {
            if
            (
                !dict.found("bottomRightRoll") ||
                !dict.isDict("bottomRightRoll")
            )
                FatalErrorIn
                (
                    "void rollingMillMesh::rollingMillGeometryHandler::"
                    "parseDictionaries(const dictionary& meshDict)"
                ) << "bottomRightRoll dictionary does not"
                  << " exist in threeRollsDict!"
                  << exit(FatalError);

            parseRollDictionary(dict.subDict("bottomRightRoll"));
        }
    }
    else
    {
        FatalError << "Unknown roll setup " << rollSetup << " in meshDict"
            << ". Please check your setup and try again" << exit(FatalError);
    }

    if( geometries_.size() == 0 )
    {
        FatalErrorIn
        (
            "void rollingMillMesh::rollingMillGeometryHandler::"
            "parseDictionaries(const dictionary& meshDict)"
        ) << "No valid geometries present in the current setup!"
          << " rollSetup " << rollSetup << " and symmetryType " << symmType
          << " do not generate a valid setup"
          << exit(FatalError);
    }

    Info << "Finished parsing roll setup from dictionary" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollingMillGeometryHandler::rollingMillGeometryHandler
(
    const Time& runTime,
    const dictionary& meshDict,
    const rollingMillPatchNamesHandler& patchHandler
)
:
    runTime_(runTime),
    patchHandler_(patchHandler),
    geometries_(),
    wireGeomPtr_(NULL),
    dieGeomPtr_(NULL),
    wireLength_(0.0),
    wireDiameter_(0.0),
    contactLength_(-1.0),
    extraContactAngle_(10.0),
    initialRollGapClearance_(-1.0),
    geometryTolerance_(1e-4),
    isRollerSetup_(true),
    areRollsClosed_(false),
    centreInAxialDirection_(false),
    singleRollerContactPatch_(false),
    symmetryType_(symmetryTypes_::NONE)
{
    parseDictionaries(meshDict);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollingMillGeometryHandler::~rollingMillGeometryHandler()
{
    deleteDemandDrivenData(wireGeomPtr_);
    deleteDemandDrivenData(dieGeomPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const PtrList<rollGeometryInfo>&
rollingMillGeometryHandler::rollPositions() const
{
    return geometries_;
}

const wireGeometryInfo& rollingMillGeometryHandler::wireGeometry() const
{
    if( !wireGeomPtr_ )
    {
        FatalErrorIn
        (
            "const rollGeometryInfo& rollingMillGeometryHandler"
            "::wireGeometry() const"
        ) << "Wire geometry does not exist" << abort(FatalError);
    }

    return *wireGeomPtr_;
}

const dieGeometryInfo& rollingMillGeometryHandler::dieGeometry() const
{
    if( !dieGeomPtr_ )
    {
        FatalErrorIn
        (
            "const dieGeometryInfo& rollingMillGeometryHandler"
            "::dieGeometry() const"
        ) << "The geometry of a die does not exist" << abort(FatalError);
    }

    return *dieGeomPtr_;
}

scalar rollingMillGeometryHandler::geometryTolerance() const
{
    return geometryTolerance_;
}

scalar rollingMillGeometryHandler::wireLength() const
{
    return wireLength_;
}

scalar rollingMillGeometryHandler::wireDiameter() const
{
    return wireDiameter_;
}

scalar rollingMillGeometryHandler::contactLength() const
{
    return contactLength_;
}

scalar rollingMillGeometryHandler::extraAngle() const
{
    return extraContactAngle_;
}

scalar rollingMillGeometryHandler::initialRollGapClearance() const
{
    return initialRollGapClearance_;
}

bool rollingMillGeometryHandler::isRollSetupValid() const
{
    return isRollerSetup_;
}

bool rollingMillGeometryHandler::areRollsClosed() const
{
    return areRollsClosed_;
}

bool rollingMillGeometryHandler::singleRollerContactPatch() const
{
    return singleRollerContactPatch_;
}

void rollingMillGeometryHandler::parseWireDictionary(const dictionary& meshDict)
{
    Info << "Starting parsing wire dictionary " << endl;

    if( !meshDict.found("wireMeshDict") || !meshDict.isDict("wireMeshDict") )
        FatalErrorIn
        (
            "void rollingMillGeometryHandler::"
            "parseWireDictionary(const dictionary&)"
        ) << "wireMeshDict does not exist in meshDict" << exit(FatalError);

    const dictionary& dict = meshDict.subDict("wireMeshDict");

    //- check/read geometry tolerance
    if( dict.found("geometryTolerance") )
    {
        geometryTolerance_ = readScalar(dict.lookup("geometryTolerance"));
    }

    if( dict.found("geometryTolerance") )
        dict.readIfPresent("geometryTolerance", geometryTolerance_);

    //- surface mesh of the cross section. Not used for circular cross section
    triSurf* surfPtr(NULL);

    //- the name of the wall contact patch
    word wallPatchName =
        patchHandler_.patchNameForWire
        (
            rollingMillPatchNamesHandler::WIRECONTACT
        );

    bool wireCoating(false);

    if( dict.found("coating") )
    {
        wireCoating = readBool(dict.lookup("coating"));
    }

    //- shall the profile be translated to the origin of the coordinate system
    bool wireProfileCentering(false);
    dict.readIfPresent("wireProfileCentering", wireProfileCentering);

    //- write the profile of the wire
    bool writeWireProfile(false);
    dict.readIfPresent("writeWireProfile", writeWireProfile);

    //- wedge angle in case of axi-symmetric setup
    scalar wedgeAngle(2.5);

    //- check if the setup is symmetric
    direction symmType = NONE;
    word symmetryCodeInDict = "none";

    if( meshDict.found("symmetryType") )
    {
        symmetryCodeInDict = word(meshDict.lookup("symmetryType"));
    }

    //- process symmetry setup into a code
    if( symmetryCodeInDict.find("posY") != word::npos )
    {
        symmType |= POSY;
    }

    if( symmetryCodeInDict.find("negY") != word::npos )
    {
        symmType |= NEGY;
    }

    if( symmetryCodeInDict.find("posZ") != word::npos )
    {
        symmType |= POSZ;
    }

    if( symmetryCodeInDict.find("negZ") != word::npos )
    {
        symmType |= NEGZ;
    }

    if( symmetryCodeInDict.find("axiSymmetric") != word::npos )
    {
        symmType |= AXISYMMETRIC;
    }

    //- check for the wedge angle
    if( symmType & AXISYMMETRIC )
    {
        //- read the wedge angle from wireMeshDict
        if( dict.found("wedgeAngle") )
            wedgeAngle = readScalar(dict.lookup("wedgeAngle"));

        //- override with wedge angle specified in dieMeshDict
        if( meshDict.found("dieMeshDict") && meshDict.isDict("dieMeshDict") )
        {
            const dictionary& dieDict = meshDict.subDict("dieMeshDict");

            if( dieDict.found("wedgeAngle") )
            {
                wedgeAngle = readScalar(dieDict.lookup("wedgeAngle"));
            }
        }
    }

    //- check the existence of rollSetup
    if( !meshDict.found("rollSetup") )
    {
        FatalErrorIn
        (
            "void rollingMillGeometryHandler::"
            "parseWireDictionary(const dictionary&)"
        ) << "rollSetup does not exist in meshDict" << exit(FatalError);
    }

    //- read wire geometry
    if( dict.found("wireDiameter") )
    {
        //- the cross section is a circle with the given diameter
        wireDiameter_ = readScalar(dict.lookup("wireDiameter"));
    }
    else if
    (
        dict.found("dxfFile") ||
        dict.found("geometryFile") ||
        dict.found("surfaceFile")
    )
    {
        //- wire profile is given as a DXF or a geometry file
        wireCrossSectionSurfaceCreator surfCreator
        (
            dict,
            patchHandler_,
            geometryTolerance_
        );

        surfPtr = new triSurf(surfCreator.surface());
    }

    if( writeWireProfile )
    {
        surfPtr->writeSurface("wireProfile.stl");
    }

    //- cut the geometry to achieve correct symmetry
    if( surfPtr && !cutGeometry(*surfPtr, "wire") )
    {
        FatalError << "Cannot cut wire geometry."
          << " Please check whether your geometry matches"
          << " you symmetry requirement." << exit(FatalError);
    }

    wireGeomPtr_ = new wireGeometryInfo(dict, surfPtr);

    //- set the name of the wall patch
    wireGeomPtr_->wallPatchName() = wallPatchName;

    //- set the existence of coating
    wireGeomPtr_->isCoatingPresent() = wireCoating;

    //- set the wire diameter
    wireGeomPtr_->wireDiameter() = wireDiameter_;

    //- set the type of symmetry
    wireGeomPtr_->typeOfSymmetry() = symmType;

    //- set the wedge angle
    wireGeomPtr_->wedgeAngle() = wedgeAngle;

    //- read the wire length if it exists
    if( dict.found("wireLength") )
    {
        wireLength_ = readScalar(dict.lookup("wireLength"));
    }

    //- read axial resolution
    if( dict.found("axialResolution") )
    {
        wireGeomPtr_->axialResolution() =
            readLabel(dict.lookup("axialResolution"));
    }

    //- read axial shift
    if( dict.found("axialShift") )
    {
        wireGeomPtr_->axialShift().x() = readScalar(dict.lookup("axialShift"));
    }

    //- read axial grading
    if( dict.found("axialGrading") )
    {
        wireGeomPtr_->axialGrading() = readScalar(dict.lookup("axialGrading"));
    }

    Info << "Finished parsing wire dictionary " << endl;
}

void rollingMillGeometryHandler::parseDieDictionary(const dictionary& meshDict)
{
    Info << "Starting parsing die's dictionary " << endl;

    if( !meshDict.found("dieMeshDict") || !meshDict.isDict("dieMeshDict") )
        FatalErrorIn
        (
            "void rollingMillGeometryHandler::"
            "parseDieDictionary(const dictionary&)"
        ) << "dieMeshDict does not exist in meshDict" << exit(FatalError);

    const dictionary& dict = meshDict.subDict("dieMeshDict");

    dieGeomPtr_ = new dieGeometryInfo(dict, *this, patchHandler_);

    Info << "Finished parsing die's dictionary " << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void rollingMillPatchNamesHandler::defaultRollerPatchNames()
{
    rollerPatchNames_[ROLLERCONTACT] = "rollerToWire";
    rollerPatchNames_[ROLLERAXIS] = "rollerAxis";
    rollerPatchNames_[ROLLERBACK] = "rollerBack";
    rollerPatchNames_[ROLLERFRONT] = "rollerFront";
    rollerPatchNames_[ROLLERSYMMY] = "rollerSymmetryY";
    rollerPatchNames_[ROLLERSYMMZ] = "rollerSymmetryZ";
    rollerPatchNames_[ROLLERTOROLLERBACK] = "rollerToRollerBack";
    rollerPatchNames_[ROLLERTOROLLERFRONT] = "rollerToRollerFront";
    rollerPatchNames_[ROLLERTOAIR] = "rollerToAir";

    if( meshDict_.found("rollingMillPatchNames") )
    {
        if( !meshDict_.isDict("rollingMillPatchNames") )
        {
            FatalError << "rollingMillPatchNames is not a dictionary"
                << exit(FatalError);
        }

        const dictionary& dict = meshDict_.subDict("rollingMillPatchNames");

        if( dict.found("rollerContactPatchName") )
        {
            rollerPatchNames_[ROLLERCONTACT] =
                word(dict.lookup("rollerContactPatchName"));
        }

        if( dict.found("rollerAxisPatchName") )
        {
            rollerPatchNames_[ROLLERAXIS] =
                word(dict.lookup("rollerAxisPatchName"));
        }

        if( dict.found("rollerBackPatchName") )
        {
            rollerPatchNames_[ROLLERBACK] =
                word(dict.lookup("rollerBackPatchName"));
        }

        if( dict.found("rollerFrontPatchName") )
        {
            rollerPatchNames_[ROLLERFRONT] =
                word(dict.lookup("rollerFrontPatchName"));
        }

        if( dict.found("rollerSymmetryYPatchName") )
        {
            rollerPatchNames_[ROLLERSYMMY] =
                word(dict.lookup("rollerSymmetryYPatchName"));
        }

        if( dict.found("rollerSymmetryZPatchName") )
        {
            rollerPatchNames_[ROLLERSYMMZ] =
                word(dict.lookup("rollerSymmetryZPatchName"));
        }

        if( dict.found("rollerToRollerBackPatchName") )
        {
            rollerPatchNames_[ROLLERTOROLLERBACK] =
                word(dict.lookup("rollerToRollerBackPatchName"));
        }

        if( dict.found("rollerToRollerFrontPatchName") )
        {
            rollerPatchNames_[ROLLERTOROLLERFRONT] =
                word(dict.lookup("rollerToRollerFrontPatchName"));
        }

        if( dict.found("rollerToAirPatchName") )
        {
            rollerPatchNames_[ROLLERTOAIR] =
                word(dict.lookup("rollerToAirPatchName"));
        }
    }
}

void rollingMillPatchNamesHandler::defaultWirePatchNames()
{
    wirePatchNames_[WIRECONTACT] = "wireContact";
    wirePatchNames_[WIRECOATING] = "wireToCoating";
    wirePatchNames_[WIREUPSTREAM] = "wireUpstream";
    wirePatchNames_[WIREDOWNSTREAM] = "wireDownstream";
    wirePatchNames_[WIRESYMMY] = "wireSymmetryY";
    wirePatchNames_[WIRESYMMZ] = "wireSymmetryZ";
    wirePatchNames_[WIREBACK] = "wireBack";
    wirePatchNames_[WIREFRONT] = "wireFront";

    if( meshDict_.found("rollingMillPatchNames") )
    {
        if( !meshDict_.isDict("rollingMillPatchNames") )
        {
            FatalError << "rollingMillPatchNames is not a dictionary"
                << exit(FatalError);
        }

        const dictionary& dict = meshDict_.subDict("rollingMillPatchNames");

        if( dict.found("wireContactPatchName") )
        {
            wirePatchNames_[WIRECONTACT] =
                word(dict.lookup("wireContactPatchName"));
        }

        if( dict.found("wireCoatingPatchName") )
        {
            wirePatchNames_[WIRECOATING] =
                word(dict.lookup("wireCoatingPatchName"));
        }

        if( dict.found("wireUpstreamPatchName") )
        {
            wirePatchNames_[WIREUPSTREAM] =
                word(dict.lookup("wireUpstreamPatchName"));
        }

        if( dict.found("wireDownstreamPatchName") )
        {
            wirePatchNames_[WIREDOWNSTREAM] =
                word(dict.lookup("wireDownstreamPatchName"));
        }

        if( dict.found("wireSymmetryYPatchName") )
        {
            wirePatchNames_[WIRESYMMY] =
                word(dict.lookup("wireSymmetryYPatchName"));
        }

        if( dict.found("wireSymmetryZPatchName") )
        {
            wirePatchNames_[WIRESYMMZ] =
                word(dict.lookup("wireSymmetryZPatchName"));
        }

        if( dict.found("wireFrontPatchName") )
        {
            wirePatchNames_[WIREFRONT] =
                word(dict.lookup("wireFrontPatchName"));
        }

        if( dict.found("wireBackPatchName") )
        {
            wirePatchNames_[WIREBACK] =
                word(dict.lookup("wireBackPatchName"));
        }
    }
}

void rollingMillPatchNamesHandler::defaultWireCoatingPatchNames()
{
    coatingPatchNames_[COATINGCONTACT] = "coatingContact";
    coatingPatchNames_[COATINGWIRE] = "coatingToWire";
    coatingPatchNames_[COATINGUPSTREAM] = "coatingUpstream";
    coatingPatchNames_[COATINGDOWNSTREAM] = "coatingDownstream";
    coatingPatchNames_[COATINGSYMMY] = "coatingSymmetryY";
    coatingPatchNames_[COATINGSYMMZ] = "coatingSymmetryZ";
    coatingPatchNames_[COATINGFRONT] = "coatingFront";
    coatingPatchNames_[COATINGBACK] = "coatingBack";

    if( meshDict_.found("rollingMillPatchNames") )
    {
        if( !meshDict_.isDict("rollingMillPatchNames") )
        {
            FatalError << "rollingMillPatchNames is not a dictionary"
                << exit(FatalError);
        }

        const dictionary& dict = meshDict_.subDict("rollingMillPatchNames");

        if( dict.found("coatingContactPatchName") )
        {
            coatingPatchNames_[COATINGCONTACT] =
                word(dict.lookup("coatingContactPatchName"));
        }

        if( dict.found("coatingToWirePatchName") )
        {
            coatingPatchNames_[COATINGWIRE] =
                word(dict.lookup("coatingToWirePatchName"));
        }

        if( dict.found("coatingUpstreamPatchName") )
        {
            coatingPatchNames_[COATINGUPSTREAM] =
                word(dict.lookup("coatingUpstreamPatchName"));
        }

        if( dict.found("coatingDownstreamPatchName") )
        {
            coatingPatchNames_[COATINGDOWNSTREAM] =
                word(dict.lookup("coatingDownstreamPatchName"));
        }

        if( dict.found("coatingSymmetryYPatchName") )
        {
            coatingPatchNames_[COATINGSYMMY] =
                word(dict.lookup("coatingSymmetryYPatchName"));
        }

        if( dict.found("coatingSymmetryZPatchName") )
        {
            coatingPatchNames_[COATINGSYMMZ] =
                word(dict.lookup("coatingSymmetryZPatchName"));
        }

        if( dict.found("coatingFrontPatchName") )
        {
            coatingPatchNames_[COATINGFRONT] =
                word(dict.lookup("coatingFrontPatchName"));
        }

        if( dict.found("coatingBackPatchName") )
        {
            coatingPatchNames_[COATINGBACK] =
                word(dict.lookup("coatingBackPatchName"));
        }
    }
}

void rollingMillPatchNamesHandler::defaultDiePatchNames()
{
    diePatchNames_[DIEWIRE] = "dieToWire";
    diePatchNames_[DIEHOUSING] = "dieToHousing";
    diePatchNames_[DIEUPSTREAM] = "dieUpstream";
    diePatchNames_[DIEDOWNSTREAM] = "dieDownstream";
    diePatchNames_[DIESYMMY] = "dieSymmetryY";
    diePatchNames_[DIESYMMZ] = "dieSymmetryZ";
    diePatchNames_[DIEFRONT] = "dieFront";
    diePatchNames_[DIEBACK] = "dieBack";
    diePatchNames_[DIEENTRANCECONE] = "dieEntranceCone";
    diePatchNames_[DIEEXITCONE] = "dieExitCone";

    if( meshDict_.found("rollingMillPatchNames") )
    {
        if( !meshDict_.isDict("rollingMillPatchNames") )
        {
            FatalError << "rollingMillPatchNames is not a dictionary"
                << exit(FatalError);
        }

        const dictionary& dict = meshDict_.subDict("rollingMillPatchNames");

        if( dict.found("dieToWirePatchName") )
        {
            diePatchNames_[DIEWIRE] =
                word(dict.lookup("dieToWirePatchName"));
        }

        if( dict.found("dieToHousingPatchName") )
        {
            diePatchNames_[DIEHOUSING] =
                word(dict.lookup("dieToHousingPatchName"));
        }

        if( dict.found("dieUpstreamPatchName") )
        {
            diePatchNames_[DIEUPSTREAM] =
                word(dict.lookup("dieUpstreamPatchName"));
        }

        if( dict.found("dieDownstreamPatchName") )
        {
            diePatchNames_[DIEDOWNSTREAM] =
                word(dict.lookup("dieDownstreamPatchName"));
        }

        if( dict.found("dieSymmetryYPatchName") )
        {
            diePatchNames_[DIESYMMY] =
                word(dict.lookup("dieSymmetryYPatchName"));
        }

        if( dict.found("dieSymmetryZPatchName") )
        {
            diePatchNames_[DIESYMMZ] =
                word(dict.lookup("dieSymmetryZPatchName"));
        }

        if( dict.found("dieFrontPatchName") )
        {
            diePatchNames_[DIEFRONT] = word(dict.lookup("dieFrontPatchName"));
        }

        if( dict.found("dieBackPatchName") )
        {
            diePatchNames_[DIEBACK] = word(dict.lookup("dieBackPatchName"));
        }

        if( dict.found("dieEntranceConePatchName") )
        {
            diePatchNames_[DIEENTRANCECONE] =
                word(dict.lookup("dieEntranceConePatchName"));
        }

        if( dict.found("dieExitConePatchName") )
        {
            diePatchNames_[DIEEXITCONE] =
                word(dict.lookup("dieExitConePatchName"));
        }
    }
}

void rollingMillPatchNamesHandler::defaultCasingPatchNames()
{
    casingPatchNames_[CASINGDOWNSTREAM] = "casingDownstream";
    casingPatchNames_[CASINGUPSTREAM] = "casingUpstream";
    casingPatchNames_[CASINGTOOUTSIDE] = "casingToOutside";
    casingPatchNames_[CASINGENTRANCE] = "casingEntrance";
    casingPatchNames_[CASINGTODIERADIAL] = "casingToDieRadial";
    casingPatchNames_[CASINGTODIEAXIAL] = "casingToDieAxial";
    casingPatchNames_[CASINGEXIT] = "casingExit";
    casingPatchNames_[CASINGSYMMY] = "casingSymmY";
    casingPatchNames_[CASINGSYMMZ] = "casingSymmZ";
    casingPatchNames_[CASINGFRONT] = "casingFront";
    casingPatchNames_[CASINGBACK] = "casingBack";

    if( meshDict_.found("rollingMillPatchNames") )
    {
        if( !meshDict_.isDict("rollingMillPatchNames") )
        {
            FatalError << "rollingMillPatchNames is not a dictionary"
                << exit(FatalError);
        }

        const dictionary& dict = meshDict_.subDict("rollingMillPatchNames");

        if( dict.found("casingDownstreamPatchName") )
        {
            casingPatchNames_[CASINGDOWNSTREAM] =
                word(dict.lookup("casingDownstreamPatchName"));
        }

        if( dict.found("casingUpstreamPatchName") )
        {
            casingPatchNames_[CASINGUPSTREAM] =
                word(dict.lookup("casingUpstreamPatchName"));
        }

        if( dict.found("casingToOutsidePatchName") )
        {
            casingPatchNames_[CASINGTOOUTSIDE] =
                word(dict.lookup("casingToOutsidePatchName"));
        }

        if( dict.found("casingEntrancePatchName") )
        {
            casingPatchNames_[CASINGENTRANCE] =
                word(dict.lookup("casingEntrancePatchName"));
        }

        if( dict.found("casingToDieRadialPatchName") )
        {
            casingPatchNames_[CASINGTODIERADIAL] =
                word(dict.lookup("casingToDieRadialPatchName"));
        }

        if( dict.found("casingToDieAxialPatchName") )
        {
            casingPatchNames_[CASINGTODIEAXIAL] =
                word(dict.lookup("casingToDieAxialPatchName"));
        }

        if( dict.found("casingExitPatchName") )
        {
            casingPatchNames_[CASINGEXIT] =
                word(dict.lookup("casingExitPatchName"));
        }

        if( dict.found("casingSymmetryYPatchName") )
        {
            casingPatchNames_[CASINGSYMMY] =
                word(dict.lookup("casingSymmetryYPatchName"));
        }

        if( dict.found("casingSymmetryZPatchName") )
        {
            casingPatchNames_[CASINGSYMMZ] =
                word(dict.lookup("casingSymmetryZPatchName"));
        }

        if( dict.found("casingFrontPatchName") )
        {
            casingPatchNames_[CASINGFRONT] =
                word(dict.lookup("casingFrontPatchName"));
        }

        if( dict.found("casingBackPatchName") )
        {
            casingPatchNames_[CASINGBACK] =
                word(dict.lookup("casingBackPatchName"));
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollingMillPatchNamesHandler::rollingMillPatchNamesHandler
(
    const dictionary& meshDict
)
:
    meshDict_(meshDict),
    rollerPatchNames_(),
    wirePatchNames_(),
    coatingPatchNames_(),
    diePatchNames_(),
    casingPatchNames_()
{
    defaultRollerPatchNames();

    defaultWirePatchNames();

    defaultWireCoatingPatchNames();

    defaultDiePatchNames();

    defaultCasingPatchNames();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollingMillPatchNamesHandler::~rollingMillPatchNamesHandler()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

word rollingMillPatchNamesHandler::patchNameForRoll
(
    const rollerPatchKeys& patchKey
) const
{
    std::map<label, word>::const_iterator it = rollerPatchNames_.find(patchKey);

    if( it == rollerPatchNames_.end() )
    {
        FatalErrorIn
        (
            "word rollingMillPatchNamesHandler::patchNameForRoll"
            "(const rollerPatchKeys&) const"
        ) << "Patch position " << label(patchKey)
          << " does not exist!" << abort(FatalError);
    }

    return it->second;
}

word rollingMillPatchNamesHandler::patchTypeForRoll
(
    const rollerPatchKeys& patchKey,
    const rollingMillGeometryHandler::symmetryTypes_& symm
) const
{
    switch( patchKey )
    {
        case ROLLERBACK:
        {
            if
            (
                (symm & rollingMillGeometryHandler::symmetryTypes_::POSY) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::POSZ) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::NEGY) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::NEGZ)
            )
            {
                return "symmetryPlane";
            }
            else if
            (
                symm & rollingMillGeometryHandler::symmetryTypes_::AXISYMMETRIC
            )
            {
                return "wedge";
            }
            else
            {
                return "patch";
            }
        } break;
        case ROLLERFRONT:
        {
            if
            (
                (symm & rollingMillGeometryHandler::symmetryTypes_::POSY) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::POSZ) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::NEGY) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::NEGZ)
            )
            {
                return "symmetryPlane";
            }
            else if
            (
                symm & rollingMillGeometryHandler::symmetryTypes_::AXISYMMETRIC
            )
            {
                return "wedge";
            }
            else
            {
                return "patch";
            }
        } break;
        default:
        {
            return "patch";
        }
    };

    return "patch";
}

word rollingMillPatchNamesHandler::patchNameForWire
(
    const wirePatchKeys& patchKey
) const
{
    std::map<label, word>::const_iterator it = wirePatchNames_.find(patchKey);

    if( it == wirePatchNames_.end() )
    {
        FatalErrorIn
        (
            "word rollingMillPatchNamesHandler::patchNameForWire"
            "(const wirePatchKeys&) const"
        ) << "Patch position " << label(patchKey)
          << " does not exist!" << abort(FatalError);
    }

    return it->second;
}

word rollingMillPatchNamesHandler::patchTypeForWire
(
    const wirePatchKeys& patchKey,
    const rollingMillGeometryHandler::symmetryTypes_& symm
) const
{
    switch( patchKey )
    {
        case WIRESYMMY: case WIRESYMMZ:
        {
            if
            (
                (symm & rollingMillGeometryHandler::symmetryTypes_::POSY) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::POSZ) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::NEGY) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::NEGZ)
            )
            {
                return "symmetryPlane";
            }
            else
            {
                return "patch";
            }
        } break;
        case WIREBACK: case WIREFRONT:
        {
            if
            (
                symm & rollingMillGeometryHandler::symmetryTypes_::AXISYMMETRIC
            )
            {
                return "wedge";
            }
            else
            {
                return "patch";
            }
        } break;
        default:
        {
            return "patch";
        }
    };

    return "patch";
}

word rollingMillPatchNamesHandler::patchNameForCoating
(
    const coatingPatchKeys& patchKey
) const
{
    std::map<label, word>::const_iterator it =
        coatingPatchNames_.find(patchKey);

    if( it == coatingPatchNames_.end() )
    {
        FatalErrorIn
        (
            "word rollingMillPatchNamesHandler::patchNameForCoating"
            "(const coatingPatchKeys&) const"
        ) << "Patch position " << label(patchKey)
          << " does not exist!" << abort(FatalError);
    }

    return it->second;
}

word rollingMillPatchNamesHandler::patchTypeForCoating
(
    const coatingPatchKeys& patchKey,
    const rollingMillGeometryHandler::symmetryTypes_& symm
) const
{
    switch( patchKey )
    {
        case COATINGSYMMY: case COATINGSYMMZ:
        {
            if
            (
                (symm & rollingMillGeometryHandler::symmetryTypes_::POSY) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::POSZ) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::NEGY) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::NEGZ)
            )
            {
                return "symmetryPlane";
            }
            else
            {
                return "patch";
            }
        } break;
        case COATINGBACK: case COATINGFRONT:
        {
            if
            (
                symm & rollingMillGeometryHandler::symmetryTypes_::AXISYMMETRIC
            )
            {
                return "wedge";
            }
            else
            {
                return "patch";
            }
        } break;
        default:
        {
            return "patch";
        }
    };

    return "patch";
}


word rollingMillPatchNamesHandler::patchNameForDie
(
    const diePatchKeys& patchKey
) const
{
    std::map<label, word>::const_iterator it = diePatchNames_.find(patchKey);

    if( it == diePatchNames_.end() )
    {
        FatalErrorIn
        (
            "word rollingMillPatchNamesHandler::patchNameForDie"
            "(const diePatchKeys&) const"
        ) << "Patch position " << label(patchKey)
          << " does not exist!" << abort(FatalError);
    }

    return it->second;
}

word rollingMillPatchNamesHandler::patchTypeForDie
(
    const diePatchKeys& patchKey,
    const rollingMillGeometryHandler::symmetryTypes_& symm
) const
{
    switch( patchKey )
    {
        case DIESYMMY: case DIESYMMZ:
        {
            if
            (
                (symm & rollingMillGeometryHandler::symmetryTypes_::POSY) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::POSZ) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::NEGY) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::NEGZ)
            )
            {
                return "symmetryPlane";
            }
            else
            {
                return "patch";
            }
        } break;
        case DIEFRONT: case DIEBACK:
        {
            if
            (
                symm & rollingMillGeometryHandler::symmetryTypes_::AXISYMMETRIC
            )
            {
                return "wedge";
            }
            else
            {
                return "patch";
            }
        } break;
        default:
        {
            return "patch";
        }
    };

    return "patch";
}

word rollingMillPatchNamesHandler::patchNameForCasing
(
    const casingPatchKeys& patchKey
) const
{
    std::map<label, word>::const_iterator it = casingPatchNames_.find(patchKey);

    if( it == casingPatchNames_.end() )
    {
        FatalErrorIn
        (
            "word rollingMillPatchNamesHandler::patchNameForCasing"
            "(const casingPatchKeys&) const"
        ) << "Patch position " << label(patchKey)
          << " does not exist!" << abort(FatalError);
    }

    return it->second;
}

word rollingMillPatchNamesHandler::patchTypeForCasing
(
    const casingPatchKeys& patchKey,
    const rollingMillGeometryHandler::symmetryTypes_& symm
) const
{
    switch( patchKey )
    {
        case CASINGSYMMY: case CASINGSYMMZ:
        {
            if
            (
                (symm & rollingMillGeometryHandler::symmetryTypes_::POSY) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::POSZ) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::NEGY) ||
                (symm & rollingMillGeometryHandler::symmetryTypes_::NEGZ)
            )
            {
                return "symmetryPlane";
            }
            else
            {
                return "patch";
            }
        } break;
        case CASINGFRONT: case CASINGBACK:
        {
            if
            (
                symm & rollingMillGeometryHandler::symmetryTypes_::AXISYMMETRIC
            )
            {
                return "wedge";
            }
            else
            {
                return "patch";
            }
        } break;
        default:
        {
            return "patch";
        }
    };

    return "patch";
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
