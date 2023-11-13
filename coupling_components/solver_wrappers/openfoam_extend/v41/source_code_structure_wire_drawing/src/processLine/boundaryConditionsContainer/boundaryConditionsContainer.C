/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "boundaryConditionsContainer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boundaryConditionsContainer, 0);
}


// * * * * * * * * * *  Private Member Functions  * * * * * * * * * * * * * * //

void Foam::boundaryConditionsContainer::printDefaultBCInfo
(
    word patchName, word BCtype
)
{
    Info<< "            " << patchName << ": set to " << string(BCtype)
        << " default boundary condition" << endl;
}


void Foam::boundaryConditionsContainer::printCopyBCInfo(word patchName)
{
    Info<< "            Boundary condition subDict for patch " << patchName
        << " will be copied" << endl;
}


bool Foam::boundaryConditionsContainer::copyDict(word patchName)
{
    // Check if a full dictionary is specified which can be copied.
    // This is assumed when there is a dictionary with keyword = patchName
    // specified and contains the keyword type.
    bool copyDict = false;
    if (dict_.found(patchName))
    {
        if (dict_.subDict(patchName).found("type"))
        {
            copyDict = true;
        }
    }

    return copyDict;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryConditionsContainer::boundaryConditionsContainer
(
    dictionary& dict,
    dataHandler* dataContainerPtr
)
:
    dict_(dict),
    dataContainerPtr_(dataContainerPtr)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::boundaryConditionsContainer::~boundaryConditionsContainer()
{}


// * * * * * * * * * * * * * * * Public Member Functions  * * *  * * * * * * //


// Momentum boundary conditions

void Foam::boundaryConditionsContainer::solidSymmetry
(
    word patchName
)
{
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "solidSymmetry");
        dictionary BCdict;
        {
            BCdict.add("type", "solidSymmetry");
            BCdict.add("patchType", "symmetryPlane");
            BCdict.add("value", "uniform ( 0 0 0 )");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::fixedTangentialNormalPressure
(
    const word& patchName,
    const word& subDirectory
)
{
    // TODO: dict_ is the BC dict and can contain updateTangentialDisplacement
    // keyword. Instead of hardcoding it, make it a lookupOrDefault from dict_
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "fixedTangentialNormalPressure");
        dictionary BCdict;
        {
            BCdict.add("type", "fixedTangentialNormalPressure");

            const fileName timeVsDownstreamPressure =
                fileName
                (
                    "../../dictFiles"/subDirectory/"timeVsDownstreamPressure"
                );

            BCdict.add("fileName", timeVsDownstreamPressure);
            BCdict.add("outOfBounds", "clamp");
            BCdict.add("updateTangentialDisplacement", "no");
            BCdict.add("updateTangentialDisplacementEndTime", "0.01");
            BCdict.add("relaxationFactor", "0.025");
            BCdict.add("normal", "(1 0 0)");
            BCdict.add("value", "uniform ( 0 0 0 )");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::fixedDisplacementZeroShear
(
    const word& patchName,
    const word& subDirectory
)
{
    if (!copyDict(patchName))
    {
        // Set default input parameters
        // Drawing speed (in m/s)
        scalar drawingSpeed = 5.0;

        // Change default values in case there are input values given
        if (dict_.found(patchName))
        {
            // Check if there are input parameters specified and read them
            // drawing speed
            if (dict_.subDict(patchName).found("drawingSpeed"))
            {
                drawingSpeed =
                    readScalar(dict_.subDict(patchName).lookup("drawingSpeed"));

                WarningIn(type() + "::fixedDisplacementZeroShear(...)")
                    << "drawingSpeed is currently not used!" << endl;
            }
            else
            {
                Info<< "            No drawing speed found. Default value "
                    << drawingSpeed << " will be used." << endl;
            }
        }

        // Write the boundary conditions
        printDefaultBCInfo(patchName, "fixedDisplacementZeroShear");
        dictionary BCdict;
        {
            BCdict.add("type", "fixedDisplacementZeroShear");
            dictionary displacementSeries;
            {
                WarningIn
                (
                    "void Foam::boundaryConditionsContainer::"
                    "fixedDisplacementZeroShear"
                )   << "Reading timeVsDisp from dictFiles" << endl;

                const fileName timeVsDisp =
                    fileName("../../dictFiles"/subDirectory/"timeVsDisp");

                displacementSeries.add("fileName", timeVsDisp);
                displacementSeries.add("outOfBounds", "clamp");
            }
            BCdict.add("displacementSeries", displacementSeries);
            BCdict.add("value", "uniform ( 0 0 0 )");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::solidContactSlave
(
    word patchName, word shadowPatchName
)
{
    if (!copyDict(patchName))
    {
        const dictionary patchNames =
                dataContainer().dataContainerDict().subDict("patchNames");

        printDefaultBCInfo(patchName, "solidContactSlave");
        dictionary BCdict;
        {
            BCdict.add("type", "solidContact");
            BCdict.add("shadowPatch", shadowPatchName);
//            if
//            (
//                shadowPatchName ==
//                word(patchNames.lookup("rollerContactPatchName"))
//            )
            {
                BCdict.add
                (
                    "regionOfInterestBottomCorner",
                    "(-0.04 -0.01 -0.015)"
                );
                BCdict.add
                (
                    "regionOfInterestTopCorner",
                    "(0.04 0.01 0.015)"
                );
            }

            BCdict.add("master", "no");
            BCdict.add("value", "uniform ( 0 0 0 )");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::solidContactMaster
(
    const word& patchName,
    const word& shadowPatchName,
    const scalar frictionCoefficientAxial,
    const scalar targetAveragePenetration,
    const scalar minPenaltyScale,
    const scalar maxPenaltyScale
)
{
    if (!copyDict(patchName))
    {
        const dictionary patchNames =
                dataContainer().dataContainerDict().subDict("patchNames");

        printDefaultBCInfo(patchName, "solidContactMaster");

        dictionary BCdict;
        {
            BCdict.add("type", "solidContact");
            BCdict.add("useNewPointDistanceMethod", "no");
            BCdict.add("usePrevCandidateMasterNeighbors", "no");
            BCdict.add("master", "yes");
            BCdict.add("rigidMaster", "no");
            BCdict.add("shadowPatch", shadowPatchName);
//            BCdict.add("interpolationMethod", "ggi");

//            if
//            (
//                patchName ==
//                word(patchNames.lookup("rollerContactPatchName"))
//            )
            {
                BCdict.add
                (
                    "regionOfInterestBottomCorner",
                    "(-0.04 -0.01 -0.015)"
                );
                BCdict.add
                (
                    "regionOfInterestTopCorner",
                    "(0.04 0.01 0.015)"
                );
            }

            BCdict.add("normalContactModel", "standardPenalty");
            dictionary normalModelDict;
            {
                normalModelDict.add("penaltyScale", "1");
                normalModelDict.add("relaxationFactor", "0.05");
                normalModelDict.add("infoFrequency", "10");
                normalModelDict.add("adjustPenaltyScale", "yes");
                dictionary adjustPSDict;
                {
                    adjustPSDict.add("maxPenaltyScale", maxPenaltyScale);
                    adjustPSDict.add("minPenaltyScale", minPenaltyScale);
                    adjustPSDict.add
                    (
                        "targetAveragePenetration", targetAveragePenetration
                    );
                }
                normalModelDict.add("adjustPenaltyScaleDict", adjustPSDict);
            }
            BCdict.add("standardPenaltyNormalModelDict", normalModelDict);

            BCdict.add("frictionContactModel", "standardPenalty");
            dictionary frictionModelDict;
            {
                frictionModelDict.add("penaltyScale", "1");
                frictionModelDict.add("relaxationFactor", "0.02");
                frictionModelDict.add("infoFrequency", "10");
//                frictionModelDict.add("frictionLaw", "anisotropicCoulomb");
                frictionModelDict.add("frictionLaw", "coulombOrowan");
//                frictionModelDict.add("frictionLaw", "stribeckCoulomb");
//                frictionModelDict.add("frictionLaw", "coulomb");
                dictionary frictionLawDict;
                {
                    frictionLawDict.add
                    (
                        "frictionCoeff",
                        frictionCoefficientAxial
                    );
                    frictionLawDict.add
                    (
                        "primaryFrictionCoeff",
                        0.12
                    );
                    frictionLawDict.add
                    (
                        "secondaryFrictionCoeff",
                        0.02
                    );
                    frictionLawDict.add
                    (
                        "primaryDirection",
                        "(1 0 0)"
                    );
                    frictionLawDict.add
                    (
                        "shearStressLimit",
                        0.3e9
                    );
                }
                frictionModelDict.add("frictionLawDict", frictionLawDict);
            }
            BCdict.add("standardPenaltyFrictionModelDict", frictionModelDict);
            BCdict.add("value", "uniform (0 0 0)");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}

void Foam::boundaryConditionsContainer::solidContactMaster
(
    word patchName,
    wordList shadowPatchNames,
    scalarList frictionCoefficientsAxial
)
{
    if (!copyDict(patchName))
    {
        const dictionary patchNames =
                dataContainer().dataContainerDict().subDict("patchNames");
        // TODO: thermalContactMaster 1-many!!
        printDefaultBCInfo(patchName, "solidContactMaster");
        dictionary BCdict;
        {
            BCdict.add("type", "solidContact");
            BCdict.add("master", "yes");
            BCdict.add("rigidMaster", "no");
            BCdict.add("shadowPatches", shadowPatchNames);
//            BCdict.add("interpolationMethod", "ggi");

            BCdict.add
            (
                "regionOfInterestBottomCorner",
                "(-0.04 -0.01 -0.015)"
            );
            BCdict.add
            (
                "regionOfInterestTopCorner",
                "(0.04 0.01 0.015)"
            );

            dictionary contactModelDict;
            forAll (shadowPatchNames, patchI)
            {
                Info<< "friction coefficient: "
                    << frictionCoefficientsAxial[patchI] << endl;

                contactModelDict.clear();
                contactModelDict.add("normalContactModel", "standardPenalty");
                dictionary normalModelDict;
                {
                    normalModelDict.add("penaltyScale", "1");
                    normalModelDict.add("relaxationFactor", "0.05");
                    normalModelDict.add("infoFrequency", "10");
                    normalModelDict.add("adjustPenaltyScale", "yes");
                    dictionary adjustPSDict;
                    {
                        adjustPSDict.add("maxPenaltyScale", "20");
                        adjustPSDict.add("minPenaltyScale", "0.5");
                        adjustPSDict.add("targetAveragePenetration", "2e-6");
                    }
                    normalModelDict.add("adjustPenaltyScaleDict", adjustPSDict);
                }
                contactModelDict.add
                (
                    "standardPenaltyNormalModelDict", normalModelDict
                );
                contactModelDict.add("frictionContactModel", "standardPenalty");
                dictionary frictionModelDict;
                {
                    frictionModelDict.add("penaltyScale", "0.2");
                    frictionModelDict.add("relaxationFactor", "0.02");
                    frictionModelDict.add("infoFrequency", "10");
                    frictionModelDict.add("frictionLaw", "coulomb");
                    dictionary frictionLawDict;
                    {
                        frictionLawDict.add
                        (
                            "frictionCoeff",
                            frictionCoefficientsAxial[patchI]
                        );
                    }
                    frictionModelDict.add("frictionLawDict", frictionLawDict);
                }
                contactModelDict.add
                (
                    "standardPenaltyFrictionModelDict", frictionModelDict
                );

                BCdict.add
                (
                    word
                    (
                        patchName + "_to_" + shadowPatchNames[patchI] + "_dict"
                    ),
                    contactModelDict
                );
            }
            BCdict.add("value", "uniform (0 0 0)");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::solidRigidContact
(
    word patchName,
    dictionary triSurfacesDict,
//    word frictionLaw,
    scalarList frictionCoefficientsAxial
)
{
    if (!copyDict(patchName))
    {
        const dictionary patchNames =
                dataContainer().dataContainerDict().subDict("patchNames");

        // TODO: thermalContactMaster 1-many!!
        printDefaultBCInfo(patchName, "solidRigidContact");

        wordList triSurfDictToc = triSurfacesDict.toc();

        if (triSurfDictToc.size() == 0)
        {
            FatalErrorIn
            (
                "void Foam::boundaryConditionsContainer::solidRigidContact"
            )   << "No triSurfaces are defined!" << abort(FatalError);
        }

        wordList shadowPatchNames;
        shadowPatchNames.setSize(triSurfDictToc.size());
        shadowPatchNames = triSurfDictToc;

        dictionary BCdict;
        {
            BCdict.add("triSurfaces", triSurfacesDict);

            BCdict.add("type", "solidRigidContact");

            BCdict.add
            (
                "regionOfInterestBottomCorner",
                "(-0.04 -0.01 -0.015)"
            );
            BCdict.add
            (
                "regionOfInterestTopCorner",
                "(0.04 0.01 0.015)"
            );

            dictionary contactModelDict;
            forAll (shadowPatchNames, patchI)
            {
                const fileName stlFile =
                        word(shadowPatchNames[patchI]) + ".stl";

                word systemCommand;

                if
                (
                    word(shadowPatchNames[patchI]) == "topRoll" ||
                    word(shadowPatchNames[patchI]) == "bottomRoll" ||
                    word(shadowPatchNames[patchI]) == "rightRoll" ||
                    word(shadowPatchNames[patchI]) == "leftRoll"
                )
                {
                    systemCommand =
                        "surfaceMeshTriangulate "
                      + stlFile + " -patches \"("
                      + word(patchNames.lookup("rollerContactPatchName"))
                      + ")\" -region " + word(shadowPatchNames[patchI]);
                }
                else
                {
                    systemCommand =
                        "surfaceMeshTriangulate "
                      + stlFile + " -patches \"("
                      + word(patchNames.lookup("dieToWirePatchName"))
                      + ")\" -region " + word(shadowPatchNames[patchI]);
                }

                dataContainer().runSystemCommand
                (
                    systemCommand,
                    "log.surfaceMeshTriangulate",
                    "void Foam::boundaryConditionsContainer::solidRigidContact"
                );

                cp(stlFile, "constant/triSurfaces/" + stlFile);

                Info<< "friction coefficient: "
                    << frictionCoefficientsAxial[patchI] << endl;

                contactModelDict.clear();
                contactModelDict.add("normalContactModel", "standardPenalty");
                dictionary normalModelDict;
                {
                    normalModelDict.add("penaltyScale", "1");
                    normalModelDict.add("relaxationFactor", "0.05");
                    normalModelDict.add("infoFrequency", "10");
                    normalModelDict.add("adjustPenaltyScale", "yes");
                    dictionary adjustPSDict;
                    {
                        adjustPSDict.add("maxPenaltyScale", "20");
                        adjustPSDict.add("minPenaltyScale", "0.2");
                        adjustPSDict.add("targetAveragePenetration", "2e-6");
                    }
                    normalModelDict.add("adjustPenaltyScaleDict", adjustPSDict);
                }
                contactModelDict.add
                (
                    "standardPenaltyNormalModelDict", normalModelDict
                );
                contactModelDict.add("frictionContactModel", "standardPenalty");
                dictionary frictionModelDict;
                {
                    frictionModelDict.add("penaltyScale", "1");
                    frictionModelDict.add("relaxationFactor", "0.02");
                    frictionModelDict.add("infoFrequency", "10");
//                    frictionModelDict.add("frictionLaw", "anisotropicCoulomb");
                    frictionModelDict.add("frictionLaw", "coulombOrowan");
//                    frictionModelDict.add("frictionLaw", "stribeckCoulomb");
//                    frictionModelDict.add("frictionLaw", "coulomb");
                    dictionary frictionLawDict;
                    {
                        frictionLawDict.add
                        (
                            "frictionCoeff",
                            frictionCoefficientsAxial[patchI]
                        );
                        frictionLawDict.add
                        (
                            "primaryFrictionCoeff",
                            0.12
                        );
                        frictionLawDict.add
                        (
                            "secondaryFrictionCoeff",
                            0.02
                        );
                        frictionLawDict.add
                        (
                            "primaryDirection",
                            "(1 0 0)"
                        );
                        frictionLawDict.add
                        (
                            "shearStressLimit",
                            0.125e9
                        );
                    }
                    frictionModelDict.add("frictionLawDict", frictionLawDict);
                }
                contactModelDict.add
                (
                    "standardPenaltyFrictionModelDict", frictionModelDict
                );

                BCdict.add
                (
                    word
                    (
                        patchName + "_to_" + shadowPatchNames[patchI] + "_dict"
                    ),
                    contactModelDict
                );
            }
            BCdict.add("value", "uniform (0 0 0)");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}

void Foam::boundaryConditionsContainer::solidTraction(word patchName)
{
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "solidTraction");
        dictionary BCdict;
        {
            BCdict.add("type", "solidTraction");
            BCdict.add("traction", "uniform ( 0 0 0 )");
            BCdict.add("pressure", "uniform 0");
            BCdict.add("value", "uniform ( 0 0 0 )");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::fixedDisplacement(word patchName)
{
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "fixedDisplacement");
        dictionary BCdict;
        {
            BCdict.add("type", "fixedDisplacement");
            BCdict.add("value", "uniform ( 0 0 0 )");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::solidWedge
(
    word patchName
)
{
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "wedge");
        dictionary BCdict;
        {
            BCdict.add("type", "solidWedge");
            BCdict.add("patchType", "wedge");
            BCdict.add("value", "uniform ( 0 0 0 )");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::roller
(
    word rollerType, word rollerPatchName, vector rollerPosition, scalar rollGap
)
{
    Info << "dict: " << dict_ << endl;
    vector defaultVector(vector::zero);
    vector rotationAxis(vector::zero);
    if (rollerType == "topRoll")
    {
        if (rollGap == 0)
        {
            defaultVector = vector(0, -0.001, 0);
        }
        else
        {
            defaultVector = vector(0, -(rollerPosition.y() - rollGap/2), 0);
        }
        rotationAxis = vector(0, 0, 1);
    }
    else if (rollerType == "bottomRoll")
    {
        if (rollGap == 0)
        {
            defaultVector = vector(0, 0.001, 0);
        }
        else
        {
            defaultVector = vector(0, -(rollerPosition.y() + rollGap/2), 0);
        }
        rotationAxis = vector(0, 0, -1);
    }
    else if (rollerType == "rightRoll")
    {
        if (rollGap == 0)
        {
            defaultVector = vector(0, 0, -0.001);
        }
        else
        {
            defaultVector = vector(0, 0, -(rollerPosition.z() - rollGap/2));
        }
        rotationAxis = vector(0, -1, 0);
    }
    else if (rollerType == "leftRoll")
    {
        if (rollGap == 0)
        {
            defaultVector = vector(0, 0, 0.001);
        }
        else
        {
            defaultVector = vector(0, 0, -(rollerPosition.z() + rollGap/2));
        }
        rotationAxis = vector(0, 1, 0);
    }
    else
    {
        FatalErrorIn
        (
        "void Foam::boundaryConditionsContainer::roller"
        "(word rollerType, word rollerPatchName)"
        )
            << "Unknown roller type " << rollerType << ". Possible options are "
            << "topRoll, bottomRoll, rightRoll, leftRoll." << abort(FatalError);
    }
    const vector displacement =
        dict_.lookupOrAddDefault<vector>("displacement", defaultVector);

    // Lookup RPM value or default to 100 rpm
    const scalar RPM = dict_.lookupOrAddDefault<scalar>("rpm", 100);

//    word patchName = word(rollerType + rollerPatchName);

    if (!copyDict(rollerPatchName))
    {
        printDefaultBCInfo(rollerPatchName, "roller");
        dictionary BCdict;
        {
            BCdict.add("type", "roller");
            BCdict.add("displacement", displacement);
            BCdict.add("dispRampTime", "0.01");
            BCdict.add("rpm", RPM);
            BCdict.add("rotationRampTime", "0.01");
            BCdict.add("rotationAxis", rotationAxis);
            BCdict.add("value", "uniform (0 0 0)");
        }
        dict_.set(rollerPatchName, BCdict);
    }
    else printCopyBCInfo(rollerPatchName);
}

void Foam::boundaryConditionsContainer::freeRoller
(
    word rollerType, word rollerPatchName, vector rollerPosition, scalar rollGap
)
{
    vector defaultVector(vector::zero);
    vector rotationAxis(vector::zero);
    if (rollerType == "topRoll")
    {
        if (rollGap == 0)
        {
            defaultVector = vector(0, -0.001, 0);
        }
        else
        {
            defaultVector = vector(0, -(rollerPosition.y() - rollGap/2), 0);
        }
        rotationAxis = vector(0, 0, 1);
    }
    else if (rollerType == "bottomRoll")
    {
        if (rollGap == 0)
        {
            defaultVector = vector(0, 0.001, 0);
        }
        else
        {
            defaultVector = vector(0, -(rollerPosition.y() + rollGap/2), 0);
        }
        rotationAxis = vector(0, 0, -1);
    }
    else if (rollerType == "rightRoll")
    {
        if (rollGap == 0)
        {
            defaultVector = vector(0, 0, -0.001);
        }
        else
        {
            defaultVector = vector(0, 0, -(rollerPosition.z() - rollGap/2));
        }
        rotationAxis = vector(0, -1, 0);
    }
    else if (rollerType == "leftRoll")
    {
        if (rollGap == 0)
        {
            defaultVector = vector(0, 0, 0.001);
        }
        else
        {
            defaultVector = vector(0, 0, -(rollerPosition.z() + rollGap/2));
        }
        rotationAxis = vector(0, 1, 0);
    }
    else
    {
        FatalErrorIn
        (
        "void Foam::boundaryConditionsContainer::roller"
        "(word rollerType, word rollerPatchName)"
        )
            << "Unknown roller type " << rollerType << ". Possible options are "
            << "topRoll, bottomRoll, rightRoll, leftRoll." << abort(FatalError);
    }
    const vector displacement =
        dict_.lookupOrAddDefault<vector>("displacement", defaultVector);

    if (!copyDict(rollerPatchName))
    {
        printDefaultBCInfo(rollerPatchName, "freeRoller");
        dictionary BCdict;
        {
            BCdict.add("type", "freeRollerNew");
            BCdict.add("displacement", displacement);
            BCdict.add("dispRampTime", "0.01");
            BCdict.add("rotationAxis", rotationAxis);
            BCdict.add("accelerationTime", "0.01");
            // TODO: add roller data to the dataContainer and use it here
            // instead of hard coding it
            BCdict.add("zoneName", "rightRoll");
            BCdict.add("rpm", "30");
            BCdict.add("value", "uniform (0 0 0)");
        }
        dict_.set(rollerPatchName, BCdict);
    }
    else printCopyBCInfo(rollerPatchName);
}


// Thermal boundary conditions

void Foam::boundaryConditionsContainer::fixedTemperatureGradient(word patchName)
{
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "fixedTemperatureGradient");
        dictionary BCdict;
        {
            BCdict.add("type", "fixedTemperatureGradient");
            BCdict.add("gradient", "uniform 0");
            BCdict.add("value", "uniform 20");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::fixedTemperature(word patchName)
{
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "fixedTemperature");
        dictionary BCdict;
        {
            BCdict.add("type", "fixedTemperature");
            BCdict.add("value", "uniform 20");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::thermalContactSlave
(
    word patchName, word shadowPatchName
)
{
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "thermalContactSlave");
        dictionary BCdict;
        {
            BCdict.add("type", "thermalContact");
            BCdict.add("shadowPatch", shadowPatchName);
            BCdict.add("master", "no");
            BCdict.add("underRelaxation", "0.01");
            BCdict.add("alpha", "50");
            BCdict.add
            (
                "Tinf",
                dataContainer().returnKey<scalar>("ambientTemperature")
            );
            BCdict.add("Rc", "0.000005");
            BCdict.add("value", "uniform 20");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::thermalContactMaster
(
    word patchName, word shadowPatchName
)
{
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "thermalContactMaster");
        dictionary BCdict;
        {
            BCdict.add("type", "thermalContact");
            BCdict.add("master", "yes");
            BCdict.add("shadowPatch", shadowPatchName);
            BCdict.add("underRelaxation", "0.01");
            BCdict.add("thermalConductivityName", "k");
            BCdict.add("alpha", "uniform 300");
            BCdict.add
            (
                "Tinf",
                dataContainer().returnKey<scalar>("ambientTemperature")
            );
            BCdict.add("Rc", "0.000005");
            BCdict.add("value", "uniform 20");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::thermalContactMaster
(
    word patchName, wordList shadowPatchNames
)
{
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "thermalContactMaster");
        dictionary BCdict;
        {
            BCdict.add("type", "thermalContact");
            BCdict.add("master", "yes");
            BCdict.add("shadowPatches", shadowPatchNames);
            BCdict.add("underRelaxation", "0.01");
            BCdict.add("thermalConductivityName", "k");
            BCdict.add("alpha", "uniform 300");
            BCdict.add
            (
                "Tinf",
                dataContainer().returnKey<scalar>("ambientTemperature")
            );
            BCdict.add("Rc", "0.000005");
            BCdict.add("value", "uniform 20");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::thermalSymmetry(word patchName)
{
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "thermalSymmetry");
        dictionary BCdict;
        {
            BCdict.add("type", "thermalSymmetry");
            BCdict.add("patchType", "symmetryPlane");
            BCdict.add("value", "uniform 20");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::thermalConvection(word patchName)
{
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "thermalConvection");
        dictionary BCdict;
        {
            BCdict.add("type", "thermalConvection");
            BCdict.add("alpha", "uniform 77");
            BCdict.add
            (
                "Tinf",
                dataContainer().returnKey<scalar>("ambientTemperature")
            );
            BCdict.add("value", "uniform 20");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}


void Foam::boundaryConditionsContainer::thermalWedge
(
    word patchName
)
{
    if (!copyDict(patchName))
    {
        printDefaultBCInfo(patchName, "wedge");
        dictionary BCdict;
        {
            BCdict.add("type", "solidWedge");
            BCdict.add("patchType", "wedge");
            BCdict.add("value", "uniform 0");
        }
        dict_.set(patchName, BCdict);
    }
    else printCopyBCInfo(patchName);
}

// ************************************************************************* //
