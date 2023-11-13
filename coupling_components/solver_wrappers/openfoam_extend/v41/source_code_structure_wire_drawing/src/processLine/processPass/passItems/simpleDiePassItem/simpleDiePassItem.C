/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "simpleDiePassItem.H"
#include "addToRunTimeSelectionTable.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleDiePassItem, 0);
    addToRunTimeSelectionTable
    (
        passItem, simpleDiePassItem, dictionary
    );


// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

void Foam::simpleDiePassItem::writeMeshDict(Time& runTime)
{
    // Get reference to process line symmetry for convenience.
    word processLineSymmetry =
        dataContainer().returnKey<word>("processLineSymmetry");

    if
    (
        processLineSymmetry != "axiSymmetric" &&
        processLineSymmetry != "posYposZ" &&
        processLineSymmetry != "negYposZ" &&
        processLineSymmetry != "negYnegZ" &&
        processLineSymmetry != "posYnegZ" &&
        processLineSymmetry != "posY" &&
        processLineSymmetry != "posZ" &&
        processLineSymmetry != "negY" &&
        processLineSymmetry != "negZ" &&
        processLineSymmetry != "none"
    )
    {
        FatalErrorIn
        (
            "void Foam::simpleDiePassItem::writeMeshDict(Time& runTime)"
        )
            << "Unknown processLineSymmetry " << processLineSymmetry
            << ", for passItem " << name()
            << ". The options are: posYposZ, negYposZ, negYnegZ, posYnegZ, "
            << "posY, posZ, negY, negZ, none, axiSymmetry" << abort(FatalError);
    }

    // Get die mesh input dictionary from the passItem dictionary
    dictionary meshInputDict = passItem::meshInputDict();

    // Creating meshDict dictionary
    IOdictionary meshDict
    (
        IOobject
        (
            //name()/"meshDict",
            "meshDict",
            "system",
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    // Get wedgeAngle in case of an axisymmetric process line
    scalar wedgeAngle = 2.5;
    if (processLineSymmetry == "axiSymmetric")
    {
        wedgeAngle = dataContainer().returnKey<scalar>("wedgeAngle");
    }

    // Logical checks
    if
    (
        meshInputDict.found("outletDiameter") &&
        meshInputDict.found("outerDiameter") &&
        meshInputDict.found("inletDiameter"))
    {
        if
        (
            readScalar(meshInputDict.lookup("outletDiameter")) >=
            readScalar(meshInputDict.lookup("outerDiameter"))
        )
        {
            FatalErrorIn
            (
                "Foam::simpleDiePassItem::writeMeshDict(Time& runTime)"
            )    << "Specified die outlet diameter is equal or greater than "
                 << "the die outer diameter!" << abort(FatalError);
        }
        if
        (
            readScalar(meshInputDict.lookup("inletDiameter")) >=
            readScalar(meshInputDict.lookup("outerDiameter"))
        )
        {
            FatalErrorIn
            (
                "Foam::simpleDiePassItem::writeMeshDict(Time& runTime)"
            )    << "Specified die inlet diameter is equal or greater than the "
                 << "die outer diameter!" << abort(FatalError);
        }
    }

    // const scalar wireAxialSpacing =
    //     passItem::wireLength()/passItem::wireAxialResolution();

    // Add settings to meshDict Dict
    // meshDict.add
    // (
    //     "enforceGeometryConstraints",
    //     meshInputDict.lookupOrDefault<int>("enforceGeometryConstraints", 0)
    // );
    // meshDict.add
    // (
    //     "keepCellsIntersectingBoundary",
    //     meshInputDict.lookupOrDefault<int>("keepCellsIntersectingBoundary", 0)
    // );
    // meshDict.add
    // (
    //     "maxCellSize",
    //     readScalar(meshInputDict.lookup("maxCellSize"))
    // );
    // meshDict.add
    // (
    //     "minContactCellSize",
    //     meshInputDict.lookupOrDefault<scalar>("minContactCellSize", wireAxialSpacing*2.5)
    // );
    // meshDict.add
    // (
    //     "geometryTolerance",
    //     meshInputDict.lookupOrDefault<scalar>("geometryTolerance", 1e-05)
    // );

    // Create boundaryLayers subDict
    // dictionary boundaryLayers;
    // boundaryLayers.add
    // (
    //     "nLayers",
    //     meshInputDict.lookupOrDefault<int>("nLayers", 1)
    // );
    // boundaryLayers.add
    // (
    //     "optimiseLayer",
    //     meshInputDict.lookupOrDefault<int>("optimiseLayer", 1)
    // );
    // boundaryLayers.add
    // (
    //     "symmetryPlaneLayerTopology",
    //     meshInputDict.lookupOrDefault<int>("symmetryPlaneLayerTopology", 1)
    // );
    // meshDict.add("boundaryLayers", boundaryLayers);

    meshDict.set("symmetryType", processLineSymmetry);

    meshDict.set("rollSetup", "none");

    // Create patchNames dictionary to set the names as specified in the
    // dataContainer
    dictionary patchNames =
        dataContainer().dataContainerDict().subDict("patchNames");

    meshDict.set("rollingMillPatchNames", patchNames);


    // Add local refinement dict
    // TODO: read cellSize or additionalRefinementLevels from input mesh dict
//     dictionary localRefinement;
//     dictionary dieContactPatch;
// //    dieContactPatch.add("additionalRefinementLevels", 3);
//     dieContactPatch.add("cellSize", wireAxialSpacing*1.5);
//     localRefinement.add(dieContactPatchName_, dieContactPatch );
    // meshDict.add("localRefinement", localRefinement);


    // Add dieMeshDict
    {
        dictionary& dieMeshDict = meshDict.subDict("dieMeshDict");

        const word dieType = meshInputDict.lookup("type");

        if (dieType == "bekaertSpecs")
        {
            dieMeshDict.set("type", "bekaertSpecs");
            dieMeshDict.set("H", readScalar(meshInputDict.lookup("H")));
            dieMeshDict.set("W", readScalar(meshInputDict.lookup("W")));
            dieMeshDict.set("Con", readScalar(meshInputDict.lookup("Con")));
            dieMeshDict.set("TC", readScalar(meshInputDict.lookup("TC")));
            dieMeshDict.set("gamma2", readScalar(meshInputDict.lookup("gamma2")));
            dieMeshDict.set("BC", readScalar(meshInputDict.lookup("BC")));
            dieMeshDict.set("gamma1", readScalar(meshInputDict.lookup("gamma1")));
            dieMeshDict.set("C_R", readScalar(meshInputDict.lookup("CR")));
            dieMeshDict.set("C_A", readScalar(meshInputDict.lookup("CR_A")));
            dieMeshDict.set("EB_R", readScalar(meshInputDict.lookup("EB_R")));
            dieMeshDict.set("EB_L", readScalar(meshInputDict.lookup("EB_L")));
            dieMeshDict.set("EC_A", readScalar(meshInputDict.lookup("EC_A")));
            dieMeshDict.set("EC_L", readScalar(meshInputDict.lookup("EC_L")));
            dieMeshDict.set("RC_A", readScalar(meshInputDict.lookup("RC_A")));
            dieMeshDict.set("RC_L", readScalar(meshInputDict.lookup("RC_L")));
            dieMeshDict.set("RC_R", readScalar(meshInputDict.lookup("RC_R")));
            dieMeshDict.set("D", readScalar(meshInputDict.lookup("D")));
            dieMeshDict.set("BL", readScalar(meshInputDict.lookup("BL")));
            dieMeshDict.set("BR_A", readScalar(meshInputDict.lookup("BR_A")));
            dieMeshDict.set("BR_L", readScalar(meshInputDict.lookup("BR_L")));
            dieMeshDict.set("BR_R", readScalar(meshInputDict.lookup("BR_R")));
            dieMeshDict.set("XC_A", readScalar(meshInputDict.lookup("XC_A")));
            dieMeshDict.set("XC_L", readScalar(meshInputDict.lookup("XC_L")));
            dieMeshDict.set("XC_R", readScalar(meshInputDict.lookup("XC_R")));
        }
        else if (dieType == "conical")
        {
            dieMeshDict.set("type", "conical");
            dieMeshDict.set("wedgeAngle", wedgeAngle);
            dieMeshDict.set
            (
                "outerDiameter",
                readScalar(meshInputDict.lookup("outerDiameter"))
            );
            dieMeshDict.set
            (
                "axialLength",
                readScalar(meshInputDict.lookup("axialLength"))
            );
            dieMeshDict.set
            (
                "inletDiameter",
                readScalar(meshInputDict.lookup("inletDiameter"))
            );
            dieMeshDict.set
            (
                "outletDiameter",
                readScalar(meshInputDict.lookup("outletDiameter"))
            );
            // TODO: Check if geometries are specified; distinguish between
            // circular and profile die
            WarningIn(type())
                << "TODO: Check if geometries are specified; distinguish "
                << "between circular and profile dies" << endl;
            if (meshInputDict.found("outletGeometryFile"))
            {
                dieMeshDict.set
                (
                    "outletGeometryFile",
                    meshInputDict.lookup("outletGeometryFile")
                );
                dieMeshDict.set
                (
                    "dieGeometryFile",
                    meshInputDict.lookup("dieGeometryFile")
                );
            }
        }
        //meshDict.add("dieMeshDict", dieMeshDict);
    }

//    Info << "meshDict: " << meshDict << endl;
    // Write meshDict dict
    meshDict.regIOobject::write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
simpleDiePassItem::simpleDiePassItem
(
    const processPass& parent,
    const word& name,
    const dictionary& dict,
    const int& passNr,
    dataHandler* dataContainerPtr
)
:
    passItem(parent, name, dict, passNr, dataContainerPtr),
    boundaryConditionsDUGenerated_(false),
    boundaryConditionsTGenerated_(false),
    dieDownstreamPatchName_(),
    dieUpstreamPatchName_(),
    dieContactPatchName_(),
    dieContactShadowPatchName_(),
    dieToHousingPatchName_(),
    dieSymmetryYPatchName_(),
    dieSymmetryZPatchName_(),
    dieFrontPatchName_(),
    dieBackPatchName_(),
    dieEntranceConePatchName_(),
    dieExitConePatchName_(),
    frictionCoefficientAxial_()
{
    const dictionary patchNames =
        dataContainer().dataContainerDict().subDict("patchNames");

    dieDownstreamPatchName_ = word(patchNames.lookup("dieDownstreamPatchName"));

    dieUpstreamPatchName_ = word(patchNames.lookup("dieUpstreamPatchName"));

    dieContactPatchName_ = word(patchNames.lookup("dieToWirePatchName"));

    dieContactShadowPatchName_ =
        word(patchNames.lookup("wireContactPatchName"));

    dieToHousingPatchName_ = word(patchNames.lookup("dieToHousingPatchName"));

    dieSymmetryYPatchName_ = word(patchNames.lookup("dieSymmetryYPatchName"));

    dieSymmetryZPatchName_ = word(patchNames.lookup("dieSymmetryZPatchName"));

    dieFrontPatchName_ = word(patchNames.lookup("dieFrontPatchName"));

    dieBackPatchName_ = word(patchNames.lookup("dieBackPatchName"));

    dieEntranceConePatchName_ =
        word(patchNames.lookup("dieEntranceConePatchName"));

    dieExitConePatchName_ = word(patchNames.lookup("dieExitConePatchName"));

    frictionCoefficientAxial_ =
        passItem::boundaryConditionsDU().lookupOrDefault<scalar>
        (
            "frictionCoefficientAxial",
            0.01
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

simpleDiePassItem::~simpleDiePassItem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void simpleDiePassItem::setupMesh(Time& runTime)
{
    Info<< "    " << name() << ": setting up die mesh" << endl;

    // cfMesh will be used for creating the mesh
    // Change to the processPass directory and write the meshDict
    chDir(runTime.path());

    // Copy in default meshDict
    WarningIn("void simpleDiePassItem::setupMesh(Time& runTime)")
        << "Why are we using runSystemCommand here for a cp?" << endl;
    dataContainer().runSystemCommand
    (
        "cp ../../dictFiles/defaultDicts/meshDict system/",
        "log.cp_simpleDiePassItem_setupMesh_meshDict",
        "void Foam::processPass::setup()"
    );

    // Edit meshDict
    writeMeshDict(runTime);

    // Make a copy of the meshDict for debugging
    cp("system/meshDict", "system/" + name() + "/meshDict");

    //cp("system/" + name() + "/meshDict","system/meshDict");

    // Run the mesh utility
    dataContainer().runSystemCommand
    (
        "export OMP_NUM_THREADS=1; "
        "generateDieMesh",
        "log.generateMesh" + name(),
        "void rollerPassItem::setupMesh"
    );

    // Move mesh from constant/polymesh to the subMesh directory and clean
    // meshDict from system directory
    mv(runTime.path()/"constant/polyMesh","constant"/name()/"polyMesh");
    rm("system/meshDict");

    // Position the die so that the maximum x coordinate is x = 0
    boundBox meshBB(mesh(runTime).points());

    const vector translation = vector(-meshBB.max().x(), 0, 0);

    positionMesh(runTime, translation);
}


Foam::dictionary& Foam::simpleDiePassItem::boundaryConditionsDU()
{
    // Get reference to the DU input boundary conditions dictionary
    dictionary& dict = passItem::boundaryConditionsDU();

    // Check if the DU boundary conditions dicitonary was already generated.
    // If no, generate it. If yes, do not generate it and return the existing
    // one.
    if (!boundaryConditionsDUGenerated_)
    {
        Info<< "        " << name() << ": setting up DU boundary condition"
            << endl;

        // Create boundary conditions container
        boundaryConditionsContainer BCCont(dict, &dataContainer());

        // Get reference to process line symmetry for convenience.
        word processLineSymm =
            dataContainer().returnKey<word>("processLineSymmetry");

        // Set default boundary conditions
        BCCont.fixedDisplacement(dieDownstreamPatchName_);
        BCCont.solidTraction(dieUpstreamPatchName_);
        BCCont.fixedDisplacement(dieToHousingPatchName_);
        BCCont.solidTraction(dieEntranceConePatchName_);
        BCCont.solidTraction(dieExitConePatchName_);

        // A roller is considered master by default, except when the rollSetup
        // is singleRoll (nrOfPassItems = 2) or a single roller contact patch
        // is specified

        const scalar nrOfPassItems =
            dataContainer().returnKey<scalar>
            (
                "nrOfPassItems",
                passItem::passNr()
            );

        const Switch singleRollerContact =
            dataContainer().returnKey<Switch>
            (
                "singleRollerContact",
                passItem::passNr()
            );

        if (nrOfPassItems > 2 && !singleRollerContact)
        {
            BCCont.solidContactSlave
            (
                dieContactPatchName_, dieContactShadowPatchName_
            );
        }
        else
        {
            BCCont.solidContactMaster
            (
                dieContactPatchName_,
                dieContactShadowPatchName_,
                frictionCoefficientAxial_,
                passItem::boundaryConditionsDU().lookupOrDefault<scalar>
                (
                    "targetAveragePenetration", 2e-06
                ),
                passItem::boundaryConditionsDU().lookupOrDefault<scalar>
                (
                    "minPenaltyScale", 0.5
                ),
                passItem::boundaryConditionsDU().lookupOrDefault<scalar>
                (
                    "maxPenaltyScale", 20
                )
            );
        }
//        BCCont.solidContactSlave
//            (
//                dieContactPatchName_, dieContactShadowPatchName_
//            );

        // Specify symmetry boundary conditions
        if (processLineSymm == "posYposZ" || processLineSymm == "negYposZ" ||
            processLineSymm == "negYnegZ" || processLineSymm == "posYnegZ" ||
            processLineSymm == "posY" || processLineSymm == "negY" ||
            processLineSymm == "posZ" || processLineSymm == "negZ")
        {
            BCCont.solidSymmetry(dieSymmetryYPatchName_);
            BCCont.solidSymmetry(dieSymmetryZPatchName_);
        }
        else if (processLineSymm == "axiSymmetric")
        {
            BCCont.solidWedge(dieFrontPatchName_);
            BCCont.solidWedge(dieBackPatchName_);
        }

        boundaryConditionsDUGenerated_ = true;

        // At this point, all entries in dict should be sub-dicts.
        // All other entries will be removed.
        const wordList patchNames = dict.toc();

        forAll(patchNames, patchI)
        {
            const word patchName = patchNames[patchI];

            if (!dict.isDict(patchName))
            {
                dict.remove(patchName);
            }
        }
    }

    return dict;
}


Foam::dictionary& Foam::simpleDiePassItem::boundaryConditionsT()
{
    // Get reference to the DU input boundary conditions dictionary
    dictionary& dict = passItem::boundaryConditionsT();

    if (!boundaryConditionsTGenerated_)
    {
        Info<< "        " << name() << ": setting up T boundary condition"
            << endl;

        // Create boundary conditions container
        boundaryConditionsContainer BCCont(dict, &dataContainer());

        // Get reference to process line symmetry for convenience.
        word processLineSymm =
            dataContainer().returnKey<word>("processLineSymmetry");

        // Set default boundary conditions
        BCCont.thermalConvection(dieDownstreamPatchName_);
        BCCont.thermalConvection(dieUpstreamPatchName_);
        BCCont.thermalConvection(dieToHousingPatchName_);
        BCCont.thermalConvection(dieEntranceConePatchName_);
        BCCont.thermalConvection(dieExitConePatchName_);

        // A roller is considered master by default, except when the rollSetup
        // is singleRoll (nrOfPassItems = 2) or a single roller contact patch
        // is specified

        const scalar nrOfPassItems =
            dataContainer().returnKey<scalar>
            (
                "nrOfPassItems",
                passItem::passNr()
            );

        const Switch singleRollerContact =
            dataContainer().returnKey<Switch>
            (
                "singleRollerContact",
                passItem::passNr()
            );

        if (nrOfPassItems > 2 && !singleRollerContact)
        {
            BCCont.thermalContactSlave
            (
                dieContactPatchName_, dieContactShadowPatchName_
            );
        }
        else
        {
            BCCont.thermalContactMaster
            (
                dieContactPatchName_,
                dieContactShadowPatchName_
            );
        }

//        BCCont.thermalContactSlave
//            (
//                dieContactPatchName_, dieContactShadowPatchName_
//            );

        if (processLineSymm == "posYposZ" || processLineSymm == "negYposZ" ||
            processLineSymm == "negYnegZ" || processLineSymm == "posYnegZ" ||
            processLineSymm == "posY" || processLineSymm == "negY" ||
            processLineSymm == "posZ" || processLineSymm == "negZ")
        {
            BCCont.thermalSymmetry(dieSymmetryYPatchName_);
            BCCont.thermalSymmetry(dieSymmetryZPatchName_);
        }
        else if (processLineSymm == "axiSymmetric")
        {
            BCCont.thermalWedge(dieSymmetryYPatchName_);
            BCCont.thermalWedge(dieSymmetryZPatchName_);
        }

        boundaryConditionsDUGenerated_ = true;

        // At this point, all entries in dict should be sub-dicts.
        // All other entries will be removed.
        const wordList patchNames = dict.toc();

        forAll(patchNames, patchI)
        {
            const word patchName = patchNames[patchI];

            if (!dict.isDict(patchName))
            {
                dict.remove(patchName);
            }
        }
    }

    return dict;
}


Foam::word Foam::simpleDiePassItem::contactPatch()
{
    return dieContactPatchName_;
}


Foam::scalar Foam::simpleDiePassItem::frictionCoefficientAxial()
{
    return frictionCoefficientAxial_;
}


List<dictionary> simpleDiePassItem::functionObjects()
{
    // By default we will return no function objects; this can be
    // overwritten by the specific pass items

    List<dictionary> dicts(2);

    {
        dictionary subDict;
        subDict.add("type", "solidForces");
        subDict.add("historyPatch", dieToHousingPatchName_);
        subDict.add("stressName", "sigmaCauchy");

        dicts[0] = subDict;
    }

    {
        dictionary subDict;
        subDict.add("type", "solidForces");
        subDict.add("historyPatch", dieContactPatchName_);
        subDict.add("stressName", "sigmaCauchy");

        dicts[1] = subDict;
    }

    return dicts;
}


// ************************************************************************* //

} // end of namespace foam
