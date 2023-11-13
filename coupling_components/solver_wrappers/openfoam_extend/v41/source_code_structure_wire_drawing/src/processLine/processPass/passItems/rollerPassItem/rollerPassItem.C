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

#include "rollerPassItem.H"
#include "addToRunTimeSelectionTable.H"
#include "standAlonePatch.H"
#include "patchToPatchInterpolationNew.H"
#include "processPass.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rollerPassItem, 0);
    addToRunTimeSelectionTable
    (
        passItem, rollerPassItem, dictionary
    );


// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //


void Foam::rollerPassItem::writeMeshDict(Time& runTime)
{
    // Get reference to process line symmetry for convenience.
    word processLineSymmetry =
        dataContainer().returnKey<word>("processLineSymmetry");

    // Read location of input data
    fileName meshInputDirectory =
        //dataContainer().returnKey<fileName>("meshInputDir");
        dataContainer().returnKey<fileName>
        (
            "processLineRootDir"
        )/"dictFiles"/parent().name();

    if
    (
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
        FatalErrorIn("void Foam::rollerPassItem::writeMeshDict(Time& runTime)")
            << "Unknown processLineSymmetry " << processLineSymmetry
            << ", for passItem " << name()
            << ". The options are: posYposZ, negYposZ, negYnegZ, posYnegZ, "
            << "posY, posZ, negY, negZ, none" << abort(FatalError);
    }

    // TODO: check if e.g. processLineSymmetry = "posY" and a bottomRoll is
    // asked, to give a FatalError since cfMesh won't create the bottomRoller
    // in this condition

    // Get wire mesh input dictionary from the passItem dictionary
    dictionary meshInputDict = passItem::meshInputDict();

    // Read and do basic check on roller diameters
    scalar innerDiameter = readScalar(meshInputDict.lookup("innerDiameter"));
    scalar outerDiameter = readScalar(meshInputDict.lookup("outerDiameter"));
    if (innerDiameter > 0.9*outerDiameter)
    {
        FatalErrorIn("void Foam::rollerPassItem::writeMeshDict(Time& runTime)")
            << "Inner roll diameter (" << innerDiameter << ") should be less "
            << " than 0.9 times the outer roll diameter ("
            << outerDiameter << "), " << " for passItem " << name()
            << ". Either specify a geometryFile or a dxfFile."
            << abort(FatalError);
    }

    // Read geometry tolerance
    scalar geometryTolerance =
        meshInputDict.lookupOrDefault<scalar>("geometryTolerance", 5e-5);

    // We force the roll setup to be singleRoll as all rollers then are
    // a separate region which later can be easily positioned
    word rollSetup = "singleRoll";

    // Creating meshDict dictionary
    IOdictionary meshDict
    (
        IOobject
        (
            name()/"meshDict",
            "system",
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dictionary("meshDict")
    );

    // Add settings to meshDict Dict
    meshDict.add
    (
        "enforceGeometryConstraints",
        meshInputDict.lookupOrDefault<int>("enforceGeometryConstraints", 0)
    );
    meshDict.add
    (
        "keepCellsIntersectingBoundary",
        meshInputDict.lookupOrDefault<int>("keepCellsIntersectingBoundary", 0)
    );
    meshDict.add
    (
        "maxCellSize",
        meshInputDict.lookupOrDefault<scalar>("maxCellSize", 0.001)
    );
    meshDict.add
    (
        "minContactCellSize",
        meshInputDict.lookupOrDefault<scalar>("minContactCellSize", 10e-6)
    );
    meshDict.add
    (
        "minNumFacesBetweenFeatures",
        meshInputDict.lookupOrDefault<scalar>("minNumFacesBetweenFeatures", 3)
    );
    meshDict.add
    (
        "geometryTolerance", geometryTolerance
    );


    // Create boundaryLayers subDict
    {
        dictionary boundaryLayers;
        boundaryLayers.add
        (
            "nLayers",
            meshInputDict.lookupOrDefault<int>("nLayers", 1)
        );
        boundaryLayers.add
        (
            "optimiseLayer",
            meshInputDict.lookupOrDefault<int>("optimiseLayer", 1)
        );
        boundaryLayers.add
        (
            "symmetryPlaneLayerTopology",
            meshInputDict.lookupOrDefault<int>("symmetryPlaneLayerTopology", 1)
        );
        meshDict.add("boundaryLayers", boundaryLayers);
    }


//    dictionary anisotropicSources;
//    dictionary plane;
//    plane.add("type", "plane");
//    plane.add("origin", "(0 0 0)");
//    plane.add("normal", "(1 0 0)");
//    plane.add("scalingDistance", 0.1);
//    plane.add("scalingFactor", 2);
//    anisotropicSources.add("plane", plane);
//    meshDict.add("anisotropicSources", anisotropicSources);

    // anisotropic sources to get different mesh density in axial versus radial
    // direction
//anisotropicSources
//{
//     plane
//     {
//        type plane;
//       origin (0 0 0);
//       normal (1 0 0);
//       scalingDistance 0.05; // distance up to which you wish to scale
                                //the cells
//       scalingFactor 2; // larger than one makes the cells larger
//    }
//}
// This is a standard setting in cfMesh, and you have to be careful how you
// position the normal in your case. For top and bottom rolls the axial
// direction is x, and for left and right roll the axial direction is y during
// meshing. 2D mesher runs in the x-y coordinate system and it transforms the
// geometry from the y-z system.


//    // Create, set and add automatic refinement settings
//    dictionary automaticRefinement;
//    // Allow max 10 refinement levels on top of maxCellSize
//    automaticRefinement.add("maxAdditionalRefinementLevels", "10");
//    // Increasing next number to get more cells over the gap you wish to mesh
//    automaticRefinement.add("nCellsOverFeature", "1");
//    // This flag causes refinement of cells that contain parts of the geometry
//    // that are not connected via a feature edge
//    automaticRefinement.add("distinctPartsRefinement", "1");
//    meshDict.add("automaticRefinement", automaticRefinement);


    // Create wireMeshDict subDict
    {
        dictionary wireMeshDict;
        // wireDiameter for roller meshing is used to position the rollers,
        // i.e., to indicate the initial roll gap. It is defaulted to 50mm as we
        // assume no wire with a diameter larger than 50mm will be processed.
        wireMeshDict.add
        (
             "wireDiameter",
             meshInputDict.lookupOrDefault<scalar>("wireDiameter", 0.05)
        );
        wireMeshDict.add
        (
            "wireLength",
            readScalar(meshInputDict.lookup("wireLength"))
        );
        meshDict.add("wireMeshDict", wireMeshDict);

        meshDict.add("symmetryType", processLineSymmetry);
        meshDict.add("rollSetup", rollSetup);
    }

    // Create patchNames dictionary to set the names as specified in the
    // dataContainer
    dictionary patchNames =
        dataContainer().dataContainerDict().subDict("patchNames");

    meshDict.add("rollingMillPatchNames", patchNames);

    // Create, set and add rollSetupDict subDict
    {
        dictionary rollSetupDict;
        dictionary rollerDict;

        // Get roll profile
        fileName geometryFileName;
        fileName dxfFileName;
        bool geometryFileFound = false;
        bool dxfFileFound = false;
        if (meshInputDict.found("geometryFile"))
        {
            geometryFileName = fileName(meshInputDict.lookup("geometryFile"));
            geometryFileFound = true;
        }
        if (meshInputDict.found("dxfFile"))
        {
            dxfFileName = fileName(meshInputDict.lookup("dxfFile"));
            dxfFileFound = true;
        }
        fileName geometryFile = meshInputDirectory + "/" + geometryFileName;
        fileName dxfFile = meshInputDirectory + "/" + dxfFileName;

        if (!geometryFileFound & !dxfFileFound)
        {
            Info<< "    Flat roller" << endl;
            rollerDict.add("type", "flatRoller");
        }
        else if (dxfFileFound)
        {
            rollerDict.add("dxfFile", dxfFile);
            rollerDict.add("type", "dxfFileRoller");
        }
        else
        {
            rollerDict.add("geometryFile", geometryFile);
            rollerDict.add("type", "geometryFileRoller");
        }

        rollerDict.add
        (
            "rollWidth",
            meshInputDict.lookupOrDefault<scalar>
            (
                "rollWidth", 0.04
            )
        );

        // Number of cells in circumferential direction, outside of refined
        // region which is in contact with the wire. Note that the tolerance is
        // kind of master meshing parameter.
        rollerDict.add
        (
            "circumferentialResolution",
            meshInputDict.lookupOrDefault<int>("circumferentialResolution", 20)
        );

        // Flip the roll about the z-axis
        rollerDict.add
        (
            "flipRoll",
            meshInputDict.lookupOrDefault<bool>("flipRoll", 0)
        );

        // Grading of cells in circumferential direction, starting from te
        // refined region which is in contact with the wire
        rollerDict.add
        (
            "circumGrading",
            meshInputDict.lookupOrDefault<scalar>("circumGrading", 1.2)
        );
        // Grading in radial roll direction. Cell at inner roll diameter will be
        // "radialGrading" times larger than cell near contact surface.
        rollerDict.add
        (
            "radialGrading",
            meshInputDict.lookupOrDefault<scalar>("radialGrading", 1)
        );
        // Grading in axial roll direction. A cell near the front or back face
        // will be "axialGrading" times larger than cell in the mid plane of
        // the roll
        rollerDict.add
        (
            "axialGrading",
            meshInputDict.lookupOrDefault<scalar>("axialGrading", 1)
        );

        // Determine axial spacing of the wire mesh. The roller mesh in
        // cricumferential direction is adapted to it.
        scalar wireAxialSpacing =
            passItem::wireLength()/passItem::wireAxialResolution();

        scalar circumScaling = 1.1*wireAxialSpacing/
            (
               0.5*outerDiameter
              *acos(1 - geometryTolerance/(0.5*outerDiameter))
            );

        // Scaling of cells in circumferential direction
        rollerDict.add
        (
            "circumScaling",
            meshInputDict.lookupOrDefault<scalar>
            (
                "circumScaling", circumScaling
            )
        );
        // Scaling of cells in radial direction
        rollerDict.add
        (
            "radialScaling",
            meshInputDict.lookupOrDefault<scalar>("radialScaling", 1.0)
        );
        // Scaling of cells in axial direction
        rollerDict.add
        (
            "axialScaling",
            meshInputDict.lookupOrDefault<scalar>("axialScaling", 1.0)
        );
        rollerDict.add
        (
            "innerDiameter",
            meshInputDict.lookupOrDefault<scalar>("innerDiameter", 0.2)
        );
        rollerDict.add
        (
            "outerDiameter", outerDiameter
        );
        rollerDict.add
        (
            "contactAreaLength",
            readScalar(meshInputDict.lookup("contactAreaLength"))
        );
        rollSetupDict.add(rollerType_, rollerDict);
        word rollSetupDictName=rollSetup+"Dict";
        meshDict.add(rollSetupDictName, rollSetupDict);
    }

    // Write meshDict dict
    meshDict.regIOobject::write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
rollerPassItem::rollerPassItem
(
    const processPass& parent,
    const word& name,
    const dictionary& dict,
    const int& passNr,
    dataHandler* dataContainerPtr
)
:
    passItem(parent, name, dict, passNr, dataContainerPtr),
    rollerType_(),
    boundaryConditionsDUGenerated_(false),
    boundaryConditionsTGenerated_(false),
    rollerAxisPatchName_(),
    rollerContactPatchName_(),
    rollerContactShadowPatchName_(),
    rollerBackPatchName_(),
    rollerFrontPatchName_(),
    rollerSymmetryYPatchName_(),
    rollerSymmetryZPatchName_(),
    rollerToRollerFrontPatchName_(),
    rollerToRollerBackPatchName_(),
    rollerToAirPatchName_(),
    frictionCoefficientAxial_(),
    rollerPosition_(),
    motionEnabler_
    (
        dict.subDict("boundaryConditions").lookupOrDefault<word>
        (
            "motionEnabler", "drivenRoller"
        )
    )
{
    const word rollSetup = dataContainer().returnKey<word>("rollSetup", passNr);

    if (rollSetup == "none")
    {
        FatalErrorIn("rollerPassItem::rollerPassItem(const word& name...)")
            << "Incorrect rollSetup specified for the process pass with "
            << "passItem " << passItem::name()
            << abort(FatalError);
    }

    if (dict.found("rollerType"))
    {
        rollerType_ = word(dict.lookup("rollerType"));
        if
        (
            rollerType_ != "topRoll" &&
            rollerType_ != "bottomRoll" &&
            rollerType_ != "leftRoll" &&
            rollerType_ != "rightRoll"
        )
        {
            FatalErrorIn("rollerPassItem::rollerPassItem(const word& name...)")
                << "Unknown rollerType " << rollerType_ << " specified for "
                << "passItem " << passItem::name()
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn("rollerPassItem::rollerPassItem(const word& name...)")
            << "rollerType not specified for passItem " << passItem::name()
            << abort(FatalError);
    }

    const dictionary patchNames =
        dataContainer().dataContainerDict().subDict("patchNames");

    rollerAxisPatchName_ = rollerType_ +
        word(patchNames.lookup("rollerAxisPatchName"));

    const Switch singleRollerContact =
        dataContainer().returnKey<Switch>
        (
            "singleRollerContact",
            passItem::passNr()
        );

    if (singleRollerContact)
    {
        rollerContactPatchName_ =
            word(patchNames.lookup("rollerContactPatchName"));
    }
    else
    {
        rollerContactPatchName_ = rollerType_ +
            word(patchNames.lookup("rollerContactPatchName"));
    }

    rollerContactShadowPatchName_ =
        word(patchNames.lookup("wireContactPatchName"));

    rollerBackPatchName_ = rollerType_ +
       word(patchNames.lookup("rollerBackPatchName"));

    rollerFrontPatchName_ = rollerType_ +
        word(patchNames.lookup("rollerFrontPatchName"));

    rollerSymmetryYPatchName_ = rollerType_ +
        word(patchNames.lookup("rollerSymmetryYPatchName"));

    rollerSymmetryZPatchName_ = rollerType_ +
        word(patchNames.lookup("rollerSymmetryZPatchName"));

    rollerToRollerFrontPatchName_ = rollerType_ +
        word(patchNames.lookup("rollerToRollerFrontPatchName"));

    rollerToRollerBackPatchName_ = rollerType_ +
        word(patchNames.lookup("rollerToRollerBackPatchName"));

    rollerToAirPatchName_ = rollerType_ +
        word(patchNames.lookup("rollerToAirPatchName"));

    frictionCoefficientAxial_ =
        passItem::boundaryConditionsDU().lookupOrDefault<scalar>
        (
            "frictionCoefficientAxial",
            0.15
        );
    rollerPosition_ = vector::zero;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rollerPassItem::~rollerPassItem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::rollerPassItem::positionMesh
(
    const Time& runTime,
    const fvMesh& wireMesh
)
{
    // Note: we do not need to refCast the passItem because the mesh function is
    // defined in the passItem, but this will give the same behaviour as if we
    // did refCast to wirePassItem
    //const fvMesh& wireMesh = parent().passItemObject("wire").mesh(runTime);
    //const fvMesh& wireMesh =
    //    refCast<wirePassItem>(parent().passItemObject("wire")).mesh(runTime);

    if (meshInputDict().found("positioningVector"))
    {
        const vector translation = meshInputDict().lookup("positioningVector");

        Info<< "            " << name() << ": applying translation of "
            << translation << ". Roll will not be aligned with wire!"
            << endl;

        passItem::positionMesh(runTime, translation);
    }
    else
    {
        // Lookup wire contact patch index
        const label wirePatchID =
            wireMesh.boundaryMesh().findPatchID(rollerContactShadowPatchName_);

        if (wirePatchID == -1)
        {
            FatalErrorIn
            (
                "void rollerPassItem::positionMesh("
                "const Time& runTime, const fvMesh& wireMesh)"
            )   << "wire contact patch not found!" << abort(FatalError);
        }

        // Create a standAlonePatch for the wire surface
        standAlonePatch wirePatch
        (
            wireMesh.boundaryMesh()[wirePatchID].localFaces(),
            wireMesh.boundaryMesh()[wirePatchID].localPoints()
        );

        // Define the positioning tolerance for moving the roller
        // We will define it as a fraction of the smallest face on the
        // wire patch
        const scalar positioningTol =
            1e-6*Foam::sqrt
            (
                min(mag(wireMesh.boundaryMesh()[wirePatchID].faceAreas()))
            );

        // We will iteratively move the roller until it is "just touching"
        // the wire surface.
        // We will do this by calculting the distances from the roller surface
        // to the wire surface and then moving the roller away from the wire
        // (positive y direction for the topRoll, negative y for the bottom
        // roll, etc.)

        scalar rollerMinY = GREAT;
        scalar rollerMaxY = -GREAT;
        scalar rollerMinZ = GREAT;
        scalar rollerMaxZ = -GREAT;

        vector translation = vector::zero;
        const int maxPositioningCorrectors = 10;
        int i = 0;

        do
        {
            Info<< "            Iteration " << (i + 1) << endl;

            // Reset translation to zero
            translation = vector::zero;

            // Read roller mesh from dict
            // This will re-read the mesh from disk if it has changed
            const fvMesh& rollerMesh = mesh(runTime);

            // Lookup roller contact patch index
            const label rollerPatchID =
                rollerMesh.boundaryMesh().findPatchID
                (
                    word(rollerContactPatchName_)
                );

            if (rollerPatchID == -1)
            {
                Info<< "Looking for patch: " << word(rollerContactPatchName_)
                    << endl;

                FatalErrorIn
                (
                    "void rollerPassItem::positionMesh("
                    "const Time& runTime, const fvMesh& wireMesh)"
                )   << "roller contact patch not found!" << abort(FatalError);
            }

            // Create a standAlonePatch for the roller surface
            standAlonePatch rollerPatch
            (
                rollerMesh.boundaryMesh()[rollerPatchID].localFaces(),
                rollerMesh.boundaryMesh()[rollerPatchID].localPoints()
            );

            // Create patchToPatch object to calculate point distances between
            // the wire and roller patches
            PatchToPatchInterpolationNew<standAlonePatch, standAlonePatch>
                patchToPatch
                (
                    rollerPatch,
                    wirePatch,
                    intersection::VISIBLE
                );

            // Store all rollerPatch boundary points for calculating the roller
            // position. This is needed to determine the displacement needed to
            // reach the required roll gap.
            const vectorField boundaryPoints = rollerPatch.localPoints();

            // Calculation of point distances from the roller patch to the wire
            // patch and convert to vector distances by multiply by the point
            // normals
            //PDJ: patchToPatch.pointDistanceToIntersection() gives an issue
            // in case of symmetry!!
            const vectorField wirePointDist =
                patchToPatch.pointDistanceToIntersection()
               *wirePatch.pointNormals();

            scalar translationY = 0;
            scalar translationZ = 0;

            // Find the most negative distance (when there is no crossover this
            // is the least positive distance)
            if (rollerType_ == "topRoll")
            {
                scalar minTranslationY = GREAT;

                // Find minimum distance between wire and roller
                forAll(wirePointDist, pI)
                {
                    if (mag(wirePointDist[pI]) < 1e14)
                    {
                        minTranslationY =
                            min(minTranslationY, wirePointDist[pI].y());
                    }
                }

                translationY = -minTranslationY;

                // Get the absolute position of the roller
                forAll(boundaryPoints, pI)
                {
                    if (boundaryPoints[pI].y() < rollerMinY)
                    {
                        rollerMinY = boundaryPoints[pI].y();
                    }
                }
                rollerPosition_ = vector(0, rollerMinY, 0);

            }
            else if (rollerType_ == "bottomRoll")
            {
                scalar maxTranslationY = -GREAT;

                forAll(wirePointDist, pI)
                {
                    if (mag(wirePointDist[pI]) < 1e14)
                    {
                        maxTranslationY =
                            max(maxTranslationY, wirePointDist[pI].y());
                    }
                }

                translationY = -maxTranslationY;

                forAll(boundaryPoints, pI)
                {
                    if (boundaryPoints[pI].y() > rollerMaxY)
                    {
                        rollerMaxY = boundaryPoints[pI].y();
                    }
                }
                rollerPosition_ = vector(0, rollerMaxY, 0);
            }
            else if (rollerType_ == "rightRoll")
            {
                scalar minTranslationZ = GREAT;

                forAll(wirePointDist, pI)
                {
                    if (mag(wirePointDist[pI]) < 1e14)
                    {
                        minTranslationZ =
                            min(minTranslationZ, wirePointDist[pI].z());
                    }
                }

                translationZ = -minTranslationZ;

                forAll(boundaryPoints, pI)
                {
                    if (boundaryPoints[pI].z() < rollerMinZ)
                    {
                        rollerMinZ = boundaryPoints[pI].z();
                    }
                }
                rollerPosition_ = vector(0, 0, rollerMinZ);
            }
            else if (rollerType_ == "leftRoll")
            {
                scalar maxTranslationZ = -GREAT;

                forAll(wirePointDist, pI)
                {
                    if (mag(wirePointDist[pI]) < 1e14)
                    {
                        maxTranslationZ =
                            max(maxTranslationZ, wirePointDist[pI].z());
                    }
                }

                translationZ = -maxTranslationZ;

                forAll(boundaryPoints, pI)
                {
                    if (boundaryPoints[pI].z() > rollerMaxZ)
                    {
                        rollerMaxZ = boundaryPoints[pI].z();
                    }
                }
                rollerPosition_ = vector(0, 0, rollerMaxZ);
            }
            else
            {
                FatalErrorIn
                (
                    "void rollerPassItem::positionMesh(const Time& runTime)"
                )   << "Positioning not implemented for roller type"
                    << rollerType_
                    << abort(FatalError);
            }

            // Calculate translation vector
            translation = vector(0, translationY, translationZ);

            Info<< "            Translating by: " << translation << endl;

            // Move the roller mesh
            passItem::positionMesh(runTime, translation);
        }
        while
        (
            i++ < maxPositioningCorrectors &&
            mag(translation) > positioningTol
        );

        Info<< "            Roller positioned to a tolerance of "
            << translation << endl;
    }
    // Check if positionVector is specified and add it on to the translation
//    translation =
//        meshInputDict().lookupOrDefault<vector>
//        (
//            "positioningVector",
//            vector::zero
//        );

//    rollerPosition_ += translation;

//    Info<< "    " << name() << ": superimposing translation from "
//        "positionVector of " << translation << nl << endl;

//    passItem::positionMesh(runTime, translation);
}


void rollerPassItem::setupMesh(Time& runTime)
{
    Info<< "    " << name() << ": setting up mesh" << endl;

    const Switch singleRollerContact =
        dataContainer().returnKey<Switch>
        (
            "singleRollerContact",
            passItem::passNr()
        );

    if (singleRollerContact)
    {
        Info<< "        Roller " << name() << " mesh already exists." << endl;
    }
    else
    {
        // Change to the processPass directory and write the meshDict
        chDir(runTime.path());

        writeMeshDict(runTime);

        // Copy meshDict to system directory for mesh generation utility
        cp("system/" + name() + "/meshDict","system/meshDict");

        // Run the mesh utility
        dataContainer().runSystemCommand
        (
            "export OMP_NUM_THREADS=1; "
            "generate3DRollerMesh",
            "log.generateMesh" + name(),
            "void rollerPassItem::setupMesh"
        );

        // Move mesh from constant/polymesh to the subMesh directory and clean
        // meshDict from system directory
        mv(runTime.path()/"constant/polyMesh","constant"/name()/"polyMesh");
        rm("system/meshDict");
    }
}


Foam::dictionary& Foam::rollerPassItem::boundaryConditionsDU()
{
    // Get reference to the DU input boundary conditions dictionary
    dictionary& dict = passItem::boundaryConditionsDU();

    // Check if the DU boundary conditions dicitonary was already generated.
    // If not, generate it. If yes, do not generate it and return the existing
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

        // Lookup settings which can be changed
        word motionEnabler =
            passItem::boundaryConditions().lookupOrDefault<word>
            (
                "motionEnabler", "drivenRoller"
            );

        if (motionEnabler != "drivenRoller" && motionEnabler != "freeRoller")
        {
            WarningIn
                (
                    "Foam::dictionary& Foam::rollerPassItem::"
                    "boundaryConditionsDU()"
                )
                << "Unknown method of enabling roller motion: " << motionEnabler
                << ". Setting to default drivenRoller method, i.e., roller "
                << " will be driven to the given RPM value" << endl;
            motionEnabler = "drivenRoller";
        }

        scalar rollGap = 0;
        if (dict.found("rollGap"))
        {
             rollGap =
                 readScalar(dict.lookup("rollGap"));
        }

        if (motionEnabler == "drivenRoller")
        {
            BCCont.roller
            (
                rollerType_,
                rollerAxisPatchName_,
                rollerPosition_,
                rollGap
            );
        }
        else if (motionEnabler == "freeRoller")
        {
            BCCont.freeRoller
            (
                rollerType_,
                rollerAxisPatchName_,
                rollerPosition_,
                rollGap
            );
        }
        else
        {
            FatalErrorIn
            (
            "Foam::dictionary& Foam::rollerPassItem::boundaryConditionsDU()"
            )
                << "Unknown motion enabler for roller type " << rollerType_
                << ", for passItem "
                << name() << "!" << abort(FatalError);
        }

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
                rollerContactPatchName_, rollerContactShadowPatchName_
            );
        }
        else
        {
            BCCont.solidContactMaster
            (
                rollerContactPatchName_,
                rollerContactShadowPatchName_,
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

        // Specify symmetry boundary conditions
        if (rollerType_ == "topRoll")
        {
            if(processLineSymm == "posY" || processLineSymm == "none")
            {
                BCCont.solidTraction(rollerBackPatchName_);
                BCCont.solidTraction(rollerFrontPatchName_);
            }
            else if
            (
                processLineSymm == "posZ" || processLineSymm == "posYposZ"
            )
            {
                BCCont.solidTraction(rollerFrontPatchName_);
                BCCont.solidSymmetry(rollerSymmetryZPatchName_);
            }
            else if
            (
                processLineSymm == "negZ" || processLineSymm == "posYnegZ"
            )
            {
                BCCont.solidSymmetry(rollerSymmetryZPatchName_);
                BCCont.solidTraction(rollerBackPatchName_);
            }
            else
            {
                FatalErrorIn
                (
                "Foam::dictionary& Foam::rollerPassItem::boundaryConditionsDU()"
                )
                    << "Symmetry " << processLineSymm << " is not possible"
                    << " with roll setup " << rollerType_ << " for passItem "
                    << name() << "!" << abort(FatalError);
            }
        }
        else if (rollerType_ == "bottomRoll")
        {
            if(processLineSymm == "negY" || processLineSymm == "none")
            {
                BCCont.solidTraction(rollerBackPatchName_);
                BCCont.solidTraction(rollerFrontPatchName_);
            }
            else if
            (
                processLineSymm == "posZ" || processLineSymm == "negYposZ"
            )
            {
                BCCont.solidTraction(rollerFrontPatchName_);
                BCCont.solidSymmetry(rollerSymmetryZPatchName_);
            }
            else if
            (
                processLineSymm == "negZ" || processLineSymm == "negYnegZ"
            )
            {
                BCCont.solidSymmetry(rollerSymmetryZPatchName_);
                BCCont.solidTraction(rollerBackPatchName_);
            }
            else
            {
                FatalErrorIn
                (
                "Foam::dictionary& Foam::rollerPassItem::boundaryConditionsDU()"
                )
                    << "Symmetry " << processLineSymm << " is not possible"
                    << "with roll setup " << rollerType_ << " for passItem "
                    << name() << "!" << abort(FatalError);
            }
         }
        else if (rollerType_ == "rightRoll")
        {
            if(processLineSymm == "posZ" || processLineSymm == "none")
            {
                BCCont.solidTraction(rollerBackPatchName_);
                BCCont.solidTraction(rollerFrontPatchName_);
            }
            else if
            (
                processLineSymm == "posY" || processLineSymm == "posYposZ"
            )
            {
                BCCont.solidTraction(rollerFrontPatchName_);
                BCCont.solidSymmetry(rollerSymmetryYPatchName_);
            }
            else if
            (
                processLineSymm == "negY" || processLineSymm == "negYposZ"
            )
            {
                BCCont.solidSymmetry(rollerSymmetryYPatchName_);
                BCCont.solidTraction(rollerBackPatchName_);
            }
            else
            {
                FatalErrorIn
                (
                "Foam::dictionary& Foam::rollerPassItem::boundaryConditionsDU()"
                )
                    << "Symmetry " << processLineSymm << " is not possible"
                    << "with roll setup " << rollerType_ << " for passItem "
                    << name() << "!" << abort(FatalError);
            }
        }
        else if (rollerType_ == "leftRoll")
        {
            if(processLineSymm == "negZ" || processLineSymm == "none")
            {
                BCCont.solidTraction(rollerBackPatchName_);
                BCCont.solidTraction(rollerFrontPatchName_);
            }
            else if
            (
                processLineSymm == "posY" || processLineSymm == "posYnegZ"
            )
            {
                BCCont.solidTraction( rollerFrontPatchName_);
                BCCont.solidSymmetry( rollerSymmetryYPatchName_);
            }
            else if
            (
                processLineSymm == "negY" || processLineSymm == "negYnegZ"
            )
            {
                BCCont.solidSymmetry( rollerSymmetryYPatchName_);
                BCCont.solidTraction( rollerBackPatchName_);
            }
            else
            {
                FatalErrorIn
                (
                "Foam::dictionary& Foam::rollerPassItem::boundaryConditionsDU()"
                )
                    << "Symmetry " << processLineSymm << " is not possible"
                    << "with roll setup " << rollerType_ << " for passItem "
                    << name() << "!" << abort(FatalError);
            }
        }

        BCCont.solidTraction(rollerToAirPatchName_);

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


Foam::dictionary& Foam::rollerPassItem::boundaryConditionsT()
{
    // Get reference to the DU input boundary conditions dictionary
    dictionary& dict = passItem::boundaryConditionsT();

    // Check if the T boundary conditions dicitonary was already generated.
    // If not, generate it. If yes, do not generate it and return the existing
    // one.
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
        BCCont.fixedTemperature( rollerAxisPatchName_);

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
                rollerContactPatchName_, rollerContactShadowPatchName_
            );
        }
        else
        {
            BCCont.thermalContactMaster
            (
                rollerContactPatchName_,
                rollerContactShadowPatchName_
            );
        }

//        BCCont.thermalContactSlave
//        (
//             rollerContactPatchName_, rollerContactShadowPatchName_
//        );

        // Specify symmetry boundary conditions
        if (rollerType_ == "topRoll")
        {
            if (processLineSymm == "posY" || processLineSymm == "none")
            {
                BCCont.fixedTemperatureGradient( rollerBackPatchName_);
                BCCont.fixedTemperatureGradient( rollerFrontPatchName_);
            }
            else if
            (
                processLineSymm == "posZ" || processLineSymm == "posYposZ"
            )
            {
                BCCont.fixedTemperatureGradient( rollerFrontPatchName_);
                BCCont.thermalSymmetry( rollerSymmetryZPatchName_);
            }
            else if
            (
                processLineSymm == "negZ" || processLineSymm == "posYnegZ"
            )
            {
                BCCont.thermalSymmetry( rollerSymmetryZPatchName_);
                BCCont.fixedTemperatureGradient( rollerBackPatchName_);
            }
            else
            {
                FatalErrorIn
                (
                "Foam::dictionary& Foam::rollerPassItem::boundaryConditionsDU()"
                )
                    << "Symmetry " << processLineSymm << " is not possible"
                    << "with roll setup " << rollerType_ << " for passItem "
                    << name() << "!" << abort(FatalError);
            }
        }
        else if (rollerType_ == "bottomRoll")
        {
            if(processLineSymm == "negY" || processLineSymm == "none")
            {
                BCCont.fixedTemperatureGradient( rollerBackPatchName_);
                BCCont.fixedTemperatureGradient( rollerFrontPatchName_);
            }
            else if
            (
                processLineSymm == "posZ" || processLineSymm == "negYposZ"
            )
            {
                BCCont.fixedTemperatureGradient( rollerFrontPatchName_);
                BCCont.thermalSymmetry( rollerSymmetryZPatchName_);
            }
            else if
            (
                processLineSymm == "negZ" || processLineSymm == "negYnegZ"
            )
            {
                BCCont.thermalSymmetry( rollerSymmetryZPatchName_);
                BCCont.fixedTemperatureGradient( rollerBackPatchName_);
            }
            else
            {
                FatalErrorIn
                (
                "Foam::dictionary& Foam::rollerPassItem::boundaryConditionsDU()"
                )
                    << "Symmetry " << processLineSymm << " is not possible"
                    << "with roll setup " << rollerType_ << " for passItem "
                    << name() << "!" << abort(FatalError);
            }
         }
        else if (rollerType_ == "rightRoll")
        {
            if(processLineSymm == "posZ" || processLineSymm == "none")
            {
                BCCont.fixedTemperatureGradient( rollerBackPatchName_);
                BCCont.fixedTemperatureGradient( rollerFrontPatchName_);
            }
            else if
            (
                processLineSymm == "posY" || processLineSymm == "posYposZ"
            )
            {
                BCCont.fixedTemperatureGradient( rollerFrontPatchName_);
                BCCont.thermalSymmetry( rollerSymmetryYPatchName_);
            }
            else if
            (
                processLineSymm == "negY" || processLineSymm == "negYposZ"
            )
            {
                BCCont.thermalSymmetry( rollerSymmetryYPatchName_);
                BCCont.fixedTemperatureGradient( rollerBackPatchName_);
            }
            else
            {
                FatalErrorIn
                (
                "Foam::dictionary& Foam::rollerPassItem::boundaryConditionsDU()"
                )
                    << "Symmetry " << processLineSymm << " is not possible"
                    << "with roll setup " << rollerType_ << " for passItem "
                    << name() << "!" << abort(FatalError);
            }
        }
        else if (rollerType_ == "leftRoll")
        {
            if(processLineSymm == "negZ" || processLineSymm == "none")
            {
                BCCont.fixedTemperatureGradient( rollerBackPatchName_);
                BCCont.fixedTemperatureGradient( rollerFrontPatchName_);
            }
            else if
            (
                processLineSymm == "posY" || processLineSymm == "posYnegZ"
            )
            {
                BCCont.fixedTemperatureGradient( rollerFrontPatchName_);
                BCCont.thermalSymmetry( rollerSymmetryYPatchName_);
            }
            else if
            (
                processLineSymm == "negY" || processLineSymm == "negYnegZ"
            )
            {
                BCCont.thermalSymmetry( rollerSymmetryYPatchName_);
                BCCont.fixedTemperatureGradient( rollerBackPatchName_);
            }
            else
            {
                FatalErrorIn
                (
                "Foam::dictionary& Foam::rollerPassItem::boundaryConditionsDU()"
                )
                    << "Symmetry " << processLineSymm << " is not possible"
                    << "with roll setup " << rollerType_ << " for passItem "
                    << name() << "!" << abort(FatalError);
            }
        }

        BCCont.thermalConvection(rollerToAirPatchName_);

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


void Foam::rollerPassItem::rollerToRollerBoundaryConditions()
{

}


Foam::word Foam::rollerPassItem::contactPatch()
{
    return word( rollerContactPatchName_);
}


Foam::scalar Foam::rollerPassItem::frictionCoefficientAxial()
{
    return frictionCoefficientAxial_;
}


List<dictionary> rollerPassItem::functionObjects()
{
    // By default we will return no function objects; this can be
    // overwritten by the specific pass items

    // TODO: there can be rollers in different positions. Make functionObjects
    // accordingly.

    int nrOfFunctionObjects = 3;
    List<dictionary> dicts(nrOfFunctionObjects);

    {
        dictionary subDict;
        subDict.add("type", "solidForces");
        subDict.add("historyPatch",  rollerAxisPatchName_);

        dicts[0] = subDict;
    }

    if (motionEnabler_ == "drivenRoller")
    {
        {
            dictionary subDict;
            subDict.add("type", "solidPower");
            subDict.add("historyPatch",  rollerAxisPatchName_);
            dicts[1] = subDict;
        }
        {
            dictionary subDict;
            subDict.add("type", "solidTorque");
            subDict.add("historyPatch",  rollerAxisPatchName_);
            dicts[2] = subDict;
        }
    }
    else
    {
        dicts.resize(nrOfFunctionObjects - 2);
    }

    return dicts;
}


// ************************************************************************* //

} // end of namespace foam
