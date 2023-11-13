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

#include "wirePassItem.H"
#include "addToRunTimeSelectionTable.H"
#include "processPass.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wirePassItem, 0);
    addToRunTimeSelectionTable(passItem, wirePassItem, dictionary);
}

// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wirePassItem::wirePassItem
(
    const processPass& parent,
    const word& name,
    const dictionary& dict,
    const int& passNr,
    dataHandler* dataContainerPtr
)
:
    passItem(parent, name, dict, passNr, dataContainerPtr),
    createMesh_(dict.subDict("mesh").lookup("generateNewWireMesh")),
    boundaryConditionsDUGenerated_(false),
    boundaryConditionsTGenerated_(false),
    wireDownstreamPatchName_(),
    wireUpstreamPatchName_(),
    wireContactPatchName_(),
    wireContactShadowPatchName_("notSet"),
    frictionCoefficientsAxial_(),
    shadowPatchNames_(),
    wireSymmetryYPatchName_(),
    wireSymmetryZPatchName_(),
    wireFrontPatchName_(),
    wireBackPatchName_()
{
    const dictionary patchNames =
        dataContainer().dataContainerDict().subDict("patchNames");

    wireUpstreamPatchName_ =
        word(patchNames.lookup("wireUpstreamPatchName"));

    wireDownstreamPatchName_ =
        word(patchNames.lookup("wireDownstreamPatchName"));

    wireContactPatchName_ =
        word(patchNames.lookup("wireContactPatchName"));

    wireSymmetryYPatchName_ =
        word(patchNames.lookup("wireSymmetryYPatchName"));

    wireSymmetryZPatchName_ =
        word(patchNames.lookup("wireSymmetryZPatchName"));

    wireFrontPatchName_ =
        word(patchNames.lookup("wireFrontPatchName"));

    wireBackPatchName_ =
        word(patchNames.lookup("wireBackPatchName"));

    this->passItem::setWireLength
    (
        readScalar(dict.subDict("mesh").lookup("steadyStateWireLength"))
    );
    this->passItem::setWireAxialResolution
    (
        readScalar(dict.subDict("mesh").lookup("axialResolution"))
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //


Foam::wirePassItem::~wirePassItem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::wirePassItem::writeMeshDict(Time& runTime, fileName geometryFile)
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
        FatalErrorIn("void Foam::wirePassItem::writeMeshDict(Time& runTime)")
            << "Unknown processLineSymmetry " << processLineSymmetry
            << ", for passItem " << name()
            << ". The options are: posYposZ, negYposZ, negYnegZ, posYnegZ, "
            << "posY, posZ, negY, negZ, none, axiSymmetric" << abort(FatalError);
    }

    // Get wire mesh input dictionary from the passItem dictionary
    dictionary meshInputDict = passItem::meshInputDict();

    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            "system",
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dictionary("meshDict")
    );

    // Get wedgeAngle in case of an axisymmetric process line
    scalar wedgeAngle = 0;
    if (processLineSymmetry == "axiSymmetric")
    {
        wedgeAngle = dataContainer().returnKey<scalar>("wedgeAngle");
    }

    // Add settings to meshDict Dict
    meshDict.add
    (
        "enforceGeometryConstraints",
        meshInputDict.lookupOrDefault<int>("enforceGeometryConstraints", 0)
    );
    meshDict.add
    (
        "keepCellsIntersectingBoundary",
        meshInputDict.lookupOrDefault<int>("keepCellsIntersectingBoundary", 1)
    );
    meshDict.add
    (
        "maxCellSize",
        readScalar(meshInputDict.lookup("maxCellSize"))
    );
    meshDict.add
    (
        "geometryTolerance",
        meshInputDict.lookupOrDefault<scalar>("geometryTolerance", 1e-5)
    );
    meshDict.add
    (
        "minContactCellSize",
        meshInputDict.lookupOrDefault<scalar>("minContactCellSize", 25e-5)
    );

    // Create boundaryLayers subDict
    dictionary boundaryLayers;
    boundaryLayers.add
    (
        "nLayers",
        meshInputDict.lookupOrDefault<int>("nLayers", 0)
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

    meshDict.add("symmetryType", processLineSymmetry);
    meshDict.add("rollSetup", "none");

    // Create patchNames dictionary to set the names as specified in the
    // dataContainer
    dictionary patchNames =
        dataContainer().dataContainerDict().subDict("patchNames");

    meshDict.add("rollingMillPatchNames", patchNames);

    // Create, set and add wireMeshDict subDict
    {
        dictionary wireMeshDict;

        if (geometryFile == "none")
        {
            wireMeshDict.add
            (
                "wireDiameter",
                readScalar(meshInputDict.lookup("wireDiameter"))
            );
            wireMeshDict.add
            (
                "wireLength",
                readScalar(meshInputDict.lookup("initialWireLength"))
            );
        }
        else if (geometryFile.ext() == "fms")
        {
            wireMeshDict.add("surfaceFile", geometryFile);
            wireMeshDict.add
            (
                "wireLength",
                readScalar(meshInputDict.lookup("newWireLength"))
            );
        }
        else
        {
            wireMeshDict.add("geometryFile", geometryFile);
            wireMeshDict.add
            (
                "wireLength",
                readScalar(meshInputDict.lookup("newWireLength"))
            );
        }
//        wireMeshDict.add
//        (
//            "wireLength",
//            readScalar(meshInputDict.lookup("initialWireLength"))
//        );
        wireMeshDict.add
        (
            "axialResolution",
            readScalar(meshInputDict.lookup("axialResolution"))
        );
        wireMeshDict.add
        (
            "axialShift",
            meshInputDict.lookupOrDefault<scalar>("axialShift", 0)
        );
        wireMeshDict.add
        (
            "axialGrading",
            meshInputDict.lookupOrDefault<scalar>("axialGrading", 1.0)
        );
        wireMeshDict.add
        (
            "maxCellSize",
            meshInputDict.lookupOrDefault<scalar>("maxCellSize", 0.002)
        );
        meshDict.add("wireMeshDict", wireMeshDict);
    }

    // Add dieMeshDict
    dictionary dieMeshDict;
    dieMeshDict.add("wedgeAngle", wedgeAngle);
    meshDict.add("dieMeshDict", dieMeshDict);

    // Write meshDict dict
    meshDict.regIOobject::write();
}


void Foam::wirePassItem::setupMesh(Time& runTime)
{
    if (createMesh_)
    {
        Info<< "    " << name() << ": setting up mesh" << endl;

        // Change to the processPass directory and write the meshDict
        chDir(runTime.path());

        fileName geometryFile("none");
        if (dict().subDict("mesh").found("geometryFile"))
        {
             geometryFile =
                 fileName(dict().subDict("mesh").lookup("geometryFile"));

             Info<< "    geometryFile found: " << geometryFile << endl;

             if (geometryFile == "wireProfile.prof")
             {
                 FatalErrorIn
                     (
                         "Foam::dictionary& "
                         "Foam::wirePassItem::setupMesh(Time& runTime)"
                     )
                     << nl << "    wireProfile.prof is a preserved fileName "
                     "for remeshing the wire between passes. Please specify "
                     "another geometryFile name"
                     << abort(FatalError);
             }
        }

        // Copy in default meshDict
        WarningIn("void Foam::wirePassItem::setupMesh(Time& runTime)")
            << "Why are we using runSystemCommand here for a cp?" << endl;
        dataContainer().runSystemCommand
        (
            "cp ../../dictFiles/defaultDicts/meshDict system/",
            "log.wirePassItem_setupMesh_meshDict",
            "void Foam::processPass::setup()"
        );

        writeMeshDict(runTime, geometryFile);

        // Copy meshDict to system directory for mesh generation utility
        //cp("system/" + name() + "/meshDict","system/meshDict");
        // Make a copy of the meshDict for debugging
        cp("system/meshDict", "system/" + name() + "/meshDict");

        // Run the mesh utility
        dataContainer().runSystemCommand
        (
            "export OMP_NUM_THREADS=1; "
            "generateWireMesh",
            "log.generateMesh" + name(),
            "void rollerPassItem::setupMesh"
        );

        // Move mesh from constant/polymesh to the subMesh directory and clean
        // meshDict from system directory
        mv(runTime.path()/"constant/polyMesh","constant"/name()/"polyMesh");
        rm("system/meshDict");
    }
    else
    {
        Info<< "    " << name() << ": mesh is not constructed."
            << " It will to be copied"
            << " from previous pass prior to merging." << endl;
    }
}


Foam::dictionary& Foam::wirePassItem::boundaryConditionsDU()
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

        // Specify downstream boundary conditions. This depends on how the wire
        // motion is applied, i.e., via a roller or a capstan.
        word motionEnabler =
            passItem::boundaryConditions().lookupOrDefault<word>
            (
                "motionEnabler", "capstan"
            );

        if (motionEnabler != "roller" && motionEnabler != "capstan")
        {
            WarningIn
                ("Foam::dictionary& Foam::wirePassItem::boundaryConditionsDU()")
                << "Unknown method of enabling wire motion: " << motionEnabler
                <<". Setting to default capstan method, i.e., wire will be "
                << "pulled" << endl;
            motionEnabler = "capstan";
        }

        if (motionEnabler == "roller")
        {
            // PC, 13Apr18: looks for the timeVsPressure in the
            // "dictFiles/<passName>" sub-directory
            BCCont.fixedTangentialNormalPressure
            (
                wireDownstreamPatchName_,
                parent().name()
            );
        }
        else if (motionEnabler == "capstan")
        {
            // PC, 13Apr18: looks for the timeVsDisp in the
            // "dictFiles/<passName>" sub-directory
            BCCont.fixedDisplacementZeroShear
            (
                wireDownstreamPatchName_,
                parent().name()
            );
        }

        // The upstream boundary of the wire always gets a backPull tension
        // PC, 13Apr18: looks for the timeVsPressure in the
        // "dictFiles/<passName>" sub-directory
        BCCont.fixedTangentialNormalPressure
        (
            wireUpstreamPatchName_,
            parent().name()
        );

        // Wire is considered slave by default, except when there are multiple
        // shadow patches and it is specified that the tool contact patches
        // form no single contact patch

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

        const Switch rigidTool =
            passItem::boundaryConditions().lookupOrDefault<Switch>
            (
                "rigidTool", "no"
            );

        if (!rigidTool)
        {
            if (nrOfPassItems > 2 && !singleRollerContact)
            {
                BCCont.solidContactMaster
                (
                    wireContactPatchName_,
                    shadowPatchNames_,
                    frictionCoefficientsAxial_
                );
            }
            else
            {
                BCCont.solidContactSlave
                (
                    wireContactPatchName_,
                    shadowPatchNames_[0]
                );
            }
        }
        else
        {
            const dictionary triSurfacesDict =
                    passItem::boundaryConditions().subDict("triSurfaces");

            chDir
            (
                dataContainer().returnKey<fileName>("processLineDir")/
                dataContainer().returnKey<word>
                (
                    "passName",
                    passItem::passNr()
                )
            );

            mkDir("constant/triSurfaces");

            rmDir("constant/polyMesh");

            cp("constant"/passItem::name()/"polyMesh", "constant/polyMesh");

            Info << "wire? "
                 << fileName("constant"/passItem::name()/"polyMesh")
                 << endl;

            BCCont.solidRigidContact
            (
                wireContactPatchName_,
                triSurfacesDict,
                frictionCoefficientsAxial_
            );
        }

        // Specify symmetry boundary conditions
        if (processLineSymm == "posYposZ" || processLineSymm == "negYposZ" ||
            processLineSymm == "negYnegZ" || processLineSymm == "posYnegZ")
        {
            BCCont.solidSymmetry(wireSymmetryYPatchName_);
            BCCont.solidSymmetry(wireSymmetryZPatchName_);
        }
        else if (processLineSymm == "posY" || processLineSymm == "negY")
        {
            BCCont.solidSymmetry(wireSymmetryYPatchName_);
        }
        else if (processLineSymm == "posZ" || processLineSymm == "negZ")
        {
            BCCont.solidSymmetry(wireSymmetryZPatchName_);
        }
        else if (processLineSymm == "axiSymmetric")
        {
            BCCont.solidWedge(wireFrontPatchName_);
            BCCont.solidWedge(wireBackPatchName_);
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


Foam::dictionary& Foam::wirePassItem::boundaryConditionsT()
{
    // Get reference to the DU input boundary conditions dictionary
    dictionary& dict = passItem::boundaryConditionsT();

    // Check if the T boundary conditions dicitonary was already generated.
    // If no, generate it. If yes, do not generate it and return the existing
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
        BCCont.fixedTemperatureGradient(wireDownstreamPatchName_);
        BCCont.fixedTemperature(wireUpstreamPatchName_);


        // Wire is considered slave by default, except when there are multiple
        // shadow patches and it is specified that the tool contact patches
        // form no single contact patch

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
            BCCont.thermalContactMaster
            (
                wireContactPatchName_,
                shadowPatchNames_
            );
        }
        else
        {
            BCCont.thermalContactSlave
            (
                wireContactPatchName_,
                shadowPatchNames_[0]
            );
        }

//        // Get shadowPatch name if specified in boundaryConditions dictionary
//        if (passItem::boundaryConditions().found("shadowPatch"))
//        {
//            wireContactShadowPatchName_ =
//                word(passItem::boundaryConditions().lookup("shadowPatch"));
//        }

//        if (shadowPatchNames_.size() > 1)
//        {
//            BCCont.thermalContactMaster
//            (
//                wireContactPatchName_,
//                shadowPatchNames_
//            );
//        }
//        else if
//        (
//                (shadowPatchNames_.size() == 1 &&
//                wireContactShadowPatchName_ == "notSet") ||
//                wireContactShadowPatchName_ != "notSet"
//        )
//        {
//            if
//            (
//                shadowPatchNames_.size() == 1 &&
//                wireContactShadowPatchName_ == "notSet"
//            )
//            {
//                wireContactShadowPatchName_ = shadowPatchNames_[0];
//            }
//            BCCont.thermalContactMaster
//            (
//                wireContactPatchName_,
//                wireContactShadowPatchName_
//            );
//        }
//        BCCont.thermalContactMaster
//        (
//            wireContactPatchName_, wireContactShadowPatchName_
//        );

        // Specify symmetry boundary conditions
        if (processLineSymm == "posYposZ" || processLineSymm == "negYposZ" ||
            processLineSymm == "negYnegZ" || processLineSymm == "posYnegZ")
        {
            BCCont.thermalSymmetry(wireSymmetryYPatchName_);
            BCCont.thermalSymmetry(wireSymmetryZPatchName_);
        }
        else if (processLineSymm == "posY" || processLineSymm == "negY")
        {
            BCCont.thermalSymmetry(wireSymmetryYPatchName_);
        }
        else if (processLineSymm == "posZ" || processLineSymm == "negZ")
        {
            BCCont.thermalSymmetry(wireSymmetryZPatchName_);
        }
        else if (processLineSymm == "axiSymmetric")
        {
            BCCont.thermalWedge(wireFrontPatchName_);
            BCCont.thermalWedge(wireBackPatchName_);
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


Foam::word Foam::wirePassItem::contactPatch()
{
    return wireContactPatchName_;
}


Foam::scalar Foam::wirePassItem::frictionCoefficientAxial()
{
    WarningIn
        (
            "Foam::scalar Foam::wirePassItem::"
            "frictionCoefficientAxial()"
        )
        << "Returning zero as axial friction coefficient of the wire!" << endl;
    return 0;
}


void Foam::wirePassItem::addShadowPatchName(word shadowPatchName)
{
    // Add shadowPatchName to shadowPatchNames_ list
    shadowPatchNames_.resize(shadowPatchNames_.size() + 1);
    shadowPatchNames_[shadowPatchNames_.size() - 1] = shadowPatchName;
}


void Foam::wirePassItem::addFrictionCoefficientAxial
(
    scalar frictionCoefficientAxial
)
{
    // Add friction coefficient to frictionCoefficientsAxial_ list
    frictionCoefficientsAxial_.resize(frictionCoefficientsAxial_.size() + 1);
    frictionCoefficientsAxial_[frictionCoefficientsAxial_.size() - 1] =
        frictionCoefficientAxial;
}


Foam::wordList& Foam::wirePassItem::shadowPatchNames()
{
    return shadowPatchNames_;
}


void Foam::wirePassItem::positionMesh
(
    const Time& runTime,
    const fvMesh& wireMesh
)
{
    vector translation = vector::zero;

    // Read positioning vector if found, else position the wire so is
    // symmetrically position about the y-z plane (x = 0); this means that half
    // the wire has negative x coordinates and half the wire has positive x
    // coordinates
    vector positioningVector = vector::zero;
    if (meshInputDict().found("positioningVector"))
    {
        positioningVector = meshInputDict().lookup("positioningVector");
    }
    boundBox wireBoundBox(wireMesh.points());

//    Info << "z coordinate: " << wireBoundBox.min().z() << " ; " << wireBoundBox.max().z() << endl;

    translation = vector(-wireBoundBox.midpoint().x(), 0, 0)
                + positioningVector;

    passItem::positionMesh(runTime, translation);
}


Foam::List<Foam::dictionary> Foam::wirePassItem::functionObjects()
{
    // By default we will return no function objects; this can be
    // overwritten by the specific pass items

    List<dictionary> dicts(5);

    {
        dictionary subDict;
        subDict.add("type", "solidForces");
        subDict.add("historyPatch", wireDownstreamPatchName_);

        dicts[0] = subDict;
    }

    {
        dictionary subDict;
        subDict.add("type", "solidForces");
        subDict.add("historyPatch", wireUpstreamPatchName_);

        dicts[1] = subDict;
    }

    {
        dictionary subDict;
        subDict.add("type", "solidPower");
        subDict.add("historyPatch", wireDownstreamPatchName_);

        dicts[2] = subDict;
    }

    {
        dictionary subDict;
        subDict.add("type", "solidPower");
        subDict.add("historyPatch", wireUpstreamPatchName_);

        dicts[3] = subDict;
    }

    {
        dictionary subDict;
        subDict.add("type", "solidForces");
        subDict.add("historyPatch", wireContactPatchName_);

        dicts[4] = subDict;
    }

    return dicts;
}


void Foam::wirePassItem::writeCurrentProfile(const Time& runTime)
{
    Info<< "    " << type() << " : writing current wire profile" << endl;

    // Read the current wire mesh
    const fvMesh& wireMesh = mesh(runTime);

    // Get the patch index for the contact patch
    const label contactPatchID =
        wireMesh.boundaryMesh().findPatchID(wireContactPatchName_);

    if (contactPatchID == -1)
    {
        FatalErrorIn
        (
            "void Foam::wirePassItem::writeCurrentProfile(const Time& runTime)"
        )   << "wire contact patch not found!" << abort(FatalError);
    }

    // Take a reference to the patch
    const polyPatch& ppatchC = wireMesh.boundaryMesh()[contactPatchID];

    // Get an ordered list of edges around the boundary of the patch
    // We expect that there is only one wire (for now)
    if (ppatchC.edgeLoops().size() != 1)
    {
        if (contactPatchID == -1)
        {
            FatalErrorIn
            (
                "void Foam::wirePassItem::"
                "writeCurrentProfile(const Time& runTime)"
            )   << "There is more than one edgeLoop on the contact patch!"
                << " Is there more than 1 wire?!"
                << abort(FatalError);
        }
    }

    const labelList& edgeLoopC = ppatchC.edgeLoops()[0];


    // Get the patch index for the downstream patch
    const label downstreamPatchID =
        wireMesh.boundaryMesh().findPatchID(wireDownstreamPatchName_);

    if (downstreamPatchID == -1)
    {
        FatalErrorIn
        (
            "void Foam::wirePassItem::writeCurrentProfile(const Time& runTime)"
        )   << "wire downstream patch not found!" << abort(FatalError);
    }

    // Take a reference to the patch
    const polyPatch& ppatch = wireMesh.boundaryMesh()[downstreamPatchID];

    // Get a list of the points on the boundary of the downstream patch
    //const labelList& boundaryPoints = ppatch.boundaryPoints();

    // Get an ordered list of edges around the boundary of the patch
    // We expect that there is only one wire (for now)
    if (ppatch.edgeLoops().size() != 1)
    {
        if (downstreamPatchID == -1)
        {
            FatalErrorIn
            (
                "void Foam::wirePassItem::"
                "writeCurrentProfile(const Time& runTime)"
            )   << "There is more than one edgeLoop on the downstream patch!"
                << " Is there more than 1 wire?!"
                << abort(FatalError);
        }
    }

    const labelList& edgeLoop = ppatch.edgeLoops()[0];

    // Get a list of the point coordinates on the downstream patch
    const pointField& localPoints = ppatch.localPoints();

    // Get a list of the point coordinates on the contact patch
    const pointField& localPointsC = ppatchC.localPoints();

    if (debug)
    {
        OFstream downStreamProfileFile("downstreamPatch.points");
        forAll(localPoints, wpI)
        {
            if (localPoints[wpI].x() == 0.04)
            {
                downStreamProfileFile << localPoints[wpI] * 1000.0 << endl;
            }
        }
        OFstream contactProfileFile("contact.points");
        forAll(localPointsC, wpI)
        {
            if (localPointsC[wpI].x() == 0.04)
            {
                contactProfileFile << localPointsC[wpI] * 1000.0 << endl;
            }
        }
    }


    DynamicList<vector> wireProfilePoints;

    // Find mutual point labels of contact and downstream patches
//    forAll(localPoints, pI)
//    {
//        forAll(localPointsC, pCI)
//        {
    // Loop over edgeloops as the points are sorted.
    forAll(edgeLoop, eI)
    {
        forAll(edgeLoopC, eCI)
        {
            const label pI = edgeLoop[eI];
            const label pCI = edgeLoopC[eCI];
            if
            (
                    mag(localPoints[pI].x() - localPointsC[pCI].x()) <= SMALL &&
                    mag(localPoints[pI].y() - localPointsC[pCI].y()) <= SMALL &&
                    mag(localPoints[pI].z() - localPointsC[pCI].z()) <= SMALL
            )
            {
                vector wirePoint =
                    vector
                    (
                        localPoints[pI].x()*1000.0,
                        localPoints[pI].y()*1000.0,
                        localPoints[pI].z()*1000.0
                    );

//                if (debug)
//                {
//                    vector wirePointC =
//                        vector
//                        (
//                            localPointsC[pCI].x()*1000.0,
//                            localPointsC[pCI].y()*1000.0,
//                            localPointsC[pCI].z()*1000.0
//                        );
//                    Info<< "Mutual point found. " << endl;

//                    Info<< "    wirePoint nr "<< pI << " with pos: "
//                        << wirePoint << " or components: " << wirePoint.x()
//                        << ", " << wirePoint.y() << ", " << wirePoint.z()
//                        << endl;

//                    Info<< "    wirePointC nr "<< pCI << " with pos: "
//                        << wirePointC << " or components: " << wirePointC.x()
//                        << ", " << wirePointC.y() << ", " << wirePointC.z()
//                        << endl;

//                    Info<< "    Checks: "
//                        << mag(wirePoint.x() - wirePointC.x()) << "; "
//                        << mag(wirePoint.y() - wirePointC.y()) << "; "
//                        << mag(wirePoint.z() - wirePointC.z()) << " - "
//                        << Switch(abs(wirePoint.x()-wirePointC.x()) < SMALL)
//                        << "; "
//                        << Switch(abs(wirePoint.y()-wirePointC.y()) < SMALL)
//                        << "; "
//                        << Switch(abs(wirePoint.z()-wirePointC.z()) < SMALL)
//                        << endl;
//                }

                wirePoint.x() = 0.0;

                wireProfilePoints.append(wirePoint);
            }
        }
    }

    if (wireProfilePoints.capacity() == 0)
    {
        FatalErrorIn
        (
            "void Foam::wirePassItem::"
            "writeCurrentProfile(const Time& runTime)"
        )   << "No mutual points for downstream and contact patch found. "
            << "No wire profile generated."
            << abort(FatalError);
    }

    // Add first point to ensure the loop is closed
    wireProfilePoints.append(wireProfilePoints[0]);

    // Write the points to file
    const fileName profileFileName("wireProfile.prof");

    OFstream profileFile(profileFileName);

    forAll(wireProfilePoints, wpI)
    {
        profileFile << wireProfilePoints[wpI] << endl;
    }

/*
    List<vector> wirePoints;
    wirePoints.setSize(edgeLoop.size() + 1);

    forAll(edgeLoop, bpI)
    {
        const label pointID = edgeLoop[bpI];

        vector wirePoint =
            vector
            (
                0,
                localPoints[pointID].y()*1000.0,
                localPoints[pointID].z()*1000.0
            );

        wirePoints[bpI] = wirePoint;
    }
    wirePoints[wirePoints.size()-1] = wirePoints[0];

    forAll(wirePoints, wpI)
    {
        profileFile << wirePoints[wpI] << endl;
    }
*/


    // TODO: write mode while writing the re-meshing meshDict?
    // Mode is needed for cfMesh to re-mesh the wire
    profileFile << "mode b-spline" << endl;

    Info << "wire profile name: "
         << "wireProfile_" + runTime.caseName() + ".fms" << endl;


}


// ************************************************************************* //
