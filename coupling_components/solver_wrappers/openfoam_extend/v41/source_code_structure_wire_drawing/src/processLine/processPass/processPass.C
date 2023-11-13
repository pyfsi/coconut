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

#include "processPass.H"
#include "polyMesh.H"
#include "IOobjectList.H"
#include "boundBox.H"
#include "IFstream.H"
#include "gnuplotGraph.H"
#include "rollerPassItem.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processPass, 0);
}

using namespace Foam::mathematicalConstant;

// * * * * * * * * * *  Private Member Functions  * * * * * * * * * * * * * * //

void Foam::processPass::readMesh() const
{
    if (!meshPtr_.empty())
    {
        FatalErrorIn("void Foam::processPass::readMesh()")
            << "pointer already set" << abort(FatalError);
    }

    if(debug)
    {
        Info<< "        Reading the mesh for time = " << runTime().timeName()
            << endl;
    }

    meshPtr_.set
    (
        new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTime().timeName(),
                runTime(),
                IOobject::MUST_READ
            )
        )
    );
}


Foam::PtrList<Foam::passItem>& Foam::processPass::passItems()
{
    if (passItems_.size() == 0)
    {
        makePassItems();
    }

    return passItems_;
}


const Foam::PtrList<Foam::passItem>& Foam::processPass::passItems() const
{
    if (passItems_.size() == 0)
    {
        makePassItems();
    }

    return passItems_;
}


void Foam::processPass::makePassItems() const
{
    if (passItems_.size() != 0)
    {
        FatalErrorIn("void processPass::makePassItems()")
            << "pointer list already set!" << abort(FatalError);
    }

    // Read the passItems dict
    const dictionary& passItemsDict = this->passItemsDict();

    // Get the list of names of the passItems
    const wordList passItemsNames = passItemsDict.toc();

    // Check all passItems are specified as subDicts
    forAll(passItemsNames, passItemI)
    {
        if (!passItemsDict.isDict(passItemsNames[passItemI]))
        {
            FatalErrorIn
            (
                "processPass::processPass\n"
                "(\n"
                "    const word& processLineDir,\n"
                "    const word& passName,\n"
                "    const dictionary& dict\n"
                ")"
             )   << "All entries in the passItems list should be "
                 << "subDicts, but passItem number " << (passItemI + 1)
                 << " is not!" << abort(FatalError);
        }
    }

    Info<< "    Number of entries in passItems dictionary: "
        << passItemsNames.size() << endl;

    // Add number of passItems as it will trigger the contact setup.
    // This number is always number of tools +1 (i.e., the wire)
    dataContainer().addKey<scalar>
    (
        "nrOfPassItems",
        passItemsNames.size(),
        passNr()
    );

    // Create a list passItem objects
    passItems_.setSize(passItemsNames.size());

    int passItemI = 0;
    forAll(passItemsNames, passItemNameI)
    {
        const word& passItemName = passItemsNames[passItemNameI];

        const dictionary& passItemDict =
            passItemsDict.subDict(passItemsNames[passItemNameI]);

        const word passItemType = word(passItemDict.lookup("type"));

        if (passItemType == "roller")
        {
            const word rollSetup = dataContainer().returnKey<word>
                (
                    "rollSetup",
                    passNr_
                );

            const word processLineSymm =
                dataContainer().returnKey<word>("processLineSymmetry");

            const word passItemRollType = passItemDict.lookup("rollerType");

            // Due to the choice of symmetry, some rollers may be inactive; only
            // make rollers that are active
            if (rollValid(rollSetup, processLineSymm, passItemRollType))
            {
                // Create roller object
                passItems_.set
                (
                    passItemI,
                    passItem::New
                    (
                        *this,
                        passItemName,
                        passItemsDict.subDict(passItemName),
                        passNr_,
                        &dataContainer()
                    )
                );
                passItemI++;
            }
            else
            {
                int passItemsSize = passItems_.size();
                passItems_.resize(passItemsSize -1);
            }
        }
        else
        {
            // Create die or wire object
            passItems_.set
            (
                passItemI,
                passItem::New
                (
                    *this,
                    passItemName,
                    passItemsDict.subDict(passItemName),
                    passNr_,
                    &dataContainer()
                )
            );
            passItemI++;
        }
    }

    // Add number of passItems as it will trigger the contact setup.
    // This number is always number of tools +1 (i.e., the wire)
    dataContainer().addKey<scalar>
    (
        "nrOfPassItems",
        passItems_.size(),
        passNr()
    );

    Info<< "    Number of passItems created: " << passItems_.size()
         << ": ";
    forAll(passItems_, passItemI)
    {
        Info<< passItems_[passItemI].name();
        if(passItemI < passItems_.size()-1)
        {
            Info<< ", ";
        }
    }
    Info<< endl;
}


void Foam::processPass::mergePassItemMeshes()
{
    Info<< "    Merging meshes" << endl;

    // Change to pass directory
    chDir(runTime().path());

    // Construct command, which specifies subMesh region names
    word mergeSubMeshesCommand = "mergeSubMeshes";
    forAll(passItems(), itemI)
    {
        mergeSubMeshesCommand += " " + passItems()[itemI].name();
    }

    dataContainer().runSystemCommand
    (
        mergeSubMeshesCommand,
        "log.mergeSubMeshes",
        "void Foam::processPass::mergePassItemMeshes()"
    );

    // Remove face zones
    rm("constant/polyMesh/faceZones");

    if(dataContainer().returnKey<word>("processLineSymmetry") == "axiSymmetric")
    {
        // Make two patches; one consisting the front die and wire patch and
        // one consisting of the back die and wire patch. This because the
        // solver insists on having exactly 2 wedge patches in case of
        // axi-symmetry.
        // Write createPatch dictionary
        writeCreateWedgePatchDict();

        dataContainer().runSystemCommand
        (
            "createPatch -overwrite",
            "log.createPatchOverwrite",
            "void Foam::processPass::mergePassItemMeshes()"
        );
    }
}


void Foam::processPass::mergePassItemMaterialProperties()
{
    Info<< "    Merging material properties" << endl;

    // Mechanical properties
    {
        // Creating mechanicalProperties dictionary
        IOdictionary mechanicalProperties
        (
            IOobject
            (
                "mechanicalProperties",
                runTime().constant(),
                runTime(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            dictionary("mechanicalProperties")
        );

        // Add settings to mechanicalProperties dict
        mechanicalProperties.add("planeStress", "no");

        // Create mechanical subDict
        dictionary mechanical;

        // Add settings to mechanical subDict
        mechanical.add("type", "multiMaterial");

        // Create empty list of laws
        dictionary laws;

        forAll(passItems(), itemI)
        {
            // Get mechanical law from the pass item and add it to the laws list
            laws.add
            (
                passItems()[itemI].name(),
                passItems()[itemI].mechanicalLawDict()
            );
        }

        // Update the laws list in the mechanical subDict
        mechanical.add("laws", laws);

        // Add mechanical subDict to mechanicalProperties dict
        mechanicalProperties.add("mechanical", mechanical);

        // Write out the mechanicalProperties dictionary
        mechanicalProperties.regIOobject::write();
    }


    // Thermal properties
    {
        // Creating thermalProperties dictionary
        IOdictionary thermalProperties
        (
            IOobject
            (
                "thermalProperties",
                runTime().constant(),
                runTime(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            dictionary("thermalProperties")
        );

        // Add settings to thermalProperties dict
        thermalProperties.add("planeStress", "no");

        // Create thermal subDict
        dictionary thermal;

        // Add settings to thermal subDict
        thermal.add("type", "multiMaterial");

        // Create empty list of laws
        PtrList<entry> laws;

        laws.resize(passItems().size());

        forAll(passItems(), itemI)
        {
            // Get thermal law from the pass item and add it to the law list
            laws.set
                (
                    itemI,
                    new dictionaryEntry
                    (
                        passItems()[itemI].name(),
                        thermal,
                        passItems()[itemI].thermalLawDict()
                    )
                );
        }

        // Update the laws list in the thermal subDict
        thermal.add("laws", laws);

        // Add thermal subDict to thermalProperties dict
        thermalProperties.add("thermal", thermal);

        // Write out the thermalProperties dictionary
        thermalProperties.regIOobject::write();
    }
}


void Foam::processPass::mergePassItemBoundaryConditions()
{
    Info<< "    Merging boundary conditions" << endl;

    const dictionary patchNames =
        dataContainer().dataContainerDict().subDict("patchNames");

    const word rollerToRollerFrontPatchName =
            patchNames.lookup("rollerToRollerFrontPatchName");
    const word rollerToRollerBackPatchName =
            patchNames.lookup("rollerToRollerBackPatchName");
    const word rollerToAirPatchName =
            patchNames.lookup("rollerToAirPatchName");

    const word topRollFrontRTRPatchName =
            word("topRoll") + rollerToRollerFrontPatchName;
    const word topRollBackRTRPatchName =
            word("topRoll") + rollerToRollerBackPatchName;
    const word bottomRollFrontRTRPatchName =
            word("bottomRoll") + rollerToRollerFrontPatchName;
    const word bottomRollBackRTRPatchName =
            word("bottomRoll") + rollerToRollerBackPatchName;
    const word rightRollFrontRTRPatchName =
            word("rightRoll") + rollerToRollerFrontPatchName;
    const word rightRollBackRTRPatchName =
            word("rightRoll") + rollerToRollerBackPatchName;
    const word leftRollFrontRTRPatchName =
            word("leftRoll") + rollerToRollerFrontPatchName;
    const word leftRollBackRTRPatchName =
            word("leftRoll") + rollerToRollerBackPatchName;

    const word rollSetup =
            dataContainer().returnKey<word>("rollSetup", passNr_);

    // DU field
    {
        // Creating DU dictionary
        IOdictionary DU
        (
            IOobject
            (
                "DU",
                "0",
                runTime(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            dictionary("DU")
        );

        // Add settings to DU field
        DU.add("dimensions", "[0 1 0 0 0 0 0]");
        DU.add("internalField", "uniform ( 0 0 0 )");

        // Create boundaryField subDict
        dictionary boundaryField;

        forAll(passItems(), itemI)
        {
            // Read passItem subMesh
            polyMesh subMesh
            (
                IOobject
                (
                    passItems()[itemI].name(),
                    runTime().timeName(),
                    runTime(),
                    IOobject::MUST_READ
                )
            );

            // Read boundary conditions for this passItem
            dictionary& itemIBCs = passItems()[itemI].boundaryConditionsDU();

            // Add boundary condition for each patch on the subMesh boundary
            forAll(subMesh.boundaryMesh(), patchI)
            {
                const word& patchName = subMesh.boundaryMesh()[patchI].name();
                if (itemIBCs.isDict(patchName))
                {
                    const dictionary& DUpatchDict = itemIBCs.subDict(patchName);
                    if (!boundaryField.found(patchName))
                    {
                        boundaryField.add(patchName, DUpatchDict);
                    }
                    else
                    {
                        Info<< "            -->" << patchName << " already exists "
                            << "in DU field. Skipping passItem "
                            << passItems()[itemI].name() << " entry" << endl;
                    }
                }
            }
        }

        // TO DO: take symmetries in account and read friction coefficient
        // instead of hardcoding it to 0.16
        if
        (
            rollSetup == "twoHighRolls" ||
            rollSetup == "sideRolls" ||
            rollSetup == "edgeRolls"
        )
        {
            if
            (
                dataContainer().returnKey<Switch>
                (
                    "rollerToRollerContact", passNr_
                )
            )
            {
                Info<< "        Adding roller-to-roller contact boundary "
                    << "condition" << endl;

                // Create boundary conditions container
                boundaryConditionsContainer BCCont(boundaryField, &dataContainer());

                BCCont.solidTraction(rollerToAirPatchName);

                if (rollSetup == "twohighRolls" || rollSetup == "edgeRolls")
                {
                    BCCont.solidContactMaster
                    (
                        topRollFrontRTRPatchName,
                        bottomRollFrontRTRPatchName,
                        0.16,
                        2e-06,
                        0.5,
                        20.0
                     );
                    BCCont.solidContactSlave
                    (
                        bottomRollFrontRTRPatchName,
                        topRollFrontRTRPatchName
                    );
                    BCCont.solidContactMaster
                    (
                        topRollBackRTRPatchName,
                        bottomRollBackRTRPatchName,
                        0.16,
                        2e-06,
                        0.5,
                        20.0
                     );
                    BCCont.solidContactSlave
                    (
                        bottomRollBackRTRPatchName,
                        topRollBackRTRPatchName
                    );
                }
                if (rollSetup == "sideRolls")
                {
                    BCCont.solidContactMaster
                    (
                        topRollFrontRTRPatchName,
                        rightRollFrontRTRPatchName,
                        0.16,
                        2e-06,
                        0.5,
                        20.0
                    );
                    BCCont.solidContactSlave
                    (
                        rightRollFrontRTRPatchName,
                        topRollFrontRTRPatchName
                    );
                    BCCont.solidContactMaster
                    (
                        bottomRollFrontRTRPatchName,
                        rightRollBackRTRPatchName,
                        0.16,
                        2e-06,
                        0.5,
                        20.0
                    );
                    BCCont.solidContactSlave
                    (
                        rightRollBackRTRPatchName,
                        bottomRollFrontRTRPatchName
                    );
                    BCCont.solidContactMaster
                    (
                        bottomRollBackRTRPatchName,
                        leftRollBackRTRPatchName,
                        0.16,
                        2e-06,
                        0.5,
                        20.0
                    );
                    BCCont.solidContactSlave
                    (
                        leftRollBackRTRPatchName,
                        bottomRollBackRTRPatchName
                    );
                    BCCont.solidContactMaster
                    (
                        topRollBackRTRPatchName,
                        leftRollFrontRTRPatchName,
                        0.16,
                        2e-06,
                        0.5,
                        20.0
                    );
                    BCCont.solidContactSlave
                    (
                        leftRollFrontRTRPatchName,
                        topRollBackRTRPatchName
                    );
                }
            }
        }

        // Update boundaryField subDict in the DU field
        DU.add("boundaryField", boundaryField);

        // Write out the DU field
        const_cast<word&>(DU.type()) = "volVectorField";
        DU.regIOobject::write();
        const_cast<word&>(DU.type()) = "dictionary";
    }


    // T field
    {
        // Creating T dictionary
        IOdictionary T
            (
                IOobject
                (
                    "T",
                    "0",
                    runTime(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                dictionary("T")
            );

        // Add settings to T field
        T.add("dimensions", "[0 0 0 1 0 0 0]");
        T.add
        (
            "internalField", "uniform 20"
//            word
//            (
//                "uniform " +
//                Foam::name(readScalar(dict().lookup("Tinf")))
//            )
        );

        // Create boundaryField subDict
        dictionary boundaryField;

        forAll(passItems(), itemI)
        {
            // Read passItem subMesh
            polyMesh subMesh
            (
                IOobject
                (
                    passItems()[itemI].name(),
                    runTime().timeName(),
                    runTime(),
                    IOobject::MUST_READ
                )
            );

            // Read boundary conditions for this passItem
            dictionary& itemIBCs = passItems()[itemI].boundaryConditionsT();

            // Add boundary condition for each patch on the subMesh boundary
            forAll(subMesh.boundaryMesh(), patchI)
            {
                const word& patchName = subMesh.boundaryMesh()[patchI].name();
                if (itemIBCs.isDict(patchName))
                {
                    const dictionary& TpatchDict = itemIBCs.subDict(patchName);
                    if (!boundaryField.found(patchName))
                    {
                        boundaryField.add(patchName, TpatchDict);
                    }
                    else
                    {
                        Info<< "            -->" << patchName << " already exists "
                            << "in T field. Skipping passItem "
                            << passItems()[itemI].name() << " entry" << endl;
                    }
                }
            }
        }

        if
        (
            rollSetup == "twoHighRolls" ||
            rollSetup == "sideRolls" ||
            rollSetup == "edgeRolls"
        )
        {
            if
            (
                dataContainer().returnKey<Switch>
                (
                    "rollerToRollerContact", passNr_
                )
            )
            {
                Info<< "        Adding roller-to-roller contact boundary condition"
                    << endl;

                // Create boundary conditions container
                boundaryConditionsContainer BCCont(boundaryField, &dataContainer());

                BCCont.thermalConvection(rollerToAirPatchName);

                if (rollSetup == "twohighRolls" || rollSetup == "edgeRolls")
                {
                    BCCont.thermalContactMaster
                    (
                        topRollFrontRTRPatchName,
                        bottomRollFrontRTRPatchName
                     );
                    BCCont.thermalContactSlave
                    (
                        bottomRollFrontRTRPatchName,
                        topRollFrontRTRPatchName
                    );
                    BCCont.thermalContactMaster
                    (
                        topRollBackRTRPatchName,
                        bottomRollBackRTRPatchName
                     );
                    BCCont.thermalContactSlave
                    (
                        bottomRollBackRTRPatchName,
                        topRollBackRTRPatchName
                    );
                }
                if (rollSetup == "sideRolls")
                {
                    BCCont.thermalContactMaster
                    (
                        topRollFrontRTRPatchName,
                        rightRollFrontRTRPatchName
                    );
                    BCCont.thermalContactSlave
                    (
                        rightRollFrontRTRPatchName,
                        topRollFrontRTRPatchName
                    );
                    BCCont.thermalContactMaster
                    (
                        bottomRollFrontRTRPatchName,
                        rightRollBackRTRPatchName
                    );
                    BCCont.thermalContactSlave
                    (
                        rightRollBackRTRPatchName,
                        bottomRollFrontRTRPatchName
                    );
                    BCCont.thermalContactMaster
                    (
                        bottomRollBackRTRPatchName,
                        leftRollBackRTRPatchName
                    );
                    BCCont.thermalContactSlave
                    (
                        leftRollBackRTRPatchName,
                        bottomRollBackRTRPatchName
                    );
                    BCCont.thermalContactMaster
                    (
                        topRollBackRTRPatchName,
                        leftRollFrontRTRPatchName
                    );
                    BCCont.thermalContactSlave
                    (
                        leftRollFrontRTRPatchName,
                        topRollBackRTRPatchName
                    );
                }
            }
        }

        // Update boundaryField subDict in the T field
        T.add("boundaryField", boundaryField);

        // Write out the T field
        const_cast<word&>(T.type()) = "volScalarField";
        T.regIOobject::write();
        const_cast<word&>(T.type()) = "dictionary";
    }
}


void Foam::processPass::mergePassItemFields(const label passItemID)
{
    Info<< "    Merge fields from: "
        << passItems()[passItemID].name() << endl;

    // Take a reference to the meshes
    const fvMesh& baseMesh = mesh();
    const fvMesh& passItemMesh = passItems()[passItemID].mesh(runTime());

    // Get a list of all fields in the base mesh
    IOobjectList baseObjects(baseMesh, baseMesh.time().timeName());

    // Get a list of all fields in the passItem region mesh
    IOobjectList passItemObjects
    (
        passItemMesh, passItemMesh.time().timeName()
    );

    // Merge volFields of passItemMesh into the baseMesh
    MergeVolFields<scalar>
    (
        baseObjects, passItemObjects, baseMesh, passItemMesh
    );

    MergeVolFields<vector>
    (
        baseObjects, passItemObjects, baseMesh, passItemMesh
    );

    MergeVolFields<tensor>
    (
        baseObjects, passItemObjects, baseMesh, passItemMesh
    );

    MergeVolFields<symmTensor>
    (
        baseObjects, passItemObjects, baseMesh, passItemMesh
    );

    MergeVolFields<diagTensor>
    (
        baseObjects, passItemObjects, baseMesh, passItemMesh
    );
}


void Foam::processPass::writeMaterialIndicatorField()
{
    chDir(runTime().path());

    dataContainer().runSystemCommand
    (
        "setMatFromCellZones",
        "log.setMatFromCellZones",
        "void Foam::processPass::writeMaterialIndicatorField()"
    );
}


void Foam::processPass::createMeshModifierZones()
{
    chDir(runTime().path());

    // Remove sets sub-directory from the base mesh
    if (isDir("constant/polyMesh/sets"))
    {
        rmDir("constant/polyMesh/sets");
    }

    // Take reference to patchNames dictionary
    const dictionary& patchNames =
        dataContainer().dataContainerDict().subDict("patchNames");

    // Create faceZones for topoSolver
    OFstream batchFileTopo("topo.setSet");
    batchFileTopo
        << "faceSet wireInletFaces new patchToFace "
        << word(patchNames.lookup("wireUpstreamPatchName")) << endl;
    batchFileTopo
        << "cellSet wireInletFacesMasterCells new zoneToCell wire" << endl;
    batchFileTopo
        << "faceSet wireOutletFaces new patchToFace "
        << word(patchNames.lookup("wireDownstreamPatchName")) << endl;
    batchFileTopo
        << "cellSet wireOutletFacesMasterCells new zoneToCell wire" << endl;

    dataContainer().runSystemCommand
    (
        "setSet -batch topo.setSet",
        "log.setSetCellZones",
        "void Foam::processPass::createMeshModifierZones()"
    );
}


void Foam::processPass::writeCreateWedgePatchDict()
{
    Info<< "    Writing createPatchDict to generate wedge patches" << endl;

    // Creating createPatchDict dictionary
    IOdictionary createPatchDict
    (
        IOobject
        (
            "createPatchDict",
            "system",
            runTime(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dictionary("createPatchDict")
    );

    // Add default settings
    createPatchDict.add("matchTolerance", "1e-6");
    createPatchDict.add("pointSync", "true");

    // Take reference to patchNames dictionary
    const dictionary& patchNames =
        dataContainer().dataContainerDict().subDict("patchNames");

    // Make frontWedge patchInfo dictionary
    word FP1 = "\"" + word(patchNames.lookup("wireFrontPatchName")) + "\"";
    word FP2 = "\"" + word(patchNames.lookup("dieFrontPatchName")) + "\"";
    dictionary frontWedge;
    frontWedge.add("name", patchNames.lookup("wireFrontPatchName"));
    frontWedge.add("constructFrom", "patches");
    frontWedge.add("patches", word( "(" + FP1 + FP2 + ")"));
    dictionary newPatchTypeDictF;
    newPatchTypeDictF.add("type", "wedge");
    frontWedge.add("dictionary", newPatchTypeDictF);

    // Make backWedge patchInfo dictionary
    word BP1 = "\"" + word(patchNames.lookup("wireBackPatchName")) + "\"";
    word BP2 = "\"" + word(patchNames.lookup("dieBackPatchName")) + "\"";
    dictionary backWedge;
    backWedge.add("name", patchNames.lookup("wireBackPatchName"));
    backWedge.add("constructFrom", "patches");
    backWedge.add("patches", word( "(" + BP1 + BP2 + ")"));
    dictionary newPatchTypeDict;
    newPatchTypeDict.add("type", "wedge");
    backWedge.add("dictionary", newPatchTypeDict);

    // Add patchInfo dictionaries to a list
    // Create empty list of patchInfos
    PtrList<entry> patchInfo;
    patchInfo.resize(2);

    patchInfo.set
    (
        0,
        new dictionaryEntry
        (
            "",
            createPatchDict,
            frontWedge
        )
    );

    patchInfo.set
    (
        1,
        new dictionaryEntry
        (
            "",
            createPatchDict,
            backWedge
        )
    );

    createPatchDict.add("patchInfo", patchInfo);

    // Write createPatchDict dict
    createPatchDict.regIOobject::write();

}


void Foam::processPass::addFunctionsObjectsToControlDict()
{
    // Read default dict in the dictFiles directory
    IOdictionary controlDict
    (
        IOobject
        (
            "controlDict",
            "system",
            runTime(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    // Add custom function objects tot he controlDict

    const dictionary wireInputDict =
        dict().subDict("passItems").subDict(wire().name());

    Switch rigidTool(0);

    if (wireInputDict.subDict("boundaryConditions").found("rigidTool"))
    {
        rigidTool = Switch
            (
                wireInputDict.subDict("boundaryConditions").lookup("rigidTool")
            );
    }

    if (!rigidTool)
    {
        // Create empty list of function objects
        PtrList<dictionaryEntry> functionObjects;

        // Count the total number of function objects
        int nFunctionObjects = 0;
        forAll(passItems(), itemI)
        {
            nFunctionObjects += passItems()[itemI].functionObjects().size();
        }

        // Info statements are put here to keep logical output to screen, i.e.,
        // first give the pass items (from passItems() function) and then
        // go to dictionaries for the processPass.
        Info<< nl << ">> WRITING DICTIONARIES" << endl;
        Info<< "    Writing controlDict" << endl;
        Info<< "        Adding functionObjects:" << endl;

        functionObjects.resize(nFunctionObjects);

        // Add function objects from each pass item
        int functionObjectI = 0;
        forAll(passItems(), itemI)
        {
            const List<dictionary> passItemFunctionObjects =
                passItems()[itemI].functionObjects();

            forAll(passItemFunctionObjects, funcI)
            {
                Info<< "            "
                    << word(passItems()[itemI].name() + "_" + Foam::name(funcI))
                    << endl;

                functionObjects.set
                (
                    functionObjectI,
                    new dictionaryEntry
                    (
                        word
                        (
                            passItems()[itemI].name() + "_" + Foam::name(funcI)
                        ),
                        controlDict,
                        passItemFunctionObjects[funcI]
                    )
                );

                functionObjectI++;
            }
        }

        // Add functionObjects to the controlDict
        controlDict.add("functions", functionObjects);
    }

    // Write controlDict
    controlDict.regIOobject::write();
}


void Foam::processPass::writeGnuplotScript
(
    word const patchName,
    word const dataType
)
{
    if (dataType == "force" || dataType == "power" || dataType == "torque")
    {
        const fileName scriptDir("history/gnuplotScripts");

        fileName dataFileName;
        fileName scriptFileName;
        word unit;
        if (dataType == "force")
        {
            dataFileName = fileName("solidForces" + patchName + ".dat");
            scriptFileName = "forces_" + patchName + ".gnuplot";
            unit = "kN";
        }
        if (dataType == "power")
        {
            dataFileName = "solidPower" + patchName + ".dat";
            scriptFileName = "power_" + patchName + ".gnuplot";
            unit = "kW";
        }
        if (dataType == "torque")
        {
            dataFileName = "solidTorque" + patchName + ".dat";
            scriptFileName = "torque_" + patchName + ".gnuplot";
            unit = "kNm";
        }

        Info<< "        Writing " << scriptFileName << " script" << endl;

        IFstream forceFile
        (
            runTime().path()/"history/0"/dataFileName
        );

        if (forceFile.opened())
        {
            // Create a gnuplot script for a graph of force versus time
            OFstream gnuplotScript
            (
                runTime().path()/scriptDir/scriptFileName
            );

            gnuplotScript
                << "set terminal pdf" << endl
                << "set output \"../PDFs/" << dataType << "_" << patchName
                << ".pdf\"" << endl
                << "set style line 80 lt 1 lw 2" << endl
                << "set style line 81 lt 0 lw 2" << endl
                << "set grid back linestyle 81" << endl
                << "set border 3 back linestyle 80" << endl
                << "set xtics nomirror" << endl
                << "set ytics nomirror" << endl
                << "set encoding iso_8859_2" << endl
                << "set title \"" << dataType << ": " << patchName << "\""
                << endl
                << "set xlabel \"Time [in s]\"" << endl
                << "set ylabel \"" << dataType << " [in " << unit << "]\""
                << endl
                << "plot \\" << endl;

            if (dataType == "force")
            {
                gnuplotScript
                    << "\"../0/" << word(dataFileName)
                    << "\" u 1:(0.001*$2) w lp t \"X component\", \\" << endl
                    << "\"\" u 1:(0.001*$3) w lp t \"Y component\", \\" << endl
                    << "\"\" u 1:(0.001*$4) w lp t \"Z component\"" << endl;
            }
            else if (dataType == "power")
            {
                gnuplotScript
                    << "\"../0/" << word(dataFileName)
                    << "\" u 1:(0.001*$2) w lp t \"External Power\"" << endl;
            }
            else if (dataType == "torque")
            {
                gnuplotScript
                    << "\"../0/" << word(dataFileName)
                    << "\" u 1:(0.001*$2) w lp t \"External Power\"" << endl;
            }

            // Run gnuplot to create the graph
            chDir(runTime().path()/scriptDir);

            dataContainer().runSystemCommand
            (
                "gnuplot " + word(scriptFileName),
                "log." + word(scriptFileName),
                "void Foam::processPass::writeGnuplotScript",
                false // error from system command is ignored
            );

            chDir(runTime().path());
        }
        else
        {
            Info<< "    Data file not found!" << endl;
        }
    }
    else
    {
        Info<< "    Unknown data type to process!" << endl;
    }
}


void Foam::processPass::overwritePassDefaults()
{
    // Overwrite with custom pass settings, if found
    if (!dict().found("overwritePassDefaults"))
    {
        Info<< "        No overwritePassDefaults dict found for custom "
            << "settings: using defaults."
            << endl;

        return;
    }

    Info<< "        Writing custom pass settings from the "
        << "overwritePassDefaults dict" << endl;

    // Read in raw list of entries
    List<token>& rawStream = dict().lookup("overwritePassDefaults");

    DynamicList< DynamicList<word> > address;
    label counter = -1;
    bool flag = true;
    forAll(rawStream, tokI)
    {
        const token& curToken = rawStream[tokI];

        if
        (
            curToken.pToken() == token::BEGIN_LIST
         || curToken.pToken() == token::END_LIST
        )
        {
            continue;
        }

        if (flag)
        {
            flag = false;
            address.append(DynamicList<word>());
            counter++;
        }

        if (curToken.isWord())
        {
            address[counter].append(curToken.wordToken());
        }
        else if (curToken.isNumber())
        {
            address[counter].append(Foam::name(curToken.number()));
        }
        else if (curToken.isPunctuation())
        {
            if (curToken.pToken() == token::END_STATEMENT)
            {
                flag = true;
            }
            else if (curToken.pToken() != token::DIVIDE)
            {
                FatalErrorIn(type())
                    << "Symbol " << curToken.pToken() << " is not recognised."
                    << " Only '/' and ';' are allowed!"
                    << abort(FatalError);
            }
        }
        else
        {
            FatalErrorIn(type())
                << "Current token not recognised! Please check the format of "
                << "entry number " << tokI << nl
                << "Example format: " << nl
                << "system/controlDict/deltaT   0.002;"
                << abort(FatalError);
        }
    }

    address.shrink();

    wordList value(address.size());

    forAll(address, i)
    {
        DynamicList<word>& curAddress = address[i];
        value[i] = curAddress.last();
        curAddress.resize(curAddress.size() - 1);
        curAddress.shrink();
    }

    forAll(address, addressI)
    {
        DynamicList<word>& curAddress = address[addressI];

        // Find file name and the file path

        fileName relativePathOfFile = "";
        fileName testRelativePathOfFile = curAddress[0];
        word nameOfFile = curAddress[0];

        label addressCompI = 0;

        while (isDir(fileName(runTime().path()/testRelativePathOfFile)))
        {
            addressCompI++;

            // Update path of file
            relativePathOfFile = testRelativePathOfFile;

            // Next component of the address
            const word& nextAddrComp = curAddress[addressCompI];

            // Update next path to check
            testRelativePathOfFile =
                fileName(testRelativePathOfFile/nextAddrComp);

            // Update potential file name
            nameOfFile = nextAddrComp;
        }

        // Open the file
        IOdictionary curDict
        (
            IOobject
            (
                nameOfFile,
                relativePathOfFile,
                runTime(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        );

        // Find subDict, if any
        dictionary* curSubDictPtr = &curDict;
        while (curSubDictPtr->isDict(curAddress[addressCompI + 1]))
        {
            curSubDictPtr =
                &(curSubDictPtr->subDict(curAddress[addressCompI + 1]));
            addressCompI++;
        }
        dictionary& curSubDict = *curSubDictPtr;

        // Set value to keyword
        const word& key = curAddress.last();
        Info<< "    setting " << key << " to " << value[addressI] << endl;
        curSubDict.set(key, value[addressI]);

        // Re-write the entire dict
        Info<< "    Re-writing dict to disk" << endl;
        curDict.regIOobject::write();
    }
}


void Foam::processPass::generateGnuplotPDFs()
{
    Info<< "    Gnuplot post-processing" << endl;

    // Create directory for gnuplot scripts
    chDir(runTime().path());
    const fileName scriptDir("history/gnuplotScripts");
    mkDir(runTime().path()/scriptDir);

    // Create directory for the PDFs
    mkDir(runTime().path()/"history/PDFs");

    const dictionary patchNames =
        dataContainer().dataContainerDict().subDict("patchNames");

    const word rollSetup = dataContainer().returnKey<word>
       (
           "rollSetup",
           passNr()
       );

    if (rollSetup == "none")
    {
        // rollSetup none assumes a drawing pass, generating the die gnuplot
        // graphs
        const word dieToHousingPatchName =
                patchNames.lookup("dieToHousingPatchName");
        writeGnuplotScript(dieToHousingPatchName, "force");

        const word dieToWirePatchName =
                patchNames.lookup("dieToHousingPatchName");
        writeGnuplotScript(dieToHousingPatchName, "force");
    }
    else
    {
        // rollSetup is assumed to be a setup of a roller pass, generating
        // the roller gnuplot graphs
        forAll(passItems_, passI)
        {
            if (passItems_[passI].type() == "roller")
            {
                const word rollerType =
                        word(passItems_[passI].dict().lookup("rollerType")) ;
                const word tempPatchName =
                        patchNames.lookup("rollerAxisPatchName");
                const word rollerAxisPatchName =
                        word(rollerType + tempPatchName);
                writeGnuplotScript(rollerAxisPatchName, "force");
                writeGnuplotScript(rollerAxisPatchName, "power");
                writeGnuplotScript(rollerAxisPatchName, "torque");
            }
        }
    }

    const word wireUpstreamPatchName =
            patchNames.lookup("wireUpstreamPatchName");
    writeGnuplotScript(wireUpstreamPatchName, "force");
    writeGnuplotScript(wireUpstreamPatchName, "power");

    const word wireDownstreamPatchName =
            patchNames.lookup("wireDownstreamPatchName");
    writeGnuplotScript(wireDownstreamPatchName, "force");
    writeGnuplotScript(wireDownstreamPatchName, "power");
}


void Foam::processPass::generateParaViewImages()
{
    Info<< "    ParaView post-processing" << endl;

    WarningIn("void Foam::processPass::generateParaViewImages()")
        << "Disabled generation of ParaView images" << endl;

    // // Change to the case directory
    // chDir(runTime().path());

    // // Lookup the paraview script
    // const fileName paraviewScript
    // (
    //     runTime().path()/"../.."/
    //         fileName("dictFiles/paraviewScripts/postProcess.py")
    // );

    // // Create a bash script to run the paraview script
    // // We could directly run the paraview command using the system command
    // // but sometimes there can be clashes between foam libraries and paraview
    // // libraries, so we will instead call it from a script
    // // We use pvpython to run paraview without the GUI, we could also use
    // // pvbatch
    // const fileName bashScriptName("runParaviewPostProcess.sh");
    // Info<< "    Writing bash script: " << word(bashScriptName) << endl;
    // {
    //     OFstream runParaviewScript
    //     (
    //         fileName(runTime().path()/bashScriptName)
    //     );

    //     runParaviewScript
    //         << "#!/bin/bash" << nl
    //         << "pvpython " << word(paraviewScript)
    //         << " >& log.pvpython" << endl;
    // }

    // // Change the permissions of the script to be executable
    // // Use OF command chMod which is equivalent to chmod 775
    // chMod(bashScriptName, 33277);

    // dataContainer().runSystemCommand
    // (
    //     "touch case.foam",
    //     "log.touchCaseDotFoam",
    //     "void Foam::processPass::generateParaViewImages()"
    // );

    // // Create directory for the paraview images
    // const fileName imageDir("history/paraviewImages");
    // Info<< "    Creating directory: " << word(imageDir) << endl;
    // mkDir(runTime().path()/imageDir);

    // // Run the script; gives an error at the moment, so not running the
    // // system
    // // command via the dataContainer function
    // Info<< "    Running script: " << word(bashScriptName)
    //     << " >& log." << word(bashScriptName) << endl;
    // system(word("./" + bashScriptName));

    // Info<< endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processPass::processPass
(
    const processLine& parent,      // parent process line
    const word& passName,           // processPass name
    const dictionary& dict,         // subDict describing processPass
    const int& passNr,              // pass number
    dataHandler* dataContainerPtr   // pocess line data container
)
:
    parentPtr_(&parent),
    name_(passName),
    dict_(dict),
    controlDict_(),
    runTimePtr_(NULL),
    meshPtr_(NULL),
    passItems_(0),
    passNr_(passNr),
    dataContainerPtr_(dataContainerPtr)
{
    // Add options to controlDict to keep time happy
    controlDict_.add("deltaT", 1.0);
    controlDict_.add("writeInterval", 1.0);

    // Add rollSetup setting to dataContainer
    const word rollSetup = word(dict_.lookup("rollSetup"));
    dataContainer().addKey<word>("rollSetup", rollSetup, passNr);

    dataContainer().addKey<word>("passName", passName, passNr);

    if (rollSetup != "none")
    {
        const Switch singleRollerContact =
            dict_.subDict("rollSetupMeshOptions").lookupOrDefault<Switch>
            (
                "singleRollerContactPatch", "yes"
            );

        dataContainer().addKey<Switch>
        (
            "singleRollerContact",
            singleRollerContact,
            passNr
        );
    }
    else
    {
        dataContainer().addKey<Switch>
        (
            "singleRollerContact",
            "no",
            passNr
        );
    }

    dataContainer().addKey<Switch>
    (
        "contactPatchGenerated",
        "no",
        passNr
    );

    // Create the time object
    runTimePtr_.set
    (
        new Time
        (
            controlDict_,       // controlDict
            dataContainer().returnKey<fileName>("processLineDir"),
            passName,           // case directory i.e. case directory
            "system",           // system directory name
            "constant",         // constant directory name
            false               // functionObjects disabled
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::processPass::~processPass()
{}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

const Foam::processLine& Foam::processPass::parent() const
{
    return *parentPtr_;
}


Foam::Time& Foam::processPass::runTime()
{
    if (runTimePtr_.empty())
    {
        FatalErrorIn("Foam::Time& Foam::processPass::runTime()")
            << "pointer not set" << abort(FatalError);
    }

    return runTimePtr_();
}


const Foam::Time& Foam::processPass::runTime() const
{
    if (runTimePtr_.empty())
    {
        FatalErrorIn("const Foam::Time& Foam::processPass::runTime() const")
            << "pointer not set" << abort(FatalError);
    }

    return runTimePtr_();
}


Foam::fvMesh& Foam::processPass::mesh()
{
    if (meshPtr_.empty())
    {
        readMesh();
    }

    return meshPtr_();
}


const Foam::fvMesh& Foam::processPass::mesh() const
{
    if (meshPtr_.empty())
    {
        readMesh();
    }

    return meshPtr_();
}


const Foam::passItem& Foam::processPass::passItemObject
(
    const word& passItemName
) const
{
    label passItemID = -1;

    forAll(passItems(), itemI)
    {
        if (passItems()[itemI].name() == passItemName)
        {
            if (passItemID != -1)
            {
                FatalErrorIn
                (
                    "const Foam::passItem& Foam::processPass::passItemObject\n"
                    "(\n"
                    "    const word& passItemName\n"
                    ")"
                )   << "There is more than 1 passItem with name '"
                    << passItemName << "': there should only be 1!"
                    << abort(FatalError);
            }

            passItemID = itemI;
        }
    }

    if (passItemID == -1)
    {
        FatalErrorIn
        (
            "const Foam::passItem& Foam::processPass::passItemObject\n"
            "(\n"
            "    const word& passItemName\n"
            ")"
        )   << "Did not find passItem with name '" << passItemName << "'"
            << abort(FatalError);
    }

    return passItems()[passItemID];
}


Foam::wirePassItem& Foam::processPass::wire()
{
    label wireID = -1;

    forAll(passItems(), itemI)
    {
        if (passItems()[itemI].type() == wirePassItem::typeName)
        {
            if (wireID != -1)
            {
                FatalErrorIn("Foam::wirePassItem& Foam::processPass::wire()")
                << "There is more than 1 wire defined in this pass: there "
                << "should only be 1!" << abort(FatalError);
            }

            wireID = itemI;
        }
    }

    if (wireID == -1)
    {
        FatalErrorIn("Foam::wirePassItem& Foam::processPass::wire()")
            << "wire passItem object not found: a wire should be defined in "
            << "each pass" << abort(FatalError);
    }

    return refCast<wirePassItem>(passItems()[wireID]);
}


void Foam::processPass::setup()
{
    Info<< nl << nl << "PASS: " << name() << endl;
    Info<< "--------------------" << endl;

    // Take a reference to the time object for convenience
    Time& runTime = processPass::runTime();


    // 1. Make directories


    // Make case directory
    mkDir(runTime.rootPath()/name_);

    // Make case constant directory
    mkDir(runTime.rootPath()/name_/"constant");

    // Make case ssytem directory
    mkDir(runTime.rootPath()/name_/"system");

    // Make case 0 directory
    mkDir(runTime.rootPath()/name_/"0");


    // 2. Copy in all dictionaries


    // Write IOdictionary for writing

    // Create default empty case dictionaries and files
    chDir(runTime.path());
    const fileName dictFiles = "../../dictFiles/defaultDicts";
    WarningIn("void Foam::processPass::setup()")
        << "Why are we using runSystemCommand here for a cp?" << endl;
    dataContainer().runSystemCommand
    (
        "cp " + dictFiles + "/controlDict system/",
        "log.cp_controlDict",
        "void Foam::processPass::setup()"
    );
    dataContainer().runSystemCommand
    (
        "cp " + dictFiles + "/fvSchemes system/",
        "log.cp_fvSchemes",
        "void Foam::processPass::setup()"
    );
    dataContainer().runSystemCommand
    (
        "cp " + dictFiles + "/fvSolution system/",
        "log.cp_fvSolution",
        "void Foam::processPass::setup()"
    );
    dataContainer().runSystemCommand
    (
        "cp " + dictFiles + "/dynamicMeshDict constant/",
        "log.cp_dynamicMeshDict",
        "void Foam::processPass::setup()"
    );
    dataContainer().runSystemCommand
    (
        "cp " + dictFiles + "/decomposeParDict system/",
        "log.cp_decomposeParDict",
        "void Foam::processPass::setup()"
    );
    dataContainer().runSystemCommand
    (
        "cp " + dictFiles + "/meshDict system/",
        "log.cp_meshDict",
        "void Foam::processPass::setup()"
    );



    // 3. Modify dicts in the case
    addFunctionsObjectsToControlDict();
    // Overwrite pass defaults in now called just before the solver so all case
    // files exist
    //overwritePassDefaults();


    Info<< ">> SETUP PASSITEMS" << endl;


    // 4. Optional: create roll setup meshes (dict and creation)


    // First check if rollSetup is specified, if it is then we will create all
    // the roller mesh simulataneously, and then split the overall roller meshes
    // into the individual roller passItem meshes
    // The when the roller passItems try setup their mesh below, they will find
    // that it already exists and need to do nothing.
    const word rollSetup =
        dataContainer().returnKey<word>("rollSetup", passNr());

    const Switch singleRollerContact =
        dataContainer().returnKey<Switch>("singleRollerContact", passNr());

    if (rollSetup != "none" && singleRollerContact)
    {
        Info<< "    rollSetup: " << rollSetup << endl;

        // First, generate meshDict in system
        generateRollSetupMeshDict(rollSetup);

        // Second, run generate3DRollerMesh
        chDir(runTime.path());

        dataContainer().runSystemCommand
        (
            "export OMP_NUM_THREADS=1; "
            "generate3DRollerMesh",
            "log.generate3DRollerMesh.rollMeshSetup",
            "void Foam::processPass::setup()"
        );

        // Get a reference to current polyMesh. Note this does not include a
        // wire mesh yet. Hence, the meshPtr_ needs to be cleared later to let
        // it reference to the complete processPass mesh.
        const fvMesh& meshRef = mesh();

        if(meshRef.cellZones().size()>1)
        {
            dataContainer().runSystemCommand
            (
                "splitMeshRegionsLastTimeStep -overwrite",
                "log.splitMeshRegionsLastTimeStep.rollMeshSetup",
                "void Foam::processPass::setup()"
            );


            // Move all region meshes from 0 to the constant directory
            Info
            << "            Moving region meshes from 0 "
            << "to the constant directory"
            << endl;
            forAll(passItems(), passItemI)
            {
                if
                (
                    word(passItems()[passItemI].dict().lookup("type")) != "wire"
                )
                {
                    system
                    (
                        "mv 0/" +
                        passItems()[passItemI].name() +
                        " constant/ >& log.mvRollSetupRegion" +
                        passItems()[passItemI].name()
                    );
                }
            }
        }
        else if (meshRef.cellZones().size() == 1)
        {
            Info<< "            One cellZone found: "
                << word(meshRef.cellZones().names()[0]) << endl;

            mkDir("constant"/word(meshRef.cellZones().names()[0]));
            mv
            (
                "constant/polyMesh",
                "constant"/word(meshRef.cellZones().names()[0])/"polyMesh"
            );
        }
        else
        {
            FatalErrorIn
            (
                "void Foam::processPass::setup()"
            )   << nl << "No cellZones found!"
                << abort(FatalError);
        }

        // Finally, clean up: remove meshDict and remove global mesh
        Info<< "            Cleaning up constant/polyMesh and system/meshDict"
            << endl;
        chDir(runTime.path());
        mv
        (
            "system/meshDict",
            "system/meshDict.rollSetup"
        );
        rmDir("constant/polyMesh");

        // Clearing meshPtr_ so next time readMesh is called it will read the
        // current mesh
        meshPtr_.clear();
    }


    // 5.

    // Tell every pass item that is not a wire, what the wire length is
    // TODO: we may be able to do this with the parent functionality
    // And then setup each pass item e.g. passItems or wire
    label nrOfPassItems = 0;
    forAll(passItems(), passItemI)
    {
        if (word(passItems()[passItemI].dict().lookup("type")) != "wire" )
        {
            passItems()[passItemI].setWireLength
            (
                wire().passItem::wireLength()
            );
            passItems()[passItemI].setWireAxialResolution
            (
                wire().passItem::wireAxialResolution()
            );
            nrOfPassItems ++;
        }

        // Call setup for each passItem
        // Note: each roller/die now knows the wire length from above
        passItems()[passItemI].setup(runTime);
    }


    // 6. Check for error in defining single roller
    // Todo: remove this check as we will not allow single roller


    if
    (
        rollSetup == "singleRoll" &&
        nrOfPassItems > 1
    )
    {
        Info<< nl << nl << nrOfPassItems << nl << endl;
        FatalErrorIn
        (
            "void Foam::processPass::setup()"
        )   << nl << "More than 1 roller specified for rollSetup singleRoll"
            << abort(FatalError);
    }


    // 7. Pass shadowPatch names between wire/die/roller
    // Todo: we need a more general way to deal with inter-passItem and
    // inter-pass communication

    // Loop through all passItems and add all contactPatchNames to the
    // shadow patch name list of the wire passItem
    forAll(passItems(), passItemI)
    {
        if (passItems()[passItemI].type() != "wire" )
        {
            wire().addShadowPatchName
            (
                passItems()[passItemI].contactPatch()
            );

            wire().addFrictionCoefficientAxial
            (
                passItems()[passItemI].frictionCoefficientAxial()
            );
        }
    }
}


void Foam::processPass::run()
{
    // Take a reference to the time object for convenience
    Time& runTime = processPass::runTime();

    chDir(runTime.path());


    // 8. Clean-up files in the new case


    rm("0/cellToRegion");
    rm("0/cellToRegion.gz");


    // 9. Position each pass item

    Info<< "    Positioning pass items" << endl;
    // Position the pass items before merging
    const fvMesh& wireMesh = wire().mesh(runTime);
    forAll(passItems(), itemI)
    {
        Info<< "        Positioning " << passItems()[itemI].name() << endl;
        passItems()[itemI].positionMesh(runTime, wireMesh);
    }


    // 10. Merge pass item meshes
    mergePassItemMeshes();


    // 11. Merge pass item material properties dicts
    mergePassItemMaterialProperties();


    // 12. Merge pass item boundary condition files
    mergePassItemBoundaryConditions();


    // 13. Merge pass item fields, if found

    // Check if the pass items have fields; if found, then merge them with the
    // base mesh fields (or create them if not found in the base mesh)
    forAll(passItems(), itemI)
    {
        // Only merge fields from a wire pass item
        if (passItems()[itemI].type() == wirePassItem::typeName) //// not needed
        {
            mergePassItemFields(itemI);
        }
    }

    // 14. Write materials indicator field

    Info<< "    Final preparation" << endl;
    // Create the material indicator field
    writeMaterialIndicatorField();


    // 15. Write zones

    // Create the zones for the topological mesh modifiers
    // Not needed anymore. The required face and cellZones will be written by
    // the wire mesh generator
    createMeshModifierZones();

    // setsToZones
    dataContainer().runSystemCommand
    (
        "setsToZones",
        "log.setsToZonesCellZones",
        "void Foam::processPass::run()"
    );


    // 16. Renumber mesh for efficiency

    // renumberMesh
    dataContainer().runSystemCommand
    (
        "renumberMesh -overwrite",
        "log.renumberMesh",
        "void Foam::processPass::run()"
    );

    // 16.5: call overwrite pass defaults a seccond time to set e.g. settigns in
    // the boundary conditions etc.
    overwritePassDefaults();

    // 17. Run the solver

    Info<< nl << ">> CASE RUNNING " << name() << endl;

    chDir(runTime.path());

    // Read the number of CPU cores to run the solver
    const int nCores = dict().lookupOrDefault<int>("nCores", 1);

    const word solverName =
        dict().lookupOrDefault<word>("solverName", "plasticNonLinSolidFoam");

    if (nCores == 1)
    {
        // solver
        dataContainer().runSystemCommand
        (
            solverName,
            "log."+ solverName,
            "void Foam::processPass::run()"
        );
    }
    else
    {
        // Overwrite nCores in the decomposeParDict
        IOdictionary decomposeParDict
        (
            IOobject
            (
                "decomposeParDict",
                "system",
                runTime,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        );
        Info<< "            Setting nCores to " << nCores << endl;
        decomposeParDict.set("numberOfSubdomains", nCores);
        decomposeParDict.regIOobject::write();

        // decomposePar
        dataContainer().runSystemCommand
        (
            "decomposePar -cellDist",
            "log.decomposePar",
            "void Foam::processPass::run()"
        );

        // solver
        dataContainer().runSystemCommand
        (
            word
            (
                "mpirun -np " + Foam::name(nCores)
                + " " + solverName + " -parallel"
            ),
            "log.plasticNonLinSolidFoam",
            "void Foam::processPass::run()"
        );

        // reconstructPar
        dataContainer().runSystemCommand
        (
            "reconstructParMeshZones",
            "log.reconstructParMeshZones",
            "void Foam::processPass::run()"
        );
    }


    // 18. preparing the wire for the next pass

    // prepare wire for next pass
    Info<< nl << "    PREPARE WIRE FOR NEXT PASS" << endl;

    chDir(runTime.path());

    // splitMeshRegionsLastTimeStep
    dataContainer().runSystemCommand
    (
        "splitMeshRegionsLastTimeStep",
        "log.splitMeshRegionsLastTimeStep",
        "void Foam::processPass::run()"
    );

    const instantList& times = runTime.times();
    const label latestTimeID = times.size() -1;
    runTime.setTime(times[latestTimeID], latestTimeID);

    if (!isDir(runTime.path()/runTime.timeName()/wire().name()))
    {
        cp
        (
            runTime.path()/runTime.timeName(),
            runTime.path()/wire().name()
        );

        mv
        (
            runTime.path()/wire().name(),
            runTime.path()/runTime.timeName()/wire().name()
        );
    }

    if (wire().meshInputDict().found("newWireLength"))
    {
        // Lookup the new wire length
        const scalar newWireLength =
            readScalar
            (
                wire().meshInputDict().lookup("newWireLength")
            );

        // Smooth the wire mesh and fields in the flow direction
        dataContainer().runSystemCommand
        (
            word
            (
                "smoothWireMeshFields " + Foam::name(newWireLength)
            ),
            "log.smoothWireMeshFields",
            "void Foam::processPass::run()"
        );
    }


    // In case of remeshAfterwards is set, re-mesh the wire and map fields for
    // next pass
    // THIS IS NOW CALLED IN THE PROCESSLINE.C
    //remeshWireFromPreviousPass();


    // 19: post processing

    // Perform post processing e.g. generating graphs
    postProcess();
}


void Foam::processPass::generateRollSetupMeshDict(const word& rollSetup)
{
    Info<< "        generateRollSetupMeshDict" << endl;

    // Creating meshDict dictionary
    //dictionary meshDict;
    // Read template mesh dict: we will then edit it below
    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            "system",
            runTime(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );


    // Get reference to process line symmetry for convenience.
    const word processLineSymmetry =
        dataContainer().returnKey<word>("processLineSymmetry");

    // Set symmetry
    meshDict.set("symmetryType", processLineSymmetry);

    // Read location of input data
    const fileName meshInputDirectory =
        dataContainer().returnKey<fileName>
        (
            "processLineRootDir"
        )/"dictFiles"/name();

    // Set roll setup variable in the meshDict
    meshDict.set("rollSetup", rollSetup);

    // Lookup rollSetupMeshOptions subDict and add to the meshDict

    // Add rollSetup mesh optioins
    const dictionary& rollSetupMeshOptions =
        dict_.subDict("rollSetupMeshOptions");

    scalar geometryTolerance = 0.0;
    {
        Info<< "            Adding general settings" << endl;

        // Add settings to meshDict Dict
        // meshDict.add
        // (
        //     "enforceGeometryConstraints",
        //     rollSetupMeshOptions.lookupOrDefault<int>
        //     (
        //         "enforceGeometryConstraints", 0
        //     )
        // );
        // meshDict.add
        // (
        //     "keepCellsIntersectingBoundary",
        //     rollSetupMeshOptions.lookupOrDefault<int>
        //     (
        //         "keepCellsIntersectingBoundary", 0
        //     )
        // );
        // meshDict.add
        // (
        //     "maxCellSize",
        //     rollSetupMeshOptions.lookupOrDefault<scalar>("maxCellSize", 0.001)
        // );
        // meshDict.add
        // (
        //     "minContactCellSize",
        //     rollSetupMeshOptions.lookupOrDefault<scalar>
        //     (
        //         "minContactCellSize", 10e-6
        //     )
        // );
        // meshDict.add
        // (
        //     "minNumFacesBetweenFeatures",
        //     rollSetupMeshOptions.lookupOrDefault<scalar>
        //     (
        //         "minNumFacesBetweenFeatures", 3
        //     )
        // );

        // Store geometry tolerance as we use it later
        geometryTolerance =
            rollSetupMeshOptions.lookupOrDefault<scalar>
            (
                "geometryTolerance", 5e-5
            );
        meshDict.set
        (
            "geometryTolerance", geometryTolerance
        );

        // meshDict.add
        // (
        //     "singleRollerContactPatch",
        //     rollSetupMeshOptions.lookupOrDefault<Switch>
        //     (
        //         "singleRollerContactPatch", "yes"
        //     )
        // );

        // Create boundaryLayers subDict
        // {
        //     dictionary boundaryLayers;
        //     boundaryLayers.add
        //     (
        //         "nLayers",
        //         rollSetupMeshOptions.lookupOrDefault<int>("nLayers", 1)
        //     );
        //     boundaryLayers.add
        //     (
        //         "optimiseLayer",
        //         rollSetupMeshOptions.lookupOrDefault<int>("optimiseLayer", 1)
        //     );
        //     boundaryLayers.add
        //     (
        //         "symmetryPlaneLayerTopology",
        //         rollSetupMeshOptions.lookupOrDefault<int>
        //         (
        //             "symmetryPlaneLayerTopology", 1
        //         )
        //     );
        //     meshDict.add("boundaryLayers", boundaryLayers);
        // }

        // Create wireMeshDict subDict
        // {
        //     dictionary wireMeshDict;
        //     // wireDiameter for roller meshing is used to position the rollers,
        //     // i.e.,
        //     // to indicate the initial roll gap. It is defaulted to 50mm as we
        //     // assume no wire with a diameter larger than 50mm will be processed
        //     wireMeshDict.add("wireDiameter", 0.05);
        //     wireMeshDict.add
        //     (
        //         "wireLength",
        //         readScalar(rollSetupMeshOptions.lookup("wireLength"))
        //     );
        //     meshDict.add("wireMeshDict", wireMeshDict);
        //     meshDict.add("symmetryType", processLineSymmetry);
        //     meshDict.add("rollSetup", rollSetup);
        // }

        // meshDict.add
        // (
        //     "extraContactAngle",
        //     rollSetupMeshOptions.lookupOrDefault<scalar>
        //             (
        //                 "extraContactAngle", 20
        //             )
        // );

        dictionary patchNames =
            dataContainer().dataContainerDict().subDict("patchNames");

//        meshDict.add("centreInAxialDirection", 1);

        meshDict.set("rollingMillPatchNames", patchNames);
    }


    // // Create, set and add rollSetupDict subDict
    // {
    //     Info<< "            Adding settings for each roller" << endl;

    //     // Create rollSetup dict e.g. twoHighRollDict, etc.
    //     dictionary rollSetupDict;
    // Read and edit rollSetupDict
    const word rollSetupDictName = rollSetup + "Dict";
    dictionary& rollSetupDict = meshDict.subDict(rollSetupDictName);

    //     // Switch indicating of roller-to-roller patches are generated
    //     Switch rollerToRollerPatchesCreated("no");
    //     // Check if vertical or horizontal roll gap is specified as
    //     // this implies that rollerToRoller contact patches will be
    //     // generated.
    //     if (rollSetupMeshOptions.found("verticalRollGap"))
    //     {
    //         scalar verticalRollGap =
    //             readScalar
    //             (
    //                 rollSetupMeshOptions.lookup("verticalRollGap")
    //             );
    //         rollSetupDict.add("verticalRollGap", verticalRollGap);
    //         rollerToRollerPatchesCreated = "yes";
    //     }

    //     if (rollSetupMeshOptions.found("horizontalRollGap"))
    //     {
    //         scalar horizontalRollGap =
    //             readScalar
    //             (
    //                 rollSetupMeshOptions.lookup("horizontalRollGap")
    //             );
    //         rollSetupDict.add("horizontalRollGap", horizontalRollGap);
    //         rollerToRollerPatchesCreated = "yes";
    //     }

    //     rollSetupDict.add
    //     (
    //         "gapOffset",
    //         rollSetupMeshOptions.lookupOrDefault<scalar>("gapOffset", 0.0)
    //     );
/*
        // Check if vertical or horizontal spindle distance is specified as
        // this implies that rollerToRoller contact patches will be
        // generated for the sideRolls roll setup.
        if (rollSetupMeshOptions.found("verticalSpindleDistance"))
        {
            scalar vertSpindDist =
                readScalar
                (
                    rollSetupMeshOptions.lookup("verticalSpindleDistance")
                );
            rollSetupDict.add("verticalSpindleDistance", vertSpindDist);
            rollerToRollerPatchesCreated = "yes";
        }

        if (rollSetupMeshOptions.found("horizontalSpindleDistance"))
        {
            scalar horSpindDist =
                readScalar
                (
                    rollSetupMeshOptions.lookup("horizontalSpindleDistance")
                );
            rollSetupDict.add("horizontalSpindleDistance", horSpindDist);
            rollerToRollerPatchesCreated = "yes";
        }
*/

    if (rollSetupMeshOptions.lookupOrDefault<Switch>("contact", Switch(true)))
        {
            rollSetupDict.add("writeRollProfiles", "1");
            dataContainer().addKey<Switch>
            (
                "rollerToRollerContact",
                "yes",
                passNr_
            );
        }
        else
        {
            dataContainer().addKey<Switch>
            (
                "rollerToRollerContact",
                "no",
                passNr_
            );
        }

        // Add contact length (length of refined mesh along the roller surface)
        // const scalar contactLength =
        // rollSetupMeshOptions.lookupOrDefault<scalar>("contactLength", 0.02);
        // rollSetupDict.add("contactLength", contactLength);
        // Info<< "                contactLength: " << contactLength << endl;

        // Loop through all passItems
        forAll(passItems(), passItemI)
        {
            // We assume all passItems are rollers except for the wire
            if (word(passItems()[passItemI].dict().lookup("type")) == "roller")
            {
                // Roller mesh input dict
                const dictionary& rollerInputDict =
                    passItems()[passItemI].dict().subDict("mesh");

                // Roller name e.g. topRoll, bottomRoll, etc.
                const word name = passItems()[passItemI].name();

                Info<< "                roller: " << name << endl;

                // Mesh dict for this roller to be added to rollSetupDict
                //dictionary rollerDict;
                // Read rollerDict
                dictionary& rollerDict = rollSetupDict.subDict(name);

                // Read and do basic check on roller diameters
                const scalar innerDiameter =
                    readScalar(rollerInputDict.lookup("innerDiameter"));
                const scalar outerDiameter =
                    readScalar(rollerInputDict.lookup("outerDiameter"));

                if (innerDiameter > 0.9*outerDiameter)
                {
                    FatalErrorIn
                    (
                        "void Foam::processPass::generateRollSetupMeshDict"
                        "(const word& rollSetup)"
                    )   << "Inner roll diameter (" << innerDiameter
                        << ") should be less "
                        << " than 0.9 times the outer roll diameter ("
                        << outerDiameter << "), " << " for passItem " << name
                        << ". Either specify a geometryFile or a dxfFile."
                        << abort(FatalError);
                }

                // Get roll profile
                fileName geometryFileName;
                fileName dxfFileName;
                bool geometryFileFound = false;
                bool dxfFileFound = false;
                if (rollerInputDict.found("geometryFile"))
                {
                    geometryFileName =
                        fileName(rollerInputDict.lookup("geometryFile"));
                    geometryFileFound = true;
                }
                if (rollerInputDict.found("dxfFile"))
                {
                    dxfFileName = fileName(rollerInputDict.lookup("dxfFile"));
                    dxfFileFound = true;
                }
                fileName geometryFile =
                    meshInputDirectory + "/" + geometryFileName;
                fileName dxfFile = meshInputDirectory + "/" + dxfFileName;

                if (!geometryFileFound & !dxfFileFound)
                {
                    Info<< "    Flat roller" << endl;
                    rollerDict.set("type", "flatRoller");
                    rollerDict.remove("dxfFile");
                }
                else if (dxfFileFound)
                {
                    rollerDict.set("dxfFile", dxfFile);
                    rollerDict.set("type", "dxfFileRoller");
                }
                else
                {
                    rollerDict.set("geometryFile", geometryFile);
                    rollerDict.set("type", "geometryFileRoller");
                }

//                if (!geometryFileFound & !dxfFileFound)
//                {
//                    Info<< "                "
//                        << "No profile specified. Flat roller created."
//                         << endl;
//                }
//                else if (dxfFileFound)
//                {
//                    rollerDict.set("dxfFile", dxfFile);
//                }
//                else
//                {
//                    rollerDict.set("geometryFile", geometryFile);
//                }

                rollerDict.set
                (
                    "rollWidth",
                    rollerInputDict.lookupOrDefault<scalar>
                    (
                        "rollWidth", 0.04
                    )
                );

                // Number of cells in circumferential direction, outside of
                // refined region which is in contact with the wire. Note that
                // the tolerance is kind of master meshing parameter.
                // rollerDict.set
                // (
                //     "circumferentialResolution",
                //     rollerInputDict.lookupOrDefault<int>
                //     (
                //         "circumferentialResolution", 20
                //     )
                // );

                // // Add option to read the radial resolution
                // if (rollerInputDict.found("radialResolution"))
                // {
                //     rollerDict.set
                //     (
                //         "radialResolution",
                //         readInt(rollerInputDict.lookup("radialResolution"))
                //     );
                // }

                // // Add option to read the axial resolution
                // if (rollerInputDict.found("axialResolution"))
                // {
                //     rollerDict.set
                //     (
                //         "axialResolution",
                //         readInt(rollerInputDict.lookup("axialResolution"))
                //     );
                // }

                // // Flip the roll about the z-axis
                // rollerDict.set
                // (
                //     "flipRoll",
                //     rollerInputDict.lookupOrDefault<bool>("flipRoll", 0)
                // );

                // // Grading of cells in circumferential direction, starting from
                // // the refined region which is in contact with the wire
                // rollerDict.set
                // (
                //     "circumGrading",
                //     rollerInputDict.lookupOrDefault<scalar>
                //     (
                //         "circumGrading", 1.2
                //     )
                // );
                // // Grading in radial roll direction. Cell at inner roll diameter
                // // will "radialGrading" times larger than cell near contact
                // // surface
                // rollerDict.set
                // (
                //     "radialGrading",
                //     rollerInputDict.lookupOrDefault<scalar>("radialGrading", 1)
                // );

                // Grading in axial roll direction. A cell near the front or
                // back face will be "axialGrading" times larger than cell in
                // the mid plane of the roll
                // rollerDict.set
                // (
                //     "axialGrading",
                //     rollerInputDict.lookupOrDefault<scalar>("axialGrading", 1)
                // );

                // Determine axial spacing of the wire mesh. The roller mesh in
                // cricumferential direction is adapted to it.
                scalar wireAxialSpacing =
                    wire().wireLength()/wire().wireAxialResolution();

                scalar circumScaling = 1.5*wireAxialSpacing/
                    (
                       0.5*outerDiameter
                      *acos(1 - geometryTolerance/(0.5*outerDiameter))
                    );

                // Scaling of cells in circumferential direction
                rollerDict.set
                (
                    "circumScaling",
                    rollerInputDict.lookupOrDefault<scalar>
                    (
                        "circumScaling", circumScaling
                    )
                );

                // Scaling of cells in radial direction
                // rollerDict.set
                // (
                //     "radialScaling",
                //     rollerInputDict.lookupOrDefault<scalar>
                //     (
                //         "radialScaling", 1.0
                //     )
                // );

                // // Scaling of cells in axial direction
                // rollerDict.set
                // (
                //     "axialScaling",
                //     rollerInputDict.lookupOrDefault<scalar>("axialScaling", 1.0)
                // );

                rollerDict.set
                (
                    "innerDiameter",
                    rollerInputDict.lookupOrDefault<scalar>
                    (
                        "innerDiameter", 0.2
                    )
                );

                rollerDict.set("outerDiameter", outerDiameter);

                // This is no longer used: contactLength is specified instead
                //rollerDict.set
                //(
                //    "contactAreaLength",
                //    readScalar(rollerInputDict.lookup("contactAreaLength"))
                //);

                // No need as rollerDict is directly a reference to the subDict
                //rollSetupDict.add(name, rollerDict);
            }
        }

        // No need as rollSetupDict is directly a reference to the meshDict
        // const word rollSetupDictName = rollSetup + "Dict";
        // meshDict.add(rollSetupDictName, rollSetupDict);
        //}



    // IOdictionary meshDictIO
    // (
    //     IOobject
    //     (
    //         "meshDict",
    //         "system",
    //         runTime(),
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     meshDict
    // );

    // Write meshDict
    Info<< "            Writing meshDict" << endl;
    meshDict.regIOobject::write();
}


Foam::Switch Foam::processPass::rollValid
(
    const word& rollSetup,
    const word& processLineSymm,
    const word& passItemRollType
) const
{
    Info<< "        Verifying " << passItemRollType
        << " validity for rollSetup '" << rollSetup
        << "'' and process line symmetry '" << processLineSymm << "': ";

    Switch validRoll(false);

    List<word> validRollers;
    validRollers.setSize(4);
    forAll (validRollers, rollerI)
    {
        validRollers[rollerI] = "none";
    }

    // TODO: implement checks to set Switch validRoll.
    if (processLineSymm == "none")
    {
        if (rollSetup == "twoHighRolls")
        {
            validRollers[0] = "topRoll";
            validRollers[1] = "bottomRoll";
        }
        else if (rollSetup == "edgeRolls")
        {
            validRollers[0] = "leftRoll";
            validRollers[1] = "rightRoll";
        }
        else if (rollSetup == "sideRolls")
        {
            validRollers[0] = "topRoll";
            validRollers[1] = "bottomRoll";
            validRollers[2] = "leftRoll";
            validRollers[3] = "rightRoll";
        }
    }
    else if (processLineSymm == "posY")
    {
        if (rollSetup == "twoHighRolls")
        {
            validRollers[0] = "topRoll";
        }
        else if (rollSetup == "edgeRolls")
        {
            validRollers[0] = "leftRoll";
            validRollers[1] = "rightRoll";
        }
        else if (rollSetup == "sideRolls")
        {
            validRollers[0] = "topRoll";
            validRollers[1] = "leftRoll";
            validRollers[2] = "rightRoll";
        }
    }
    else if (processLineSymm == "negY")
    {
        if (rollSetup == "twoHighRolls")
        {
            validRollers[0] = "bottomRoll";
        }
        else if (rollSetup == "edgeRolls")
        {
            validRollers[0] = "leftRoll";
            validRollers[1] = "rightRoll";
        }
        else if (rollSetup == "sideRolls")
        {
            validRollers[0] = "bottomRoll";
            validRollers[1] = "leftRoll";
            validRollers[2] = "rightRoll";
        }
    }
    else if (processLineSymm == "posZ")
    {
        if (rollSetup == "twoHighRolls")
        {
            validRollers[0] = "topRoll";
            validRollers[1] = "bottomRoll";
        }
        else if (rollSetup == "edgeRolls")
        {
            validRollers[0] = "rightRoll";
        }
        else if (rollSetup == "sideRolls")
        {
            validRollers[0] = "topRoll";
            validRollers[1] = "bottomRoll";
            validRollers[2] = "rightRoll";
        }
    }
    else if (processLineSymm == "negZ")
    {
        if (rollSetup == "twoHighRolls")
        {
            validRollers[0] = "topRoll";
            validRollers[1] = "bottomRoll";
        }
        else if (rollSetup == "edgeRolls")
        {
            validRollers[0] = "leftRoll";
        }
        else if (rollSetup == "sideRolls")
        {
            validRollers[0] = "topRoll";
            validRollers[1] = "bottomRoll";
            validRollers[2] = "leftRoll";
        }
    }
    else if (processLineSymm == "posYposZ")
    {
        if (rollSetup == "twoHighRolls")
        {
            validRollers[0] = "topRoll";
        }
        else if (rollSetup == "edgeRolls")
        {
            validRollers[0] = "rightRoll";
        }
        else if (rollSetup == "sideRolls")
        {
            validRollers[0] = "topRoll";
            validRollers[1] = "rightRoll";
        }
    }
    else if (processLineSymm == "negYposZ")
    {
        if (rollSetup == "twoHighRolls")
        {
            validRollers[0] = "bottomRoll";
        }
        else if (rollSetup == "edgeRolls")
        {
            validRollers[0] = "rightRoll";
        }
        else if (rollSetup == "sideRolls")
        {
            validRollers[0] = "bottomRoll";
            validRollers[1] = "rightRoll";
        }
    }
    else if (processLineSymm == "negYnegZ")
    {
        if (rollSetup == "twoHighRolls")
        {
            validRollers[0] = "bottomRoll";
        }
        else if (rollSetup == "edgeRolls")
        {
            validRollers[0] = "leftRoll";
        }
        else if (rollSetup == "sideRolls")
        {
            validRollers[0] = "bottomRoll";
            validRollers[1] = "leftRoll";
        }
    }
    else if (processLineSymm == "posYnegZ")
    {
        if (rollSetup == "twoHighRolls")
        {
            validRollers[0] = "topRoll";
        }
        else if (rollSetup == "edgeRolls")
        {
            validRollers[0] = "rightRoll";
        }
        else if (rollSetup == "sideRolls")
        {
            validRollers[0] = "topRoll";
            validRollers[1] = "rightRoll";
        }
    }
    else
    {
        FatalErrorIn
            (
                "Foam::Switch Foam::processPass::rollValid"
                 "(const word& rollSetup, const word& processLineSymm, "
                 "const word& passItemRollType)"
            )
            << nl << "Roll validity issue. Aborting."
            << abort(FatalError);
    }

    forAll (validRollers, rollerI)
    {
        if (word(validRollers[rollerI]) == passItemRollType)
        {
            validRoll = true;
        }
    }

    Info<< validRoll << endl;

    if (!validRoll)
    {
        Info<< "        Not creating " << passItemRollType << endl;
    }

    return validRoll;
}


void Foam::processPass::postProcess()
{
    Info<< nl << ">> POST PROCESSING" << endl;

    generateGnuplotPDFs();
    generateParaViewImages();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
