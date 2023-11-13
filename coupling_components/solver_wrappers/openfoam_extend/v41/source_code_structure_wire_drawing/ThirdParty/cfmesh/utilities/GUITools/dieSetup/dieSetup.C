/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    deleteSubMeshes

Description
    Die setup (Bekaert specs).

Author
    Mario Baburic, Creative Fields

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "clock.H"
#include "fileName.H"
#include "IOdictionary.H"
#include "foamTime.H"
#include "fvMeshSubset.H"
#include "polyMesh.H"
#include "dictionaryEntry.H"
#include "IFstream.H"
#include "IOobjectList.H"
#include "fileStat.H"
#include "helperFunctionsAux.H"

#include <cmath>

#include <cstdlib>

using namespace Foam;

void inspect(const scalar v, const word var, bool range, bool equality, bool positive, const scalar min = 0., const scalar max = 90., const scalar amount = 0.)
{

    // Check for min/max
    if ( range && min <= max)
    {
        if (v < min)
        {
            Info << "   Suggestion for variable " << var << ": Value given = " << v << " Suggested value >= " << min << endl;
        }
        if (v > max)
        {
            Info << "   Suggestion for variable " << var << ": Value given = " << v << " Suggested value <= " << max << "." << endl;
        }
    }

    // Check for equality
    if ( equality )
    {
        scalar tolerance = 1.e-6;
        scalar difference = v - amount;
        if ( difference < 0.)
            difference *= -1.;
        if ( difference > tolerance )
        {
            Info << "   Suggestion for variable " << var << ": Value given = " << v << " Suggested value = " << amount << "." << endl;

            if (positive && amount < 0.)
                Info << "     WARNING: Negative value suggested for a parameter that should be positive." << nl
                     << "              Please inspect the following equality: H = EB_L + EC_L + RC_L + BR_L + BL + XC_L" << endl;
        }
    }

    // Check for zero
    if (positive && v <= 0.)
    {
        Info << "   Suggestion for variable " << var << ": Value given = " << v << " Suggested value --> a positive number." << endl;
    }

}

void checkBekaertSpecs
(
    const scalar H,
    const scalar W,
    const scalar Con,
    const scalar TC,
    const scalar gamma2,
    const scalar BC,
    const scalar gamma1,
    const scalar CR,
    const scalar CRA,
    const scalar EBR,
    const scalar EBL,
    const scalar ECA,
    const scalar ECL,
    const scalar RCA,
    const scalar RCL,
    const scalar RCR,
    const scalar D,
    const scalar BL,
    const scalar BRA,
    const scalar BRL,
    const scalar BRR,
    const scalar XCA,
    const scalar XCL,
    const scalar XCR
)
{
    Info << nl << "Suggestions (if any) ..." << endl;

    // Checking
    inspect( H, "H", false, false, true);
    inspect( W, "W", false, false, true);
    inspect( Con,"Con", true, false, false);
    inspect( TC, "TC", true, false, false);
    inspect( gamma2, "gamma2", true, false, false);
    inspect( BC, "BC", true, false, false);
    inspect( gamma1, "gamma1", true, false, false);
    inspect( CR, "C_R", false, false, true);
    inspect( CRA, "CR_A", true, false, false);
    inspect( EBR, "EB_R", false, true, true, 0., 90., (EBL-CR*std::sin(CRA*M_PI/180.))/(2.*std::sin(0.5*(90.-CRA-ECA/2.)*M_PI/180.)*std::cos(0.5*(90.-CRA+ECA/2.)*M_PI/180.)));
    inspect( EBL, "EB_L", false, false, true);
    inspect( ECA, "EC_A", true, false, false);
    inspect( ECL, "EC_L", false, false, true);
    inspect( RCA, "RC_A", true, false, false);
    inspect( RCL, "RC_L", false, false, true);
    inspect( RCR, "RC_R", false, false, true);
    inspect( D, "D", false, false, true);
    inspect( BL, "BL", false, false, true);
    inspect( BRA, "BR_A", true, false, false);
    inspect( BRL, "BR_L", false, true, true, 0., 90., H-EBL-ECL-RCL-BL-XCL);
    inspect( BRR, "BR_R", false, false, true);
    inspect( XCA, "XC_A", true, false, false);
    inspect( XCL, "XC_L", false, false, true);
    inspect( XCR, "XC_R", false, false, true);

    Info << "Done" << endl;
}

void readBekaertSpecsFromFile
(
    scalar & H,
    scalar & W,
    scalar & Con,
    scalar & TC,
    scalar & gamma2,
    scalar & BC,
    scalar & gamma1,
    scalar & CR,
    scalar & CRA,
    scalar & EBR,
    scalar & EBL,
    scalar & ECA,
    scalar & ECL,
    scalar & RCA,
    scalar & RCL,
    scalar & RCR,
    scalar & D,
    scalar & BL,
    scalar & BRA,
    scalar & BRL,
    scalar & BRR,
    scalar & XCA,
    scalar & XCL,
    scalar & XCR,
    const word subDictionary,
    const fileName file,
    const fileName folder,
    bool abortOnError
)
{
    Info << nl << "Reading parameters from file ..." << endl;

    if(fileStat(cwd()/folder/file).isValid() == 0)
    {
        Info << "   WARNING: Cannot read from the file. The file is missing." << endl;
        Info << "   MISSING FILE: " << cwd()/folder/file << endl;
        Info << "Done" << endl;
        return;
    }

    Info << "   FILE: " << cwd()/folder/file << endl;
    Info << "   SUBDICTIONARY: " << subDictionary << endl;

    H = help::readScalarFromDict("H",subDictionary,abortOnError,file,folder);
    Info << "     H:      " << H << endl;

    W = help::readScalarFromDict("W",subDictionary,abortOnError,file,folder);
    Info << "     W:      " << W << endl;

    Con = help::readScalarFromDict("Con",subDictionary,abortOnError,file,folder);
    Info << "     Con:    " << Con << endl;

    TC = help::readScalarFromDict("TC",subDictionary,abortOnError,file,folder);
    Info << "     TC:     " << TC << endl;

    gamma2 = help::readScalarFromDict("gamma2",subDictionary,abortOnError,file,folder);
    Info << "     gamma2: " << gamma2 << endl;

    BC = help::readScalarFromDict("BC",subDictionary,abortOnError,file,folder);
    Info << "     BC:     " << BC << endl;

    gamma1 = help::readScalarFromDict("gamma1",subDictionary,abortOnError,file,folder);
    Info << "     gamma1: " << gamma1 << endl;

    CR = help::readScalarFromDict("C_R",subDictionary,abortOnError,file,folder);
    Info << "     C_R:    " << CR << endl;

    CRA = help::readScalarFromDict("CR_A",subDictionary,abortOnError,file,folder);
    Info << "     CR_A:   " << CRA << endl;

    EBR = help::readScalarFromDict("EB_R",subDictionary,abortOnError,file,folder);
    Info << "     EB_R:   " << EBR << endl;

    EBL = help::readScalarFromDict("EB_L",subDictionary,abortOnError,file,folder);
    Info << "     EB_L:   " << EBL << endl;

    ECA = help::readScalarFromDict("EC_A",subDictionary,abortOnError,file,folder);
    Info << "     EC_A:   " << ECA << endl;

    ECL = help::readScalarFromDict("EC_L",subDictionary,abortOnError,file,folder);
    Info << "     EC_L:   " << ECL << endl;

    RCA = help::readScalarFromDict("RC_A",subDictionary,abortOnError,file,folder);
    Info << "     RC_A:   " << RCA << endl;

    RCL = help::readScalarFromDict("RC_L",subDictionary,abortOnError,file,folder);
    Info << "     RC_L:   " << RCL << endl;

    RCR = help::readScalarFromDict("RC_R",subDictionary,abortOnError,file,folder);
    Info << "     RC_R:   " << RCR << endl;

    D = help::readScalarFromDict("D",subDictionary,abortOnError,file,folder);
    Info << "     D:      " << D << endl;

    BL = help::readScalarFromDict("BL",subDictionary,abortOnError,file,folder);
    Info << "     BL:     " << BL << endl;

    BRA = help::readScalarFromDict("BR_A",subDictionary,abortOnError,file,folder);
    Info << "     BR_A:   " << BRA << endl;

    BRL = help::readScalarFromDict("BR_L",subDictionary,abortOnError,file,folder);
    Info << "     BR_L:   " << BRL << endl;

    BRR = help::readScalarFromDict("BR_R",subDictionary,abortOnError,file,folder);
    Info << "     BR_R:   " << BRR << endl;

    XCA = help::readScalarFromDict("XC_A",subDictionary,abortOnError,file,folder);
    Info << "     XC_A:   " << XCA << endl;

    XCL = help::readScalarFromDict("XC_L",subDictionary,abortOnError,file,folder);
    Info << "     XC_L:   " << XCL << endl;

    XCR = help::readScalarFromDict("XC_R",subDictionary,abortOnError,file,folder);
    Info << "     XC_R:   " << XCR << endl;

    Info << "Done" << endl;
}

void saveBekaertSpecs
(
    const scalar H,
    const scalar W,
    const scalar Con,
    const scalar TC,
    const scalar gamma2,
    const scalar BC,
    const scalar gamma1,
    const scalar CR,
    const scalar CRA,
    const scalar EBR,
    const scalar EBL,
    const scalar ECA,
    const scalar ECL,
    const scalar RCA,
    const scalar RCL,
    const scalar RCR,
    const scalar D,
    const scalar BL,
    const scalar BRA,
    const scalar BRL,
    const scalar BRR,
    const scalar XCA,
    const scalar XCL,
    const scalar XCR,
    dictionary & dict
)
{
    Info << nl << "Saving parameters to meshDict ..." << endl;

    dict.set("H",H);
    dict.set("W",W);
    dict.set("Con",Con);
    dict.set("TC",TC);
    dict.set("gamma2",gamma2);
    dict.set("BC",BC);
    dict.set("gamma1",gamma1);
    dict.set("C_R",CR);
    dict.set("CR_A",CRA);
    dict.set("EB_R",EBR);
    dict.set("EB_L",EBL);
    dict.set("EC_A",ECA);
    dict.set("EC_L",ECL);
    dict.set("RC_A",RCA);
    dict.set("RC_L",RCL);
    dict.set("RC_R",RCR);
    dict.set("D",D);
    dict.set("BL",BL);
    dict.set("BR_A",BRA);
    dict.set("BR_L",BRL);
    dict.set("BR_R",BRR);
    dict.set("XC_A",XCA);
    dict.set("XC_L",XCL);
    dict.set("XC_R",XCR);

    Info << "Done" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Optional argument
    // The utility cannot be run in parallel
    argList::noParallel();
    argList::validOptions.insert("check", "");
    argList::validOptions.insert("file", "word");
    argList::validOptions.insert("folder", "word");
    argList::validOptions.insert("dict", "word");
    argList::validOptions.insert("abort", "");
    argList::validOptions.insert("save", "");
    argList args(argc, argv);

    // Flags
    bool isFile = args.optionFound("file");
    bool isFolder = args.optionFound("folder");
    bool isDict = args.optionFound("dict");
    bool isCheck = args.optionFound("check");
    bool isAbort = args.optionFound("abort");
    bool isSave = args.optionFound("save");

    // Check for file/dict entries
    fileName file="";
    fileName folder="";
    word subDictionary="";
    if ( isFile )
    {
        file = args.option("file");
        if ( isFolder )
            folder = args.option("folder");
        if ( isDict )
        {
            subDictionary = args.option("dict");
        }
        else
        {
            Info << "WARNING: The dictionary must be specified." << endl;
            Info << "The program stops without doing any actions." << endl;
            return 0;
        }
    }

    // Checking case folder
    fileName passDir = cwd();
    help::checkFile(passDir/"system","Missing system directory");
    help::checkFile(passDir/"system/controlDict","Missing controlDict file");
    help::checkFile(passDir/"system/meshDict","Missing meshDict file");
    help::checkFile(passDir/"constant","Missing constant directory");

    // Creating time object
    Info << nl << "Creating time object" << endl;
    Time runTime
    (
        Time::controlDictName,
        passDir.path(),
        passDir.name(),
        "system",
        "constant",
        false
    );

    // Checking pass
    IOdictionary controlDict
    (
        IOobject
        (
            "controlDict",
            "system",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    word passType = controlDict.lookup("application");
    if (passType != "drawingPass")
    {
        FatalErrorIn("dieSetup")
            << nl << "    >>>> Invalid (pass) type. Aborting. "
            << nl << "    Type: " << passType
            << abort(FatalError);
    }
    else
    {
        Info << nl << "Pass type: " << passType << endl;
    }

    // Open meshDict
    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            "system",
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    // Check for dieMeshDict
    if(!meshDict.found("dieMeshDict") || !meshDict.isDict("dieMeshDict"))
    {
        Info << nl << "Dictionary dieMeshDict not found in meshDict. Please check the meshDict for input validity." << endl;
        Info << "The program stops without doing any actions." << endl;
        return 0;
    }
    dictionary& dict = meshDict.subDict("dieMeshDict");

    // Check for die type
    word type = "-";
    if( dict.found("type") )
    {
        type = word(dict.lookup("type"));
        if(type != "bekaertSpecs")
        {
            Info << nl << "The die type is not a Bekaert die." << endl;
            Info << "The program stops without doing any actions." << endl;
            return 0;
        }
        else
        {
            Info << nl << "Die type: " << type << endl;
        }
    }
    else
    {
        Info << nl << "The die type not set in the meshDict. Please check the meshDict for input validity." << endl;
        Info << "The program stops without doing any actions." << endl;
        return 0;
    }

    // Die parameters according to Bekaert specification
    scalar H, W, Con, TC, gamma2, BC, gamma1, CR, CRA, EBR, EBL, ECA, ECL, RCA, RCL, RCR, D, BL, BRA, BRL, BRR, XCA, XCL, XCR;

    if ( isCheck && !isFile )
    {
        readBekaertSpecsFromFile
        (
            H,
            W,
            Con,
            TC,
            gamma2,
            BC,
            gamma1,
            CR,
            CRA,
            EBR,
            EBL,
            ECA,
            ECL,
            RCA,
            RCL,
            RCR,
            D,
            BL,
            BRA,
            BRL,
            BRR,
            XCA,
            XCL,
            XCR,
            "dieMeshDict",
            "meshDict",
            "system",
            isAbort
        );
        checkBekaertSpecs
        (
            H,
            W,
            Con,
            TC,
            gamma2,
            BC,
            gamma1,
            CR,
            CRA,
            EBR,
            EBL,
            ECA,
            ECL,
            RCA,
            RCL,
            RCR,
            D,
            BL,
            BRA,
            BRL,
            BRR,
            XCA,
            XCL,
            XCR
        );
    }

    if ( isFile )
    {
        readBekaertSpecsFromFile
        (
            H,
            W,
            Con,
            TC,
            gamma2,
            BC,
            gamma1,
            CR,
            CRA,
            EBR,
            EBL,
            ECA,
            ECL,
            RCA,
            RCL,
            RCR,
            D,
            BL,
            BRA,
            BRL,
            BRR,
            XCA,
            XCL,
            XCR,
            subDictionary,
            file,
            folder,
            isAbort
        );
        checkBekaertSpecs
        (
            H,
            W,
            Con,
            TC,
            gamma2,
            BC,
            gamma1,
            CR,
            CRA,
            EBR,
            EBL,
            ECA,
            ECL,
            RCA,
            RCL,
            RCR,
            D,
            BL,
            BRA,
            BRL,
            BRR,
            XCA,
            XCL,
            XCR
        );
        if ( isSave )
        {
            saveBekaertSpecs
            (
                H,
                W,
                Con,
                TC,
                gamma2,
                BC,
                gamma1,
                CR,
                CRA,
                EBR,
                EBL,
                ECA,
                ECL,
                RCA,
                RCL,
                RCR,
                D,
                BL,
                BRA,
                BRL,
                BRR,
                XCA,
                XCL,
                XCR,
                dict
            );
            meshDict.regIOobject::writeObject(meshDict.time().writeFormat(),IOstream::currentVersion,IOstream::compressionType());
        }
    }

    // Print out date and time
    Info << nl << nl << nl << "Die setup completed on " << clock::date() << " at "
         << clock::clockTime() << endl;

    Info << "END" << nl << endl;

    return 0;
}

// ************************************************************************* //
