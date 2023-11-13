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

#include "dataFriction.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "IFstream.H"
#include "interpolateXY.H"
#include "frictionContactModel.H"
#include "interpolationTable.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dataFriction, 0);
    addToRunTimeSelectionTable(frictionLaw, dataFriction, dictionary);



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dataFriction::calcFrictionCoeff
(
    const scalar entrSpeed,
    const scalar SRR,
    const scalar pressure,
    const scalar temperature
)
{
    // Put current parameters into list based on the defined priorities
    scalarList curParameters(3, 0.0);

    curParameters[SRRPriority_] = SRR;
    curParameters[pressurePriority_] = pressure;
    curParameters[temperaturePriority_] = temperature;

    // Interpolate friction coefficient from the data
    List<Tuple2<label, scalar> > primaryHighLowValues
    (
        highLowValues
        (
            curParameters[0],
            primaryParameters_,
            useNearest_[0]
        )
    );

    scalarList primaryFrictionCoeffs(primaryHighLowValues.size(), 0.0);
    scalarList primaryParameters(primaryHighLowValues.size(), 0.0);

    forAll(primaryHighLowValues, valueI)
    {
        label primaryIndex = primaryHighLowValues[valueI].first();
        primaryParameters[valueI] = primaryHighLowValues[valueI].second();

        List<Tuple2<label, scalar> > secondaryHighLowValues
        (
            highLowValues
            (
                curParameters[1],
                secondaryParameters_[primaryIndex],
                useNearest_[1]
            )
        );

        scalarList secondaryFrictionCoeffs(secondaryHighLowValues.size(), 0.0);
        scalarList secondaryParameters(secondaryHighLowValues.size(), 0.0);

        forAll(secondaryHighLowValues, valueJ)
        {
            label secondaryIndex =
                secondaryHighLowValues[valueJ].first();
            secondaryParameters[valueJ] =
                secondaryHighLowValues[valueJ].second();

            List<Tuple2<label, scalar> > tertiaryHighLowValues
            (
                highLowValues
                (
                    curParameters[2],
                    tertiaryParameters_[primaryIndex][secondaryIndex],
                    useNearest_[2]
                )
            );

            scalarList tertiaryFrictionCoeffs
            (
                tertiaryHighLowValues.size(),
                0.0
            );
            scalarList tertiaryParameters(tertiaryHighLowValues.size(), 0.0);

            forAll(tertiaryHighLowValues, valueK)
            {
                label tertiaryIndex = tertiaryHighLowValues[valueK].first();
                tertiaryParameters[valueK] =
                    tertiaryHighLowValues[valueK].second();
                label dataIndex =
                    dataIndices_[primaryIndex][secondaryIndex][tertiaryIndex];

                interpolationTable<scalar> dataInterpolationTable
                (
                    frictionCoeffList_[dataIndex],
                    interpolationTable<scalar>().wordToBoundsHandling("clamp"),
                    "dataInterpolationTable"
                );

                tertiaryFrictionCoeffs[valueK] =
                    dataInterpolationTable(entrSpeed);
            }

            if (tertiaryFrictionCoeffs.size() == 1)
            {
                secondaryFrictionCoeffs[valueJ] = tertiaryFrictionCoeffs[0];
            }
            else
            {
                secondaryFrictionCoeffs[valueJ] =
                    interpolateXY
                    (
                        curParameters[2],
                        scalarField(tertiaryParameters),
                        scalarField(tertiaryFrictionCoeffs)
                    );
            }
        }

        if (secondaryFrictionCoeffs.size() == 1)
        {
            primaryFrictionCoeffs[valueI] = secondaryFrictionCoeffs[0];
        }
        else
        {
            primaryFrictionCoeffs[valueI] =
                interpolateXY
                (
                    curParameters[1],
                    scalarField(secondaryParameters),
                    scalarField(secondaryFrictionCoeffs)
                );
        }
    }

    if (primaryFrictionCoeffs.size() == 1)
    {
        frictionCoeff_ = primaryFrictionCoeffs[0];
    }
    else
    {
        frictionCoeff_ =
            interpolateXY
            (
                curParameters[0],
                scalarField(primaryParameters),
                scalarField(primaryFrictionCoeffs)
            );
    }
}


Foam::List<Tuple2<label, scalar> > dataFriction::highLowValues
(
    const scalar value,
    const scalarList& valueList,
    const bool useNearest
)
{
    if (value <= valueList[0])
    {
        return List<Tuple2<label, scalar> >
        (
            1, Tuple2<label, scalar>(0, valueList[0])
        );
    }

    if (value >= valueList[valueList.size()-1])
    {
        return List<Tuple2<label, scalar> >
        (
            1,
            Tuple2<label, scalar>
            (
                valueList.size() - 1,
                valueList[valueList.size() - 1]
            )
        );
    }

    forAll(valueList, i)
    {
        if (value <= valueList[i+1])
        {
            // If we are interpolating using current property value
            if (!useNearest)
            {
                List<Tuple2<label, scalar> > highLowList
                (
                    2, Tuple2<label, scalar>(-1, 0.0)
                );

                highLowList[0] = Tuple2<label, scalar>(i, valueList[i]);
                highLowList[1] = Tuple2<label, scalar>(i + 1, valueList[i+1]);

                return highLowList;
            }

            if (mag(value - valueList[i]) <= mag(value - valueList[i+1]))
            {
                return List<Tuple2<label, scalar> >
                (
                    1, Tuple2<label, scalar>(i, valueList[i])
                );
            }
            else
            {
                return List<Tuple2<label, scalar> >
                (
                    1, Tuple2<label, scalar>(i + 1, valueList[i+1])
                );
            }
        }
    }

    FatalErrorIn
    (
        "dataFriction::highLowValues(const scalar value, "
        "const scalarList& valueList, const bool useNearest)"
    )   << "Something went terribly wrong."
        << abort(FatalError);

    return List<Tuple2<label, scalar> >
    (
        1, Tuple2<label, scalar>(-1, 0.0)
    );
}


void dataFriction::lookupTemperatureField()
{
    if (TPtr_)
    {
        FatalErrorIn("void dataFriction::lookupTemperatureField()")
            << "temperature field pointer already set" << abort(FatalError);
    }

    const volScalarField& TRef =
        mesh_.objectRegistry::lookupObject<volScalarField>("T");
    TPtr_ = &TRef;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
dataFriction::dataFriction
(
    const word& name,
    const frictionContactModel& fricModel,
    const dictionary& dict
)
:
    frictionLaw(name, fricModel, dict),
    mesh_(fricModel.mesh()),
    frictionLawDict_(dict.subDict("frictionLawDict")),
    frictionCoeff_(0.0),
    parameterList_(1, scalarList(3, 0.0)),
    primaryParameters_(),
    secondaryParameters_(),
    tertiaryParameters_(),
    dataIndices_(),
    frictionCoeffList_(1, List<Tuple2<scalar, scalar> >(1)),
    useNearest_(3, false),
    SRRPriority_(-1),
    pressurePriority_(-1),
    temperaturePriority_(-1),
    fixedSRR_(false, 0.0),
    fixedPressure_(false, 0.0),
    fixedTemperature_(false, 0.0),
    TPtr_(NULL)
{
    const fvMesh& mesh = mesh_;

    // Data properties (SRR, pressure, temperature) priority
    // Interpolation of the friction coefficient based on the data properties
    // is carried out in the order specified in this dictionary
    // Allowed values are: 0, 1, 2
    SRRPriority_ =
        frictionLawDict_.subDict("SRR").lookupOrDefault<label>("priority", 0);

    pressurePriority_ =
        frictionLawDict_.subDict("pressure").lookupOrDefault<label>
        (
            "priority", 1
        );

    temperaturePriority_ =
        frictionLawDict_.subDict("temperature").lookupOrDefault<label>
        (
            "priority", 2
        );

    // Check if all parameter priorities are equal to 0, 1 or 2 and that they
    // are all different (their sum must be equalt to 3)
    {
        labelList priorityList(3, 0);

        priorityList[0] = SRRPriority_;
        priorityList[1] = pressurePriority_;
        priorityList[2] = temperaturePriority_;

        label sum = 0;
        bool badPriorities = 0;

        forAll(priorityList, i)
        {
            if (priorityList[i] < 0 || priorityList[i] > 2)
            {
                badPriorities = 1;
            }

            sum += priorityList[i];
        }

        if (sum != 3)
        {
            badPriorities = 1;
        }

        if (badPriorities)
        {
            FatalErrorIn("dataFriction::dataFriction()")
                << "Bad property priorities. Priorities must be 0, 1 or 2, and "
                << "they have to be different. "
                << abort(FatalError);
        }
    }

    // Check if properties are set as fixed
    fixedSRR_.first() =
        Switch(frictionLawDict_.subDict("SRR").lookup("fixedValue"));
    if (fixedSRR_.first())
    {
        fixedSRR_.second() =
            readScalar(frictionLawDict_.subDict("SRR").lookup("value"));
    }

    fixedPressure_.first() =
        Switch(frictionLawDict_.subDict("pressure").lookup("fixedValue"));
    if (fixedPressure_.first())
    {
        fixedPressure_.second() =
            readScalar(frictionLawDict_.subDict("pressure").lookup("value"));
    }

    fixedTemperature_.first() =
        Switch(frictionLawDict_.subDict("temperature").lookup("fixedValue"));
    if (fixedTemperature_.first())
    {
        fixedTemperature_.second() =
            readScalar(frictionLawDict_.subDict("temperature").lookup("value"));
    }

    // Check if thermalStress is ON.
    Switch thermalStress
    (
        mesh.solutionDict().subDict("solidMechanics").lookupOrDefault<Switch>
        (
            "thermalStress",
            false
        )
    );

    if ((!thermalStress) && (!fixedTemperature_.first()))
    {
        FatalErrorIn
        (
            "dataFriction::dataFriction()"
        )   << "Thermal stress is OFF. Please, set temperature property as "
               "fixedValue and specify it."
            << abort(FatalError);
    }

    // If useNearest_ is OFF for a property, then friction coefficient for a
    // needed priority is linearly interpolated between two friction
    // coefficients with closest property values to the needed property value.
    // If switch is ON, then friction coefficient is equal to a friction
    // coefficient with property value closest to the needed property value
    useNearest_[SRRPriority_] =
        frictionLawDict_.subDict("SRR").lookupOrDefault<Switch>
        (
            "useNearest",
            false
        );
    useNearest_[pressurePriority_] =
        frictionLawDict_.subDict("pressure").lookupOrDefault<Switch>
        (
            "useNearest",
            false
        );
    useNearest_[temperaturePriority_] =
        frictionLawDict_.subDict("temperature").lookupOrDefault<Switch>
        (
            "useNearest",
            false
        );

    // Read interpolation data file names
    wordList dataFilenames
    (
        frictionLawDict_.lookup("dataFilenames")
    );

    frictionCoeffList_.resize(dataFilenames.size());
    parameterList_.resize(dataFilenames.size());

    // Create empty dynamic lists with one member used for initialization
    DynamicList<scalar> emptyDynamicPrimaryList;
    emptyDynamicPrimaryList.append(0.0);

    DynamicList<DynamicList<scalar> > emptyDynamicSecondaryList;
    emptyDynamicSecondaryList.append(emptyDynamicPrimaryList);

    DynamicList<DynamicList<DynamicList<scalar> > > emptyDynamicTertiaryList;
    emptyDynamicTertiaryList.append(emptyDynamicSecondaryList);

    // Creat dynamic parameter lists
    primaryParameters_ = emptyDynamicPrimaryList;
    secondaryParameters_ = emptyDynamicSecondaryList;
    tertiaryParameters_ = emptyDynamicTertiaryList;
    dataIndices_ = emptyDynamicTertiaryList;

    // Loop through all data files (which are dictionaries) and read experiment
    // parameters (temperature, pressure, SRR) and measured data (entrainment
    // speed VS friction coefficients).
    fileName constantDir(mesh.time().caseConstant());

    forAll(dataFilenames, nameI)
    {
        IOdictionary dataFile
        (
            IOobject
            (
                dataFilenames[nameI],
                constantDir/"entrainmentSpeedVsFrictionCoeff",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Read parameters
        scalar SRR(readScalar(dataFile.lookup("SRR")));
        scalar pressure(readScalar(dataFile.lookup("pressure")));
        scalar temperature(readScalar(dataFile.lookup("temperature")));

        scalarList parameterListI(3, 0.0);

        parameterListI[SRRPriority_] = SRR;
        parameterListI[pressurePriority_] = pressure;
        parameterListI[temperaturePriority_] = temperature;

        if (nameI == 0)
        {
            // Append primary parameter
            primaryParameters_[nameI] = parameterListI[0];
            // Append secondary parameter
            secondaryParameters_[nameI][nameI] = parameterListI[1];
            // Append tertiary parameter
            tertiaryParameters_[nameI][nameI][nameI] = parameterListI[2];
            // Append data index
            dataIndices_[nameI][nameI][nameI] = nameI;
        }
        else
        {
            bool found0 = 0;
            forAll(primaryParameters_, paramI)
            {
                if (primaryParameters_[paramI] == parameterListI[0])
                {
                    bool found1 = 0;
                    forAll(secondaryParameters_[paramI], paramJ)
                    {
                        if
                        (
                            secondaryParameters_[paramI][paramJ]
                         == parameterListI[1]
                        )
                        {
                            forAll(tertiaryParameters_[paramI][paramJ], paramK)
                            {
                                if
                                (
                                    tertiaryParameters_[paramI][paramJ][paramK]
                                 == parameterListI[2]
                                )
                                {
                                    FatalErrorIn("dataFriction::dataFriction()")
                                        << "Multiple friction data files with "
                                        << "the same parameters exist: " << nl
                                        << "SRR = " << SRR << nl
                                        << "Pressure = " << pressure << nl
                                        << "Temperature = " << temperature
                                        << abort(FatalError);
                                }
                            }

                            // Append tertiary parameter
                            tertiaryParameters_[paramI][paramJ].append
                            (
                                parameterListI[2]
                            );
                            // Append data index
                            dataIndices_[paramI][paramJ].append(nameI);

                            found1 = 1;
                            break;
                        }
                    }

                    if (!found1)
                    {
                        // Append secondary parameter
                        secondaryParameters_[paramI].append(parameterListI[1]);
                        label size = secondaryParameters_[paramI].size();
                        // Append tertiary parameter
                        tertiaryParameters_[paramI].append
                        (
                            emptyDynamicPrimaryList
                        );
                        tertiaryParameters_[paramI][size-1][0] =
                            parameterListI[2];
                        // Append data index
                        dataIndices_[paramI].append(emptyDynamicPrimaryList);
                        dataIndices_[paramI][size-1][0] = nameI;
                    }

                    found0 = 1;
                    break;
                }

            }

            if (!found0)
            {
                // Append primary parameter
                primaryParameters_.append(parameterListI[0]);
                label size = primaryParameters_.size();
                // Append secondary parameter
                secondaryParameters_.append(emptyDynamicPrimaryList);
                secondaryParameters_[size-1][0] = parameterListI[1];
                // Append tertiary parameter
                tertiaryParameters_.append(emptyDynamicSecondaryList);
                tertiaryParameters_[size-1][0][0] = parameterListI[2];
                // Append data index
                dataIndices_.append(emptyDynamicSecondaryList);
                dataIndices_[size-1][0][0] = nameI;
            }

        }
        // Read experimental measurements (entrainment speed VS friction coeff)
        List<Tuple2<scalar, scalar> > frictionCoeffs
        (
            dataFile.lookup("entrainmentSpeedVsFrictionCoeff")
        );
        frictionCoeffList_[nameI].resize(frictionCoeffs.size());
        frictionCoeffList_[nameI] = frictionCoeffs;
    }

    // Sort dynamic lists

    // Create copy of dynamic lists
    DynamicList<scalar> primaryParametersOld(primaryParameters_);
    DynamicList<DynamicList<scalar> > secondaryParametersOld
    (
        secondaryParameters_
    );
    DynamicList<DynamicList<DynamicList<scalar> > > tertiaryParametersOld
    (
        tertiaryParameters_
    );
    DynamicList<DynamicList<DynamicList<scalar> > > dataIndicesOld
    (
        dataIndices_
    );

    // Sort tertiary parameters
    forAll(primaryParameters_, paramI)
    {
        forAll(secondaryParameters_[paramI], paramJ)
        {
            SortableList<scalar> sortedListTertiary
            (
                tertiaryParameters_[paramI][paramJ]
            );
            tertiaryParameters_[paramI][paramJ] = sortedListTertiary;
            const labelList originalIndices(sortedListTertiary.indices());

            forAll(originalIndices, paramK)
            {
                label originalIndex = originalIndices[paramK];

                dataIndices_[paramI][paramJ][paramK] =
                    dataIndicesOld[paramI][paramJ][originalIndex];
            }
        }
    }

    tertiaryParametersOld = tertiaryParameters_;
    dataIndicesOld = dataIndices_;

    // Sort secundary parameters
    forAll(primaryParameters_, paramI)
    {
        SortableList<scalar> sortedListSecundary(secondaryParameters_[paramI]);
        secondaryParameters_[paramI] = sortedListSecundary;
        const labelList originalIndices(sortedListSecundary.indices());

        forAll(originalIndices, paramJ)
        {
            label originalIndex = originalIndices[paramJ];

            tertiaryParameters_[paramI][paramJ] =
                tertiaryParametersOld[paramI][originalIndex];

            dataIndices_[paramI][paramJ] =
                dataIndicesOld[paramI][originalIndex];
        }
    }

    secondaryParametersOld = secondaryParameters_;
    tertiaryParametersOld = tertiaryParameters_;
    dataIndicesOld = dataIndices_;

    // Sort primary parameters
    SortableList<scalar> sortedListPrimary(primaryParameters_);
    primaryParameters_ = sortedListPrimary;
    const labelList originalIndices(sortedListPrimary.indices());

    forAll(originalIndices, paramI)
    {
        label originalIndex = originalIndices[paramI];

        secondaryParameters_[paramI] =
            secondaryParametersOld[originalIndex];

        tertiaryParameters_[paramI] =
            tertiaryParametersOld[originalIndex];

        dataIndices_[paramI] =
            dataIndicesOld[originalIndex];
    }
}


// Construct as a copy
dataFriction::dataFriction
(
    const dataFriction& fricLaw
)
:
    frictionLaw(fricLaw),
    mesh_(fricLaw.mesh_),
    frictionLawDict_(fricLaw.frictionLawDict_),
    frictionCoeff_(fricLaw.frictionCoeff_),
    parameterList_(fricLaw.parameterList_),
    primaryParameters_(),
    secondaryParameters_(),
    tertiaryParameters_(),
    dataIndices_(),
    frictionCoeffList_(fricLaw.frictionCoeffList_),
    useNearest_(fricLaw.useNearest_),
    SRRPriority_(fricLaw.SRRPriority_),
    pressurePriority_(fricLaw.pressurePriority_),
    temperaturePriority_(fricLaw.temperaturePriority_),
    fixedSRR_(fricLaw.fixedSRR_),
    fixedPressure_(fricLaw.fixedPressure_),
    fixedTemperature_(fricLaw.fixedTemperature_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dataFriction::~dataFriction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar dataFriction::slipTraction
(
    const scalar pressure,                // Contact pressure
    const vector& faceSlip,               // Slip vector
    const vector& slaveVelocity,          // Velocity of slave face
    const vector& masterVelocity,         // Velocity of master face
    const label patchIndex,               // Slave patch index
    const label faceIndex                 // Local slave face ID
)
{
    // Calculate entrainment speed
    const scalar curEntrSpeed = mag(0.5*(masterVelocity + slaveVelocity));

    // Calculate slide-roll-ratio (SRR)
    scalar curSRR =
        mag(masterVelocity - slaveVelocity)
       /(mag(0.5*(masterVelocity + slaveVelocity)) + VSMALL);

    scalar curPressure = pressure;
    scalar curTemperature = 0;

    // Check if SRR is set as fixed
    if (fixedSRR_.first())
    {
        curSRR = fixedSRR_.second();
    }

    // Check if pressure is set as fixed
    if (fixedPressure_.first())
    {
        curPressure = fixedPressure_.second();
    }

    // Check if temperature is set as fixed
    if (fixedTemperature_.first())
    {
        curTemperature = fixedTemperature_.second();
    }
    else
    {
        curTemperature = T().boundaryField()[patchIndex][faceIndex];
    }

    // Calculate friction coefficient
    calcFrictionCoeff(curEntrSpeed, curSRR, curPressure, curTemperature);

    return frictionCoeff_*curPressure;
}


void dataFriction::writeDict(Ostream& os) const
{
    word keyword("frictionLawDict");
    os.writeKeyword(keyword)
        << frictionLawDict_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
