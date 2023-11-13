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

InClass
    lubricatedContact

\*---------------------------------------------------------------------------*/

#include "lubricatedContact.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"
#include "fixedValueFaPatchFields.H"
#include "interpolateXY.H"
#include "coordinateSystem.H"
#include "lubricatedContactIntegration/lubricatedContactIntegration.H"
#include "RodriguesRotation.H"
#include "DynamicList.H"
#include <gsl/gsl_integration.h>
#include "flowFactors/flowFactors.H"
#include "SortableList.H"

using namespace Foam::mathematicalConstant;
using namespace Foam::lubricatedContactIntegration;
using namespace Foam::flowFactors;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(lubricatedContact, 0);
  addToRunTimeSelectionTable(normalContactModel, lubricatedContact, dictionary);


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

void lubricatedContact::createFiniteArea()
{
    // Create faMesh
    aMeshPtr_ = new faMesh(mesh_);

    // Create Finite Area fields
    const faMesh& aMesh = *aMeshPtr_;

    surfaceSeparationPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "surfaceSeparation",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar("surfaceSeparationInitial", dimLength, 0.0),
            zeroGradientFaPatchScalarField::typeName
        );

    filmBulkModulusPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "filmBulkModulus",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar
            (
                "filmBulkModulus", dimPressure, 0.069e09
            )
        );

    filmHydrodynamicPressurePtr_ =
        new areaScalarField
        (
            IOobject
            (
                "filmHydrodynamicPressure",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh
        );

    filmDensityPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "filmDensity",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            refDensityL_,
            filmHydrodynamicPressurePtr_->boundaryField().types()
    );
    calcFilmDensity();
    calcFilmBulkModulus();

    filmViscosityPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "filmViscosity",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            refViscosityL_
        );

    filmMeanThicknessPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "filmMeanThickness",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            minFilmMeanThickness_,
            zeroGradientFaPatchScalarField::typeName
        );

    filmPressureFlowFactorPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "filmPressureFlowFactor",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar("filmPressureFlowFactor", dimless, 0.0),
            zeroGradientFaPatchScalarField::typeName
        );

    filmShearFlowFactorPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "filmShearFlowFactor",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar("filmShearFlowFactor", dimless, 1.0),
            zeroGradientFaPatchScalarField::typeName
        );

    filmShearStressFactorPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "filmShearStressFactor",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar("filmShearStressFactor", dimless, 1.0),
            zeroGradientFaPatchScalarField::typeName
        );

    filmShearStressSlidingFactorPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "filmShearStressSlidingFactor",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar("filmShearStressSlidingFactor", dimless, 1.0),
            zeroGradientFaPatchScalarField::typeName
        );

    filmShearStressPtr_ =
        new areaVectorField
        (
            IOobject
            (
                "filmShearStress",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedVector("filmShearStress", dimPressure, vector::zero),
            zeroGradientFaPatchVectorField::typeName
        );

    alphaPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "alpha",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar("alpha", dimless, 1.0),
            filmDensityPtr_->boundaryField().types()
        );

    filmTemperaturePtr_ =
        new areaScalarField
        (
            IOobject
            (
                "filmTemperature",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            refTemperatureL_,
            zeroGradientFaPatchScalarField::typeName
        );

    asperityPressurePtr_ =
        new areaVectorField
        (
            IOobject
            (
                "asperityPressure",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedVector("asperityPressure", dimPressure, vector::zero)
        );

    asperityContactAreaPtr_ =
        new areaScalarField
        (
            IOobject
            (
                "asperityContactArea",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh,
            dimensionedScalar("asperityContactArea", dimless, 0.0)
        );

    // Check for processor patches
    DynamicList<label> processorPatchId;

    forAll(aMesh.boundary(), patchI)
    {
        if (aMesh.boundary()[patchI].type() == processorFaPatch::typeName)
        {
            processorPatchId.append(patchI);
        }
    }

    if (processorPatchId.size())
    {
        processorPatchId_.setSize(processorPatchId.size());
        processorPatchId_ = processorPatchId;
    }
}


void lubricatedContact::createInterpolationTables()
{
    // Asperity contact tables

    // Read separation, asperity and contact area tables
    fileName constantDir(mesh_.time().constant());

    // Create dictionary which will be read by interpolationTable
    // Consists of fileName and outOfBounds option
    dictionary contactTablesDict("contactTablesDict");
    contactTablesDict.add("outOfBounds", "clamp");

    // Separation vs aperity pressure table
    word tableName = "separationVsAsperityContactLoad";
    contactTablesDict.add("fileName", constantDir/tableName, true);

    separationVsAsperityPressureTablePtr_ =
        new interpolationTable<scalar>(contactTablesDict);

    // Asperity pressure vs separation table
    tableName = "asperityContactLoadVsSeparation";
    contactTablesDict.add("fileName", constantDir/tableName, true);

    asperityPressureVsSeparationTablePtr_ =
        new interpolationTable<scalar>(contactTablesDict);

    // Separation vs aperity contact area table
    tableName = "separationVsAsperityContactArea";
    contactTablesDict.add("fileName", constantDir/tableName, true);

    separationVsAsperityContactAreaTablePtr_ =
        new interpolationTable<scalar>(contactTablesDict);
}


void lubricatedContact::calcLocalRoughnessRotationTensor()
{
    const fvMesh& mesh = mesh_;
    const faMesh& aMesh = *aMeshPtr_;
    const scalar previousTime =
        (mesh.time().value() - mesh.time().deltaT0().value());

    IOobject localRoughnessRotationTensorHeader
    (
        "localRoughnessRotationTensor",
         mesh.time().timeName(previousTime),
         mesh,
         IOobject::MUST_READ
    );

    // If it is already calculated read it from the file (needed for restart)
    if (localRoughnessRotationTensorHeader.headerOk())
    {
        localRoughnessRotationTensorPtr_ =
            new areaTensorField
            (
                IOobject
                (
                    "localRoughnessRotationTensor",
                    mesh.time().timeName(previousTime),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                aMesh
            );
    }
    else
    {
        localRoughnessRotationTensorPtr_ =
            new areaTensorField
            (
                IOobject
                (
                    "localRoughnessRotationTensor",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                aMesh,
                dimensionedTensor
                (
                    "zeroTensor", dimless, tensor::zero
                ),
                zeroGradientFaPatchTensorField::typeName
            );

        areaTensorField& localRoughnessRotationTensor =
            *localRoughnessRotationTensorPtr_;

        const faceList& slaveFaces = slaveFaceZonePatch_.localFaces();
        const pointField& slavePoints = slaveFaceZonePatch_.localPoints();
        const vector axialRoughnessDirection =
            axialRoughnessDirection_/mag(axialRoughnessDirection_);

        forAll(localRoughnessRotationTensor, faceI)
        {
            const face& curFace = slaveFaces[faceI];
            const pointField facePoints = curFace.points(slavePoints);
            const point faceCentre = curFace.centre(slavePoints);
            vector faceNormal = curFace.normal(slavePoints);
            faceNormal /= mag(faceNormal);

            // Create local coordinate system
            const vector& e3 = faceNormal;
            vector e1 = facePoints[0] - faceCentre;
            e1 /= mag(e1);
            coordinateSystem faceCs("faceCs", faceCentre, e3, e1);

            // Calculate local rotation tensor for rotating vector of local
            // coordinate system e1 axis to axialRoughnessDirection around e3
            // axis
            const tensor rotationTensor =
                RodriguesRotation
                (
                    faceCs.localVector(e3),
                    faceCs.localVector(e1),
                    faceCs.localVector(axialRoughnessDirection)
                );

            localRoughnessRotationTensor[faceI] = rotationTensor;
        }
        localRoughnessRotationTensor.correctBoundaryConditions();
    }
}


void lubricatedContact::calcFilmMeanThickness()
{
    // Average flow thickness of the lubricant film is the average of the
    // distance between asperity and roll surface.
    areaScalarField& surfaceSeparation = *surfaceSeparationPtr_;
    areaScalarField& filmMeanThickness = *filmMeanThicknessPtr_;
    dimensionedScalar smallLength("smallLength", dimLength, SMALL);

    filmMeanThickness =
        0.5*surfaceSeparation
       *(
            erf
            (
                sqrt(0.5)*surfaceSeparation
               /(compositeSurfaceRoughness_ + smallLength)
            )
          + 1.0
        )
      + (
            compositeSurfaceRoughness_
           *exp
            (
                -0.5
               *sqr
                (
                    surfaceSeparation/(compositeSurfaceRoughness_ + smallLength)
                )
            )
           /sqrt(2*pi)
        );

    filmMeanThickness =
        max
        (
            min(filmMeanThickness, maxFilmMeanThickness_), minFilmMeanThickness_
        );
}


Foam::tmp<Foam::areaVectorField> lubricatedContact::slipMeanVelocity
(
    const word& velocityType,
    const vectorField& slavePatchFaceNormals,
    const vectorField& slaveDU,
    const vectorField& masterDUInterpToSlave
)
{
    // Calculate and return slip or mean velocity between slave and master
    const faMesh& aMesh = *aMeshPtr_;
    areaScalarField& filmMeanThickness = *filmMeanThicknessPtr_;

    tmp<areaVectorField> tVelocity
    (
        new areaVectorField
        (
            IOobject
            (
                velocityType + "Velocity",
                mesh_.time().timeName(startTime_),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh,
            dimensionedVector
            (
                velocityType + "Velocity",
                dimVelocity,
                vector::zero
            ),
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaVectorField& velocity = tVelocity();

    forAll(velocity, faceI)
    {
        if (filmMeanThickness[faceI] < maxFilmMeanThickness_.value())
        {
            if (velocityType == "slip")
            {
                velocity[faceI] =
                    slaveDU[faceI] - masterDUInterpToSlave[faceI];
            }
            else if (velocityType == "mean")
            {
                velocity[faceI] =
                    (slaveDU[faceI] + masterDUInterpToSlave[faceI])/2;
            }
            else if (velocityType == "slave")
            {
                velocity[faceI] = slaveDU[faceI];
            }
            else
            {
                FatalErrorIn
                (
                    "lubricatedContact::slipMeanVelocity()"
                )
                << "Wrong velocity type specified."
                << abort(FatalError);
            }

            const scalar deltaT = mesh_.time().deltaTValue();
            const areaVectorField& faceAreaNormals = aMesh.faceAreaNormals();

            velocity[faceI] =
                ((I - sqr(faceAreaNormals[faceI])) & velocity[faceI])
               /deltaT;
        }
    }
    velocity.correctBoundaryConditions();

    return tVelocity;
}


void lubricatedContact::calcPressureFlowFactors()
{
    // Preliminaries
    areaScalarField& filmPressureFlowFactor = *filmPressureFlowFactorPtr_;

    if (compositeSurfaceRoughness_.value() < SMALL)
    {
        filmPressureFlowFactor.internalField() = 1.0;
    }
    else
    {
        const faceList& slaveFaces = slaveFaceZonePatch_.localFaces();
        const pointField& slavePoints = slaveFaceZonePatch_.localPoints();
        const areaTensorField& localRoughnessRotationTensor =
            *localRoughnessRotationTensorPtr_;

        // Calculate pressure flow factors (Wilson, Marsault)
        const areaScalarField& filmMeanThickness = *filmMeanThicknessPtr_;
        const scalarField normFilmMeanThickness =
            filmMeanThickness.internalField()
           /(compositeSurfaceRoughness_.value() + SMALL);

        scalarField pressureFlowFactorsAxial =
            pressureFlowFactorWilsonMarsault
            (
                compositeCorrLengthRatio_,
                normFilmMeanThickness
            );

        scalarField pressureFlowFactorsCirc =
            pressureFlowFactorWilsonMarsault
            (
                1.0/compositeCorrLengthRatio_,
                normFilmMeanThickness
            );

    //    // Calculate pressure flow factors (Patir,Cheng)
    //    const areaScalarField& surfaceSeparation = *surfaceSeparationPtr_;
    //    const scalarField normSurfaceSeparation =
    //        surfaceSeparation.internalField()/compositeSurfaceRoughness_.value();
    //    scalarField pressureFlowFactorsAxial =
    //        pressureFlowPatirCheng
    //        (
    //            compositeCorrLengthRatio_,
    //            normSurfaceSeparation
    //        );
    //
    //    scalarField pressureFlowFactorsCirc =
    //        pressureFlowPatirCheng
    //        (
    //            1.0/compositeCorrLengthRatio_,
    //            normSurfaceSeparation
    //        );

        // Calculate direction of the density gradient vector
        const areaScalarField& filmDensity = *filmDensityPtr_;
        tmp<areaVectorField> pressureGradient = fac::grad(filmDensity);
        vectorField pressureGradientDir = pressureGradient->internalField();

        forAll(pressureGradientDir, faceI)
        {
            if (mag(pressureGradientDir[faceI]) > SMALL)
            {
                pressureGradientDir[faceI] /= mag(pressureGradientDir[faceI]);
            }
        }


        // Calculate composite pressure flow factors
        forAll(filmPressureFlowFactor, faceI)
        {
            const face& curFace = slaveFaces[faceI];
            const pointField facePoints = curFace.points(slavePoints);
            const point faceCentre = curFace.centre(slavePoints);
            vector faceNormal = curFace.normal(slavePoints);
            faceNormal /= mag(faceNormal);

            // Create local coordinate system
            const vector& e3 = faceNormal;
            vector e1 = facePoints[0] - faceCentre;
            e1 /= mag(e1);
            coordinateSystem faceCs("faceCs", faceCentre, e3, e1);

            // Rotate e1 to axial surface roughness direction and convert it to
            // global coordinate system
            const vector globalAxialRoughnessVector =
                faceCs.globalVector
                (
                    localRoughnessRotationTensor[faceI] & faceCs.localVector(e1)
                );

            // Create new local coordinate system where axial roughness vector is e1
            e1 = globalAxialRoughnessVector/mag(globalAxialRoughnessVector);
            faceCs = coordinateSystem("faceCs", faceCentre, e3, e1);

            // Transform global normalised density gradient vector into local CS
            const vector localGradientDir =
                faceCs.localVector(pressureGradientDir[faceI]);

            // Using the ellipse equation (x^2/a^2 + y^2/b^2 = 1) where 'a' equals
            // axial flow factor and 'b' equals circular flow factor, we are
            // calculating the intersection point between the ellipse and
            // localGradientDir line. Magnitude of that point (vector) equals the
            // composite pressure flow factor.
            const scalar a = pressureFlowFactorsAxial[faceI];
            const scalar b = pressureFlowFactorsCirc[faceI];
            const scalar k = localGradientDir.y()/(localGradientDir.x() + SMALL);

            const scalar x = sqrt(sqr(a*b)/(sqr(b) + sqr(a*k) + SMALL));
            const scalar y = k*x;

            filmPressureFlowFactor[faceI] = mag(vector(x,y,0));
        }
    }
    filmPressureFlowFactor.correctBoundaryConditions();
}


void lubricatedContact::calcShearFlowFactors()
{
    // Preliminaries
    areaScalarField& filmShearFlowFactor = *filmShearFlowFactorPtr_;

    // Calculate shear flow factors (Wilson, Marsault)
    const areaScalarField& filmMeanThickness = *filmMeanThicknessPtr_;
    const scalarField normFilmMeanThickness =
        filmMeanThickness.internalField()
       /(compositeSurfaceRoughness_.value() + SMALL);

    filmShearFlowFactor.internalField() =
        sqr
        (
            slaveSurfaceRoughness_.value()
           /(compositeSurfaceRoughness_.value() + SMALL)
        )
       *shearFlowFactorsWilsonMarsault
        (
            slaveCorrLengthRatio_,
            normFilmMeanThickness
        )
      - sqr
        (
            masterSurfaceRoughness_.value()
           /(compositeSurfaceRoughness_.value() + SMALL)
        )
       *shearFlowFactorsWilsonMarsault
        (
            masterCorrLengthRatio_,
            normFilmMeanThickness
        );

//    // Calculate shear flow factors (Patir, Cheng)
//    const areaScalarField& surfaceSeparation = *surfaceSeparationPtr_;
//    const scalarField normSurfaceSeparation =
//        surfaceSeparation.internalField()/compositeSurfaceRoughness_.value();
//
//    filmShearFlowFactor.internalField() =
//        sqr(slaveSurfaceRoughness_/compositeSurfaceRoughness_).value()
//       *shearFlowFactorsPatirCheng
//        (
//            slaveCorrLengthRatio_,
//            normSurfaceSeparation
//        )
//      - sqr(masterSurfaceRoughness_/compositeSurfaceRoughness_).value()
//       *shearFlowFactorsPatirCheng
//        (
//            masterCorrLengthRatio_,
//            normSurfaceSeparation
//        );

    filmShearFlowFactor.correctBoundaryConditions();
}


void lubricatedContact::calcShearStressFactors()
{
    // Preliminaries
    areaScalarField& surfaceSeparation = *surfaceSeparationPtr_;
    areaScalarField& filmShearStressFactor = *filmShearStressFactorPtr_;

    // Calculate shear stress factors (Patir, Cheng)
    filmShearStressFactor.internalField() =
        sqr
        (
            slaveSurfaceRoughness_.value()
           /(compositeSurfaceRoughness_.value() + SMALL)
        )
       *shearStressFactorsPatirCheng
        (
            slaveCorrLengthRatio_,
            surfaceSeparation.internalField()
           /(compositeSurfaceRoughness_.value() + SMALL)
        )
      - sqr
        (
            masterSurfaceRoughness_.value()
           /(compositeSurfaceRoughness_.value() + SMALL)
        )
       *shearStressFactorsPatirCheng
        (
            masterCorrLengthRatio_,
            surfaceSeparation.internalField()
           /(compositeSurfaceRoughness_.value() + SMALL)
        );

    filmShearStressFactor.correctBoundaryConditions();
}


void lubricatedContact::calcShearStressSlidingFactors()
{
    // Preliminaries
    areaScalarField& surfaceSeparation = *surfaceSeparationPtr_;
    areaScalarField& filmShearStressSlidingFactor =
        *filmShearStressSlidingFactorPtr_;

    // Calculate shear stress sliding factors (Patir, Cheng)
    forAll(surfaceSeparation, faceI)
    {
        interpolationTable<scalar>& shearStressSlidingFactorsTable =
            *shearStressSlidingFactorsTablePtr_;

        filmShearStressSlidingFactor[faceI] =
            shearStressSlidingFactorsTable(surfaceSeparation[faceI]);
    }
    filmShearStressSlidingFactor.correctBoundaryConditions();
}


void lubricatedContact::calcShearStressSlidingFactorsTable()
{
    // Number of points in the interpolation table
    const label noPoints = 1000;

    // Create list
    List<Tuple2<scalar, scalar> > surfaceSeparationVsSssFactorList
    (
        noPoints, Tuple2<scalar, scalar>(0.0, 0.0)
    );

    // Create Gauss-Legendre table
    const label gaussLegendreOrder = 1024;
    gsl_integration_glfixed_table* gaussLegendreTablePtr =
        gsl_integration_glfixed_table_alloc(gaussLegendreOrder);

    // Loop through film thicknesses
    // NOTE: from h = 0.1*sigma to h = 100*sigma
    for (int i = 0; i < noPoints; ++i)
    {
        const scalar deltaH = 0.1*compositeSurfaceRoughness_.value();
        const scalar h = deltaH + i*deltaH;

        surfaceSeparationVsSssFactorList[i].first() = h;

        surfaceSeparationVsSssFactorList[i].second() =
            shearStressSlidingFactor
            (
                h,
                compositeSurfaceRoughness_.value(),
                gaussLegendreTablePtr
            );
    }

    // Interpolation table
    shearStressSlidingFactorsTablePtr_ =
        new interpolationTable<scalar>
        (
            surfaceSeparationVsSssFactorList,
            interpolationTable<scalar>().wordToBoundsHandling("clamp"),
            "shearStressSlidingFactorsTable"
        );
}


void lubricatedContact::calcFilmViscosity
(
    const areaVectorField& filmSlipVelocity
)
{
    // Calculate film viscosity depending on the fact that lubricant is either a
    // Newtonian or a Non-Newtonian fluid
    const areaScalarField& filmMeanThickness = *filmMeanThicknessPtr_;
    const areaScalarField& filmHydrodynamicPressure =
        *filmHydrodynamicPressurePtr_;
    const areaScalarField& filmTemperature = *filmTemperaturePtr_;
    areaScalarField& filmViscosity = *filmViscosityPtr_;

    // Roelands viscosity equation
    // Book: Gohar - Elastohydrodynamics (Imperial College Press)
    const scalar S0 =
        temperatureExpL_*(refTemperatureL_.value() - 138)
       /(log(refViscosityL_.value()) + 9.67);
    const scalar Z =
        pressureExpL_/(5.1e-09*(log(refViscosityL_.value()) + 9.67));

    // Pressure-viscosity coefficient
    const dimensionedScalar oneKelvin("oneKelvin", dimTemperature, 1.0);
    const dimensionedScalar onePascal("onePascal", dimPressure, 1.0);

    areaScalarField alphaStarP =
        (log(refViscosityL_.value()) + 9.67)
       *pow
        (
            (refTemperatureL_.value() - 138)/(filmTemperature/oneKelvin - 138),
            S0
        )
       *(
            pow(filmHydrodynamicPressure/onePascal/1.98e08 + 1, Z)
          - 1
        );

    // Temperature-viscosity coefficient
    areaScalarField betaStar =
        (
            (log(refViscosityL_.value()) + 9.67)
           *pow(5.1e-9*filmHydrodynamicPressure/onePascal + 1, Z)
           *(S0/(refTemperatureL_.value() - 138))
        );

    // Calculate viscosity
    filmViscosity =
        refViscosityL_
       *exp
        (
            alphaStarP
          - betaStar*(filmTemperature/oneKelvin - refTemperatureL_.value())
        );

    // If lubricant is a non-Newtonian fluid, Ree-Eyring model is used
    // PhD Thesis: Hartinger - CFD Modelling of Elastohydrodynamic
    //                         Lubrication
    if (nonNewtonianL_)
    {
        const dimensionedScalar oneSecond("oneSecond", dimTime, 1.0);

        const areaScalarField filmShearRate =
            max(mag(filmSlipVelocity)/filmMeanThickness, SMALL/oneSecond);

        filmViscosity =
            eyringStressL_/filmShearRate
           *asinh(filmViscosity*filmShearRate/eyringStressL_);
    }
}


void lubricatedContact::calcFilmDensity()
{
    // Calculating lubricant density dependant on pressure
    areaScalarField& filmDensity = *filmDensityPtr_;
    areaScalarField& filmHydrodynamicPressure = *filmHydrodynamicPressurePtr_;

    // Density-pressure relation for lubricating oil was developed by Dowson
    // PhD Thesis: Hartinger - CFD Modelling of Elastohydrodynamic Lubrication
    const dimensionedScalar onePascal("onePascal", dimPressure, 1.0);

    filmDensity =
        refDensityL_
       *(
            (C1_*onePascal + C2_*filmHydrodynamicPressure)
           /(C1_*onePascal + filmHydrodynamicPressure)
        );

    // Force assignment to filmDensity boundary field
    filmDensity.boundaryField() ==
        refDensityL_.value()
       *(
            (C1_ + C2_*filmHydrodynamicPressure.boundaryField())
           /(C1_ + filmHydrodynamicPressure.boundaryField())
        );
}


Foam::areaScalarField lubricatedContact::filmPressure
(
    const areaScalarField& filmDensity
)
{
    // Calculating lubricant pressure dependant on film density

    // Density-pressure relation for lubricating oil was developed by Dowson
    // PhD Thesis: Hartinger - CFD Modelling of Elastohydrodynamic Lubrication
    const dimensionedScalar onePascal("onePascal", dimPressure, 1.0);
    const areaScalarField limitedFilmDensity =
        min(max(filmDensity, cavDensityL_), 0.99*refDensityL_*C2_);

    return
    (
        (refDensityL_ - limitedFilmDensity)*(C1_*onePascal)
       /(limitedFilmDensity - refDensityL_*C2_)
    );
}


void lubricatedContact::calcFilmBulkModulus()
{
    // Calculating lubricant bulk modulus  = rho*d(p)/d(rho)
    // Density-pressure relation is taken from Dowson
    areaScalarField& filmDensity = *filmDensityPtr_;
    areaScalarField& filmBulkModulus = *filmBulkModulusPtr_;

    const dimensionedScalar onePascal("onePascal", dimPressure, 1.0);
    const areaScalarField limitedFilmDensity =
        min(max(filmDensity, cavDensityL_), 0.99*refDensityL_*C2_);

    filmBulkModulus =
        limitedFilmDensity
       *(C1_*(C2_ - 1)*refDensityL_/sqr(limitedFilmDensity - C2_*refDensityL_))
       *onePascal;
}


void lubricatedContact::relaxPointPenetration
(
    const scalarField& slavePointPenetration
)
{
    // Use fixed relaxation factor since variable factor calculated by the
    // Aitken procedure does not improve solution convergence
    const scalarField residual =
        slavePointPenetration - fluidSlavePointPenetration_;

    forAll(slavePointPenetration, pointI)
    {
        if
        (
            slavePointPenetration[pointI] > 100
         || fluidSlavePointPenetration_[pointI] > 100
        )
        {
            fluidSlavePointPenetration_[pointI] = slavePointPenetration[pointI];
        }
        else
        {
            fluidSlavePointPenetration_[pointI] += relaxFac_*residual[pointI];
        }
    }
}


void lubricatedContact::calcContactBoundaryVelocities
(
    const dynamicLabelList& nonContactFaces,
    const Field<label>& contact,
    areaVectorField& filmMeanVelocity,
    areaVectorField& filmSlipVelocity
)
{
    // Set velocity into faces which are disregarded by GGI due to quick search
    // algorithm, but have active faces (with measured distances by GGI) as
    // their neighbours. The inserted velocities are calculated as mean velocity
    // of the neighbour faces.
    const faMesh& aMesh = *aMeshPtr_;

    forAll(nonContactFaces, faceI)
    {
        const label curFace = nonContactFaces[faceI];
        const labelList& curFaceFaces =
                    aMesh.patch().faceFaces()[curFace];

        vector curFaceMeanMeanVelocity = vector::zero;
        vector curFaceMeanSlipVelocity = vector::zero;
        label nActiveNeighbours = 0;

        forAll(curFaceFaces, faceJ)
        {
            if (contact[curFaceFaces[faceJ]] == 1)
            {
                curFaceMeanMeanVelocity +=
                    filmMeanVelocity[curFaceFaces[faceJ]];
                curFaceMeanSlipVelocity +=
                    filmSlipVelocity[curFaceFaces[faceJ]];

                nActiveNeighbours ++;
            }
        }

        if (nActiveNeighbours)
        {
            curFaceMeanMeanVelocity /= nActiveNeighbours;
            curFaceMeanSlipVelocity /= nActiveNeighbours;

            filmMeanVelocity[curFace] = curFaceMeanMeanVelocity;
            filmSlipVelocity[curFace] = curFaceMeanSlipVelocity;
        }
    }

    const areaVectorField& faceAreaNormals = aMesh.faceAreaNormals();

    filmSlipVelocity.internalField() =
        (
            (I - sqr(faceAreaNormals.internalField()))
          & filmSlipVelocity.internalField()
        );
    filmMeanVelocity.internalField() =
        (
            (I - sqr(faceAreaNormals.internalField()))
          & filmMeanVelocity.internalField()
        );

    filmSlipVelocity.correctBoundaryConditions();
    filmMeanVelocity.correctBoundaryConditions();
}


void lubricatedContact::calcCavitationBoundary
(
    edgeScalarField& alphaEdge,
    const edgeScalarField& gammaEdge,
    scalarField& cavBoundaryDiagonal,
    scalarField& cavBoundarySource,
    const areaVectorField& filmMeanVelocity
)
{
    // Preliminaries
    const faMesh& aMesh = *aMeshPtr_;
    areaScalarField& alpha = *alphaPtr_;
    areaScalarField& filmDensity = *filmDensityPtr_;

    // Get internal edge owner and neighbour
    const unallocLabelList& owner = aMesh.owner();
    const unallocLabelList& neighbour = aMesh.neighbour();

    // Get delta coeffs and edge length magnitudes
    const edgeScalarField& deltaCoeffs =
        aMesh.edgeInterpolation::deltaCoeffs();
    const edgeScalarField& magLe = aMesh.magLe();

    // Calculate mean film velocity flux
    edgeScalarField meanFlux =
        fac::interpolate(filmMeanVelocity, "interpolate(filmMeanVelocity)")
      & aMesh.Le();

    // Cavitation density
    const scalar cavDensity = cavDensityL_.value();


    // Loop through all mesh internal edges and calculate cavitation boundary
    // contribution
    forAll(meanFlux, edgeI)
    {
        const label curOwner = owner[edgeI];
        const label curNeighbour = neighbour[edgeI];

        // If edge is between active and non-active zone
        if (alpha[curOwner] != alpha[curNeighbour])
        {
            // Set current alpha edge to 0
            alphaEdge[edgeI] = 0.0;

            // Note active and non-active face of this edge
            label curEdgeActiveFace = -1;
            label curEdgeNonActiveFace = -1;

            if (alpha[curOwner] == 1.0 && alpha[curNeighbour] == 0.0)
            {
                curEdgeActiveFace = curOwner;
                curEdgeNonActiveFace = curNeighbour;
            }
            else if (alpha[curNeighbour] == 1.0 && alpha[curOwner] == 0.0)
            {
                curEdgeActiveFace = curNeighbour;
                curEdgeNonActiveFace = curOwner;
            }
            else
            {
                FatalErrorIn("boundary.H")
                    << "Strange internal edge." << abort(FatalError);
            }

            // Determine edge type and calculate cavitation boundary
            // contribution, i.e. fluxes

            // Laplacian upper coefficient
            const scalar upperCoeff =
                deltaCoeffs[edgeI]*gammaEdge[edgeI]*magLe[edgeI];

            // Formation or rupture edge
            if (mag(meanFlux[edgeI]) > SMALL)
            {
                cavBoundaryDiagonal[curEdgeActiveFace] += upperCoeff;

                cavBoundarySource[curEdgeActiveFace] += upperCoeff*cavDensity;
                cavBoundarySource[curEdgeNonActiveFace] +=
                    upperCoeff*(filmDensity[curEdgeActiveFace] - cavDensity);
            }
            // Zero flux cavitation boundary edge
            else
            {
                cavBoundaryDiagonal[curEdgeActiveFace] += upperCoeff;
                cavBoundarySource[curEdgeActiveFace] += upperCoeff*cavDensity;
            }

        }
    }


    // Loop through all processor edges and calculate cavitation boundary
    // contribution
    forAll(processorPatchId_, patchI)
    {
        const label curPatchId = processorPatchId_[patchI];
        // Edge owners of patch edges
        const unallocLabelList& edgeFaces =
            aMesh.boundary()[curPatchId].edgeFaces();
        // Patch neighbout field of filmDensity
        const scalarField& filmDensityPnf =
            filmDensity.boundaryField()[curPatchId].patchNeighbourField()();

        forAll(meanFlux.boundaryField()[curPatchId], edgeI)
        {
            const scalar alphaOwner =
                alpha.boundaryField()[curPatchId].patchInternalField()()[edgeI];
            const scalar alphaNeighbour =
                alpha.boundaryField()[curPatchId].patchNeighbourField()()[edgeI];

            // If edge is between active and non-active zone
            if (alphaOwner != alphaNeighbour)
            {
                // Set current alpha edge to 0
                alphaEdge.boundaryField()[curPatchId][edgeI] = 0.0;

                // Determine edge type and calculate cavitation boundary
                // attribution, i.e. fluxes
                const scalar upperCoeff =
                    deltaCoeffs.boundaryField()[curPatchId][edgeI]
                   *gammaEdge.boundaryField()[curPatchId][edgeI]
                   *magLe.boundaryField()[curPatchId][edgeI];

                if (mag(meanFlux.boundaryField()[curPatchId][edgeI]) > SMALL)
                {
                    if (alphaOwner == 1.0 && alphaNeighbour == 0.0)
                    {
                        cavBoundaryDiagonal[edgeFaces[edgeI]] += upperCoeff;

                        cavBoundarySource[edgeFaces[edgeI]] +=
                            upperCoeff*cavDensity;
                    }
                    else if (alphaOwner == 0.0 && alphaNeighbour == 1.0)
                    {
                        cavBoundarySource[edgeFaces[edgeI]] +=
                            upperCoeff*(filmDensityPnf[edgeI] - cavDensity);
                    }
                }
                else
                {
                    if (alphaOwner == 1.0 && alphaNeighbour == 0.0)
                    {
                        cavBoundaryDiagonal[edgeFaces[edgeI]] += upperCoeff;

                        cavBoundarySource[edgeFaces[edgeI]] +=
                            upperCoeff*cavDensity;
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lubricatedContact::lubricatedContact
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const standAlonePatch& masterFaceZonePatch,
    const standAlonePatch& slaveFaceZonePatch
)
:
    normalContactModel
    (
        name,
        patch,
        dict,
        masterPatchID,
        slavePatchID,
        masterFaceZonePatch,
        slaveFaceZonePatch
    ),
    normalContactModelDict_(dict.subDict(name + "NormalModelDict")),
    mesh_(patch.boundaryMesh().mesh()),
    slaveFaceZonePatch_(slaveFaceZonePatch),
    startTime_(mesh_.time().startTime().value()),
    curTime_(0.0),
    curCorrector_(0),
    slavePressure_(mesh_.boundaryMesh()[slavePatchID].size(), vector::zero),
    relaxFac_
    (
        normalContactModelDict_.lookupOrDefault<scalar>
        (
            "relaxationFactor", 0.02
        )
    ),
    standardPenaltyPtr_(NULL),
    lubricationDict_(normalContactModelDict_.subDict("lubricationDict")),
    maxLubricationCorr_
    (
        lubricationDict_.lookupOrDefault<label>("nCorr", 200)
    ),
    lubricationStartTime_
    (
        readScalar(lubricationDict_.lookup("lubricationStartTime"))
    ),
    surfaceRoughnessDict_
    (
        IOobject
        (
            "surfaceRoughnessProperties",
            mesh_.time().constant(),
            mesh_.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    aMeshPtr_(NULL),
    surfaceSeparationPtr_(NULL),
    filmHydrodynamicPressurePtr_(NULL),
    filmDensityPtr_(NULL),
    filmViscosityPtr_(NULL),
    filmBulkModulusPtr_(NULL),
    filmMeanThicknessPtr_(NULL),
    compositeSurfaceRoughness_
    (
        surfaceRoughnessDict_.lookup("compositeSurfaceRoughness")
    ),
    masterSurfaceRoughness_
    (
        surfaceRoughnessDict_.lookup("masterSurfaceRoughness")
    ),
    slaveSurfaceRoughness_
    (
        surfaceRoughnessDict_.lookup("slaveSurfaceRoughness")
    ),
    maxFilmMeanThickness_
    (
        "maxFilmMeanThicknessDefault",
        dimLength,
        lubricationDict_.lookupOrDefault("maxFilmMeanThickness", 1e-03)
    ),
    minFilmMeanThickness_
    (
        "minFilmMeanThicknessDefault",
        dimLength,
        lubricationDict_.lookupOrDefault("minFilmMeanThickness", 1e-08)
    ),
    refTemperatureL_
    (
        "refTemperature",
        dimTemperature,
        readScalar(lubricationDict_.lookup("referenceTemperature"))
    ),
    refViscosityL_
    (
        "refViscosity",
        dimPressure*dimTime,
        readScalar(lubricationDict_.lookup("referenceViscosity"))
    ),
    refDensityL_
    (
        "refDensity",
        dimDensity,
        readScalar(lubricationDict_.lookup("referenceDensity"))
    ),
    cavDensityL_(refDensityL_),
    cavPressureL_
    (
        "cavPressure",
        dimPressure,
        lubricationDict_.lookupOrDefault<scalar>("cavPressure", 0.0)
    ),
    pressureExpL_
    (
        readScalar(lubricationDict_.lookup("pressureExponent"))
    ),
    temperatureExpL_
    (
        readScalar(lubricationDict_.lookup("temperatureExponent"))
    ),
    C1_(2.22e09),
    C2_(1.66),
    calcViscosityL_
    (
        readBool(lubricationDict_.lookup("calculateViscosity"))
    ),
    nonNewtonianL_
    (
        readBool(lubricationDict_.lookup("nonNewtonian"))
    ),
    eyringStressL_("eyringStress", dimPressure, 0.0),
    separationVsAsperityPressureTablePtr_(NULL),
    asperityPressureVsSeparationTablePtr_(NULL),
    separationVsAsperityContactAreaTablePtr_(NULL),
    compositeCorrLengthRatio_
    (
        surfaceRoughnessDict_.lookupOrDefault<scalar>
        (
            "compositeCorrelationLengthRatio", 1.0
        )
    ),
    masterCorrLengthRatio_
    (
        surfaceRoughnessDict_.lookupOrDefault<scalar>
        (
            "masterCorrelationLengthRatio", 1.0
        )
    ),
    slaveCorrLengthRatio_
    (
        surfaceRoughnessDict_.lookupOrDefault<scalar>
        (
            "slaveCorrelationLengthRatio", 1.0
        )
    ),
    filmPressureFlowFactorPtr_(NULL),
    filmShearFlowFactorPtr_(NULL),
    filmShearStressFactorPtr_(NULL),
    filmShearStressSlidingFactorPtr_(NULL),
    filmShearStressPtr_(NULL),
    alphaPtr_(NULL),
    filmTemperaturePtr_(NULL),
    asperityPressurePtr_(NULL),
    asperityContactAreaPtr_(NULL),
    axialRoughnessDirection_
    (
        surfaceRoughnessDict_.lookup("axialRoughnessDirection")
    ),
    localRoughnessRotationTensorPtr_(NULL),
    surfaceRoughnessDirections_
    (
        2, vectorField(mesh_.boundaryMesh()[slavePatchID].size(), vector::zero)
    ),
    shearStressSlidingFactorsTablePtr_(NULL),
    maxFilmPressure_(readScalar(lubricationDict_.lookup("maxFilmPressure"))),
    debugL_(lubricationDict_.lookupOrDefault<bool>("debugMode", false)),
    rigidContact_
    (
        lubricationDict_.lookupOrDefault<bool>("rigidContact", false)
    ),
    fluidSlavePointPenetration_
    (
        mesh_.boundaryMesh()[slavePatchID].nPoints(), 0.0
    ),
    processorPatchId_(),
    boxesOfInterest_
    (
        lubricationDict_.lookupOrDefault("boxesOfInterest", List<boundBox>())
    )
{
    // Calculate cavitation density
    const dimensionedScalar onePascal("onePascal", dimPressure, 1.0);

    cavDensityL_ =
        refDensityL_
       *((C1_*onePascal + C2_*cavPressureL_)/(C1_*onePascal + cavPressureL_));

    // If lubricant is non-Newtonian, read Eyring stress
    if (nonNewtonianL_)
    {
        eyringStressL_.value() =
            readScalar(lubricationDict_.lookup("EyringStress"));
    }

    // Create interpolation tables
    createInterpolationTables();

    // The new point distance method is not working at the moment, so we need to
    // use the old one.
    if (dict.lookupOrDefault<Switch>("useNewPointDistanceMethod", true))
    {
        FatalErrorIn
        (
            "lubricatedContact::lubricatedContact()"
        )
        << "The new point distance method does not work with lubricatedContact"
        << "\ncontact model. Please set \"useNewPointDistanceMethod\" to "
        << "\"false\". "
        << abort(FatalError);
    }
}


lubricatedContact::lubricatedContact(const lubricatedContact& nm)
:
    normalContactModel(nm),
    normalContactModelDict_(nm.normalContactModelDict_),
    mesh_(nm.mesh_),
    slaveFaceZonePatch_(nm.slaveFaceZonePatch_),
    startTime_(nm.startTime_),
    curTime_(nm.curTime_),
    curCorrector_(nm.curCorrector_),
    slavePressure_(nm.slavePressure_),
    relaxFac_(nm.relaxFac_),
    standardPenaltyPtr_(NULL),
    lubricationDict_(nm.lubricationDict_),
    maxLubricationCorr_(nm.maxLubricationCorr_),
    lubricationStartTime_(nm.lubricationStartTime_),
    surfaceRoughnessDict_(nm.surfaceRoughnessDict_),
    aMeshPtr_(nm.aMeshPtr_),
    surfaceSeparationPtr_(NULL),
    filmHydrodynamicPressurePtr_(NULL),
    filmDensityPtr_(NULL),
    filmViscosityPtr_(NULL),
    filmBulkModulusPtr_(NULL),
    filmMeanThicknessPtr_(NULL),
    compositeSurfaceRoughness_(nm.compositeSurfaceRoughness_),
    masterSurfaceRoughness_(nm.masterSurfaceRoughness_),
    slaveSurfaceRoughness_(nm.slaveSurfaceRoughness_),
    maxFilmMeanThickness_(nm.maxFilmMeanThickness_),
    minFilmMeanThickness_(nm.minFilmMeanThickness_),
    refTemperatureL_(nm.refTemperatureL_),
    refViscosityL_(nm.refViscosityL_),
    refDensityL_(nm.refDensityL_),
    cavDensityL_(nm.cavDensityL_),
    cavPressureL_(nm.cavPressureL_),
    pressureExpL_(nm.pressureExpL_),
    temperatureExpL_(nm.temperatureExpL_),
    C1_(nm.C1_),
    C2_(nm.C2_),
    calcViscosityL_(nm.calcViscosityL_),
    nonNewtonianL_(nm.nonNewtonianL_),
    eyringStressL_(nm.eyringStressL_),
    separationVsAsperityPressureTablePtr_(NULL),
    asperityPressureVsSeparationTablePtr_(NULL),
    separationVsAsperityContactAreaTablePtr_(NULL),
    compositeCorrLengthRatio_(nm.compositeCorrLengthRatio_),
    masterCorrLengthRatio_(nm.masterCorrLengthRatio_),
    slaveCorrLengthRatio_(nm.slaveCorrLengthRatio_),
    filmPressureFlowFactorPtr_(NULL),
    filmShearFlowFactorPtr_(NULL),
    filmShearStressFactorPtr_(NULL),
    filmShearStressSlidingFactorPtr_(NULL),
    filmShearStressPtr_(NULL),
    alphaPtr_(NULL),
    filmTemperaturePtr_(NULL),
    asperityPressurePtr_(NULL),
    asperityContactAreaPtr_(NULL),
    axialRoughnessDirection_(nm.axialRoughnessDirection_),
    localRoughnessRotationTensorPtr_(NULL),
    surfaceRoughnessDirections_(nm.surfaceRoughnessDirections_),
    shearStressSlidingFactorsTablePtr_(NULL),
    maxFilmPressure_(nm.maxFilmPressure_),
    debugL_(nm.debugL_),
    rigidContact_(nm.rigidContact_),
    fluidSlavePointPenetration_(nm.fluidSlavePointPenetration_),
    processorPatchId_(nm.processorPatchId_),
    boxesOfInterest_(nm.boxesOfInterest_)
{
    if (nm.standardPenaltyPtr_)
    {
        standardPenaltyPtr_ = nm.standardPenaltyPtr_->clone().ptr();
    }

    if (nm.surfaceSeparationPtr_)
    {
        surfaceSeparationPtr_ = new areaScalarField(*nm.surfaceSeparationPtr_);
    }

    if (nm.filmHydrodynamicPressurePtr_)
    {
        filmHydrodynamicPressurePtr_ =
            new areaScalarField(*nm.filmHydrodynamicPressurePtr_);
    }

    if (nm.filmDensityPtr_)
    {
        filmDensityPtr_ =
            new areaScalarField(*nm.filmDensityPtr_);
    }

    if (nm.filmViscosityPtr_)
    {
        filmViscosityPtr_ =
            new areaScalarField(*nm.filmViscosityPtr_);
    }

    if (nm.filmBulkModulusPtr_)
    {
        filmBulkModulusPtr_ =
            new areaScalarField(*nm.filmBulkModulusPtr_);
    }

    if (nm.filmMeanThicknessPtr_)
    {
        filmMeanThicknessPtr_ =
            new areaScalarField(*nm.filmMeanThicknessPtr_);
    }

    if (nm.separationVsAsperityPressureTablePtr_)
    {
        separationVsAsperityPressureTablePtr_ =
            new interpolationTable<scalar>
            (
                *nm.separationVsAsperityPressureTablePtr_
            );
    }

    if (nm.asperityPressureVsSeparationTablePtr_)
    {
        asperityPressureVsSeparationTablePtr_ =
            new interpolationTable<scalar>
            (
                *nm.asperityPressureVsSeparationTablePtr_
            );
    }

    if (nm.separationVsAsperityContactAreaTablePtr_)
    {
        separationVsAsperityContactAreaTablePtr_ =
            new interpolationTable<scalar>
            (
                *nm.separationVsAsperityContactAreaTablePtr_
            );
    }

    if (nm.filmPressureFlowFactorPtr_)
    {
        filmPressureFlowFactorPtr_ =
            new areaScalarField(*nm.filmPressureFlowFactorPtr_);
    }

    if (nm.filmShearFlowFactorPtr_)
    {
        filmShearFlowFactorPtr_ =
            new areaScalarField(*nm.filmShearFlowFactorPtr_);
    }

    if (nm.filmShearStressFactorPtr_)
    {
        filmShearStressFactorPtr_ =
            new areaScalarField(*nm.filmShearStressFactorPtr_);
    }

    if (nm.localRoughnessRotationTensorPtr_)
    {
        localRoughnessRotationTensorPtr_ =
            new areaTensorField(*nm.localRoughnessRotationTensorPtr_);
    }

    if (nm.shearStressSlidingFactorsTablePtr_)
    {
        shearStressSlidingFactorsTablePtr_ =
            new interpolationTable<scalar>
            (
                *nm.shearStressSlidingFactorsTablePtr_
            );
    }

    if (nm.filmShearStressSlidingFactorPtr_)
    {
        filmShearStressSlidingFactorPtr_ =
            new areaScalarField(*nm.filmShearStressSlidingFactorPtr_);
    }

    if (nm.filmShearStressPtr_)
    {
        filmShearStressPtr_ =
            new areaVectorField(*nm.filmShearStressPtr_);
    }

    if (nm.alphaPtr_)
    {
        alphaPtr_ = new areaScalarField(*nm.alphaPtr_);
    }

    if (nm.filmTemperaturePtr_)
    {
        filmTemperaturePtr_ = new areaScalarField(*nm.filmTemperaturePtr_);
    }

    if (nm.asperityPressurePtr_)
    {
        asperityPressurePtr_ = new areaVectorField(*nm.asperityPressurePtr_);
    }

    if (nm.asperityContactAreaPtr_)
    {
        asperityContactAreaPtr_ = new areaScalarField(*nm.asperityContactAreaPtr_);
    }
}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //


void lubricatedContact::correct
(
    const vectorField& slavePatchFaceNormals,
    const scalarField& slavePointPenetration,
    const vectorField& slaveDU,
    const vectorField& masterDUInterpToSlave
)
{
    // Preliminaries
    const fvMesh& mesh = mesh_;

    // Create Finite Area mesh and fields
    if (!aMeshPtr_)
    {
        createFiniteArea();
        calcLocalRoughnessRotationTensor();
//        calcShearStressSlidingFactorsTable();
    }

    // Check if current time step and skip the first corrector
    if (curTime_ != mesh.time().value())
    {
        curTime_ = mesh.time().value();
        curCorrector_ = 0;

        return;
    }
    else
    {
        curCorrector_++;
    }

    // Make references to Finite Area pointers
    const faMesh& aMesh = *aMeshPtr_;

    areaScalarField& filmHydrodynamicPressure = *filmHydrodynamicPressurePtr_;
    areaScalarField& filmDensity = *filmDensityPtr_;
    areaScalarField& filmViscosity = *filmViscosityPtr_;
    areaScalarField& filmBulkModulus = *filmBulkModulusPtr_;
    areaScalarField& alpha = *alphaPtr_;
    areaScalarField& filmMeanThickness = *filmMeanThicknessPtr_;
    areaScalarField& filmPressureFlowFactor = *filmPressureFlowFactorPtr_;
    areaScalarField& filmShearFlowFactor = *filmShearFlowFactorPtr_;
    areaScalarField& surfaceSeparation = *surfaceSeparationPtr_;
    areaVectorField& asperityPressure = *asperityPressurePtr_;
    areaScalarField& asperityContactArea = *asperityContactAreaPtr_;


    // Make references to interpolation table pointers
    interpolationTable<scalar>& separationVsAsperityPressureTable =
        *separationVsAsperityPressureTablePtr_;
    interpolationTable<scalar>& separationVsAsperityContactAreaTable =
        *separationVsAsperityContactAreaTablePtr_;

    // Relax point penetration
    relaxPointPenetration(slavePointPenetration);

    // Interpolate fluid slave separation (penetration) from point to face
    // values
    primitivePatchInterpolation localSlaveInterpolator
    (
        mesh.boundaryMesh()[slavePatchID()]
    );

    // NOTE:
    // Surface separation represents distance between mean plane of surface
    // heights:
    // Surface separation --> h = d + yS
    surfaceSeparation.internalField() =
        localSlaveInterpolator.pointToFaceInterpolate<scalar>
        (
            fluidSlavePointPenetration_
        );
    surfaceSeparation.correctBoundaryConditions();

    // Calculate point distances for the slave zone and interpolate to face
    // centres. These distances (separations) are used to interpolate asperity
    // contact pressures and areas.
    forAll(asperityPressure, faceI)
    {
        asperityPressure[faceI] =
            separationVsAsperityPressureTable(surfaceSeparation[faceI])
           *(-slavePatchFaceNormals[faceI]);

        asperityContactArea[faceI] =
            separationVsAsperityContactAreaTable(surfaceSeparation[faceI]);
    }

    // Skip lubrication if it is not turned on yet
    if ((mesh.time().value() + SMALL) < lubricationStartTime_)
    {
        slavePressure_ = asperityPressure.internalField();

        return;
    }

    // Calculate average film flow thickness
    calcFilmMeanThickness();

    // Locate faces out of contact
    dynamicLabelList nonContactFaces;
    dynamicLabelList contactFaces;

    forAll(filmMeanThickness, faceI)
    {
        if (filmMeanThickness[faceI] >= maxFilmMeanThickness_.value())
        {
            nonContactFaces.append(faceI);
        }
        else
        {
            contactFaces.append(faceI);
        }
    }

    // Locate contact faces which are not of interest
    dynamicLabelList outOfInterestFaces;
    const areaVectorField& areaCentres = aMesh.areaCentres();

    if (boxesOfInterest_.size())
    {
        forAll(contactFaces, faceI)
        {
            const label curFace = contactFaces[faceI];
            bool faceOfInterest = false;

            forAll(boxesOfInterest_, boxI)
            {
                const boundBox& curBoundBox = boxesOfInterest_[boxI];

                if (curBoundBox.contains(areaCentres[curFace]))
                {
                    faceOfInterest = true;
                    break;
                }
            }

            if (!faceOfInterest)
            {
                outOfInterestFaces.append(faceI);
            }
        }
    }

    // Create film contact field (0 for faces not in contact with lubricant, 1
    // for in contact with lubricant)
    Field<label> contact(filmMeanThickness.size(), 1);

    forAll(nonContactFaces, faceI)
    {
        const label curFace = nonContactFaces[faceI];

        contact[curFace] = 0;
    }

    // Calculate mean and slip surface velocities
    tmp<areaVectorField> tFilmSlipVelocity =
        slipMeanVelocity
        (
            "slip",
            slavePatchFaceNormals,
            slaveDU,
            masterDUInterpToSlave
        );
    areaVectorField& filmSlipVelocity = tFilmSlipVelocity();

    tmp<areaVectorField> tFilmMeanVelocity =
        slipMeanVelocity
        (
            "mean",
            slavePatchFaceNormals,
            slaveDU,
            masterDUInterpToSlave
        );
    areaVectorField& filmMeanVelocity = tFilmMeanVelocity();

    tmp<areaVectorField> tSlaveVelocity =
        slipMeanVelocity
        (
            "slave",
            slavePatchFaceNormals,
            slaveDU,
            masterDUInterpToSlave
        );
    areaVectorField& slaveVelocity = tSlaveVelocity();

    // Set velocity into faces which are disregarded by GGI due to quick search
    // algorithm, but have active faces (with measured distances by GGI) as
    // their neighbours. The inserted velocities are calculated as mean velocity
    // of the neighbour faces.
    calcContactBoundaryVelocities
    (
        nonContactFaces, contact, filmMeanVelocity, filmSlipVelocity
    );

    // Calculate flow factors
    calcPressureFlowFactors();
    calcShearFlowFactors();
//    calcShearStressFactors();


    // Outer lubrication loop
    for (int lCorr = 0; lCorr <= maxLubricationCorr_; ++lCorr)
    {
        if (debugL_)
        {
            Info << nl << nl << "Outer corrector " << lCorr << endl;
        }

        // Calculate alpha field
        const scalar cavDensity = cavDensityL_.value();

        alpha.internalField() = pos(filmDensity.internalField() - cavDensity);
        alpha.boundaryField() == pos(filmDensity.boundaryField() - cavDensity);

        // Calculate film viscosity
        if (calcViscosityL_)
        {
            calcFilmViscosity(filmSlipVelocity);
        }

        // Calculate film bulk modulus
        calcFilmBulkModulus();


        // Calculate coefficients required for solving Reynolds equation

        // Calculate gamma with Wilson and Marsault pressure flow factors
        areaScalarField gamma =
            filmBulkModulus*filmPressureFlowFactor*pow3(filmMeanThickness)
           /(12*filmViscosity);

//        // Calculate gamma with Patir and Cheng pressure flow factors
//        areaScalarField gamma =
//            filmBulkModulus*filmPressureFlowFactor*pow3(surfaceSeparation)
//           /(12*filmViscosity);

        // Calculate gamma edge field
        edgeScalarField gammaEdge =
            fac::interpolate(gamma, "interpolate(gamma)");

        // Calculate cavitation boundary contribution
        edgeScalarField alphaEdge
        (
            "alphaEdge", fac::interpolate(alpha, "interpolate(alpha)")
        );
        scalarField cavBoundaryDiagonal(filmDensity.size(), 0.0);
        scalarField cavBoundarySource(filmDensity.size(), 0.0);

        calcCavitationBoundary
        (
            alphaEdge,
            gammaEdge,
            cavBoundaryDiagonal,
            cavBoundarySource,
            filmMeanVelocity
        );

        // Poiseuille coefficient
        edgeScalarField poiseuilleCoeff = alphaEdge*gammaEdge;

        // Couette coefficient
        edgeScalarField couetteCoeff =
            fac::interpolate(filmMeanVelocity*filmMeanThickness)
          & aMesh.Le();

        // Surface slip coefficient
        edgeScalarField slipCoeff =
            fac::interpolate
            (
                filmSlipVelocity/2*compositeSurfaceRoughness_
               *filmShearFlowFactor
            ) & aMesh.Le();

        // Moving mesh coefficient
        edgeScalarField movingMeshCoeff =
            fac::interpolate(slaveVelocity*filmMeanThickness)
          & aMesh.Le();


        // Solve the Reynolds equation
        if (debugL_)
        {
            blockLduMatrix::debug = 1;
        }

        // Density equation
        faScalarMatrix rhoEqn
        (
          - alpha*fam::laplacian
            (
                poiseuilleCoeff, filmDensity, "laplacian(poiseuilleCoeff,rho)"
            )
          + fam::div
            (
//                couetteCoeff + slipCoeff - movingMeshCoeff,
                couetteCoeff + slipCoeff,
                filmDensity,
                "div(couetteCoeff+slipCoeff,rho)"
            )
          + fam::ddt(filmMeanThickness, filmDensity)
        );

        // Add diagonal and source contributions from cavitation boundary
        rhoEqn.diag() += cavBoundaryDiagonal;
        rhoEqn.source() += cavBoundarySource;

        // Set values for faces not in contact and faces out of interest
        dynamicLabelList setValueFaces = nonContactFaces;

        if (outOfInterestFaces.size())
        {
            setValueFaces.append(outOfInterestFaces);
        }

        rhoEqn.setValues
        (
            setValueFaces,
            scalarField(setValueFaces.size(), refDensityL_.value())
        );

        // Solve the density equation
        rhoEqn.solve();

        if (debugL_)
        {
            blockLduMatrix::debug = 0;
        }

        // Bound density explicitly in order to stabilize convergence
        forAll(filmDensity, faceI)
        {
            if (alpha[faceI] == 0.0)
            {
                filmDensity[faceI] = min(filmDensity[faceI], cavDensity);
            }

            filmDensity = max(filmDensity, dimensionedScalar("zeroDensity", dimDensity, 0.0));
        }

        // Store old hydrodynamic pressure
        areaScalarField oldHydrodynamicPressure = filmHydrodynamicPressure;

        // Calculate film pressure
        filmHydrodynamicPressure = filmPressure(filmDensity);

        // Limit film pressure due to stability
        filmHydrodynamicPressure =
            min
            (
                filmHydrodynamicPressure,
                dimensionedScalar("maxPressure", dimPressure, maxFilmPressure_)
            );

        // Calculate pressure relative residual
        scalar pResidual =
            max
            (
                mag
                (
                    (
                        oldHydrodynamicPressure.internalField()
                      - filmHydrodynamicPressure.internalField()
                    )/max(mag(oldHydrodynamicPressure.internalField()) + SMALL)
                )
            );

        reduce(pResidual, maxOp<scalar>());

        if (debugL_)
        {
            Info << "Max residual = " << pResidual << endl;
            Info << "Max pressure = " << gMax(filmHydrodynamicPressure)
                 << endl;
            Info << "Min thickness = " << gMin(filmMeanThickness)
                 << endl;
            Info << "sum(alpha) = " << gSum(alpha) << endl;
            Info << "min(filmDensity) = " << gMin(filmDensity) << endl;
            Info << "max(filmDensity) = " << gMax(filmDensity) << endl;
            Info << "min(filmPressureFlowFactor) = "
                 << gMin(filmPressureFlowFactor) << endl;
        }

        // Exit if pressure residual is below treshold
        if (pResidual < 1e-08)
        {
            break;
        }
    }

    // Calculate slave pressure
    if (rigidContact_)
    {
        slavePressure_ = vector::zero;
    }
    else
    {
        slavePressure_ =
            asperityPressure.internalField()
          + (
                (1.0 - asperityContactArea.internalField())
               *filmHydrodynamicPressure.internalField()
            )*(-slavePatchFaceNormals);
    }

    // Calculate film shear stress without shear stress factors
    areaVectorField& filmShearStress = *filmShearStressPtr_;

    filmShearStress =
        filmViscosity*(-filmSlipVelocity)/filmMeanThickness;
    filmShearStress.internalField() *=
        1.0 - asperityContactArea.internalField();

    // Set shear stress to 0.0 for faces which are not in contact with the
    // lubricant
    forAll(nonContactFaces, faceI)
    {
        const label curFace = nonContactFaces[faceI];
        filmShearStress.internalField()[curFace] = vector::zero;
    }
    filmShearStress.correctBoundaryConditions();

//    // Calculate film shear stress with shear stress factors
//    areaVectorField& filmShearStress = *filmShearStressPtr_;
//    areaScalarField& filmShearStressFactor = *filmShearStressFactorPtr_;
//    areaScalarField& filmShearStressSlidingFactor =
//        *filmShearStressSlidingFactorPtr_;
//
//    filmShearStress =
//        filmViscosity*(-filmSlipVelocity)/surfaceSeparation
//       *(
//            filmShearStressSlidingFactor
//          + filmShearStressFactor
//           *(2*sqr(slaveSurfaceRoughness_/compositeSurfaceRoughness_) - 1.0)
//        );
//    filmShearStress.internalField() *=
//        1.0 - asperityContactArea.internalField();
//    filmShearStress.correctBoundaryConditions();
}


void lubricatedContact::writeDict(Ostream& os) const
{
    word keyword(name() + "NormalModelDict");
    os.writeKeyword(keyword)
        << normalContactModelDict_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
