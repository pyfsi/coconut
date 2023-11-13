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

Class
    lubricatedContact

Description
    Normal contact model which takes into account asperity contact pressure and
    lubrication.
    Asperity contact pressure depends on distance between master and slave, and
    is interpolated from tables.
    Lubricant film pressure is calculated using the Reynolds lubrication
    equation.

SourceFiles
    lubricatedContact.C

Author
    Vanja Skuric, UNIZAG FSB

\*---------------------------------------------------------------------------*/

#ifndef lubricatedContact_H
#define lubricatedContact_H

#include "normalContactModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "faCFD.H"
#include "interpolationTable.H"
#include "FIFOStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class lubricatedContact Declaration
\*---------------------------------------------------------------------------*/

class lubricatedContact
:
    public normalContactModel
{
    // Private data

        //- Contact law dictionary
        dictionary normalContactModelDict_;

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Slave faces zone patch
        const standAlonePatch& slaveFaceZonePatch_;

        //- Start time
        const scalar startTime_;

        //- Current time
        scalar curTime_;

        //- Current displacement corrector
        label curCorrector_;

        //- Slave pressure
        vectorField slavePressure_;

        //- Under-relaxation factor
        const scalar relaxFac_;

        //- Standard penalty contact model
        normalContactModel* standardPenaltyPtr_;

        //- Lubrication dictionary
        dictionary lubricationDict_;

        //- Maximum number of lubrication correctors
        const label maxLubricationCorr_;

        //- Lubrication start time
        const scalar lubricationStartTime_;

        //- Surface rouhness dictionary
        IOdictionary surfaceRoughnessDict_;

        //- Const reference to FA mesh
        faMesh* aMeshPtr_;

        //- Surface separation represents distance between mean  plane of
        // asperity heights and opposite flat surface
        // Surface separation --> d = h - yS
        areaScalarField* surfaceSeparationPtr_;

        //- Film hydrodynamic pressure (gauge pressure)
        areaScalarField* filmHydrodynamicPressurePtr_;

        //- Film density (liquid)
        areaScalarField* filmDensityPtr_;

        //- Film viscosity (liquid)
        areaScalarField* filmViscosityPtr_;

        //- Film bulk modulus
        areaScalarField* filmBulkModulusPtr_;

        //- Film mean thickness
        areaScalarField* filmMeanThicknessPtr_;

        //- Composite surface roughness Rq (standard deviation of asperity
        //- heights, sigma)
        const dimensionedScalar compositeSurfaceRoughness_;

        //- Master surface roughness
        const dimensionedScalar masterSurfaceRoughness_;

        //- Slave surface roughness
        const dimensionedScalar slaveSurfaceRoughness_;

        //- Maximum film mean thickness
        const dimensionedScalar maxFilmMeanThickness_;

        //- Minimum film mean thickness
        const dimensionedScalar minFilmMeanThickness_;

        //- Lubricant properties reference temperature
        const dimensionedScalar refTemperatureL_;

        //- Lubricant viscosity at reference temperature and at hydrodynamic
        //- pressure equal to zero
        const dimensionedScalar refViscosityL_;

        //- Lubricant density at reference temperature and at hydrodynamic
        //- pressure equal to zero
        const dimensionedScalar refDensityL_;

        //- Lubricant density at cavitating pressure and at reference
        //  temperature
        dimensionedScalar cavDensityL_;

        //- Lubricant cavitation pressure
        const dimensionedScalar cavPressureL_;

        //- Lubricant pressure exponent
        const scalar pressureExpL_;

        //- Lubricant temperature exponent
        const scalar temperatureExpL_;

        //- Dowson-Higginson contsant C1
        const scalar C1_;

        //- Dowson-Higginson contsant C2
        const scalar C2_;

        //- Should I calculate viscosity or not?
        const bool calcViscosityL_;

        //- Non-Newtonian lubricant switch
        const bool nonNewtonianL_;

        //- Lubricant Eyring stress
        dimensionedScalar eyringStressL_;

        //- Separation vs aperity pressure table
        interpolationTable<scalar>* separationVsAsperityPressureTablePtr_;

        //- Asperity pressure vs separation table
        interpolationTable<scalar>* asperityPressureVsSeparationTablePtr_;

        //- Separation vs aperity contact area table
        interpolationTable<scalar>* separationVsAsperityContactAreaTablePtr_;

        //- Ratio of composite surface roughness correlation lengths (gamma)
        //- Patir and Chang, 1978
        const scalar compositeCorrLengthRatio_;

        //- Ratio of master surface roughness correlation lengths
        const scalar masterCorrLengthRatio_;

        //- Ratio of slave surface roughness correlation lengths
        const scalar slaveCorrLengthRatio_;

        //- Film pressure flow factor tensor
        areaScalarField* filmPressureFlowFactorPtr_;

        //- Film shear flow factor
        areaScalarField* filmShearFlowFactorPtr_;

        //- Film shear stress factor
        areaScalarField* filmShearStressFactorPtr_;

        //- Film shear stress sliding factor (phi_f)
        areaScalarField* filmShearStressSlidingFactorPtr_;

        //- Film shear stress
        areaVectorField* filmShearStressPtr_;

        //- Binary field used for distinguishing active (1.0) from cavitating
        //  (0.0) film area
        areaScalarField* alphaPtr_;

        //- Film temperature
        areaScalarField* filmTemperaturePtr_;

        //- Asperity contact pressure
        areaVectorField* asperityPressurePtr_;

        //- Asperity contact area ratio
        areaScalarField* asperityContactAreaPtr_;

        //- Initial roughness axial direction in global coordinate system
        const vector axialRoughnessDirection_;

        //- Local coordinate system roughness direction rotation tensor
        areaTensorField* localRoughnessRotationTensorPtr_;

        //- Two vector fields designatig surface roughness directions in the
        //- global coordinate system:
        //- [0] - axial roughness direction
        //- [1] - circumferential roughness direction
        Field<vectorField> surfaceRoughnessDirections_;

        //- Shear stress sliding factor interpolation table
        interpolationTable<scalar>* shearStressSlidingFactorsTablePtr_;

        //- Maximum film hydrodynamic pressure
        scalar maxFilmPressure_;

        //- Debug mode?
        const bool debugL_;

        //- Rigid contact (contact pressure going back to the solver is zero)
        const bool rigidContact_;

        //- Relaxed slave point penetration going into Reynolds equation
        scalarField fluidSlavePointPenetration_;

        //- Processor patch ID list
        labelList processorPatchId_;

        //- List of bounding boxes of interest
        List<boundBox> boxesOfInterest_;


    // Private Member Functions

        //- Create Finite Area mesh and fields
        void createFiniteArea();

        //- Create interpolation tables
        void createInterpolationTables();

        //- Calculate local axial roughness directions
        void calcLocalRoughnessRotationTensor();

        //- Calculate mean film thickness
        void calcFilmMeanThickness();

        //- Calculate slip or mean surface velocity
        tmp<areaVectorField> slipMeanVelocity
        (
            const word& velocityType,
            const vectorField& slavePatchFaceNormals,
            const vectorField& slaveDU,
            const vectorField& masterDUInterpToSlave
        );

        //- Calculate pressure flow factors
        void calcPressureFlowFactors();

        //- Calculate shear flow factors
        void calcShearFlowFactors();

        //- Calculate shear stress factors
        void calcShearStressFactors();

        //- Calculate shear stress sliding factors
        void calcShearStressSlidingFactors();

        //- Calculate shear stress factors interpolation table due to averaging
        //  the sliding velocity
        void calcShearStressSlidingFactorsTable();

        //- Calculate lubricant viscosity
        void calcFilmViscosity
        (
            const areaVectorField& filmSlipVelocity
        );

        //- Calculate film density
        void calcFilmDensity();

        //- Calculate lubricant gauge pressure
        areaScalarField filmPressure(const areaScalarField& filmDensity);

        //- Calculate lubricant bulk modulus
        void calcFilmBulkModulus();

        //- Relax point penetration field
        void relaxPointPenetration
        (
            const scalarField& slavePointPenetration
        );

        //- Calculate valocities in non-contact faces which are located right
        //- next to the the contact boundary
        void calcContactBoundaryVelocities
        (
            const dynamicLabelList& nonContactFaces,
            const Field<label>& contact,
            areaVectorField& filmMeanVelocity,
            areaVectorField& filmSlipVelocity
        );

        //- Calculate cavitation boundary contribution
        void calcCavitationBoundary
        (
            edgeScalarField& alphaEdge,
            const edgeScalarField& gammaEdge,
            scalarField& cavBoundaryDiagonal,
            scalarField& cavBoundarySource,
            const areaVectorField& filmMeanVelocity
        );

        //- Disallow copy construct
        //lubricatedContact(const lubricatedContact&);

        //- Disallow default bitwise assignment
        void operator=(const lubricatedContact&);



public:

    //- Runtime type information
    TypeName("lubricatedContact");


    // Constructors

        //- Construct from dictionary
        lubricatedContact
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID,
            const standAlonePatch& masterFaceZonePatch,
            const standAlonePatch& slaveFaceZonePatch
        );

        //- Construct as copy
        lubricatedContact(const lubricatedContact&);

        //- Construct and return a clone
        virtual autoPtr<normalContactModel> clone() const
        {
            return autoPtr<normalContactModel>(new lubricatedContact(*this));
        }


    // Destructor

        virtual ~lubricatedContact()
        {}


    // Member Functions


        //- Correct contact model
        virtual void correct
        (
            const vectorField& slavePatchFaceNormals,
            const scalarField& slavePointPenetration,
            const vectorField& slaveDU,
            const vectorField& masterDUInterpToSlave
        );

        //- Return slave pressure
        virtual const vectorField& slavePressure() const
        {
            return slavePressure_;
        }

        //- Return slave pressure
        virtual vectorField& slavePressure()
        {
            return slavePressure_;
        }

        //- Return area in contact of slave patch
        virtual const scalarField& areaInContact() const
        {
            return asperityContactAreaPtr_->internalField();
        };

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
