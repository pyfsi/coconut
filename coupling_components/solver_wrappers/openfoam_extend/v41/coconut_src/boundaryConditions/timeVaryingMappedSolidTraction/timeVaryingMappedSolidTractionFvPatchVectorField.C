/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "timeVaryingMappedSolidTractionFvPatchVectorField.H"
#include "foamTime.H"
#include "triSurfaceTools.H"
#include "triSurface.H"
#include "vector2D.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"
#include "fvc.H"
#include "surfaceInterpolationScheme.H"
#include "fvBlockMatrix.H"
#include "IOReferencer.H"
#include "vectorIOField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

timeVaryingMappedSolidTractionFvPatchVectorField::
timeVaryingMappedSolidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    pressureFieldTableName_("pressure"),
    tractionFieldTableName_("traction"),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    referenceCS_(nullptr),
    nearestVertex_(0),
    nearestVertexWeight_(0),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledPressureValues_(0),
    startSampledTractionValues_(0),
    endSampleTime_(-1),
    endSampledPressureValues_(0),
    endSampledTractionValues_(0)
{
   fvPatchVectorField::operator=(patchInternalField());
   gradient() = vector::zero;
}


timeVaryingMappedSolidTractionFvPatchVectorField::
timeVaryingMappedSolidTractionFvPatchVectorField
(
    const timeVaryingMappedSolidTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
    pressureFieldTableName_("pressure"),
    tractionFieldTableName_("traction"),
    traction_(p.size(), vector::zero),
    pressure_(stpvf.pressure_, mapper),
    referenceCS_(nullptr),
    nearestVertex_(0),
    nearestVertexWeight_(0),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledPressureValues_(0),
    startSampledTractionValues_(0),
    endSampleTime_(-1),
    endSampledPressureValues_(0),
    endSampledTractionValues_(0)
{}


timeVaryingMappedSolidTractionFvPatchVectorField::
timeVaryingMappedSolidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    pressureFieldTableName_("pressure"),
    tractionFieldTableName_("traction"),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    referenceCS_(nullptr),
    nearestVertex_(0),
    nearestVertexWeight_(0),
    sampleTimes_(0),
    startSampleTime_(-1),
    startSampledPressureValues_(0),
    startSampledTractionValues_(0),
    endSampleTime_(-1),
    endSampledPressureValues_(0),
    endSampledTractionValues_(0)
{
    	if (debug)
    	{
		Info<< patch().name() << ":" << type() << endl;
    	}
	
        //Info << "GRADIENT ########### " << dict.found("gradient") << nl << endl;

	if (dict.found("gradient"))
	{
		gradient() = vectorField("gradient", dict, p.size());
	}
	else
	{
		gradient() = vector::zero;
	}

	if (dict.found("value"))
	{
		Field<vector>::operator==(vectorField("value", dict, p.size()));
	}
	else
	{
		fvPatchVectorField::operator=(patchInternalField());
	}

	if(dict.found("findSampleTimeMultiple"))
	{
		//Info << "TRUE###########" << endl;
		dict.lookup("findSampleTimeMultiple") >> findSampleTimeMultiple_;
		//Info << "sampledTime " << findSampleTimeMultiple_ << endl;
	}
}


timeVaryingMappedSolidTractionFvPatchVectorField::
timeVaryingMappedSolidTractionFvPatchVectorField
(
const timeVaryingMappedSolidTractionFvPatchVectorField& stpvf
)
:
fixedGradientFvPatchVectorField(stpvf),
pressureFieldTableName_("pressure"),
tractionFieldTableName_("traction"),
traction_(stpvf.traction_),
pressure_(stpvf.pressure_),
referenceCS_(stpvf.referenceCS_),
nearestVertex_(stpvf.nearestVertex_),
nearestVertexWeight_(stpvf.nearestVertexWeight_),
sampleTimes_(stpvf.sampleTimes_),
startSampleTime_(stpvf.startSampleTime_),
startSampledPressureValues_(stpvf.startSampledPressureValues_),
startSampledTractionValues_(stpvf.startSampledTractionValues_),
endSampleTime_(stpvf.endSampleTime_),
endSampledPressureValues_(stpvf.endSampledPressureValues_),
endSampledTractionValues_(stpvf.endSampledTractionValues_)
{}

timeVaryingMappedSolidTractionFvPatchVectorField::
timeVaryingMappedSolidTractionFvPatchVectorField
(
const timeVaryingMappedSolidTractionFvPatchVectorField& stpvf,
const DimensionedField<vector, volMesh>& iF
)
:
fixedGradientFvPatchVectorField(stpvf, iF),
pressureFieldTableName_("pressure"),
tractionFieldTableName_("traction"),
traction_(stpvf.traction_),
pressure_(stpvf.pressure_),
referenceCS_(stpvf.referenceCS_),
nearestVertex_(stpvf.nearestVertex_),
nearestVertexWeight_(stpvf.nearestVertexWeight_),
sampleTimes_(stpvf.sampleTimes_),
startSampleTime_(stpvf.startSampleTime_),
startSampledPressureValues_(stpvf.startSampledPressureValues_),
startSampledTractionValues_(stpvf.startSampledTractionValues_),
endSampleTime_(stpvf.endSampleTime_),
endSampledPressureValues_(stpvf.endSampledPressureValues_),
endSampledTractionValues_(stpvf.endSampledTractionValues_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void timeVaryingMappedSolidTractionFvPatchVectorField::autoMap
(
const fvPatchFieldMapper& m
)
{
fixedGradientFvPatchVectorField::autoMap(m);
if (startSampledPressureValues_.size() > 0)
    {
      startSampledPressureValues_.autoMap(m);
      endSampledPressureValues_.autoMap(m);
    }
if (startSampledTractionValues_.size() > 0)
    {
      startSampledTractionValues_.autoMap(m);
      endSampledTractionValues_.autoMap(m);
    }
}


void timeVaryingMappedSolidTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const timeVaryingMappedSolidTractionFvPatchVectorField& dmptf =
        refCast<const timeVaryingMappedSolidTractionFvPatchVectorField>(ptf);

    startSampledPressureValues_.rmap(dmptf.startSampledPressureValues_, addr);
    endSampledPressureValues_.rmap(dmptf.endSampledPressureValues_, addr);
    startSampledTractionValues_.rmap(dmptf.startSampledTractionValues_, addr);
    endSampledTractionValues_.rmap(dmptf.endSampledTractionValues_, addr);
}


void timeVaryingMappedSolidTractionFvPatchVectorField::readSamplePoints()
{
    // Read the sample points

    pointIOField samplePoints
    (
        IOobject
        (
            "points",
            this->db().time().constant(),
            "boundaryData"/this->patch().name(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            false
        )
    );

    const fileName samplePointsFile = samplePoints.filePath();

    if (debug)
    {
        Info<< "timeVaryingMappedSolidTractionFvPatchVectorField :"
            << " Read " << samplePoints.size() << " sample points from "
            << samplePointsFile << endl;
    }

    // Determine coordinate system from samplePoints

    if (samplePoints.size() < 3)
    {
        FatalErrorIn
        (
            "timeVaryingMappedSolidTractionFvPatchVectorField<Type>::readSamplePoints()"
        )   << "Only " << samplePoints.size() << " points read from file "
            << samplePoints.objectPath() << nl
            << "Need at least three non-colinear samplePoints"
            << " to be able to interpolate."
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    const point& p0 = samplePoints[0];

    // Find point separate from p0
    vector e1;
    label index1 = -1;

    for (label i = 1; i < samplePoints.size(); i++)
    {
        e1 = samplePoints[i] - p0;

        scalar magE1 = mag(e1);

        if (magE1 > SMALL)
        {
            e1 /= magE1;
            index1 = i;
            break;
        }
    }

    // Find point that makes angle with n1
    label index2 = -1;
    vector e2;
    vector n;

    for (label i = index1+1; i < samplePoints.size(); i++)
    {
        e2 = samplePoints[i] - p0;

        scalar magE2 = mag(e2);

        if (magE2 > SMALL)
        {
            e2 /= magE2;

            n = e1^e2;

            scalar magN = mag(n);

            if (magN > SMALL)
            {
                index2 = i;
                n /= magN;
                break;
            }
        }
    }

    if (index2 == -1)
    {
        FatalErrorIn
        (
            "timeVaryingMappedSolidTractionFvPatchVectorField<Type>::readSamplePoints()"
        )   << "Cannot find points that make valid normal." << nl
            << "Need at least three sample points which are not in a line."
            << "\n    on patch " << this->patch().name()
            << " of points " << samplePoints.name()
            << " in file " << samplePoints.objectPath()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "timeVaryingMappedSolidTractionFvPatchVectorField :"
            << " Used points " << p0 << ' ' << samplePoints[index1]
            << ' ' << samplePoints[index2]
            << " to define coordinate system with normal " << n << endl;
    }

    referenceCS_.reset
    (
        new coordinateSystem
        (
            "reference",
            p0,  // origin
            n,   // normal
            e1   // 0-axis
        )
    );

    tmp<vectorField> tlocalVertices
    (
        referenceCS_().localPosition(samplePoints)
    );
    const vectorField& localVertices = tlocalVertices();

    // Determine triangulation
    List<vector2D> localVertices2D(localVertices.size());
    forAll(localVertices, i)
    {
        localVertices2D[i][0] = localVertices[i][0];
        localVertices2D[i][1] = localVertices[i][1];
    }

    triSurface s(triSurfaceTools::delaunay2D(localVertices2D));

    tmp<pointField> localFaceCentres
    (
        referenceCS_().localPosition
        (
            this->patch().patch().faceCentres()
        )
    );

    if (debug)
    {
        Pout<< "readSamplePoints :"
            <<" Dumping triangulated surface to triangulation.stl" << endl;
        s.write(this->db().time().path()/"triangulation.stl");

        OFstream str(this->db().time().path()/"localFaceCentres.obj");
        Pout<< "readSamplePoints :"
            << " Dumping face centres to " << str.name() << endl;

        forAll(localFaceCentres(), i)
        {
            const point& p = localFaceCentres()[i];
            str<< "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
        }
    }

    // Determine interpolation onto face centres.
    triSurfaceTools::calcInterpolationWeights
    (
        s,
        localFaceCentres,   // points to interpolate to
        nearestVertex_,
        nearestVertexWeight_
    );



    // Read the times for which data is available

    const fileName samplePointsDir = samplePointsFile.path();

    sampleTimes_ = Time::findTimes(samplePointsDir);

    if (debug)
    {
        Info<< "timeVaryingMappedSolidTractionFvPatchVectorField : In directory "
            << samplePointsDir << " found times " << timeNames(sampleTimes_)
            << endl;
    }
}


wordList timeVaryingMappedSolidTractionFvPatchVectorField::timeNames
(
    const instantList& times
)
{
    wordList names(times.size());

    forAll(times, i)
    {
        names[i] = times[i].name();
    }
    return names;
}


void timeVaryingMappedSolidTractionFvPatchVectorField::findTime
(
    const fileName& instance,
    const fileName& local,
    const scalar timeVal,
    label& lo,
    label& hi
) const
{
    lo = startSampleTime_;
    hi = -1;

    for (label i = startSampleTime_+1; i < sampleTimes_.size(); i++)
    {
        if (sampleTimes_[i].value() > timeVal)
        {
            break;
        }
        else
        {
            lo = i;
        }
    }

    if (lo == -1)
    {
        FatalErrorIn("findTime")
            << "Cannot find starting sampling values for current time "
            << timeVal << nl
            << "Have sampling values for times "
            << timeNames(sampleTimes_) << nl
            << "In directory "
            <<  this->db().time().constant()/"boundaryData"/this->patch().name()
            << "\n    on patch " << this->patch().name()
            << " of field " << pressureFieldTableName_ << " and " << tractionFieldTableName_
            << exit(FatalError);
    }

    if (lo < sampleTimes_.size()-1)
    {
        hi = lo+1;
    }


    if (debug)
    {
        if (hi == -1)
        {
            Pout<< "findTime : Found time " << timeVal << " after"
                << " index:" << lo << " time:" << sampleTimes_[lo].value()
                << endl;
        }
        else
        {
            Pout<< "findTime : Found time " << timeVal << " inbetween"
                << " index:" << lo << " time:" << sampleTimes_[lo].value()
                << " and index:" << hi << " time:" << sampleTimes_[hi].value()
                << endl;
        }
    }
}


void timeVaryingMappedSolidTractionFvPatchVectorField::checkTable()
{
    // Initialise
    readSamplePoints();
    
    // Find current time in sampleTimes
    label lo = -1;
    label hi = -1;


    findTime
    (
        this->db().time().constant(),
        "boundaryData"/this->patch().name(),
        this->db().time().value(),
        lo,
        hi
    );
 
    // Update sampled data fields.
      // if (true)// if ((lo != startSampleTime_)||(findSampleTimeMultiple_))
    //{   
        startSampleTime_ = lo; 

         //if (false) //(startSampleTime_ == endSampleTime_)
        //{

            // No need to reread since are end values
           // if (true)
           // {
               // Pout<< "checkTable : Setting startValues to (already read) "
                  //  <<   "boundaryData"
                       // /this->patch().name()
                       // /sampleTimes_[startSampleTime_].name()
                   // << endl;
            //}
            //startSampledPressureValues_ = endSampledPressureValues_;
	    //startSampledTractionValues_ = endSampledTractionValues_;
        //}
       // else
        //{
            if (debug)
            {
                Pout<< "checkTable : Reading startValues from "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[lo].name()
                    << endl;
            }

       // Reread scalar values and interpolate
	  IOField<scalar> pressureValsStart
             (
                 IOobject
                 (
                     pressureFieldTableName_,
                     this->db().time().constant(),
                     "boundaryData"
                    /this->patch().name()
                    /sampleTimes_[startSampleTime_].name(),
                     this->db(),
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE,
                     false
                 )
             );

            startSampledPressureValues_ = interpolate(pressureValsStart);
           
	// Reread vector values and interpolate	    
          IOField<vector> tractionValsStart
             (
                 IOobject
                 (
                     tractionFieldTableName_,
                     this->db().time().constant(),
                     "boundaryData"
                    /this->patch().name()
                    /sampleTimes_[startSampleTime_].name(),
                     this->db(),
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE,
                     false
                 )
             );

       	    startSampledTractionValues_ = interpolate(tractionValsStart);
	//}
    //}

       // if (true) //((hi != endSampleTime_)||(findSampleTimeMultiple_))
    //{
        endSampleTime_ = hi;

        //if (endSampleTime_ == -1)
       // {
            // endTime no longer valid. Might as well clear endValues.
            //if (true)
            //{
              //  Pout<< "checkTable : Clearing endValues" << endl;
            //}
            //endSampledPressureValues_.clear();
        //}
        //else
        //{
            if (true)
            {
                Pout<< "checkTable : Reading endValues from "
                    <<   "boundaryData"
                        /this->patch().name()
                        /sampleTimes_[endSampleTime_].name()
                    << endl;
            }
          // Reread scalar values and interpolate
     	  IOField<scalar> pressureValsEnd
             (
                 IOobject
                 (
                     pressureFieldTableName_,
                     this->db().time().constant(),
                     "boundaryData"
                    /this->patch().name()
                    /sampleTimes_[endSampleTime_].name(),
                     this->db(),
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE,
                     false
                 )
             );
 
            endSampledPressureValues_ = interpolate(pressureValsEnd);

	  // Reread vector values and interpolate
	      IOField<vector> tractionValsEnd
             (
                 IOobject
                 (
                     tractionFieldTableName_,
                     this->db().time().constant(),
                     "boundaryData"
                    /this->patch().name()
                    /sampleTimes_[endSampleTime_].name(),
                     this->db(),
                     IOobject::MUST_READ,
                     IOobject::AUTO_WRITE,
                     false
                 )
             );
	
	   endSampledTractionValues_ = interpolate(tractionValsEnd);
  
        //}
    //}
}

template<class Type>
tmp<Field<Type> > timeVaryingMappedSolidTractionFvPatchVectorField::interpolate
(
    const Field<Type>& sourceFld
) const
{
    tmp<Field<Type> > tfld(new Field<Type>(nearestVertex_.size()));
    Field<Type>& fld = tfld();

    forAll(fld, i)
    {
        const FixedList<label, 3>& verts = nearestVertex_[i];
        const FixedList<scalar, 3>& w = nearestVertexWeight_[i];

        if (verts[2] == -1)
        {
            if (verts[1] == -1)
            {
                // Use vertex0 only
                fld[i] = sourceFld[verts[0]];
            }
            else
            {
                // Use vertex 0,1
                fld[i] =
                    w[0]*sourceFld[verts[0]]
                  + w[1]*sourceFld[verts[1]];
            }
        }
        else
        {
            fld[i] =
                w[0]*sourceFld[verts[0]]
              + w[1]*sourceFld[verts[1]]
              + w[2]*sourceFld[verts[2]];
        }
    }
    return tfld;
}


void timeVaryingMappedSolidTractionFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

   //Pout<<"updateCoeffs"<<endl;

    checkTable();

    // Interpolate between the sampled data

    if (endSampleTime_ == -1)
    {
        // only start value
        if (debug)
        {
      		Pout<< "updateCoeffs : Sampled, non-interpolated values"
          	<< " from start time:"
          	<< sampleTimes_[startSampleTime_].name() << nl;
        }

        pressure_=(startSampledPressureValues_);
 	traction_=(startSampledTractionValues_);
    }
    else
    {
        scalar start = sampleTimes_[startSampleTime_].value();
        scalar end = sampleTimes_[endSampleTime_].value();

        scalar s = (this->db().time().value()-start)/(end-start);

        if (true)
        {
            Pout<< "updateCoeffs : Sampled, interpolated values"
                << " between start time:"
                << sampleTimes_[startSampleTime_].name()
                << " and end time:" << sampleTimes_[endSampleTime_].name()
                << " with weight:" << s << endl;
        }

        pressure_=((1-s)*startSampledPressureValues_ + s*endSampledPressureValues_);
        traction_=((1-s)*startSampledTractionValues_ + s*endSampledTractionValues_);

    }
  
	const word fieldName = dimensionedInternalField().name();

    gradient() =
            tractionBoundaryGradient().snGrad
            (
                traction_,                  // surface traction
                pressure_,                  // surface pressure
                fieldName,                  // working field name 
                "U",                        // total field name
                patch(),                    // polyPatch
                bool(fieldName == "DU")     // incremental
	    );

    // If a block coupled hybrid pressure-displacement approach is used then we 
    // will include the pressure implicity (From solidTraction boundary condition)
    if
    (
	patch().boundaryMesh().mesh().foundObject
	<
	    IOReferencer< fvBlockMatrix<vector4> >
	>
	(
	   "DUpEqn"
	)
    )
    {
	// Lookup the block matrix
	fvBlockMatrix<vector4>& DUpEqn = 
   	    const_cast<fvBlockMatrix<vector4>&>
	    (
		patch().boundaryMesh().mesh().lookupObject
		<
		IOReferencer< fvBlockMatrix<vector4> >
		> 
		(
			"DUpEqn"
		) ()
	    );
	
	Field<tensor4>& diag = DUpEqn.diag().asSquare();
	Field<vector4>& source = DUpEqn.source();

	const unallocLabelList& faceCells = patch().faceCells();
	const vectorField& Sf = patch().Sf();

	// Lookup the pressure
	const volScalarField* pPtr = NULL;
	if (patch().boundaryMesh().mesh().foundObject<volScalarField>("Dp"))
	{
	    pPtr = 
		&patch().boundaryMesh().mesh().lookupObject<volScalarField>
		(
			"Dp"
		);
	}
	else
	{
	    pPtr = 
		&patch().boundaryMesh().mesh().lookupObject<volScalarField>
		(
			"p"
		);
	}
	const volScalarField& p = *pPtr;

	// Add coeffs to the matrix
	forAll(faceCells, faceI)
	{
	    const label cellID = faceCells[faceI];	
	
	    diag[cellID](0, 3) -= Sf[faceI][vector::X];
	    diag[cellID](1, 3) -= Sf[faceI][vector::Y];
	    diag[cellID](2, 3) -= Sf[faceI][vector::Z];

	    source[cellID][0] -= p[cellID]*Sf[faceI][vector::X];
	    source[cellID][1] -= p[cellID]*Sf[faceI][vector::Y];
	    source[cellID][2] -= p[cellID]*Sf[faceI][vector::Z];
	}
    }

    fixedGradientFvPatchVectorField::updateCoeffs();
}

void timeVaryingMappedSolidTractionFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
	this->updateCoeffs();
    }
  
    if
    (
 	db().foundObject<volTensorField>
	(
	    "grad(" + dimensionedInternalField().name() + ")"
	)
    )
    {
	// Lookup the gradient field
	const fvPatchField<tensor>& gradField = 
	    patch().lookupPatchField<volTensorField, tensor>
	    (
		"grad(" + dimensionedInternalField().name() + ")"
	    );
     
        // Patch unit normals
        const vectorField n = patch().nf();

        // Patch delta vectors
        const vectorField delta = patch().delta();

        // Non-orthogonal correction vectors
        const vectorField k = (I - sqr(n)) & delta;

	//Set value with non-orthogonal correction
	Field<vector>::operator=
	(
	    patchInternalField()
          + (k & gradField.patchInternalField())
          + gradient()/patch().deltaCoeffs()
        );
    }
    else
    {
	// Set value without non-orthogonal correction
	Field<vector>::operator= 
	(
	    patchInternalField()
          + gradient()/patch().deltaCoeffs()
	);
    }
 
    fvPatchField<vector>::evaluate();
 
}

void timeVaryingMappedSolidTractionFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);

    if (pressureFieldTableName_ != this->dimensionedInternalField().name())
    {
        os.writeKeyword("pressureFieldTableName") << pressureFieldTableName_ << token::END_STATEMENT << nl;
        pressure_.writeEntry("pressure", os);
    }
    if (tractionFieldTableName_ != this->dimensionedInternalField().name())
    {
        os.writeKeyword("tractionFieldTableName") << tractionFieldTableName_ << token::END_STATEMENT << nl;
        traction_.writeEntry("traction", os);
    }

    this->writeEntry("value", os);
  
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, timeVaryingMappedSolidTractionFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
