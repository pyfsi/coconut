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

\*----------------------------------------------------------------------------*/

#include "solidPointStress.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidPointStress, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        solidPointStress,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solidPointStress::writeData()
{
    if (pointAndFieldFound_)
    {
        const fvMesh& mesh = time_.lookupObject<fvMesh>("region0");

        if (mesh.foundObject<volSymmTensorField>(fieldName_))
        {
            const volSymmTensorField& field =
                mesh.lookupObject<volSymmTensorField>(fieldName_);

            pointMesh pMesh(mesh);

            // pointSymmTensorField pointField
            //     (
            //         IOobject
            //         (
            //             "point" + fieldName_,
            //             mesh.time().timeName(),
            //             mesh,
            //             IOobject::NO_READ,
            //             IOobject::AUTO_WRITE
            //         ),
            //         pMesh,
            //         dimensionedSymmTensor
            //         (
            //             "zero", field.dimensions(), symmTensor::zero
            //         )
            //     );

            // Interpolate
            //interpPtr_->interpolate(field, pointField);

            // Values are wrong on a symmetryPlane because transform sets wrong
            // components to zero
            // We can do it component-wise

            volScalarField fieldyy
                (
                    IOobject
                    (
                        fieldName_+ "yy",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("zero", field.dimensions(), 0.0)
                );

            fieldyy = field.component(symmTensor::YY);

            pointScalarField pointFieldyy
                (
                    IOobject
                    (
                        "point" + fieldName_ + "yy",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    pMesh,
                    dimensionedScalar("zero", field.dimensions(), 0.0)
                );

            interpPtr_->interpolate(fieldyy, pointFieldyy);

            if (Pstream::master())
            {
                historyFilePtr_()
                    << time_.time().value()
                    << " " << pointFieldyy[pointID_] << endl;

            }
        }
        else
        {
            InfoIn(this->name() + " function object constructor")
                << "pointScalarField " << fieldName_ << " not found" << endl;
        }
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidPointStress::solidPointStress
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    fieldName_("unspecified"),
    pointID_(-1),
    pointAndFieldFound_(true),
    interpPtr_(NULL),
    historyFilePtr_(NULL)
{
    Info<< "Creating " << this->name() << " function object" << endl;

    vector point = vector::zero;
    if (dict.found("point"))
    {
        point = vector(dict.lookup("point"));
    }
    else
    {
        pointAndFieldFound_ = false;
        WarningIn(this->name() + " function object constructor")
            << "solidPointStress: point not specified" << endl;
    }

    if (dict.found("field"))
    {
        dict.lookup("field") >> fieldName_;
    }
    else
    {
        pointAndFieldFound_ = false;
        WarningIn(this->name() + " function object constructor")
            << "solidPointStress: fieldName not specified" << endl;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>("region0");

    // Create history file if not already created
    if (historyFilePtr_.empty() && pointAndFieldFound_)
    {
        // Find the point
        scalar minDist = GREAT;
        forAll(mesh.points(), pI)
        {
            scalar dist = mag(mesh.points()[pI] - point);
            if (dist < minDist)
            {
                minDist = dist;
                pointID_ = pI;
            }
        }

        Info<< this->name() << ": distance from specified point is " << minDist
            << endl;

        // Create interpolator
        interpPtr_ = new newLeastSquaresVolPointInterpolation(mesh);

        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                time_.timeName(mesh.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"history"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"history"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            // Open new file at start up
            historyFilePtr_.reset
                (
                    new OFstream
                    (
                        historyDir/"solidPointStress"+fieldName_+".dat"
                    )
                );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time value" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidPointStress::start()
{
    return writeData();
}


#if FOAMEXTEND > 40
bool Foam::solidPointStress::execute(const bool forceWrite)
#else
bool Foam::solidPointStress::execute()
#endif
{
    return writeData();
}


bool Foam::solidPointStress::read(const dictionary& dict)
{
    return true;
}

// ************************************************************************* //
