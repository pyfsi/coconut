/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "magLongDeltaFa.H"
#include "edgeFields.H"
#include "areaFields.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(magLongDeltaFa, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::magLongDeltaFa::magLongDeltaFa(const faMesh& mesh)
:
    MeshObject<faMesh, magLongDeltaFa>(mesh),
    magLongDeltaFaPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::magLongDeltaFa::~magLongDeltaFa()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::magLongDeltaFa::clearOut() const
{
    deleteDemandDrivenData(magLongDeltaFaPtr_);
}


void Foam::magLongDeltaFa::makeMagLongDistance() const
{
    if (magLongDeltaFaPtr_)
    {
        FatalErrorIn("void magLongDeltaFa::makeMagLongDistance() const")
            << "Long face distances already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        InfoIn("magLongDeltaFa::makeMagLongDistance()")
            << "Constructing magnitude of long face distance"
            << endl;
    }

    magLongDeltaFaPtr_ = new edgeScalarField
    (
        IOobject
        (
            "magLongDeltaFa",
            mesh().pointsInstance(),
            mesh().db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh(),
        dimensionedScalar("zero", dimLength, 0.0)
    );
    edgeScalarField& mldp = *magLongDeltaFaPtr_;

    // Set local references to mesh data
    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const vectorField& Cf = mesh().edgeCentres();
    const vectorField& C = mesh().areaCentres();
    const vectorField& Sf = mesh().Le();
    const scalarField& magSf = mesh().magLe();

    forAll (owner, edgei)
    {
        // This must be the same as in edgeInterpolation.C
        scalar SfdOwn = mag(Sf[edgei] & (Cf[edgei] - C[owner[edgei]]));
        scalar SfdNei = mag(Sf[edgei] & (C[neighbour[edgei]] - Cf[edgei]));
        mldp[edgei] = (SfdOwn + SfdNei)/magSf[edgei];
    }

    forAll (mldp.boundaryField(), patchi)
    {
        mldp.boundaryField()[patchi] = calcMagLongDistance(patchi);
    }

    if (debug)
    {
        InfoIn("magLongDeltaFa::makeMagLongDistance()")
            << "Finished magnitude of long cell distance"
            << endl;
    }
}


Foam::tmp<Foam::scalarField> Foam::magLongDeltaFa::calcMagLongDistance
(
    const label patchi
) const
{
    const faPatch& p = mesh().boundary()[patchi];

    vectorField d = p.faPatch::delta();

    if (p.coupled())
    {
        return (mag(p.edgeLengths() & d) + mag(p.edgeLengths() & (p.delta() - d)))
               /p.magEdgeLengths();
    }
    else
    {
        return mag(p.edgeLengths() & d)/p.magEdgeLengths();
    }
}


const Foam::edgeScalarField& Foam::magLongDeltaFa::magDelta() const
{
    if (!magLongDeltaFaPtr_)
    {
        makeMagLongDistance();
    }

    return *magLongDeltaFaPtr_;
}


const Foam::scalarField& Foam::magLongDeltaFa::magDelta
(
    const label patchi
) const
{
    return magDelta().boundaryField()[patchi];
}


bool Foam::magLongDeltaFa::movePoints() const
{
    if (debug)
    {
        InfoIn("bool magLongDeltaFa::movePoints() const")
            << "Clearing long face distance data" << endl;
    }

    clearOut();

    return true;
}


bool Foam::magLongDeltaFa::updateMesh(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        InfoIn("bool magLongDeltaFa::updateMesh(const mapPolyMesh&) const")
            << "Clearing long face distance data" << endl;
    }

    clearOut();

    return true;
}


// ************************************************************************* //
