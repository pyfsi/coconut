/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

sizeType polyMeshGenAddressing::printAllocated() const
{
    sizeType totalStorage(sizeof(polyMeshGenAddressing));
    const label toMb(1048576);

    Pout << "polyMeshGenAddressing allocated :" << endl;

    if( edgesPtr_ )
    {
        Pout << "Edges " << (edgesPtr_->byteSize()/toMb) << endl;
        totalStorage += edgesPtr_->byteSize();
    }

    if( ccPtr_ )
    {
        Pout << "Cell-cells " << (ccPtr_->printAllocated()/toMb) << endl;
        totalStorage += ccPtr_->printAllocated();
    }

    if( ecPtr_ )
    {
        Pout << "Edge-cells " << (ecPtr_->printAllocated()/toMb) << endl;
        totalStorage += ecPtr_->printAllocated();
    }

    if( pcPtr_ )
    {
        Pout << "Point-cells " << (pcPtr_->printAllocated()/toMb) << endl;
        totalStorage += pcPtr_->printAllocated();
    }

    if( efPtr_ )
    {
        Pout << "Edge-faces " << (efPtr_->printAllocated()/toMb) << endl;
        totalStorage += efPtr_->printAllocated();
    }

    if( pfPtr_ )
    {
        Pout << "Point-faces " << (pfPtr_->printAllocated()/toMb) << endl;
        totalStorage += pfPtr_->printAllocated();
    }

    if( cePtr_ )
    {
        Pout << "Cell-edges " << (cePtr_->printAllocated()/toMb) << endl;
        totalStorage += cePtr_->printAllocated();
    }

    if( fePtr_ )
    {
        Pout << "Face-edges " << (fePtr_->printAllocated()/toMb) << endl;
        totalStorage += fePtr_->printAllocated();
    }

    if( pePtr_ )
    {
        Pout << "Point-edges " << (pePtr_->printAllocated()/toMb) << endl;
        totalStorage += pePtr_->printAllocated();
    }

    if( ppPtr_ )
    {
        Pout << "Point-point " << (ppPtr_->printAllocated()/toMb) << endl;
        totalStorage += ppPtr_->printAllocated();
    }

    if( cpPtr_ )
    {
        Pout << "Cell-point " << (cpPtr_->printAllocated()/toMb) << endl;
        totalStorage += cpPtr_->printAllocated();
    }

    // Geometry
    if( cellCentresPtr_ )
    {
        Pout << "Cell-centres " << (cellCentresPtr_->byteSize()/toMb) << endl;
        totalStorage += cellCentresPtr_->byteSize();
    }

    if( faceCentresPtr_ )
    {
        Pout << "Face-centres " << (faceCentresPtr_->byteSize()/toMb) << endl;
        totalStorage += faceCentresPtr_->byteSize();
    }

    if( cellVolumesPtr_ )
    {
        Pout << "Cell-volumes " << (cellVolumesPtr_->byteSize()/toMb) << endl;
        totalStorage += cellVolumesPtr_->byteSize();
    }

    if( faceAreasPtr_ )
    {
        Pout << "Face-areas " << (faceAreasPtr_->byteSize()/toMb) << endl;
        totalStorage += faceAreasPtr_->byteSize();
    }

    // Parallel
    if( globalPointLabelPtr_ )
    {
        Pout << " Global point labels "
             << (globalPointLabelPtr_->printAllocated()/toMb) << endl;
        totalStorage += globalPointLabelPtr_->printAllocated();
    }

    if( globalFaceLabelPtr_ )
    {
        Pout << " Global face labels "
             << (globalFaceLabelPtr_->printAllocated()/toMb) << endl;
        totalStorage += globalFaceLabelPtr_->printAllocated();
    }

    if( globalCellLabelPtr_ )
    {
        Pout << " Global cell labels "
             << (globalCellLabelPtr_->printAllocated()/toMb) << endl;
        totalStorage += globalCellLabelPtr_->printAllocated();
    }

    if( globalEdgeLabelPtr_ )
    {
        Pout << "Global edge labels "
             << (globalEdgeLabelPtr_->printAllocated()/toMb) << endl;
        totalStorage += globalEdgeLabelPtr_->printAllocated();
    }

    if( pProcsPtr_ )
    {
        Pout << "Point at processors "
             << (pProcsPtr_->printAllocated()/toMb) << endl;
        totalStorage += pProcsPtr_->printAllocated();
    }

    if( globalToLocalPointAddressingPtr_ )
    {
        Pout << "Global to local point addressing " << endl;
             //<< globalToLocalEdgeAddressingPtr_->byteSize() << endl;
    }

    if( pointNeiProcsPtr_ )
    {
        Pout << "Nei processors over vertices "
             << (pointNeiProcsPtr_->printAllocated()/toMb) << endl;
        totalStorage += pointNeiProcsPtr_->printAllocated();
    }

    if( eProcsPtr_ )
    {
        Pout << "Edges at processors "
             << (eProcsPtr_->printAllocated()/toMb) << endl;
        totalStorage += eProcsPtr_->printAllocated();
    }

    if( globalToLocalEdgeAddressingPtr_ )
    {
        Pout << "Global-to-local edge addressing " << endl;
             //<< globalToLocalEdgeAddressingPtr_->byteSize() << endl;
    }

    if( edgeNeiProcsPtr_ )
    {
        Pout << "Neighbour processors over edges "
             << (edgeNeiProcsPtr_->printAllocated()/toMb) << endl;
        totalStorage += edgeNeiProcsPtr_->printAllocated();
    }

    Pout << "polyMeshGenAddressing storage " << (totalStorage/toMb) << endl;

    return totalStorage;
}

void polyMeshGenAddressing::clearGeom() const
{
    if( debug )
    {
        Pout<< "polyMeshGenAddressing::clearGeom() : "
            << "clearing geometric data"
            << endl;
    }

    deleteDemandDrivenData(cellCentresPtr_);
    deleteDemandDrivenData(faceCentresPtr_);
    deleteDemandDrivenData(cellVolumesPtr_);
    deleteDemandDrivenData(faceAreasPtr_);
}

void polyMeshGenAddressing::clearAddressing() const
{
    if( debug )
    {
        Pout<< "polyMeshGenAddressing::clearAddressing() : "
            << "clearing topology"
            << endl;
    }

    clearOutEdges();

    deleteDemandDrivenData(ccPtr_);
    deleteDemandDrivenData(ecPtr_);
    deleteDemandDrivenData(pcPtr_);

    deleteDemandDrivenData(efPtr_);
    deleteDemandDrivenData(pfPtr_);

    deleteDemandDrivenData(cePtr_);
    deleteDemandDrivenData(fePtr_);
    deleteDemandDrivenData(pePtr_);
    deleteDemandDrivenData(ppPtr_);
    deleteDemandDrivenData(cpPtr_);
}

void polyMeshGenAddressing::clearParallelAddressing() const
{
    deleteDemandDrivenData(globalPointLabelPtr_);
    deleteDemandDrivenData(globalFaceLabelPtr_);
    deleteDemandDrivenData(globalCellLabelPtr_);
    deleteDemandDrivenData(globalEdgeLabelPtr_);

    deleteDemandDrivenData(pProcsPtr_);
    deleteDemandDrivenData(globalToLocalPointAddressingPtr_);
    deleteDemandDrivenData(pointNeiProcsPtr_);
    deleteDemandDrivenData(eProcsPtr_);
    deleteDemandDrivenData(globalToLocalEdgeAddressingPtr_);
    deleteDemandDrivenData(edgeNeiProcsPtr_);
}

void polyMeshGenAddressing::clearOut() const
{
    clearGeom();
    clearAddressing();
    clearParallelAddressing();
}


void polyMeshGenAddressing::clearAll() const
{
    clearGeom();
    clearAddressing();
    clearParallelAddressing();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
