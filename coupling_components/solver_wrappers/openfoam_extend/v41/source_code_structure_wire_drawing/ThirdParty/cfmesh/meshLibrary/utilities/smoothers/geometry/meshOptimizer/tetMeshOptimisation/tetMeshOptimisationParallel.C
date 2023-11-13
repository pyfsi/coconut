/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "demandDrivenData.H"
#include "tetMeshOptimisation.H"
#include "partTetMesh.H"
#include "VRWGraph.H"
#include "polyMeshGenAddressing.H"
#include "helperFunctions.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void tetMeshOptimisation::unifyNegativePoints(boolList& negativeNode) const
{
    //- make sure that processor nodes are selected on all processors
    //- where they exist
    const LongList<direction>& smoothVertex = tetMesh_.smoothVertex();
    const polyMeshGenAddressing& addr = tetMesh_.origMesh().addressingData();
    const DynList<label>& neiProcs = addr.pointNeiProcs();
    const VRWGraph& pProcs = addr.pointAtProcs();
    const Map<label>& globalToLocal = addr.globalToLocalPointAddressing();

    std::map<label, labelLongList> selectedNegativeNodes;
    forAll(neiProcs, procI)
        selectedNegativeNodes.insert
        (
            std::make_pair(neiProcs[procI], labelLongList())
        );

    forAllConstIter(Map<label>, globalToLocal, it)
    {
        const label nI = it();

        if( !negativeNode[nI] )
            continue;
        if( !(smoothVertex[nI] & partTetMesh::PARALLELBOUNDARY) )
            continue;

        forAllRow(pProcs, nI, procI)
        {
            const label neiProc = pProcs(nI, procI);

            if( neiProc == Pstream::myProcNo() )
                continue;

            selectedNegativeNodes[neiProc].append(it.key());
        }
    }

    labelLongList receivedNodes;
    help::exchangeMap(selectedNegativeNodes, receivedNodes);

    forAll(receivedNodes, i)
        negativeNode[globalToLocal[receivedNodes[i]]] = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
