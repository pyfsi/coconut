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

Class
    Foam::metisWireLayerAdditionRemovalDecomp

Description
    Metis domain decomposition, where the "wire" cellZone is read and decomposed
    only in the YZ plane; the rest of the cells are decomposed using metis.

    Default metis will be used if no "wire" cellZone exists.

SourceFiles
    metisWireLayerAdditionRemovalDecomp.C

Author
    Philp Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#ifndef metisWireLayerAdditionRemovalDecomp_H
#define metisWireLayerAdditionRemovalDecomp_H

#include "metisDecomp.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class metisWireLayerAdditionRemovalDecomp Declaration
\*---------------------------------------------------------------------------*/

class metisWireLayerAdditionRemovalDecomp
:
    public metisDecomp
{
    // Private data

        const polyMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        metisWireLayerAdditionRemovalDecomp
        (
            const metisWireLayerAdditionRemovalDecomp&
        );

        //- Disallow default bitwise assignment
        void operator=(const metisWireLayerAdditionRemovalDecomp&);

        label decompose
        (
            const List<int>& adjncy,
            const List<int>& xadj,
            const scalarField& cellWeights,
            List<int>& finalDecomp,
            int nProcs
        );


public:

    //- Runtime type information
    TypeName("metisWireLayerAdditionRemoval");


    // Constructors

        //- Construct given the decomposition dictionary and mesh
        metisWireLayerAdditionRemovalDecomp
        (
            const dictionary& decompositionDict,
            const polyMesh& mesh
        );


    // Destructor

        virtual ~metisWireLayerAdditionRemovalDecomp()
        {}


    // Member Functions

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        //  Weights get normalised so the minimum value is 1 before truncation
        //  to an integer so the weights should be multiples of the minimum
        //  value. The overall sum of weights might otherwise overflow.
        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Return for every coordinate the wanted processor number. Gets
        //  passed agglomeration map (from fine to coarse cells) and coarse cell
        //  location. Can be overridden by decomposers that provide this
        //  functionality natively.
        //  See note on weights above.
        virtual labelList decompose
        (
            const labelList& agglom,
            const pointField& regionPoints,
            const scalarField& regionWeights
        );

        //- Same but with uniform weights
        virtual labelList decompose
        (
            const labelList& agglom,
            const pointField& regionPoints
        )
        {
            return decompose
            (
                agglom,
                regionPoints,
                scalarField(regionPoints.size(), 1.0)
            );
        }

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided mesh connectivity.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        //  See note on weights above.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
