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

#include "metisWireLayerAdditionRemovalDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "floatScalar.H"
#include "foamTime.H"
#include "wireStreamlines.H"
//#include "scotchDecomp.H"

extern "C"
{
#define OMPI_SKIP_MPICXX
#   include "metis.h"
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(metisWireLayerAdditionRemovalDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        metisWireLayerAdditionRemovalDecomp,
        dictionaryMesh
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::metisWireLayerAdditionRemovalDecomp::decompose
(
    const List<int>& adjncy,
    const List<int>& xadj,
    const scalarField& cWeights,
    List<int>& finalDecomp,
    int nProcs
)
{
    // Method of decomposition
    // recursive: multi-level recursive bisection (default)
    // k-way: multi-level k-way
    word method("k-way");

    int numCells = xadj.size()-1;

    // decomposition options. 0 = use defaults
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    // processor weights initialised with no size, only used if specified in
    // a file - use the data type that metis expects here
    Field<real_t> processorWeights;

    // cell weights (so on the vertices of the dual)
    List<int> cellWeights;

    // face weights (so on the edges of the dual)
    List<int> faceWeights;

    // Check for externally provided cellweights and if so initialise weights
    scalar minWeights = gMin(cWeights);
    if (cWeights.size() > 0)
    {
        if (minWeights <= 0)
        {
            WarningIn
            (
                "metisDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Illegal minimum weight " << minWeights
                << endl;
        }

        if (cWeights.size() != numCells)
        {
            FatalErrorIn
            (
                "metisDecomp::decompose"
                "(const pointField&, const scalarField&)"
            )   << "Number of cell weights " << cWeights.size()
                << " does not equal number of cells " << numCells
                << exit(FatalError);
        }
        // Convert to integers.
        cellWeights.setSize(cWeights.size());
        forAll(cellWeights, i)
        {
            cellWeights[i] = int(cWeights[i]/minWeights);
        }
    }


    // Check for user supplied weights and decomp options
    if (decompositionDict_.found("metisCoeffs"))
    {
        const dictionary& metisCoeffs =
            decompositionDict_.subDict("metisCoeffs");
        word weightsFile;

        if (metisCoeffs.readIfPresent("method", method))
        {
            if (method != "recursive" && method != "k-way")
            {
                FatalErrorIn("metisDecomp::decompose()")
                    << "Method " << method
                    << " in metisCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be 'recursive' or 'k-way'"
                    << exit(FatalError);
            }

            Info<< "metisDecomp : Using Metis method     " << method
                << nl << endl;
        }

        List<int> mOptions;
        if (metisCoeffs.readIfPresent("options", mOptions))
        {
            if (mOptions.size() != METIS_NOPTIONS)
            {
                FatalErrorIn("metisDecomp::decompose()")
                    << "Number of options in metisCoeffs in dictionary : "
                    << decompositionDict_.name()
                    << " should be " << METIS_NOPTIONS
                    << exit(FatalError);
            }

            forAll(mOptions, i)
            {
                options[i] = mOptions[i];
            }

            Info<< "metisDecomp : Using Metis options     " << mOptions
                << nl << endl;
        }

        if (metisCoeffs.readIfPresent("processorWeights", processorWeights))
        {
            processorWeights /= sum(processorWeights);

            //if (processorWeights.size() != nProcessors_)
            if (processorWeights.size() != nProcs)
            {
                FatalErrorIn("metisDecomp::decompose(const pointField&)")
                    << "Number of processor weights "
                    << processorWeights.size()
                    //<< " does not equal number of domains " << nProcessors_
                    << " does not equal number of domains " << nProcs
                    << exit(FatalError);
            }
        }

        if (metisCoeffs.readIfPresent("cellWeightsFile", weightsFile))
        {
            Info<< "metisDecomp : Using cell-based weights." << endl;

            IOList<int> cellIOWeights
            (
                IOobject
                (
                    weightsFile,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
            cellWeights.transfer(cellIOWeights);

            if (cellWeights.size() != xadj.size()-1)
            {
                FatalErrorIn("metisDecomp::decompose(const pointField&)")
                    << "Number of cell weights " << cellWeights.size()
                    << " does not equal number of cells " << xadj.size()-1
                    << exit(FatalError);
            }
        }
    }

    //int nProcs = nProcessors_;

    // output: cell -> processor addressing
    finalDecomp.setSize(numCells);

    // output: number of cut edges
    int edgeCut = 0;

    // Vertex weight info
    int* vwgtPtr = NULL;
    int* adjwgtPtr = NULL;

    if (cellWeights.size())
    {
        vwgtPtr = cellWeights.begin();
    }
    if (faceWeights.size())
    {
        adjwgtPtr = faceWeights.begin();
    }

    int one = 1;

    if (method == "recursive")
    {
        METIS_PartGraphRecursive
        (
            &numCells,         // num vertices in graph
            &one,
            const_cast<List<int>&>(xadj).begin(),   // indexing into adjncy
            const_cast<List<int>&>(adjncy).begin(), // neighbour info
            vwgtPtr,           // vertexweights
            NULL,
            adjwgtPtr,         // no edgeweights
            &nProcs,
            processorWeights.begin(),
            NULL,
            options,
            &edgeCut,
            finalDecomp.begin()
        );
    }
    else
    {
        METIS_PartGraphKway
        (
            &numCells,         // num vertices in graph
            &one,
            const_cast<List<int>&>(xadj).begin(),   // indexing into adjncy
            const_cast<List<int>&>(adjncy).begin(), // neighbour info
            vwgtPtr,           // vertexweights
            NULL,
            adjwgtPtr,         // no edgeweights
            &nProcs,
            processorWeights.begin(),
            NULL,
            options,
            &edgeCut,
            finalDecomp.begin()
        );
    }

    return edgeCut;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::metisWireLayerAdditionRemovalDecomp::metisWireLayerAdditionRemovalDecomp
(
    const dictionary& decompositionDict,
    const polyMesh& mesh
)
:
    metisDecomp(decompositionDict, mesh),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::metisWireLayerAdditionRemovalDecomp::decompose
(
    const pointField& points,
    const scalarField& pointWeights
)
{
    Info<< "Performing initial standard metis decomposition of the entire case"
        << nl << endl;
    labelList decomp = metisDecomp::decompose(points, pointWeights);

    // We will modify the decomposition in the wire cellZone so that it is only
    // decomposed in the YZ plane
    // This is to avoid problems for the layer addition/removal
    const label cellZoneID = mesh_.cellZones().findZoneID("wire");

    if (cellZoneID == -1)
    {
        WarningIn
        (
            "Foam::labelList Foam::metisWireLayerAdditionRemovalDecomp::"
            "decompose\n"
            "(\n"
            "    const pointField& points,\n"
            "    const scalarField& pointWeights\n"
            ")"
        )   << "No 'wire' cellZone found!" << endl;
    }
    else
    {
        const labelList& wireCellIDs = mesh_.cellZones()[cellZoneID];

        // Count number of cells within the wiure assigned to each processor

        const int nProcs = nProcessors_;

        labelList cellsPerProcInWire(nProcs, 0);

        forAll(wireCellIDs, cI)
        {
            const label cellID = wireCellIDs[cI];

            cellsPerProcInWire[decomp[cellID]]++;
        }

        // Calculate the percentage breakdown with the wire
        // and set the processor weights
        // We will also remember the processors that have cells in the wire
        Field<real_t> processorWeights(nProcs, 0.0);
        SLList<label> wireProcIDsSL;

        Info<< "Checking the decomposition with the wire:" << endl;

        forAll(processorWeights, procI)
        {
            processorWeights[procI] =
                scalar(cellsPerProcInWire[procI])/scalar(wireCellIDs.size());

            Info<< "    proc " << procI << " has "
                << 100*processorWeights[procI]
                << "% of the cells in the wire" << endl;

            // Metis does not like when a processor weight is zero so we will
            // find the processors that do have cells in the wire
            if (processorWeights[procI] > SMALL)
            {
                wireProcIDsSL.append(procI);
            }
        }

        processorWeights /= sum(processorWeights);

        // These are the processors that have cells in the wire
        const labelList wireProcIDs = labelList(wireProcIDsSL);
        const label nWireProcs = wireProcIDs.size();

        Info<< nl << "Processors with at least 1 wire cell: " << endl;
        forAll(wireProcIDs, wireProcI)
        {
            Info<< "    proc " << wireProcIDs[wireProcI] << endl;
        }

        const int nProcsOnWIre = wireProcIDs.size();

        // No need to modify the decomposition if there is only one processor
        // cells on the wire
        if (nProcsOnWIre > 1)
        {
            // Create processorWeight just for the wire processors
            scalarList wireProcessorWeights(nWireProcs, 0.0);
            forAll(wireProcessorWeights, wireProcI)
            {
                wireProcessorWeights[wireProcI] =
                    processorWeights[wireProcIDs[wireProcI]];
            }

            // Add processor weights to the decomposeParDict based on these
            // percentage breakdown
            dictionary& decomposeParDict =
                const_cast<dictionary&>(decompositionDict_);

            dictionary metisCoeffs;
            metisCoeffs.add("processorWeights", wireProcessorWeights);
            //metisCoeffs.add("processorWeights", processorWeights);

            decomposeParDict.add("metisCoeffs", metisCoeffs);

            // Lookup the upstream wire patch
            const word wireUpstreamPatchName = word("wireUpstream");
            const label upstreamPatchID =
                mesh_.boundaryMesh().findPatchID(wireUpstreamPatchName);

            if (upstreamPatchID == -1)
            {
                FatalErrorIn
                (
                    "Foam::labelList Foam::"
                    "metisWireLayerAdditionRemovalDecomp::"
                    "decompose\n"
                    "(\n"
                    "    const pointField& points,\n"
                    "    const scalarField& pointWeights\n"
                    ")"
                )   << "Patch 'wireUpstream' not found!" << abort(FatalError);
            }

            const polyPatch& ppatch = mesh_.boundaryMesh()[upstreamPatchID];

            const vectorField& patchCf = ppatch.faceCentres();
            const scalarField patchCfWeights(patchCf.size(), 1.0);
            const labelListList& patchFaceFaces = ppatch.faceFaces();

            List<int> adjncy;
            List<int> xadj;
            calcCSR
            (
                patchFaceFaces,
                adjncy,
                xadj
            );

            List<int> patchDecompWire;
            Info<< nl
                << "Decomposing the wire upstream patch using metis" << endl;
            decompose
            (
                adjncy, xadj, patchCfWeights, patchDecompWire, nWireProcs
            );

            // Convert local wire procIDs to procIDs
            List<int> patchDecomp(patchDecompWire.size(), 0);
            forAll(patchDecompWire, faceI)
            {
                const label wireProcID = patchDecompWire[faceI];
                patchDecomp[faceI] = wireProcIDs[wireProcID];
            }

            // We will now propagate this patch decomposition down the length of
            // the wire
            Info<< nl
                << "Propagating the upstream patch decomposition down the wire"
                << " streamlines" << endl;

            // Create a wire streamline object
            wireStreamlines streamlines(mesh_, wireUpstreamPatchName);

            // Get the cell streamlines
            const List< SLList<label> >& cellStreamlines =
                streamlines.cellStreamlines();

            if (cellStreamlines.size() != patchDecomp.size())
            {
                FatalErrorIn
                (
                    "Foam::labelList Foam::"
                    "metisWireLayerAdditionRemovalDecomp::"
                    "decompose\n"
                    "(\n"
                    "    const pointField& points,\n"
                    "    const scalarField& pointWeights\n"
                    ")"
                )   << "Something is wrong with the creation of the "
                    << "streamlines!"
                    << " There are " << cellStreamlines.size() << " streamlines"
                    << " but the patch size is " << patchDecomp.size()
                    << abort(FatalError);
            }

            // Update cell decomposition along streamlines
            Info<< nl << "Updating the final cell decomposition" << endl;
            forAll(cellStreamlines, streamlineI)
            {
                const labelList curStreamline = cellStreamlines[streamlineI];

                // Processor ID for this streamline
                const label streamlineProcID = patchDecomp[streamlineI];

                forAll(curStreamline, cI)
                {
                    const label cellID = curStreamline[cI];

                    decomp[cellID] = streamlineProcID;
                }
            }

            // Print out the final decomposition in the wire
            {
                labelList cellsPerProcInWire(nProcs, 0);

                forAll(wireCellIDs, cI)
                {
                    const label cellID = wireCellIDs[cI];

                    cellsPerProcInWire[decomp[cellID]]++;
                }

                scalarList processorWeights(nProcs, 0.0);
                Info<< nl
                    << "The final % breakdown of cells weights for the wire "
                    << "is:" << endl;
                forAll(processorWeights, procI)
                {
                    processorWeights[procI] =
                        scalar(cellsPerProcInWire[procI])
                        /scalar(wireCellIDs.size());

                    Info<< "    proc " << procI << " has "
                        << 100*processorWeights[procI]
                        << "% of the cells in the wire" << endl;
                }
            }
            Info<< endl;
        }
        else
        {
            Info<< nl << "There is only one processor with cells in the wire"
                << nl << endl;
        }
    }

    return decomp;
}


Foam::labelList Foam::metisWireLayerAdditionRemovalDecomp::decompose
(
    const labelList& fineToCoarse,
    const pointField& coarsePoints,
    const scalarField& coarseWeights
)
{
    return metisDecomp::decompose(fineToCoarse, coarsePoints, coarseWeights);
}


Foam::labelList Foam::metisWireLayerAdditionRemovalDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cc,
    const scalarField& cWeights
)
{
    return metisDecomp::decompose(globalCellCells, cc, cWeights);
}


// ************************************************************************* //
