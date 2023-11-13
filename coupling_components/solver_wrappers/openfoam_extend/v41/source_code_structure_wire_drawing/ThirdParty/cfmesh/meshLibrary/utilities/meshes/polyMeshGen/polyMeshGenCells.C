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

#include "polyMeshGenCells.H"
#include "polyMeshGenAddressing.H"
#include "IOobjectList.H"
#include "cellSet.H"
#include "demandDrivenData.H"

#include "labelPair.H"

#include "OFstream.H"

# ifdef USE_OMP
#include <omp.h>
# endif

namespace Foam
{

// * * * * * * * * * * Private member functions * * * * * * * * * * * * * * * //

void polyMeshGenCells::calculateOwnersAndNeighbours() const
{
    if( ownerPtr_ || neighbourPtr_ )
        FatalErrorIn
        (
            "void polyMeshGenCells::calculateOwnersAndNeighbours() const"
        ) << "Owners and neighbours are already allocated" << abort(FatalError);

    //- allocate owners
    ownerPtr_ =
        new labelIOLongList
        (
            IOobject
            (
                "owner",
                instance_,
                meshDir_,
                runTime_
            ),
            faces_.size()
        );
    labelIOLongList& own = *ownerPtr_;

    //- allocate neighbours
    neighbourPtr_ =
        new labelIOLongList
        (
            IOobject
            (
                "neighbour",
                instance_,
                meshDir_,
                runTime_
            ),
            faces_.size()
        );
    labelIOLongList& nei = *neighbourPtr_;

    //- start calculating owners and neighbours
    nIntFaces_ = 0;

    label nInternalFaces(0);

    List<List<LongList<labelPair> > > dataForOtherThreads;

    # ifdef USE_OMP
    # pragma omp parallel reduction(+ : nInternalFaces)
    # endif
    {
        # ifdef USE_OMP
        const label nThreads = omp_get_num_threads();
        const label threadI = omp_get_thread_num();

        # pragma omp single
        {
            dataForOtherThreads.setSize(nThreads);
        }
        # else
        const label nThreads = 1;
        const label threadI(0);

        dataForOtherThreads.setSize(nThreads);
        # endif

        const label chunkSize = max(1, faces_.size() / nThreads);
        const label startingFace = min(threadI * chunkSize, faces_.size());
        label endFace = min(startingFace + chunkSize, faces_.size());
        if( threadI == (nThreads - 1) )
            endFace = faces_.size();

        List<LongList<labelPair> >& dot = dataForOtherThreads[threadI];
        dot.setSize(nThreads);

        for(label faceI=startingFace;faceI<endFace;++faceI)
        {
            own[faceI] = -1;
            nei[faceI] = -1;
        }

        # ifdef USE_OMP
        # pragma omp for schedule(static)
        # endif
        forAll(cells_, cellI)
        {
            const cell& c = cells_[cellI];

            forAll(c, fI)
            {
                const label faceI = c[fI];

                const label threadNo = min(faceI / chunkSize, nThreads - 1);

                if( threadNo == threadI )
                {
                    if( own[faceI] == -1 )
                    {
                        own[faceI] = cellI;
                    }
                    else if( nei[faceI] == -1 )
                    {
                        nei[faceI] = cellI;
                        ++nInternalFaces;
                    }
                    else
                    {
                        Serr << "Face " << faces_[faceI] << endl;
                        Serr << "Owner " << own[faceI] << endl;
                        Serr << "Neighbour " << nei[faceI] << endl;
                        Serr << "Current cell " << cellI << endl;
                        FatalErrorIn
                        (
                            "void polyMeshGenCells::"
                            "calculateOwnersAndNeighbours()"
                        ) << Pstream::myProcNo() << "Face " << faceI
                            << " appears in more than 2 cells!!"
                            << abort(FatalError);
                    }
                }
                else
                {
                    dot[threadNo].append(labelPair(faceI, cellI));
                }
            }
        }

        # ifdef USE_OMP
        # pragma omp critical(otherData)
        # endif
        for(label i=0;i<nThreads;++i)
        {
            const LongList<labelPair>& data =
                dataForOtherThreads[i][threadI];

            forAll(data, j)
            {
                const label faceI = data[j].first();
                const label cellI = data[j].second();

                if( own[faceI] == -1 )
                {
                    own[faceI] = cellI;
                }
                else if( own[faceI] > cellI )
                {
                    if( nei[faceI] == -1 )
                    {
                        nei[faceI] = own[faceI];
                        own[faceI] = cellI;
                        ++nInternalFaces;
                    }
                    else
                    {
                        Serr << "Face " << faces_[faceI] << endl;
                        Serr << "Owner " << own[faceI] << endl;
                        Serr << "Neighbour " << nei[faceI] << endl;
                        Serr << "Current cell " << cellI << endl;
                        FatalErrorIn
                        (
                            "void polyMeshGenCells::"
                            "calculateOwnersAndNeighbours()"
                        ) << Pstream::myProcNo() << "Face " << faceI
                            << " appears in more than 2 cells!!"
                            << abort(FatalError);
                    }
                }
                else if( nei[faceI] == -1 )
                {
                    nei[faceI] = cellI;
                    ++nInternalFaces;
                }
                else
                {
                    Serr << "Face " << faces_[faceI] << endl;
                    Serr << "Owner " << own[faceI] << endl;
                    Serr << "Neighbour " << nei[faceI] << endl;
                    Serr << "Current cell " << cellI << endl;
                    FatalErrorIn
                    (
                        "void polyMeshGenCells::"
                        "calculateOwnersAndNeighbours()"
                    ) << Pstream::myProcNo() << "Face " << faceI
                        << " appears in more than 2 cells!!"
                        << abort(FatalError);
                }
            }
        }
    }

    nIntFaces_ = nInternalFaces;
}

void polyMeshGenCells::calculateAddressingData() const
{
    if( !ownerPtr_ || !neighbourPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline label polyMeshGenCells::calculateAddressingData() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateOwnersAndNeighbours();
    }

    addressingDataPtr_ = new polyMeshGenAddressing(*this);
}

void polyMeshGenCells::clearOut() const
{
    polyMeshGenFaces::clearOut();
    deleteDemandDrivenData(addressingDataPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGenCells::polyMeshGenCells
(
    const Time& runTime,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenFaces(runTime, instance, meshDir),
    polyMeshGenCellZones(cells_),
    cells_(),
    cellSubsets_(),
    addressingDataPtr_(NULL)
{}

//- Construct from components without the boundary
polyMeshGenCells::polyMeshGenCells
(
    const Time& runTime,
    const pointField& points,
    const faceList& faces,
    const cellList& cells,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenFaces(runTime, points, faces, instance, meshDir),
    polyMeshGenCellZones(cells_),
    cells_(),
    cellSubsets_(),
    addressingDataPtr_(NULL)
{
    cells_ = cells;
}

//- Construct from components with the boundary
polyMeshGenCells::polyMeshGenCells
(
    const Time& runTime,
    const pointField& points,
    const faceList& faces,
    const cellList& cells,
    const wordList& patchNames,
    const labelList& patchStart,
    const labelList& nFacesInPatch,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenFaces
    (
        runTime,
        points,
        faces,
        patchNames,
        patchStart,
        nFacesInPatch,
        instance,
        meshDir
    ),
    polyMeshGenCellZones(cells_),
    cells_(),
    cellSubsets_(),
    addressingDataPtr_(NULL)
{
    cells_ = cells;
}

//- Construct from components with the boundary and boundary types
polyMeshGenCells::polyMeshGenCells
(
    const Time& runTime,
    const pointField& points,
    const faceList& faces,
    const cellList& cells,
    const wordList& patchNames,
    const wordList& patchTypes,
    const labelList& patchStart,
    const labelList& nFacesInPatch,
    const labelList& procPatchMyProcNo,
    const labelList& procPatchNeighbProcNo,
    const fileName instance,
    const fileName meshDir
)
:
    polyMeshGenFaces
    (
        runTime,
        points,
        faces,
        patchNames,
        patchTypes,
        patchStart,
        nFacesInPatch,
        procPatchMyProcNo,
        procPatchNeighbProcNo,
        instance,
        meshDir
    ),
    polyMeshGenCellZones(cells_),
    cells_(),
    cellSubsets_(),
    addressingDataPtr_(NULL)
{
    cells_ = cells;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Destructor
polyMeshGenCells::~polyMeshGenCells()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- return addressing which may be needed
const polyMeshGenAddressing& polyMeshGenCells::addressingData() const
{
    if( !addressingDataPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "inline label polyMeshGenCells::addressingData() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calculateAddressingData();
    }

    return *addressingDataPtr_;
}

void polyMeshGenCells::clearAddressingData() const
{
    deleteDemandDrivenData(addressingDataPtr_);
}

label polyMeshGenCells::addCellSubset(const word& selName)
{
    label id = cellSubsetIndex(selName);
    if( id >= 0 )
    {
        Warning << "Cell subset " << selName << " already exists!" << endl;
        return id;
    }

    id = 0;
    for
    (
        std::map<label, meshSubset>::const_iterator it=cellSubsets_.begin();
        it!=cellSubsets_.end();
        ++it
    )
        id = Foam::max(id, it->first+1);

    cellSubsets_.insert
    (
        std::make_pair
        (
            id,
            meshSubset(selName, meshSubset::CELLSUBSET)
        )
    );

    return id;
}

void polyMeshGenCells::removeCellSubset(const label setI)
{
    if( cellSubsets_.find(setI) == cellSubsets_.end() )
        return;

    cellSubsets_.erase(setI);
}

word polyMeshGenCells::cellSubsetName(const label setI) const
{
    std::map<label, meshSubset>::const_iterator it =
        cellSubsets_.find(setI);
    if( it == cellSubsets_.end() )
    {
        Warning << "Subset " << setI << " is not a cell subset" << endl;
        return word();
    }

    return it->second.name();
}

label polyMeshGenCells::cellSubsetIndex(const word& selName) const
{
    std::map<label, meshSubset>::const_iterator it;
    for(it=cellSubsets_.begin();it!=cellSubsets_.end();++it)
    {
        if( it->second.name() == selName )
            return it->first;
    }

    return -1;
}

void polyMeshGenCells::read()
{
    polyMeshGenFaces::read();

    Info << "Starting creating cells" << endl;
    //- count the number of cells and create the cells
    label nCells(0);
    const labelLongList& own = this->owner();
    const labelLongList& nei = this->neighbour();

    forAll(own, faceI)
    {
        if( own[faceI] >= nCells )
            nCells = own[faceI] + 1;

        if( nei[faceI] >= nCells )
            nCells = nei[faceI] + 1;
    }

    LongList<direction> nFacesInCell(nCells, direction(0));
    forAll(own, faceI)
        ++nFacesInCell[own[faceI]];

    forAll(nei, faceI)
        if( nei[faceI] != -1 )
            ++nFacesInCell[nei[faceI]];

    cells_.setSize(nCells);
    forAll(cells_, cellI)
        cells_[cellI].setSize(nFacesInCell[cellI]);

    nFacesInCell = 0;
    forAll(own, faceI)
    {
        cells_[own[faceI]][nFacesInCell[own[faceI]]++] = faceI;
        if( nei[faceI] != -1 )
            cells_[nei[faceI]][nFacesInCell[nei[faceI]]++] = faceI;
    }

    // read cell subsets
    const fileName localDir = meshDir_ / "sets";
    fileName path = instance_.path();
    word name;
    if( path == "." )
    {
        path = instance_.name();
    }
    else
    {
        name = instance_.name();
    }

    IOobjectList allSets
    (
        runTime_,
        path,
        fileName(name/localDir)
    );

    wordList setNames = allSets.names("cellSet");
    forAll(setNames, setI)
    {
        IOobject* obj = allSets.lookup(setNames[setI]);

        cellSet cSet(*obj);

        const labelList content = cSet.toc();
        const label id = addCellSubset(setNames[setI]);

        forAll(content, i)
            addCellToSubset(id, content[i]);
    }

    readCellZones(runTime_, instance_, meshDir_);
}

void polyMeshGenCells::readFromLatestTime()
{
    polyMeshGenFaces::readFromLatestTime();

    const instantList instances = runTime_.times();

    fileName origInstance = instance_;

    forAllReverse(instances, instanceI)
    {
        instance_ = instances[instanceI].name() / origInstance;

        IOobject readFaces
        (
            "faces",
            instance_,
            meshDir_,
            runTime_,
            IOobject::MUST_READ
        );

        if( !readFaces.headerOk() )
            continue;

        polyMeshGenCells::read();
        break;
    }

    instance_ = origInstance;
}

void polyMeshGenCells::write() const
{
    polyMeshGenFaces::write();

    //- write cell subsets
    std::map<label, meshSubset>::const_iterator setIt;
    for(setIt=cellSubsets_.begin();setIt!=cellSubsets_.end();++setIt)
    {
        labelLongList containedElements;
        setIt->second.containedElements(containedElements);

        //- create a cell set
        cellSet set
        (
            IOobject
            (
                setIt->second.name(),
                instance_,
                meshDir_/"sets",
                runTime_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        );

        forAll(containedElements, i)
            set.insert(containedElements[i]);

        set.write();
    }

    writeCellZones(runTime_, instance_, meshDir_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
