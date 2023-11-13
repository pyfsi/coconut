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
    Writes the mesh in fpma format readable by AVL's CfdWM

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "argList.H"
#include "Time.H"
#include "sizeType.H"
#include "LongList.H"
#include "FRWGraph.H"
#include "VRWGraph.H"
#include "IFstream.H"
#include "helperFunctions.H"

#include "meshOctreeCreator.H"

#include <climits>
#include <cstdlib>

# ifdef USE_OMP
#include <omp.h>
# endif

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    int maxIntVal = INT_MAX;
    std::cout << "INT_MAX has value" << maxIntVal << std::endl;
    unsigned long maxLongVal = LONG_MAX;
    std::cout << "ULONG_MAX has value " << maxLongVal << std::endl;

    unsigned long long maxLongLongVal = LONG_LONG_MAX;
    std::cout << "LONG_LONG_MAX has value " << maxLongLongVal << std::endl;

    sizeType maxSizeTypeVal = FOAM_SIZETYPE_MAX;
    std::cout << "Max sizeType value " << maxSizeTypeVal << std::endl;

    Info << "Size of sizeType " << label(sizeof(sizeType)) << endl;
/*
    sizeType sRead;
    IFstream file("sizeType.dat");
    file >> sRead;
    Info << "sRead " << sRead << endl;

    LongList<sizeType> ll;

    sizeType s = FOAM_LABEL_MAX;
    Info << "Current s " << s << endl;
    s += 2;
    Info << "Requested list size " << s << endl;
    for(sizeType i=0;i<s;++i)
        ll.append(i);
    Info << "Size of the list " << ll.size() << endl;
    ll.setSize(0);

    VRWGraph vg(s/10+5000, 10, -1);

    # ifdef USE_OMP
    const scalar startTime = omp_get_wtime();

    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp single nowait
        # endif
        {
            # ifdef USE_OMP
            const label nTasks = 10 * omp_get_num_threads();
            # else
            const label nTasks = 1;
            # endif

            const label chSize = max(vg.size() / nTasks, 1);

            for(label taskI=0;taskI<nTasks;++taskI)
            {
                # ifdef USE_OMP
                # pragma omp task shared(vg) firstprivate(taskI)
                # endif
                {
                    const label sr = taskI * chSize;
                    const label er = taskI!=(nTasks-1)?sr+chSize:vg.size();

                    for(label i=sr;i<er;++i)
                    {
                        forAllRow(vg, i, j)
                            vg(i, j) = j;

                        forAllRow(vg, i, j)
                            if( vg(i, j) < 0 )
                                Info << "Drekec" << endl;
                    }
                }
            }
        }
    }

    # ifdef USE_OMP
    Info << "Time for shuffling " << (omp_get_wtime()-startTime) << endl;
    # endif
    vg.setSize(0);
*/
    if( Pstream::parRun() )
    {
        std::map<label, labelLongList> exchangeData;

        for(label procI=0;procI<Pstream::nProcs();++procI)
        {
            if( procI == Pstream::myProcNo() )
                continue;

            labelLongList& dts = exchangeData[procI];
            dts.clear();

            for(label i=1;i<10;++i)
                dts.append(Pstream::myProcNo() * i);
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        for(label procI=0;procI<Pstream::nProcs();++procI)
        {
            if( procI == Pstream::myProcNo() )
            {
                Pout << "Received data " << receivedData << endl;
            }

            returnReduce(1, sumOp<label>());
        }
    }

    IOdictionary meshDict
    (
        IOobject
        (
            "meshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Pout << "REading surface file" << endl;
    fileName surfaceFile = meshDict.lookup("surfaceFile");
//    if( Pstream::parRun() )
//        surfaceFile = ".."/surfaceFile;
    triSurf surf(surfaceFile);

    Pout << "Constructing octree" << endl;
    meshOctree octree(surf);
    Pout << "Creating octree" << endl;
    meshOctreeCreator(octree, meshDict).createOctreeBoxes();

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
