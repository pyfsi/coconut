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
    Testing of LongList

\*---------------------------------------------------------------------------*/

#include "labelLongList.H"
#include "labelledPoint.H"
#include "refLabelledPoint.H"
#include "labelList.H"

#include "OFstream.H"
#include "IFstream.H"

#include <omp.h>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    labelLongList firstList;

    const label size = 100000;

    const scalar startFirst = omp_get_wtime();
    firstList.setSize(size);
    for(label i=0;i<size;++i)
    {
        firstList.append(i);
        //firstList[i] = i;
    }
    const scalar endFirst = omp_get_wtime();
    Info << "Time for first append " << (endFirst - startFirst) << endl;

    {
//        labelList aaa(20);
//        forAll(aaa, i)
//            aaa[i] = i;

        //firstList.clear();

        OFstream file("abc.dat", OFstream::BINARY, OFstream::currentVersion, OFstream::COMPRESSED);
        file << firstList;
        //file << aaa;
        file.flush();
    }

    {
        labelLongList readFirstList;
        IFstream file("abc.dat", IFstream::BINARY);

        file >> readFirstList;

        Info << "Read from file " << readFirstList <<endl;

        if( firstList.size()  != readFirstList.size() )
            Info << "Wrong size of read list " << readFirstList.size()
                 << " should be " << firstList.size() << endl;

        forAll(firstList, i)
            if( firstList[i] != readFirstList[i] )
                Info << "Something is wrong with element " << i
                     << " orig value " << firstList[i]
                     << " read value " << readFirstList << endl;
    }

    Info << "Size of labelledPoint " << label(sizeof(labelledPoint)) << endl;
    const scalar startSecond = omp_get_wtime();
    LongList<labelledPoint> secondList;
    const labelledPoint lp(0, point::zero);
    secondList.setSize(size);
    for(label i=0;i<size;++i)
    {
        //const labelledPoint lp(i, point(i, i, i));
        //secondList.append(lp);
        secondList[i] = lp;
    }
    const scalar endSecond = omp_get_wtime();
    Info << "Time for second append " << (endSecond - startSecond) << endl;
    return 0;
}

// ************************************************************************* //
