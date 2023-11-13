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

#include "error.H"
#include "polyMeshGenModifier.H"
#include "meshSurfaceDistanceFromGeometry.H"
#include "meshSurfaceEngine.H"
#include "meshOctree.H"
#include "triSurf.H"
#include "helperFunctions.H"

# ifdef USE_OMP
#include <omp.h>
# endif

//#define DEBUGEdgeExtractor

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
// Private member functions

void meshSurfaceDistanceFromGeometry::findActiveFacesAndPoints
(
    const wordList& avoidPatches
)
{
    activeFaces_.setSize(surfaceEngine_.boundaryFaces().size());
    activePoints_.setSize(surfaceEngine_.boundaryPoints().size());

    //- find the patches that shall be avoided
    const polyMeshGen& mesh = surfaceEngine_.mesh();
    std::set<label> avoidIds;
    forAll(mesh.boundaries(), patchI)
    {
        const word pName = mesh.boundaries()[patchI].patchName();

        forAll(avoidPatches, i)
            if( avoidPatches[i] == pName )
            {
                avoidIds.insert(patchI);

                break;
            }
    }

    //- mark active faces and points
    const labelLongList& facePatches = surfaceEngine_.boundaryFacePatches();
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        //- marking of active faces
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(facePatches, bfI)
        {
            if( avoidIds.find(facePatches[bfI]) == avoidIds.end() )
            {
                //- face is active
                activeFaces_[bfI] = true;
            }
            else
            {
                //- face shall not be considered
                activeFaces_[bfI] = false;
            }
        }

        //- marking of active points. A point attached to an active faces
        //- is considered active
        # ifdef USE_OMP
        # pragma omp for schedule(static, 1)
        # endif
        forAll(pointFaces, bpI)
        {
            bool anyActive(false);
            forAllRow(pointFaces, bpI, pfI)
            {
                if( activeFaces_[pointFaces(bpI, pfI)] )
                {
                    anyActive = true;
                    break;
                }
            }

            activePoints_[bpI] = anyActive;
        }
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();

        //- allocate the map for exchanging data
        std::map<label, labelLongList> exchangeData;
        forAll(surfaceEngine_.bpNeiProcs(), i)
            exchangeData[surfaceEngine_.bpNeiProcs()[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( activePoints_[bpI] )
            {
                forAllRow(bpAtProcs, bpI, i)
                {
                    const label neiProc = bpAtProcs(bpI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append(it.key());
                }
            }
        }

        labelLongList receiveData;
        help::exchangeMap(exchangeData, receiveData);

        forAll(receiveData, i)
        {
            const label bpI = globalToLocal[receiveData[i]];

            activePoints_[bpI] = true;
        }
    }
}

void meshSurfaceDistanceFromGeometry::calculatePointDistances()
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();

    boundaryPointDistance_.setSize(bPoints.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 20)
    # endif
    forAll(bPoints, bpI)
    {
        if( !activePoints_[bpI] )
        {
            boundaryPointDistance_[bpI] = 0.0;
            continue;
        }

        const point& p = points[bPoints[bpI]];

        point np;
        scalar dSq;
        label nt, patch;
        octree_.findNearestSurfacePoint(np, dSq, nt, patch, p);

        boundaryPointDistance_[bpI] = mag(np - p);
    }
}

void meshSurfaceDistanceFromGeometry::calculateFaceCentreDistances()
{
    const vectorLongList& fCentres = surfaceEngine_.faceCentres();

    faceCentreDistance_.setSize(fCentres.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 20)
    # endif
    forAll(fCentres, bfI)
    {
        if( !activeFaces_[bfI] )
        {
            faceCentreDistance_[bfI] = 0.0;
            continue;
        }

        const point& p = fCentres[bfI];

        point np;
        scalar dSq;
        label nt, patch;
        octree_.findNearestSurfacePoint(np, dSq, nt, patch, p);

        faceCentreDistance_[bfI] = mag(np - p);
    }
}

void meshSurfaceDistanceFromGeometry::calculatePointNormalDeviationAngle()
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const vectorLongList& pNormals = surfaceEngine_.pointNormals();

    const triSurf& surf = octree_.surface();
    const vectorField& tNormals = surf.facetNormals();

    pointNormalDeviationAngle_.setSize(bPoints.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 20)
    # endif
    forAll(bPoints, bpI)
    {
        if( !activePoints_[bpI] )
        {
            pointNormalDeviationAngle_[bpI] = 0.0;
            continue;
        }

        const point& p = points[bPoints[bpI]];

        point np;
        scalar dSq;
        label nt, patch;
        octree_.findNearestSurfacePoint(np, dSq, nt, patch, p);

        vector n = pNormals[bpI];
        n /= (mag(n) + VSMALL);

        vector sn = tNormals[nt];
        sn /= (mag(sn) + VSMALL);

        const scalar val = mag(acos(max(-1.0, min(n & sn, 1.0))));
        pointNormalDeviationAngle_[bpI] = val;
    }
}

void meshSurfaceDistanceFromGeometry::calculateFaceNormalDeviationAngle()
{
    const vectorLongList& fCentres = surfaceEngine_.faceCentres();
    const vectorLongList& fNormals = surfaceEngine_.faceNormals();

    const triSurf& surf = octree_.surface();
    const vectorField& tNormals = surf.facetNormals();

    faceNormalDeviationAngle_.setSize(fNormals.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(guided, 20)
    # endif
    forAll(fCentres, bfI)
    {
        if( !activeFaces_[bfI] )
        {
            faceNormalDeviationAngle_[bfI] = 0.0;
            continue;
        }

        const point& p = fCentres[bfI];

        point np;
        scalar dSq;
        label nt, patch;
        octree_.findNearestSurfacePoint(np, dSq, nt, patch, p);

        vector n = fNormals[bfI];
        n /= (mag(n) + VSMALL);

        vector sn = tNormals[nt];
        sn /= (mag(sn) + VSMALL);

        const scalar val = mag(acos(max(-1.0, min(n & sn, 1.0))));
        faceNormalDeviationAngle_[bfI] = val;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

meshSurfaceDistanceFromGeometry::meshSurfaceDistanceFromGeometry
(
    const meshSurfaceEngine& surfaceEngine,
    const meshOctree& octree,
    const wordList& avoidPatches
)
:
    surfaceEngine_(surfaceEngine),
    octree_(octree),
    activePoints_(),
    activeFaces_(),
    boundaryPointDistance_(),
    faceCentreDistance_(),
    pointNormalDeviationAngle_(),
    faceNormalDeviationAngle_()
{
    findActiveFacesAndPoints(avoidPatches);

    calculatePointDistances();

    calculateFaceCentreDistances();

    calculatePointNormalDeviationAngle();

    calculateFaceNormalDeviationAngle();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
// Destructor

meshSurfaceDistanceFromGeometry::~meshSurfaceDistanceFromGeometry()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

scalar meshSurfaceDistanceFromGeometry::boundaryPointDistance
(
    const label bpI
) const
{
    return boundaryPointDistance_[bpI];
}

scalar meshSurfaceDistanceFromGeometry::boundaryFaceCentreDistance
(
    const label bfI
) const
{
    return faceCentreDistance_[bfI];
}

scalar meshSurfaceDistanceFromGeometry::angleDeviationAtBoundaryPoint
(
    const label bpI
) const
{
    return pointNormalDeviationAngle_[bpI];
}

scalar meshSurfaceDistanceFromGeometry::angleDeviationAtBoundaryFace
(
    const label bfI
) const
{
    return faceNormalDeviationAngle_[bfI];
}

void meshSurfaceDistanceFromGeometry::writeToVTK(const fileName& fName) const
{
    OFstream file(fName);

    //- write the header
    file << "# vtk DataFile Version 3.0\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    const pointField& points = surfaceEngine_.points();
    const labelLongList& bPoints = surfaceEngine_.boundaryPoints();
    const labelLongList& bp = surfaceEngine_.bp();
    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();

    //- write points
    file << "POINTS " << bPoints.size() << " float\n";
    forAll(bPoints, bpI)
    {
        const point& p = points[bPoints[bpI]];

        file << p.x() << ' ' << p.y() << ' ' << p.z() << ' ';

        if( bpI % 5 == 0 )
            file << "\n";
    }

    //- write triangles
    label nFacePoints(0);
    forAll(bFaces, bfI)
        nFacePoints += bFaces[bfI].size() + 1;

    file << "\n";
    file << "\nPOLYGONS " << bFaces.size() << " " << nFacePoints << endl;
    forAll(bFaces, bfI)
    {
        const face& bf = bFaces[bfI];

        file << bf.size();
        forAll(bf, pI)
            file << " " << bp[bf[pI]];
        file << nl;
    }

    //- point distance
//    file << "\n";
//    file << "\nPOINT_DATA " << boundaryPointDistance_.size() << "\n";

//    file << "SCALARS PointDistances double\n";
//    file << "LOOKUP_TABLE default\n";
//    forAll(boundaryPointDistance_, bpI)
//    {
//        file << boundaryPointDistance_[bpI] << " ";

//        if( bpI % 5 == 0 )
//            file << endl;
//    }

    //- point normal deviation
//    file << "\n";
//    file << "\nPOINT_DATA " << pointNormalDeviationAngle_.size() << "\n";

//    file << "SCALARS PointNormalDeviationAngle double\n";
//    file << "LOOKUP_TABLE default\n";
//    forAll(pointNormalDeviationAngle_, bpI)
//    {
//        file << pointNormalDeviationAngle_[bpI] << " ";

//        if( bpI % 5 == 0 )
//            file << endl;
//    }

    //- face centre distance
    file << "\n";
    file << "\nCELL_DATA " << faceCentreDistance_.size() << "\n";

    file << "SCALARS faceDistances double\n";
    file << "LOOKUP_TABLE default\n";
    forAll(faceCentreDistance_, bpI)
    {
        file << faceCentreDistance_[bpI] << " ";

        if( bpI % 5 == 0 )
            file << endl;
    }

    //- face centre normal deviation
//    file << "\n";
//    file << "\nCELL_DATA " << faceNormalDeviationAngle_.size() << "\n";

//    file << "SCALARS faceNormalDeviationAngle double\n";
//    file << "LOOKUP_TABLE default\n";
//    forAll(faceNormalDeviationAngle_, bpI)
//    {
//        file << faceNormalDeviationAngle_[bpI] << " ";

//        if( bpI % 5 == 0 )
//            file << endl;
//    }

    //- flush the file on disk before exitting
    file.flush();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
