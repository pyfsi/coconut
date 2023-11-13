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

#include "rollingMillMesh.H"
#include "rollerSurfaceCreator.H"
#include "triSurf.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "triSurfaceCleanupDuplicates.H"
#include "triSurfaceCopyParts.H"
#include "triSurfaceChecks.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"
#include "helperFunctions.H"

//#define DEBUGContact

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool rollingMillContactHandler::readGeometryData()
{
    //- get the information about the roll setup
    const word rollSetup(meshDict_.lookup("rollSetup"));

    //- read the the geometry data
    bool evaluateContact(false);

    if( meshDict_.found(rollSetup+"Dict") )
    {
        const dictionary& dict = meshDict_.subDict(rollSetup+"Dict");

        Switch useContact;
        dict.readIfPresent("contact", useContact);

        if( useContact )
        {
            if( !dict.isDict("contactDict") )
            {
                FatalErrorIn
                (
                    "bool rollingMillContactHandler::readGeometryData()"
                ) << "contactDict dictionary does not exist in " << dict
                  << exit(FatalError);
            }

            const dictionary& contactDict = dict.subDict("contactDict");

            contactDict.readIfPresent("writeRollProfiles", writeRollProfile_);

            evaluateContact = true;

            //- calculate tolerance from roll width
            const wordList rollPositions = dict.toc();
            forAll(rollPositions, i)
            {
                if
                (
                    (rollPositions[i].find("Roll") != word::npos) &&
                    dict.isDict(rollPositions[i])
                )
                {
                    const dictionary& rollDict = dict.subDict(rollPositions[i]);
                    //- update the contact tolerance
                    scalar w;
                    if( rollDict.readIfPresent("rollWidth", w) )
                    {
                        contactDistanceTolerance_ =
                            max(contactDistanceTolerance_, 0.0025 * w);
                    }
                }
            }

            //- if exists, read the distance tolerance from disk
            if( dict.found("contactDistanceTolerance") )
            {
                contactDistanceTolerance_ =
                    readScalar(dict.lookup("contactDistanceTolerance"));
            }

            //- get the surfaces and move them to the correct position
            crossSectionSurfPtr_ = new triSurf();
            triSurfModifier sMod(*crossSectionSurfPtr_);

            manifold_.clear();
            rollPositionToManifold_.clear();

            forAll(rollGeometries_, geomI)
            {
                const rollGeometryInfo& geom = rollGeometries_[geomI];

                const word rollPosition = geom.rollPosition();

                //- set the roll to the correct manifold
                rollPositionToManifold_[rollPosition] = geomI;

                //- set roll displacement
                vector displacement(vector::zero);
                if( rollPosition == "topRoll" )
                {
                    displacement.y() = 0.5 * geom.outerDiameter();
                    displacement.y() += geom.radialShiftDistance();
                }
                else if( rollPosition == "bottomRoll" )
                {
                    displacement.y() = -0.5 * geom.outerDiameter();
                    displacement.y() -= geom.radialShiftDistance();
                }
                else if( rollPosition == "leftRoll" )
                {
                    displacement.x() = -0.5 * geom.outerDiameter();
                    displacement.x() -= geom.radialShiftDistance();
                }
                else if( rollPosition == "rightRoll" )
                {
                    displacement.x() = 0.5 * geom.outerDiameter();
                    displacement.x() += geom.radialShiftDistance();
                }
                else
                {
                    FatalErrorIn
                    (
                        "bool rollingMillContactHandler::readGeometryData()"
                    ) << "Unknown roll position " << rollPosition
                        << " cannot be moved to its position"
                        << exit(FatalError);
                }

                const triSurf* sPtr = geom.transformedSurface();

                //- copy and transform points
                const label nPts = crossSectionSurfPtr_->nPoints();
                sMod.pointsAccess().setSize(nPts+sPtr->nPoints());
                forAll(sPtr->points(), pI)
                {
                    sMod.pointsAccess()[nPts+pI] =
                        sPtr->points()[pI] + displacement;
                }

                //- copy patches
                const geometricSurfacePatchList& patches = sPtr->patches();
                DynList<label> patchIndex(patches.size(), -1);

                geometricSurfacePatchList& sPatches = sMod.patchesAccess();
                forAll(patches, patchI)
                {
                    forAll(sPatches, i)
                    {
                        if( sPatches[i].name() == patches[patchI].name() )
                            patchIndex[patchI] = i;
                    }

                    if( patchIndex[patchI] < 0 )
                    {
                        const label s = sPatches.size();
                        sPatches.setSize(s+1);
                        sPatches[s].name() = patches[patchI].name();
                        sPatches[s].geometricType() =
                            patches[patchI].geometricType();

                        patchIndex[patchI] = s;
                    }
                }

                forAll(sPtr->facets(), tI)
                {
                    const labelledTri& t = sPtr->facets()[tI];

                    //- store the origin of the triangle
                    manifold_.append(geomI);

                    //- append the triangle to the surface
                    //- rollerToAir region is not allow, replace with contact
                    crossSectionSurfPtr_->appendTriangle
                    (
                        labelledTri
                        (
                            t[0] + nPts,
                            t[1] + nPts,
                            t[2] + nPts,
                            patchIndex[t.region()]
                        )
                    );
                }
            }
        }

        # ifdef DEBUGContact
        crossSectionSurfPtr_->writeSurface("surfaceMeshes/crossSection.fms");
        # endif
    }

    return evaluateContact;
}

void rollingMillContactHandler::createOctree()
{
    if( octreePtr_ )
        return;

    octreePtr_ = new meshOctree(*crossSectionSurfPtr_, true);

    meshOctreeCreator(*octreePtr_).createOctreeWithRefinedBoundary(15, 50);
}

void rollingMillContactHandler::checkForContactRegions()
{
    createOctree();

    const scalar tolSq = sqr(contactDistanceTolerance_);

    Info << "Contact tolerance " << tolSq << endl;

    const triSurf& surf = *crossSectionSurfPtr_;
    const pointField& points = surf.points();
    const edgeLongList& edges = surf.edges();
    const VRWGraph& fEdges = surf.facetEdges();
    const VRWGraph& eFaces = surf.edgeFacets();

    //- detect contact patches
    const word contactName =
        patchNamesHandler_.patchNameForRoll
        (
            rollingMillPatchNamesHandler::ROLLERCONTACT
        );

    const word rollerToAirName =
        patchNamesHandler_.patchNameForRoll
        (
            rollingMillPatchNamesHandler::ROLLERTOAIR
        );

    const word backName =
        patchNamesHandler_.patchNameForRoll
        (
            rollingMillPatchNamesHandler::ROLLERBACK
        );

    const word frontName =
        patchNamesHandler_.patchNameForRoll
        (
            rollingMillPatchNamesHandler::ROLLERFRONT
        );

    const word symmYName =
        patchNamesHandler_.patchNameForRoll
        (
            rollingMillPatchNamesHandler::ROLLERSYMMY
        );
    const word symmZName =
        patchNamesHandler_.patchNameForRoll
        (
            rollingMillPatchNamesHandler::ROLLERSYMMZ
        );

    //- classify contact patches and pther patches
    labelList contactPatches(surf.patches().size(), NONE);
    forAll(contactPatches, i)
    {
        const word pName = surf.patches()[i].name();

        if
        (
            (pName.find(contactName) != word::npos) ||
            (pName.find(rollerToAirName) != word::npos)
        )
        {
            contactPatches[i] |= CONTACT;
        }
        else if( pName.find(backName) != word::npos )
        {
            contactPatches[i] |= BACK;
        }
        else if( pName.find(frontName) != word::npos )
        {
            contactPatches[i] |= FRONT;
        }
        else if
        (
            (pName.find(symmYName) != word::npos) ||
            (pName.find(symmZName) != word::npos)
        )
        {
            contactPatches[i] |= SYMMETRY;
        }
    }

    //- original patch tpe for each surface triangle
    //- it is needed in case if the whole front and back patches are in contact
    List<direction> origPatchType(surf.size(), NONE);
    forAll(surf, tI)
    {
        origPatchType[tI] |= contactPatches[surf[tI].region()];
    }

    //- find a contact patch for each manifold
    labelList contactPatchInManifold(rollPositionToManifold_.size(), -1);
    forAll(surf, tI)
    {
        if( contactPatches[surf[tI].region()] & CONTACT )
            contactPatchInManifold[manifold_[tI]] = surf[tI].region();
    }

    //- find overlapping surface parts
    labelLongList overlapTriangle(surf.size(), -1);

    bool foundContact(false);
    bool changed;
    do
    {
        //- detect overlap between triangles
        changed = false;

        forAll(overlapTriangle, tI)
        {
            if( overlapTriangle[tI] != -1 )
                continue;

            const labelledTri& t = surf[tI];

            vector n = t.normal(points);
            n /= (mag(n) + VSMALL);

            if( !(contactPatches[t.region()] & CONTACT) )
                continue;

            const point c = t.centre(points);

            DynList<point> searchPoints;
            searchPoints.append(c);
            scalar sr(contactDistanceTolerance_);
            forAll(t, pI)
            {
                const point& p = points[t[pI]];

                searchPoints.append(p);

                sr = max(sr, mag(p - c));
            }

            DynList<label> triangles;

            boundBox bb(c - point(sr, sr, sr), c + point(sr, sr, sr));
            octreePtr_->findTrianglesInBox(bb, triangles);

            const label manifoldI = manifold_[tI];

            scalar minDistSq(VGREAT);
            label nearestTri(-1);

            forAll(triangles, i)
            {
                const label ctI = triangles[i];

                if( manifold_[ctI] == manifoldI )
                    continue;
                if
                (
                    !(contactPatches[surf[ctI].region()] & (CONTACT|FRONT|BACK))
                )
                    continue;

                vector nn = surf[ctI].normal(points);
                nn /= (mag(nn) + VSMALL);

                if( mag(n & nn) < 0.99 )
                    continue;

                forAll(searchPoints, spI)
                {
                    const point& p = searchPoints[spI];
                    const point np =
                        help::nearestPointOnTheTriangle(ctI, surf, p);

                    const scalar dSq = magSqr(np - p);

                    if( dSq < minDistSq )
                    {
                        minDistSq = dSq;
                        nearestTri = ctI;

                        foundContact = true;
                    }
                }
            }

            if( minDistSq < tolSq )
            {
                overlapTriangle[tI] = nearestTri;
                overlapTriangle[nearestTri] = tI;

                if( !(contactPatches[surf[nearestTri].region()] & CONTACT) )
                {
                    //- set the patch to CONTACT
                    const label mI = manifold_[nearestTri];

                    const_cast<labelledTri&>(surf[nearestTri]).region() =
                            contactPatchInManifold[mI];

                    forAllRow(fEdges, nearestTri, teI)
                    {
                        const label eI = fEdges(nearestTri, teI);

                        if( eFaces.sizeOfRow(eI) != 2 )
                            continue;

                        vector ev = edges[eI].vec(points);
                        ev /= (mag(ev) + VSMALL);

                        if( (mag(ev.x()) > SMALL) || (mag(ev.y()) > SMALL) )
                        {
                            //- change the patch of the other triangle, too.
                            label neiTri = eFaces(eI, 0);
                            if( neiTri == nearestTri )
                                neiTri = eFaces(eI, 1);

                            overlapTriangle[neiTri] = tI;
                            const_cast<labelledTri&>(surf[neiTri]).region() =
                                    contactPatchInManifold[mI];

                            foundContact = true;
                        }
                    }
                }
            }
        }

        //- ensure continuous front and back patches
        //- in case continuity has been violated due to detected contacts
        for
        (
            std::map<word, label>::const_iterator it=
                rollPositionToManifold_.begin();
            it!=rollPositionToManifold_.end();
            ++it
        )
        {
            const word& rollPosition = it->first;
            const label manifoldI = it->second;

            vector axialDir(vector::zero);
            vector radialDir(vector::zero);

            if( rollPosition == "topRoll" )
            {
                axialDir.x() = 1.0;
                radialDir.y() = 1.0;
            }
            else if( rollPosition == "bottomRoll" )
            {
                axialDir.x() = 1.0;
                radialDir.y() = -1.0;
            }
            else if( rollPosition == "leftRoll" )
            {
                axialDir.y() = 1.0;
                radialDir.x() = -1.0;
            }
            else if( rollPosition == "rightRoll" )
            {
                axialDir.y() = 1.0;
                radialDir.x() = 1.0;
            }

            //- find min and max axial value
            scalar minAxial(VGREAT), maxAxial(-VGREAT);

            forAll(manifold_, tI)
            {
                if( manifold_[tI] == manifoldI )
                {
                    const labelledTri& tri = surf[tI];

                    forAll(tri, pI)
                    {
                        const scalar axial = points[tri[pI]] & axialDir;

                        minAxial = min(minAxial, axial);
                        maxAxial = max(maxAxial, axial);
                    }
                }
            }

            //- patches at minAxial
            scalar maxRadial(-VGREAT);
            forAll(manifold_, tI)
            {
                if( manifold_[tI] == manifoldI )
                {
                    const labelledTri& tri = surf[tI];

                    const point c = tri.centre(points);

                    if( mag((c & axialDir) - minAxial) < SMALL )
                    {
                        if( contactPatches[tri.region()] & CONTACT )
                        {
                            const scalar radial = c & radialDir;

                            maxRadial = max(radial, maxRadial);
                        }
                    }
                }
            }

            forAll(manifold_, tI)
            {
                if( manifold_[tI] == manifoldI )
                {
                    const labelledTri& tri = surf[tI];
                    const point c = tri.centre(points);

                    if( contactPatches[tri.region()] & BACK )
                    {
                        const scalar radial = c & radialDir;

                        if( radial < maxRadial )
                        {
                            changed = true;
                            const_cast<labelledTri&>(tri).region() =
                                    contactPatchInManifold[manifoldI];
                        }
                    }
                }
            }

            //- patches at maxAxial
            maxRadial = -VGREAT;
            forAll(manifold_, tI)
            {
                if( manifold_[tI] == manifoldI )
                {
                    const labelledTri& tri = surf[tI];

                    const point c = tri.centre(points);

                    if( mag((c & axialDir) - maxAxial) < SMALL )
                    {
                        if( contactPatches[tri.region()] & CONTACT )
                        {
                            const scalar radial = c & radialDir;

                            maxRadial = max(radial, maxRadial);
                        }
                    }
                }
            }

            forAll(manifold_, tI)
            {
                if( manifold_[tI] == manifoldI )
                {
                    const labelledTri& tri = surf[tI];
                    const point c = tri.centre(points);

                    if( contactPatches[tri.region()] & FRONT )
                    {
                        const scalar radial = c & radialDir;

                        if( radial < maxRadial )
                        {
                            changed = true;
                            const_cast<labelledTri&>(tri).region() =
                                    contactPatchInManifold[manifoldI];
                        }
                    }
                }
            }
        }
    }
    while( changed );

    if( !foundContact )
    {
        Warning << "No contact detected between the rolls" << endl;
        return;
    }

    //- check whether all triangles within a certain range of coordinates
    //- are marked as roller-to-roller contact
    for
    (
        std::map<word, label>::const_iterator it=
            rollPositionToManifold_.begin();
        it!=rollPositionToManifold_.end();
        ++it
    )
    {
        //- find the range in the positive the negative side
        const word& rollPosition = it->first;
        const label manifoldI = it->second;

        scalar maxPos(-VGREAT), minPos(VGREAT);
        scalar maxNeg(-VGREAT), minNeg(VGREAT);

        bool hasPos(false), hasNeg(false);

        //- set the search direction for rolls
        vector dir(vector::zero);
        if( rollPosition == "topRoll" || rollPosition == "bottomRoll" )
        {
            dir.x() = 1.0;
        }
        else
        {
            dir.y() = 1.0;
        }

        forAll(manifold_, tI)
        {
            const labelledTri& t = surf[tI];

            if
            (
                (contactPatches[t.region()] & CONTACT) &&
                (manifold_[tI] == manifoldI) &&
                (overlapTriangle[tI] != -1)
            )
            {
                const labelledTri& t = surf[tI];

                forAll(t, pI)
                {
                    const scalar d = (dir & points[t[pI]]);

                    if( d < -SMALL )
                    {
                        hasNeg = true;
                        minNeg = min(minNeg, d);
                        maxNeg = max(maxNeg, d);
                    }
                    else
                    {
                        hasPos = true;
                        minPos = min(minPos, d);
                        maxPos = max(maxPos, d);
                    }
                }
            }
        }

        # ifdef DEBUGContact
        Info << "Roll " << rollPosition << " in manifold " << manifoldI
             << " minNeg " << minNeg << " maxNeg " << maxNeg
             << " minPos " << minPos << " maxPos " << maxPos << endl;
        # endif

        //- assign all triangles within this range to the contact region
        forAll(manifold_, tI)
        {
            const labelledTri& t = surf[tI];

            if
            (
                (contactPatches[t.region()] & CONTACT) &&
                (manifold_[tI] == manifoldI) &&
                (overlapTriangle[tI] < 0)
            )
            {
                const point fc = t.centre(points);

                const scalar d = (fc & dir);

                if( hasNeg && (d >= minNeg) && (d <= maxNeg) )
                {
                    overlapTriangle[tI] = 0;
                }
                else if( hasPos && (d >= minPos) && (d <= maxPos) )
                {
                    overlapTriangle[tI] = 0;
                }
            }
        }
    }

    //- cleanup of triangles at the boundary of roller contact regions
    # ifdef DEBUGContact
    labelLongList cleanedTriangles;
    # endif

    forAll(overlapTriangle, i)
    {
        if( overlapTriangle[i] < 0 )
            continue;

        forAllRow(fEdges, i, j)
        {
            const label eI = fEdges(i, j);

            if( eFaces.sizeOfRow(eI) != 2 )
                continue;

            label neiTriangle = eFaces(eI, 0);
            if( neiTriangle == i )
                neiTriangle = eFaces(eI, 1);

            if( overlapTriangle[neiTriangle] < 0 )
            {
                vector ev = edges[eI].vec(points);

                if( mag(ev.y()) > SMALL || mag(ev.x()) > SMALL )
                {
                    overlapTriangle[i] = -1;

                    # ifdef DEBUGContact
                    Info << "Cleaning triangle " << i << endl;
                    cleanedTriangles.append(i);
                    # endif
                }
            }
        }
    }

    # ifdef DEBUGContact
    const label contId = crossSectionSurfPtr_->addFacetSubset("contactFaces");
    const label cleanedId = crossSectionSurfPtr_->addFacetSubset("cleanedFace");
    forAll(overlapTriangle, tI)
    {
        if( overlapTriangle[tI] != -1 )
            crossSectionSurfPtr_->addFacetToSubset(contId, tI);
    }
    forAll(cleanedTriangles, i)
        crossSectionSurfPtr_->addFacetToSubset(cleanedId, cleanedTriangles[i]);
    # endif

    //- detect regions of triangles
    labelLongList triRegion(surf.size(), -1);

    label regionI(0);
    forAll(triRegion, tI)
    {
        if( triRegion[tI] != -1 )
            continue;

        labelLongList front;
        triRegion[tI] = regionI;
        front.append(tI);

        while( front.size() )
        {
            const label tLabel = front.removeLastElement();

            forAllRow(fEdges, tLabel, feI)
            {
                const label eI = fEdges(tLabel, feI);

                if( eFaces.sizeOfRow(eI) != 2 )
                    continue;

                label neiTriangle = eFaces(eI, 0);
                if( neiTriangle == tLabel )
                    neiTriangle = eFaces(eI, 1);

                //- skip already visited triangles
                if( triRegion[neiTriangle] != -1 )
                    continue;

                //- do not cross borders of contact regions
                if
                (
                    (overlapTriangle[neiTriangle] != -1) ^
                    (overlapTriangle[tLabel] != -1)
                )
                    continue;

                //- triangles shall be in the same patch
                if( surf[neiTriangle].region() != surf[tLabel].region() )
                    continue;

                triRegion[neiTriangle] = regionI;
                front.append(neiTriangle);
            }
        }

        ++regionI;
    }

    # ifdef DEBUGContact
    Info << "Number of regions " << regionI << endl;
    Info << "Regions " << triRegion << endl;
    labelList regIds(regionI);
    forAll(regIds, i)
    {
        regIds[i] =
            crossSectionSurfPtr_->addFacetSubset
            (
                "region_"+help::labelToText(i)
            );
    }

    forAll(triRegion, i)
        crossSectionSurfPtr_->addFacetToSubset(regIds[triRegion[i]], i);
    # endif

    //- check adjacencies between regions
    std::map<label, std::set<label> > regionToRegions;
    labelList regionToPatch(regionI, -1);
    labelList regionType(regionI, NONE);
    forAll(eFaces, eI)
    {
        if( eFaces.sizeOfRow(eI) != 2 )
            continue;

        const label t0 = eFaces(eI, 0);
        const label r0 = triRegion[t0];
        regionToPatch[r0] = surf[t0].region();
        if( overlapTriangle[t0] != -1 )
        {
            regionType[r0] |= ROLLERTOROLLERCONTACTPATCH;

            if( origPatchType[t0] & BACK )
            {
                regionType[r0] |= ROLLERTOROLLERCONTACTBACK;
            }
            else if( origPatchType[t0] & FRONT )
            {
                regionType[r0] |= ROLLERTOROLLERCONTACTFRONT;
            }
        }

        const label t1 = eFaces(eI, 1);
        const label r1 = triRegion[t1];
        regionToPatch[r1] = surf[t1].region();
        if( overlapTriangle[t1] != -1 )
        {
            regionType[r1] |= ROLLERTOROLLERCONTACTPATCH;

            if( origPatchType[t1] & BACK )
            {
                regionType[r1] |= ROLLERTOROLLERCONTACTBACK;
            }
            else if( origPatchType[t1] & FRONT )
            {
                regionType[r1] |= ROLLERTOROLLERCONTACTFRONT;
            }
        }

        if( r0 != r1 )
        {
            regionToRegions[r0].insert(r1);
            regionToRegions[r1].insert(r0);
        }
    }

    # ifdef DEBUGContact
    Info << "Region to patch " << regionToPatch << endl;
    Info << "Region type " << regionType << endl;

    for
    (
        std::map<label, std::set<label> >::const_iterator it=
            regionToRegions.begin();
        it!=regionToRegions.end();
        ++it
    )
    {
        Info << "Region " << it->first << " neighbour regions ";
        forAllConstIter(std::set<label>, it->second, nIt)
            Info << " " << *nIt;

        Info << endl;
    }
    # endif

    //- classify region type
    DynList<DynList<label> > sortedRegions;
    boolList usedRegion(regionI, false);
    for
    (
        std::map<label, std::set<label> >::const_iterator it=
            regionToRegions.begin();
        it!=regionToRegions.end();
        ++it
    )
    {
        if( usedRegion[it->first] )
            continue;

        DynList<label> connectedRegions;
        usedRegion[it->first] = true;
        connectedRegions.append(it->first);

        bool finished;
        do
        {
            finished = true;

            const label region = connectedRegions.lastElement();

            const std::set<label>& neiRegions = regionToRegions[region];
            forAllConstIter(std::set<label>, neiRegions, nIt)
            {
                if( !connectedRegions.contains(*nIt) )
                {
                    usedRegion[*nIt] = true;
                    connectedRegions.append(*nIt);
                    finished = false;
                    break;
                }
            }
        } while( !finished );

        sortedRegions.append(connectedRegions);
    }

    # ifdef DEBUGContact
    Info << "Connected regions " << sortedRegions << endl;
    # endif

    //- assign type to all regions
    forAll(sortedRegions, i)
    {
        const DynList<label>& regionsInChain = sortedRegions[i];

        //- find valid patterns and assign patch types accordingly
        forAll(regionsInChain, j)
        {
            //- check the combinations of previous and next regions
            //- and evaluate possible combinations
            const label prevRegion = regionsInChain.rcElement(j);
            const label prevPatch = regionToPatch[prevRegion];

            const label region = regionsInChain[j];
            const label patch = regionToPatch[region];

            const label nextRegion = regionsInChain.fcElement(j);
            const label nextPatch = regionToPatch[nextRegion];

            //- symm -> rollerToWire -> rollerToRoller
            if
            (
                (contactPatches[prevPatch] & SYMMETRY) &&
                (contactPatches[patch] & CONTACT) &&
                (regionType[nextRegion] & ROLLERTOROLLERCONTACTPATCH)
            )
            {
                regionType[region] |= ROLLERTOWIRECONTACTPATCH;
            }
            else if
            (
                (contactPatches[nextPatch] & SYMMETRY) &&
                (contactPatches[patch] & CONTACT) &&
                (regionType[prevRegion] & ROLLERTOROLLERCONTACTPATCH)
            )
            {
                regionType[region] |= ROLLERTOWIRECONTACTPATCH;
            }

            //- back -> rollerToWire -> front
            if
            (
                (
                    (contactPatches[prevPatch] & BACK) &&
                    (contactPatches[patch] & CONTACT) &&
                    (contactPatches[nextPatch] & FRONT)
                ) ||
                (
                    (contactPatches[prevPatch] & FRONT) &&
                    (contactPatches[patch] & CONTACT) &&
                    (contactPatches[nextPatch] & BACK)
                )
            )
            {
                # ifdef DEBUGContact
                Info << "1.Region " << region << " rollerToWire" << endl;
                # endif

                //- this region is a rollerToWire contact
                regionType[region] |= ROLLERTOWIRECONTACTPATCH;
            }

            //- front -> air -> rollerToRoller front
            if
            (
                (contactPatches[prevPatch] & FRONT) &&
                (contactPatches[patch] & CONTACT) &&
                (regionType[nextRegion] & ROLLERTOROLLERCONTACTPATCH)
            )
            {
                # ifdef DEBUGContact
                Info << "2.Region " << region << " rollerToAir" << endl;
                # endif

                //- the part of the roller profile is not in contact
                //- with the wire or another roller
                regionType[region] |= ROLLERTOAIR;
                regionType[nextRegion] |= ROLLERTOROLLERCONTACTFRONT;
            }
            else if
            (
                (contactPatches[nextPatch] & FRONT) &&
                (contactPatches[patch] & CONTACT) &&
                (regionType[prevRegion] & ROLLERTOROLLERCONTACTPATCH)
            )
            {
                # ifdef DEBUGContact
                Info << "3. Region " << region << " rollerToAir" << endl;
                # endif

                //- the part of the roller profile is not in contact
                //- with the wire or another roller
                regionType[region] |= ROLLERTOAIR;
                regionType[prevRegion] |= ROLLERTOROLLERCONTACTFRONT;
            }

            //- back -> air -> rollerToRoller back
            if
            (
                (contactPatches[prevPatch] & BACK) &&
                (contactPatches[patch] & CONTACT) &&
                (regionType[nextRegion] & ROLLERTOROLLERCONTACTPATCH)
            )
            {
                # ifdef DEBUGContact
                Info << "4.Region " << region << " rollerToAir" << endl;
                # endif

                //- the part of the roller profile is not in contact
                //- with the wire or another roller
                regionType[region] |= ROLLERTOAIR;
                regionType[nextRegion] |= ROLLERTOROLLERCONTACTBACK;
            }
            else if
            (
                (contactPatches[nextPatch] & BACK) &&
                (contactPatches[patch] & CONTACT) &&
                (regionType[prevRegion] & ROLLERTOROLLERCONTACTPATCH)
            )
            {
                # ifdef DEBUGContact
                Info << "5.Region " << region << " rollerToAir" << endl;
                # endif

                //- the part of the roller profile is not in contact
                //- with the wire or another roller
                regionType[region] |= ROLLERTOAIR;
                regionType[prevRegion] |= ROLLERTOROLLERCONTACTBACK;
            }

            //- front -> rollerToRoller -> rollerToWire
            if
            (
                (contactPatches[prevPatch] & FRONT) &&
                (regionType[region] & ROLLERTOROLLERCONTACTPATCH) &&
                (contactPatches[nextPatch] & CONTACT)
            )
            {
                # ifdef DEBUGContact
                Info << "6.Region " << region << " rollerToRoller" << endl;
                # endif

                //- the part of the roller profile is not in contact
                //- with the wire or another roller
                regionType[region] |= ROLLERTOROLLERCONTACTFRONT;
                regionType[nextRegion] |= ROLLERTOWIRECONTACTPATCH;
            }
            else if
            (
                (contactPatches[nextPatch] & FRONT) &&
                (regionType[region] & ROLLERTOROLLERCONTACTPATCH) &&
                (contactPatches[prevPatch] & CONTACT)
            )
            {
                # ifdef DEBUGContact
                Info << "7.Region " << region << " rollerToRoller" << endl;
                # endif

                //- the part of the roller profile is not in contact
                //- with the wire or another roller
                regionType[region] |= ROLLERTOROLLERCONTACTFRONT;
                regionType[prevRegion] |= ROLLERTOWIRECONTACTPATCH;
            }

            //- back -> rollerToRoller -> rollerToWire
            if
            (
                (contactPatches[prevPatch] & BACK) &&
                (regionType[region] & ROLLERTOROLLERCONTACTPATCH) &&
                (contactPatches[nextPatch] & CONTACT)
            )
            {
                # ifdef DEBUGContact
                Info << "8.Region " << region << " rollerToRoller" << endl;
                # endif

                //- the part of the roller profile is not in contact
                //- with the wire or another roller
                regionType[region] |= ROLLERTOROLLERCONTACTBACK;
                regionType[nextRegion] |= ROLLERTOWIRECONTACTPATCH;
            }
            else if
            (
                (contactPatches[nextPatch] & BACK) &&
                (regionType[region] & ROLLERTOROLLERCONTACTPATCH) &&
                (contactPatches[prevPatch] & CONTACT)
            )
            {
                # ifdef DEBUGContact
                Info << "9.Region " << region << " rollerToRoller" << endl;
                # endif

                //- the part of the roller profile is not in contact
                //- with the wire or another roller
                regionType[region] |= ROLLERTOROLLERCONTACTBACK;
                regionType[prevRegion] |= ROLLERTOWIRECONTACTPATCH;
            }

            //- rollerToRoller -> rollerToWire -> rollerToRoller
            if
            (
                (regionType[prevRegion] & ROLLERTOROLLERCONTACTPATCH) &&
                (contactPatches[patch] & CONTACT) &&
                (regionType[nextRegion] & ROLLERTOROLLERCONTACTPATCH)
            )
            {
                # ifdef DEBUGContact
                Info << "10.Region " << region << " rollerToWire" << endl;
                # endif

                //- the part of the roller profile is in contact
                //- with the wire
                regionType[region] = ROLLERTOWIRECONTACTPATCH;
            }
        }
    }

    # ifdef DEBUGContact
    const label rrFrontId =
        crossSectionSurfPtr_->addFacetSubset("rollerToRollerFront");
    const label rrBackId =
        crossSectionSurfPtr_->addFacetSubset("rollerToRollerBack");
    const label rrId = crossSectionSurfPtr_->addFacetSubset("rollerToRoller");
    const label rtWireId =
        crossSectionSurfPtr_->addFacetSubset("rollerToWire");
    const label rtAirId =
        crossSectionSurfPtr_->addFacetSubset("rollerToAir");

    forAll(triRegion, i)
    {
        if( regionType[triRegion[i]] & ROLLERTOWIRECONTACTPATCH )
            crossSectionSurfPtr_->addFacetToSubset(rtWireId, i);
        if( regionType[triRegion[i]] & ROLLERTOROLLERCONTACTBACK )
            crossSectionSurfPtr_->addFacetToSubset(rrBackId, i);
        if( regionType[triRegion[i]] & ROLLERTOROLLERCONTACTPATCH )
            crossSectionSurfPtr_->addFacetToSubset(rrId, i);
        if( regionType[triRegion[i]] & ROLLERTOROLLERCONTACTFRONT )
            crossSectionSurfPtr_->addFacetToSubset(rrFrontId, i);
        if( regionType[triRegion[i]] & ROLLERTOAIR )
            crossSectionSurfPtr_->addFacetToSubset(rtAirId, i);
    }

    Info << "Surface patches before " << surf.patches() << endl;
    # endif

    //- create patch names
    DynList<word> newPatchNames, newPatchTypes;
    labelList newPatchForRegion(regionI);
    forAll(regionType, i)
    {
        const label patchI = regionToPatch[i];
        const word origPatchName = surf.patches()[patchI].name();
        const word origPatchType = surf.patches()[patchI].geometricType();

        # ifdef DEBUGContact
        Info << "Region " << i << " patch name " << origPatchName
             << " type " << origPatchType << endl;
        # endif

        if( regionType[i] & ROLLERTOROLLERCONTACTFRONT )
        {
            //- get patch name
            const word pName =
                patchNamesHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERTOROLLERFRONT
                );

            # ifdef DEBUGContact
            Info << "Front patch name " << pName << endl;
            Info << "Front patch type " << origPatchType << endl;
            # endif

            if( !newPatchNames.contains(pName) )
            {
                newPatchForRegion[i] = newPatchNames.size();
                newPatchNames.append(pName);
                newPatchTypes.append("patch");
            }
            else
            {
                newPatchForRegion[i] = newPatchNames.containsAtPosition(pName);
            }
        }
        else if( regionType[i] & ROLLERTOROLLERCONTACTBACK )
        {
            //- get patch name
            const word pName =
                patchNamesHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERTOROLLERBACK
                );

            # ifdef DEBUGContact
            Info << "Back patch name " << pName << endl;
            Info << "Back patch type " << origPatchType << endl;
            # endif

            //- set patch name
            if( !newPatchNames.contains(pName) )
            {
                newPatchForRegion[i] = newPatchNames.size();
                newPatchNames.append(pName);
                newPatchTypes.append("patch");
            }
            else
            {
                newPatchForRegion[i] = newPatchNames.containsAtPosition(pName);
            }
        }
        else if( regionType[i] & ROLLERTOAIR )
        {
            const word pName =
                patchNamesHandler_.patchNameForRoll
                (
                    rollingMillPatchNamesHandler::ROLLERTOAIR
                );

            if( !newPatchNames.contains(pName) )
            {
                newPatchForRegion[i] = newPatchNames.size();
                newPatchNames.append(pName);
                newPatchTypes.append("patch");
            }
            else
            {
                newPatchForRegion[i] = newPatchNames.containsAtPosition(pName);
            }
        }
        else
        {
            //- other patches do not change
            if( newPatchNames.contains(origPatchName) )
            {
                newPatchForRegion[i] =
                    newPatchNames.containsAtPosition(origPatchName);
            }
            else
            {
                newPatchForRegion[i] = newPatchNames.size();
                newPatchNames.append(origPatchName);
                newPatchTypes.append(origPatchType);
            }
        }
    }

    # ifdef DEBUGContact
    forAll(newPatchForRegion, i)
    {
        Info << "Region " << i << " new patch "
             << newPatchNames[newPatchForRegion[i]] << endl;
    }

    crossSectionSurfPtr_->writeSurface
    (
        "surfaceMeshes/crossSectionWithSubsets.fms"
    );

    Info << "New patch names " << newPatchNames << endl;
    Info << "New patch types " << newPatchTypes << endl;
    # endif

    //- assign patches back to the input surface mesh
    label startTriangle(0);
    forAll(rollGeometries_, geomI)
    {
        const rollGeometryInfo& geom = rollGeometries_[geomI];
        const triSurf& rollSurf = geom.surface();

        triSurfModifier sMod(rollSurf);
        geometricSurfacePatchList& patches = sMod.patchesAccess();

        //- get the patches
        std::map<label, label> patchIdMapping;
        forAll(rollSurf, triI)
        {
            //- get the new patch for the triangle
            const label newRegion = triRegion[startTriangle+triI];
            const label newPatch = newPatchForRegion[newRegion];

            if( patchIdMapping.find(newPatch) == patchIdMapping.end() )
            {
                //- create mapping
                const word& newName = newPatchNames[newPatch];

                //- check if the patch already exists
                label localPatchId(-1);
                forAll(patches, patchI)
                {
                    if( patches[patchI].name() == newName )
                    {
                        localPatchId = patchI;
                        break;
                    }
                }

                if( localPatchId >= 0 )
                {
                    patchIdMapping[newPatch] = localPatchId;
                }
                else
                {
                    const label nPatches = patches.size();
                    patches.setSize(nPatches+1);

                    patchIdMapping[newPatch] = nPatches;

                    patches[nPatches].name() = newPatchNames[newPatch];
                    patches[nPatches].geometricType() = newPatchTypes[newPatch];
                }
            }

            //- assign new patch index to a triangle
            sMod.facetsAccess()[triI].region() = patchIdMapping[newPatch];

            # ifdef DEBUGContact
            Info << "Surf " << geomI << " triangle " << triI
                 << " is in patch " << patches[rollSurf[triI].region()].name()
                 << " region patch " << newPatchNames[newPatch] << endl;
            # endif
        }

        //- update the starting triangle
        startTriangle += rollSurf.size();
    }

    if( writeRollProfile_ )
    {
        label startPoint(0);

        forAll(rollGeometries_, geomI)
        {
            const rollGeometryInfo& geom = rollGeometries_[geomI];

            triSurf newSurf(geom.surface());
            triSurfModifier sm(newSurf);
            pointField& pts = sm.pointsAccess();

            forAll(pts, pI)
            {
                pts[pI] = points[startPoint+pI];
                const scalar helper = pts[pI].z();
                pts[pI].z() = pts[pI].x();
                pts[pI].x() = -helper;
            }

            startPoint += pts.size();

            newSurf.writeSurface
            (
                "updatedPatches_"+geom.rollPosition()+".stl"
            );
        }
    }

    Info << "Updated surface patches" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollingMillContactHandler::rollingMillContactHandler
(
    const dictionary& meshDict,
    const PtrList<rollGeometryInfo>& rollGeometries,
    const rollingMillPatchNamesHandler& namesHandler
)
:
    meshDict_(meshDict),
    rollGeometries_(rollGeometries),
    patchNamesHandler_(namesHandler),
    writeRollProfile_(false),
    contactDistanceTolerance_(0.0),
    crossSectionEdges_(),
    pointMapping_(),
    crossSectionSurfPtr_(NULL),
    manifold_(),
    rollPositionToManifold_(),
    octreePtr_(NULL)
{
    if( !readGeometryData() )
    {
        FatalError << "Cannot read all necessary data to evaluate the contact"
            << exit(FatalError);
    }

    checkForContactRegions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

rollingMillContactHandler::~rollingMillContactHandler()
{
    deleteDemandDrivenData(octreePtr_);
    deleteDemandDrivenData(crossSectionSurfPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
