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

#include "wireBlockMeshGenerator.H"
#include "rollingMillMesh.H"
#include "blockMesh.H"
#include "blockDescriptor.H"
#include "curvedEdgeList.H"
#include "arcEdge.H"

#include "cellShape.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

blockWriter::blockWriter()
:
    type_(),
    nodes_(),
    resolution_(),
    grading_()
{}

blockWriter::blockWriter
(
    const word& shapeName,
    const DynList<label>& nodes,
    const Vector<label>& resolution,
    const vector& grading
)
:
    type_(shapeName),
    nodes_(nodes),
    resolution_(resolution),
    grading_(grading)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

blockWriter::~blockWriter()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void blockWriter::operator=(const blockWriter& bw)
{
    type_ = bw.type_;
    nodes_ = bw.nodes_;
    resolution_ = bw.resolution_;
    grading_ = bw.grading_;
}

bool blockWriter::operator!=(const blockWriter& bw) const
{
    if( type_ != bw.type_ )
        return true;

    if( nodes_ != bw.nodes_ )
        return true;

    if( resolution_ != bw.resolution_ )
        return true;

    if( grading_ != bw.grading_ )
        return true;

    return false;
}

Ostream& operator<<(Ostream& os, const blockWriter& bw)
{
    os.check("Ostream& operator<<(Ostream& os, const blockWriter& bw)");

    os << bw.type_ << token::SPACE << bw.nodes_
       << token::SPACE << bw.resolution_
       << " simpleGrading " << bw.grading_ << nl;

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

arcEdgeWriter::arcEdgeWriter()
:
    type_(),
    start_(),
    end_(),
    midPoint_()
{}

arcEdgeWriter::arcEdgeWriter
(
    const word& shapeName,
    const label start,
    const label end,
    const point& midPoint
)
:
    type_(shapeName),
    start_(start),
    end_(end),
    midPoint_(midPoint)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

arcEdgeWriter::~arcEdgeWriter()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void arcEdgeWriter::operator=(const arcEdgeWriter& ew)
{
    type_ = ew.type_;
    start_ = ew.start_;
    end_ = ew.end_;
    midPoint_ = ew.midPoint_;
}

bool arcEdgeWriter::operator!=(const arcEdgeWriter& ew) const
{
    if( type_ != ew.type_ )
        return true;

    if( start_ != ew.start_ )
        return true;

    if( end_ != ew.end_ )
        return true;

    if( midPoint_ != ew.midPoint_ )
        return true;

    return false;
}

Ostream& operator<<(Ostream& os, const arcEdgeWriter& ew)
{
    os.check("Ostream& operator<<(Ostream& os, const arcEdgeWriter& ew)");

    os << ew.type_ << token::SPACE
       << ew.start_ << token::SPACE
       << ew.end_ << token::SPACE
       << ew.midPoint_ << nl;

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namedDictionary::namedDictionary()
:
    name_(),
    dict_()
{}

namedDictionary::namedDictionary
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    dict_(dict)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namedDictionary::~namedDictionary()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void namedDictionary::operator=(const namedDictionary& nd)
{
    name_ = nd.name_;
    dict_ = nd.dict_;
}

bool namedDictionary::operator!=(const namedDictionary& nd) const
{
    if( name_ != nd.name_ )
        return true;

    if( dict_ != nd.dict_ )
        return true;

    return false;
}

Ostream& operator<<(Ostream& os, const namedDictionary& nd)
{
    os.check("Ostream& operator<<(Ostream& os, const namedDictionary& nd)");

    os << nd.name_ << token::SPACE << nd.dict_ << nl;

    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void wireBlockMeshGenerator::createAxiSymmetricMeshDictionary()
{
    const scalar boundCos =
        max(min(geometricTolerance_ / (radius_ + VSMALL), 1.0), 0.0);
    scalar alpha = Foam::acos(1.0 - boundCos);

    //- wedge angle in each direction
    const scalar wedgeAngle = wireGeom_.wedgeAngle() * M_PI / 360.0;
    if( alpha < wedgeAngle )
    {
        Warning << "Wedge angle " << wireGeom_.wedgeAngle() << " is greater "
            << "than " << wireGeom_.wedgeAngle()
            << ", needed to achieve geometry tolerance " << endl;
    }

    alpha = wedgeAngle;

    //- patch names
    const word contactPatchName =
        patchHandler_.patchNameForWire
        (
            rollingMillPatchNamesHandler::WIRECONTACT
        );
    const word frontPatchName =
        patchHandler_.patchNameForWire
        (
            rollingMillPatchNamesHandler::WIREFRONT
        );
    const word backPatchName =
        patchHandler_.patchNameForWire
        (
            rollingMillPatchNamesHandler::WIREBACK
        );

    //- create points
    pointField points(6);
    const label nPoints = 3;

    points[0] = vector::zero;
    points[1] = point(radius_ * sin(alpha), radius_ * cos(alpha), 0.);
    points[2] = point(-radius_ * sin(alpha), radius_ * cos(alpha), 0.);

    for(label i=0;i<nPoints;++i)
        points[i+nPoints] = points[i] + point(0., 0., 0.1*radius_);

    //- create edges
    DynList<arcEdgeWriter> edges;

    edges.append(arcEdgeWriter("arc", 1, 2, point(0., radius_, 0.)));
    edges.append(arcEdgeWriter("arc", 4, 5, point(0., radius_, 0.1 * radius_)));

    //- create blocks
    label nRadialDivisions = ceil(0.5 * radius_ / (mag(points[2].x())+VSMALL));

    if( wireGeom_.wireDict().found("contactCellSize") )
    {
        scalar cs;
        wireGeom_.wireDict().readIfPresent("contactCellSize", cs);

        nRadialDivisions = ceil(radius_ / (cs+VSMALL));
    }


    FixedList<blockWriter, 1> blocks;
    DynList<label> nodes(8);
    nodes[0] = 0;
    nodes[1] = 0;
    nodes[2] = 1;
    nodes[3] = 2;
    nodes[4] = 3;
    nodes[5] = 3;
    nodes[6] = 4;
    nodes[7] = 5;
    blocks[0] =
        blockWriter
        (
            "hex",
            nodes,
            Vector<label>(1, nRadialDivisions, 1),
            vector(1., 1., 1.)
        );

    //- create patches
    List<namedDictionary> patches(5);

    //- bottom empty patch
    dictionary bottomEmptyFaces;
    bottomEmptyFaces.add("type", "empty");
    face bf(4);
    bf[0] = 0;
    bf[1] = 0;
    bf[2] = 2;
    bf[3] = 1;
    bottomEmptyFaces.add("faces", List<face>(1, bf));
    patches[0] = namedDictionary("bottomEmptyFaces", bottomEmptyFaces);

    //- top empty faces
    dictionary topEmptyFaces;
    topEmptyFaces.add("type", "empty");
    bf[0] = 3;
    bf[1] = 3;
    bf[2] = 4;
    bf[3] = 5;
    topEmptyFaces.add("faces", List<face>(1, bf));
    patches[1] = namedDictionary("topEmptyFaces", topEmptyFaces);

    //- contact patch
    dictionary contactFaces;
    contactFaces.add("type", "patch");
    bf[0] = 2;
    bf[1] = 1;
    bf[2] = 4;
    bf[3] = 5;
    contactFaces.add("faces", List<face>(1, bf));
    patches[2] = namedDictionary(contactPatchName, contactFaces);

    //- wedge patch on the back side
    dictionary wedgeFaces;
    wedgeFaces.add("type", "wedge");
    faceList faces(1);
    bf[0] = 3;
    bf[1] = 5;
    bf[2] = 2;
    bf[3] = 0;
    faces[0] = bf;

    wedgeFaces.add("faces", faces);
    patches[3] = namedDictionary(backPatchName, wedgeFaces);

    //- wedge patch on the front side
    wedgeFaces.remove("faces");
    faces.setSize(2);
    bf[0] = 3;
    bf[1] = 4;
    bf[2] = 1;
    bf[3] = 0;
    faces[0] = bf;

    bf[0] = 0;
    bf[1] = 0;
    bf[2] = 3;
    bf[3] = 3;
    faces[1] = bf;
    wedgeFaces.add("faces", faces);
    patches[4] = namedDictionary(frontPatchName, wedgeFaces);

    //- insert values in dictionary
    blockMeshDict_.add("convertToMeters", 1.0);

    //- insert points
    blockMeshDict_.add("vertices", points);

    //- insert blocks
    blockMeshDict_.add("blocks", blocks);

    //- insert curved edges
    blockMeshDict_.add("edges", edges);

    //- insert patches
    blockMeshDict_.add("boundary", patches);
}

void wireBlockMeshGenerator::createMeshDictionary()
{
    const scalar pi_2 = 0.5 * M_PI;
    const scalar pi_4 = 0.5 * pi_2;
    const scalar pi_8 = 0.5 * pi_4;

    //- calculate the number of subdivisions
    scalar alpha;

    if( wireGeom_.wireDict().found("contactCellSize") )
    {
        scalar cs;
        wireGeom_.wireDict().readIfPresent("contactCellSize", cs);

        alpha = cs / (radius_ + VSMALL);
    }
    else
    {
        Warning << "contactCellSize is not given in wireMeshDict. Cell size "
                << "will be determined based on geometryTolerance" << endl;

        alpha =
            2.0 * Foam::acos(1.0 - geometricTolerance_ / (radius_ + VSMALL));
    }

    const label nSplits = ceil(pi_4 / (alpha + VSMALL));

    //- patch names
    const word contactPatchName =
        patchHandler_.patchNameForWire
        (
            rollingMillPatchNamesHandler::WIRECONTACT
        );
    const word symmYName =
        patchHandler_.patchNameForWire
        (
            rollingMillPatchNamesHandler::WIRESYMMY
        );
    const word symmZName =
        patchHandler_.patchNameForWire
        (
            rollingMillPatchNamesHandler::WIRESYMMZ
        );

    //- reference points
    pointField points(34);
    points[16] = vector::zero;

    FixedList<label, 4> radiusOnAxes(-1);
    FixedList<label, 4> edgeCentres(-1);
    FixedList<label, 4> onAxes(-1);
    FixedList<label, 4> inCentre(-1);

    for(label i=0;i<4;++i)
    {
        const scalar rx = radius_ * cos(i * pi_2);
        const scalar ry = radius_ * sin(i * pi_2);

        const scalar rcx = radius_ * cos(i * pi_2 + pi_4);
        const scalar rcy = radius_ * sin(i * pi_2 + pi_4);

        //- create a point on the axis at outer radius
        radiusOnAxes[i] = 4 * i;
        points[4*i] = point(rx, ry, 0.);

        //- create points at the half of the axis
        onAxes[i] = 4 * i + 1;
        points[4 * i + 1] = point(0.5 * rx, 0.5 * ry, 0.);

        //- create edge centres
        edgeCentres[i] = 4 * i + 2;
        points[4 * i + 2] = point(rcx, rcy, 0.);

        //- create a vertex in the centre of the quadrant
        inCentre[i] = 4 * i + 3;
        points[4 * i + 3] = point(0.6 * rcx, 0.6 * rcy, 0.);
    }

    const label nPoints = 17;

    for(label i=0;i<nPoints;++i)
    {
        point p = points[i];
        p.z() = 0.1 * radius_;
        points[i+17] = p;
    }

    //- calculate the number of active quadrants
    FixedList<bool, 4> activeQuadrant(true);

    if
    (
        !(symmetryType_ & rollingMillGeometryHandler::NEGY) &&
        (symmetryType_ & rollingMillGeometryHandler::POSY))
    {
        //- deactivate quadrant on the negative size of the y axis
        activeQuadrant[2] = activeQuadrant[3] = false;
    }

    if
    (
        !(symmetryType_ & rollingMillGeometryHandler::POSY) &&
        (symmetryType_ & rollingMillGeometryHandler::NEGY)
    )
    {
        //- deactivate quadrant on the position side of the y axis
        activeQuadrant[0] = activeQuadrant[1] = false;
    }

    if
    (
        !(symmetryType_ & rollingMillGeometryHandler::POSZ) &&
        (symmetryType_ & rollingMillGeometryHandler::NEGZ)
    )
    {
        //- deactivate quadrants on the positive side of the x axis
        activeQuadrant[0] = activeQuadrant[3] = false;
    }

    if
    (
        (symmetryType_ & rollingMillGeometryHandler::POSZ) &&
        !(symmetryType_ & rollingMillGeometryHandler::NEGZ)
    )
    {
        //- deactivate qudrants on the negative side of the x axis
        activeQuadrant[1] = activeQuadrant[2] = false;
    }

    //- create curved edges
    FixedList<arcEdgeWriter, 16> edges;

    for(label i=0;i<4;++i)
    {
        const scalar sp_8 = sin(i*pi_2+pi_8);
        const scalar cs_8 = cos(i*pi_2+pi_8);
        const scalar sp_3_8 = sin(i*pi_2+3.*pi_8);
        const scalar cs_3_8 = cos(i*pi_2+3.*pi_8);

        edges[4*i] =
            arcEdgeWriter
            (
                "arc",
                radiusOnAxes[i],
                edgeCentres[i],
                point(radius_*cs_8, radius_*sp_8, 0.)
            );

        edges[4*i+1] =
            arcEdgeWriter
            (
                "arc",
                edgeCentres[i],
                radiusOnAxes[(i+1)%4],
                point(radius_*cs_3_8, radius_*sp_3_8, 0.)
            );

        edges[4*i+2] =
            arcEdgeWriter
            (
                "arc",
                radiusOnAxes[i]+17,
                edgeCentres[i]+17,
                point(radius_*cs_8, radius_*sp_8, 0.1*radius_)
            );

        edges[4*i+3] =
            arcEdgeWriter
            (
                "arc",
                edgeCentres[i]+17,
                radiusOnAxes[(i+1)%4]+17,
                point(radius_*cs_3_8, radius_*sp_3_8, 0.1*radius_)
            );
    }

    //- create a list of blocks
    DynList<blockWriter> blocks;

    Vector<label> nDivisions(nSplits, nSplits, 1);
    vector grading(1.0, 1.0, 1.0);

    typedef std::map<word, std::pair<word, DynList<face> > > patchFacesMap;
    patchFacesMap patchFaces;

    for(label i=0;i<4;++i)
    {
        if( !activeQuadrant[i] )
            continue;

        //- first block
        DynList<label> nodes(8);
        nodes[0] = 16;
        nodes[1] = onAxes[i];
        nodes[2] = inCentre[i];
        nodes[3] = onAxes[(i+1)%4];
        for(label j=0;j<4;++j)
            nodes[j+4] = nodes[j] + 17;
        blocks.append
        (
            blockWriter
            (
                "hex",
                nodes,
                nDivisions,
                grading
            )
        );

        //- boundary faces for the first block
        face bf(4);
        bf[0] = nodes[0];
        bf[1] = nodes[3];
        bf[2] = nodes[2];
        bf[3] = nodes[1];
        patchFaces["bottomEmptyFaces"].second.append(bf);

        bf[0] = nodes[4];
        bf[1] = nodes[5];
        bf[2] = nodes[6];
        bf[3] = nodes[7];
        patchFaces["topEmptyFaces"].second.append(bf);

        if( !activeQuadrant[(i+3)%4] )
        {
            //- previous quadrant is not active. Create a symmetry plane B.C.
            bf[0] = nodes[0];
            bf[1] = nodes[1];
            bf[2] = nodes[5];
            bf[3] = nodes[4];

            if( i % 2 == 0 )
            {
                patchFaces[symmYName].second.append(bf);
                patchFaces[symmYName].first = "symmetryPlane";
            }
            else
            {
                patchFaces[symmZName].second.append(bf);
                patchFaces[symmZName].first = "symmetryPlane";
            }
        }

        if( !activeQuadrant[(i+1)%4] )
        {
            //- next quadrant is not active. Create a symmetry plane B.C.
            bf[0] = nodes[0];
            bf[1] = nodes[4];
            bf[2] = nodes[7];
            bf[3] = nodes[3];

            if( i % 2 == 0 )
            {
                patchFaces[symmZName].second.append(bf);
                patchFaces[symmZName].first = "symmetryPlane";
            }
            else
            {
                patchFaces[symmYName].second.append(bf);
                patchFaces[symmYName].first = "symmetryPlane";
            }
        }

        //- second block
        nodes[0] = inCentre[i];
        nodes[1] = onAxes[i];
        nodes[2] = radiusOnAxes[i];
        nodes[3] = edgeCentres[i];
        for(label j=0;j<4;++j)
            nodes[j+4] = nodes[j] + 17;
        blocks.append
        (
            blockWriter
            (
                "hex",
                nodes,
                nDivisions,
                grading
            )
        );

        //- boundary faces for the second block
        bf[0] = nodes[0];
        bf[1] = nodes[3];
        bf[2] = nodes[2];
        bf[3] = nodes[1];
        patchFaces["bottomEmptyFaces"].second.append(bf);

        bf[0] = nodes[4];
        bf[1] = nodes[5];
        bf[2] = nodes[6];
        bf[3] = nodes[7];
        patchFaces["topEmptyFaces"].second.append(bf);

        bf[0] = nodes[2];
        bf[1] = nodes[3];
        bf[2] = nodes[7];
        bf[3] = nodes[6];
        patchFaces[contactPatchName].second.append(bf);

        if( !activeQuadrant[(i+3)%4] )
        {
            //- previous quadrant is not active. Create a symmetry plane B.C.
            bf[0] = nodes[1];
            bf[1] = nodes[2];
            bf[2] = nodes[6];
            bf[3] = nodes[5];

            if( i % 2 == 0 )
            {
                patchFaces[symmYName].second.append(bf);
                patchFaces[symmYName].first = "symmetryPlane";
            }
            else
            {
                patchFaces[symmZName].second.append(bf);
                patchFaces[symmZName].first = "symmetryPlane";
            }
        }

        //- third block
        nodes[0] = inCentre[i];
        nodes[1] = edgeCentres[i];
        nodes[2] = radiusOnAxes[(i+1)%4];
        nodes[3] = onAxes[(i+1)%4];
        for(label j=0;j<4;++j)
            nodes[j+4] = nodes[j] + 17;
        blocks.append
        (
            blockWriter
            (
                "hex",
                nodes,
                nDivisions,
                grading
            )
        );

        //- boundary faces for the third block
        bf[0] = nodes[0];
        bf[1] = nodes[3];
        bf[2] = nodes[2];
        bf[3] = nodes[1];
        patchFaces["bottomEmptyFaces"].second.append(bf);
        patchFaces["bottomEmptyFaces"].first = "empty";

        bf[0] = nodes[4];
        bf[1] = nodes[5];
        bf[2] = nodes[6];
        bf[3] = nodes[7];
        patchFaces["topEmptyFaces"].second.append(bf);
        patchFaces["topEmptyFaces"].first = "empty";

        bf[0] = nodes[1];
        bf[1] = nodes[2];
        bf[2] = nodes[6];
        bf[3] = nodes[5];
        patchFaces[contactPatchName].second.append(bf);
        patchFaces[contactPatchName].first = "patch";

        if( !activeQuadrant[(i+1)%4] )
        {
            //- next qudrant is not active. Create a symmetry plane B.C.
            bf[0] = nodes[2];
            bf[1] = nodes[3];
            bf[2] = nodes[7];
            bf[3] = nodes[6];

            if( i % 2 == 0 )
            {
                patchFaces[symmZName].second.append(bf);
                patchFaces[symmZName].first = "symmetryPlane";
            }
            else
            {
                patchFaces[symmYName].second.append(bf);
                patchFaces[symmYName].first = "symmetryPlane";
            }
        }
    }

    //- create patch faces
    PtrList<namedDictionary> patches;

    patches.setSize(patchFaces.size());

    label patchI(0);
    forAllConstIter(patchFacesMap, patchFaces, it)
    {
        const DynList<face>& faces = it->second.second;

        dictionary patchDict(it->first);
        patchDict.add("type", it->second.first);
        patchDict.add("faces", faces);
        patchDict.name() = it->first;

        patches.set(patchI, new namedDictionary(it->first, patchDict));
        ++patchI;
    }

    //- insert values in dictionary
    blockMeshDict_.add("convertToMeters", 1.0);

    //- insert points
    blockMeshDict_.add("vertices", points);

    //- insert blocks
    blockMeshDict_.add("blocks", blocks);

    //- insert curved edges
    blockMeshDict_.add("edges", edges);

    //- insert patches
    blockMeshDict_.add("boundary", patches);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

wireBlockMeshGenerator::wireBlockMeshGenerator
(
    const Time& runTime,
    const rollingMillPatchNamesHandler& patchHandler,
    const wireGeometryInfo& wireGeom,
    const word regionName,
    const scalar tol
)
:
    runTime_(runTime),
    patchHandler_(patchHandler),
    radius_(0.5 * wireGeom.wireDiameter()),
    regionName_(regionName),
    symmetryType_(wireGeom.typeOfSymmetry()),
    geometricTolerance_(tol),
    wireGeom_(wireGeom),
    blockMeshDict_
    (
        IOobject
        (
            "blockMeshDict",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    )
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

wireBlockMeshGenerator::~wireBlockMeshGenerator()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

polyMeshGen* wireBlockMeshGenerator::generateCrossSectionMesh()
{
    Info << "Generating wire cross-section" << endl;

    //- patch names
    const word contactPatchName =
        patchHandler_.patchNameForWire
        (
            rollingMillPatchNamesHandler::WIRECONTACT
        );

    //- create a dictionary for the block mesher
    if( symmetryType_ & rollingMillGeometryHandler::AXISYMMETRIC )
    {
        createAxiSymmetricMeshDictionary();
    }
    else
    {
        createMeshDictionary();
    }

    //- create a block mesh
    blockMesh bm(blockMeshDict_, "wire");

    polyMesh mesh
    (
        IOobject
        (
            "region0",
            runTime_.constant(),
            runTime_
        ),
        xferCopy(bm.points()),
        bm.cells(),
        bm.patches(),
        bm.patchNames(),
        bm.patchDicts(),
        contactPatchName,
        "patch"
    );

    //- copy the mesh
    polyMeshGen* resultMeshPtr =
        new polyMeshGen(runTime_, "constant", regionName_/"polyMesh");
    polyMeshGenModifier(*resultMeshPtr).addMesh(mesh);

    return resultMeshPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
