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

#include "dieSurfaceCreatorDXF.H"
#include "rollingMillMesh.H"
#include "demandDrivenData.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "DxfFile.H"
#include "DxfFileParser.H"
#include "DxfFileToTriSurfConverter.H"
#include "IFstream.H"
#include "boundBox.H"
#include "addToRunTimeSelectionTable.H"

//#define DEBUGDie

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(dieSurfaceCreatorDXF, 0);
addToRunTimeSelectionTable
(
    dieSurfaceCreator,
    dieSurfaceCreatorDXF,
    dictionary
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dieSurfaceCreatorDXF::parseDXFFile()
{
    fileName fName;
    dict_.readIfPresent("dxfFile", fName);

    if( fName.ext() != "dxf" || fName.ext() != "DXF" )
    {
        FatalError << "Filename " << fName << " is not a dxf file"
            << exit(FatalError);
    }

    //
    // Opening input file
    //

    std::ifstream ifs(fName, std::ios::binary);
    if (!ifs.is_open())
    {
        Info << "ERROR: Cannot open file \"" << fName << "\"" << endl;
        std::abort();
    }

    //
    // Parsing input file
    //

    DxfFile file;

    try
    {
        DxfFileParser parser;
        parser.Parse(ifs, file);
    }
    catch (const std::exception& e)
    {
        Info << "ERROR: Error while parsing DXF file \"" << fName << "\": "
             << e.what() << endl;
        std::abort();
    }

    if (file.GetUnits() != DxfFile::Units::Millimeters)
    {
        Warning << "THE UNITS IN THE FILE \"" << fName << "\" ARE NOT mm."
                << " CONTINUING WITH THE ASSUMPTION THEY ARE PROVIDED IN mm!!"
                << endl;

        file.SetUnits(DxfFile::Units::Millimeters);
    }

    //
    // Converting DXF file to triSurf
    //

    triSurf surf;

    try
    {
        DxfFileToTriSurfConverter converter(file, 0.1);
        converter.Convert(surf);
    }
    catch (const std::exception& e)
    {
        Info << "ERROR: Error while converting DXF file "
             << "\"" << fName << "\" to triSurf: " << e.what() << endl;
        std::abort();
    }

    scalar minX(VGREAT), maxX(-VGREAT);
    label minId(-1), maxId(-1);

    forAll(surf.points(), pI)
    {
        const point& p = surf.points()[pI];

        if( mag(p.z()) > SMALL )
        {
            FatalErrorIn
            (
                "void dieSurfaceCreatorDXF::parseDXFFile()"
            ) << "Geometry is not in the x-y plane!" << exit(FatalError);
        }

        if( p.x() < minX )
        {
            minX = p.x();
            minId = pI;
            inletRadius_ = p.y();
        }

        if( p.x() > maxX )
        {
            maxX = p.x();
            maxId = pI;
            outletRadius_ = p.y();
        }
    }

    if( !dict_.readIfPresent("outerDiameter", outerRadius_) )
    {
        FatalError << "outerDiameter is not present in dictionary " << dict_
            << exit(FatalError);
    }
    outerRadius_ *= 0.5;

    //- create edge subsets that shall be extruded into patches
    DynList<label> edgeSubsetIDs;
    surf.edgeSubsetIndices(edgeSubsetIDs);
    forAll(edgeSubsetIDs, i)
        surf.removeEdgeSubset(edgeSubsetIDs[i]);
    edgeSubsetIDs.clear();

    edgeSubsetIDs.append
    (
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEWIRE
            )
        )
    );

    forAll(surf.featureEdges(), eI)
        surf.addEdgeToSubset(edgeSubsetIDs.lastElement(), eI);

    //- create downstream patch
    const label nPointsBefore = surf.nPoints();
    const label nEdgesBefore = surf.nFeatureEdges();

    //- add additional vertices
    surf.appendVertex(point(maxX, outerRadius_, 0.0));
    surf.appendVertex(point(minX, outerRadius_, 0.0));

    //- create missing edges
    surf.appendFeatureEdge(edge(maxId, nPointsBefore));
    surf.appendFeatureEdge(edge(nPointsBefore, nPointsBefore+1));
    surf.appendFeatureEdge(edge(nPointsBefore+1, minId));

    //- add new edges into subsets
    //- downstream
    edgeSubsetIDs.append
    (
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEDOWNSTREAM
            )
        )
    );

    surf.addEdgeToSubset(edgeSubsetIDs.lastElement(), nEdgesBefore);

    //- die to housing
    edgeSubsetIDs.append
    (
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEHOUSING
            )
        )
    );

    surf.addEdgeToSubset(edgeSubsetIDs.lastElement(), nEdgesBefore+1);

    //- die upstream
    edgeSubsetIDs.append
    (
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEUPSTREAM
            )
        )
    );

    surf.addEdgeToSubset(edgeSubsetIDs.lastElement(), nEdgesBefore+2);

    //- extrude feature edges to a 2D surface
    triSurfaceExtrude2DEdges extruder(surf);

    //- extruded edges into a ribbon
    surfPtr_ = new triSurf();
    extruder.extrudeSurface(*surfPtr_);

    //- create patches
    geometricSurfacePatchList& sPatches =
        triSurfModifier(*surfPtr_).patchesAccess();

    forAll(sPatches, patchI)
        sPatches[patchI].geometricType() = "patch";
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dieSurfaceCreatorDXF::dieSurfaceCreatorDXF
(
    const rollingMillPatchNamesHandler& patchHandler,
    const dictionary& dict,
    const scalar tol
)
:
    dieSurfaceCreator(patchHandler, dict, tol),
    inletRadius_(0.0),
    outletRadius_(0.0),
    outerRadius_(0.0)
{
    if( dict.found("dxfFile") )
    {
        parseDXFFile();
    }
    else
    {
        FatalError << "dxfFile is not found in ditionary " << dict
            << exit(FatalError);
    }
}

dieSurfaceCreatorDXF::dieSurfaceCreatorDXF
(
    const dieSurfaceCreatorDXF& die
)
:
    dieSurfaceCreator
    (
        die.patchHandler_,
        die.dict_,
        die.geometryTol_
    ),
    inletRadius_(die.inletRadius_),
    outletRadius_(die.outletRadius_),
    outerRadius_(die.outerRadius_)
{
    surfPtr_ = new triSurf(die.surface());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar dieSurfaceCreatorDXF::inletDiameter() const
{
    return 2. * inletRadius_;
}

scalar dieSurfaceCreatorDXF::outletDiameter() const
{
    return 2. * outletRadius_;
}

scalar dieSurfaceCreatorDXF::outerDiameter() const
{
    return 2. * outerRadius_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
