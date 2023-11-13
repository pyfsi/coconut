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
    Reads the AVL's surface mesh

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurf.H"
#include "DxfFile.H"
#include "DxfFileParser.H"
#include "DxfFileToTriSurfConverter.H"
#include "triSurfaceExtrude2DEdges.H"
#include "triSurfaceCleanupDuplicates.H"
#include "meshOctree.H"
#include "meshOctreeCreator.H"

#include <fstream>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input dxf file");
    argList::validArgs.append("output surface file");

    argList::validOptions.insert("exportToVTK", "");
    argList::validOptions.insert("discretisationAngle", "scalar");
    argList::validOptions.insert("discretisationLength", "scalar");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);

    scalar angleTol(1.0);
    if( args.options().found("discretisationAngle") )
    {
        const scalar ang
        (
            readScalar(IStringStream(args.options()["discretisationAngle"])())
        );

        angleTol = ang;
    }
    else
    {
        Info << "Using the default 1 deg angle!" << endl;
    }

    scalar discretisationLength(VGREAT);
    if( args.options().found("discretisationLength") )
    {
        discretisationLength =
            readScalar(IStringStream(args.options()["discretisationLength"])());
    }
    else
    {
        Info << "Using default VGREAT discretisation length!" << endl;
    }

    if( inFileName.ext() != "dxf" && inFileName.ext() != "DXF"  )
    {
        Info << "Input file does not have a dxf extension!! Exitting!!" << endl;
        return 0;
    }

    //
    // Opening input file
    //

    std::ifstream ifs(inFileName, std::ios::binary);
    if (!ifs.is_open())
    {
        Info << "ERROR: Cannot open file \"" << inFileName << "\"" << endl;
        return 1;
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
        Info << "ERROR: Error while parsing: " << e.what() << endl;
        return 1;
    }

    if (file.GetUnits() != DxfFile::Units::Millimeters)
    {
        Warning << "THE UNITS IN THE FILE \"" << inFileName << "\" ARE NOT mm."
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
        DxfFileToTriSurfConverter converter(file, angleTol, discretisationLength);
        converter.Convert(surf);
    }
    catch (const std::exception& e)
    {
        Info << "ERROR: Error while converting DXF file to triSurf: "
             << e.what() << endl;
        return 1;
    }

    if( args.options().found("exportToVTK") )
    {
        const pointField& pts = surf.points();
        const edgeLongList& featureEdges = surf.featureEdges();

        OFstream file("profileEdges.vtk");

        //- write the header
        file << "# vtk DataFile Version 3.0\n";
        file << "vtk output\n";
        file << "ASCII\n";
        file << "DATASET POLYDATA\n";

        //- write points
        file << "\nPOINTS " << pts.size() << " float\n";
        forAll(pts, pI)
        {
            const point& p = pts[pI];

            file << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
        }

        //- write lines
        file << "\nLINES " << featureEdges.size()
             << " " << 3*featureEdges.size() << nl;
        forAll(featureEdges, eI)
        {
            const edge& e = featureEdges[eI];

            file << 2 << " " << e.start() << " " << e.end() << nl;
        }

        file << "\n";
        file.flush();
    }

    //- extrude lines into a surface triangulation
    triSurfaceExtrude2DEdges extruder(surf);

    const triSurf* extrudedSurfPtr = extruder.extrudeSurface();
    if( !extrudedSurfPtr )
    {
        FatalError << "Failed extruding a 2D mesh from the file "
            << inFileName << endl;

        return 1;
    }

    //- merge duplicate points
    meshOctree mo(*extrudedSurfPtr, true);
    meshOctreeCreator(mo).createOctreeWithRefinedBoundary(20, 20);
    triSurfaceCleanupDuplicates cleaner(mo);
    cleaner.mergeIdentities();

    //- write the surface
    extrudedSurfPtr->writeSurface(outFileName);
    delete extrudedSurfPtr;

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //
