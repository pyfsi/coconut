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

Class
    profiledDieGeometryInterpolator

Description
    Creates a set of cross-section of a die needed to reconstruct the whole
    3D mesh

SourceFiles
    profiledDieGeometryInterpolator.C

\*---------------------------------------------------------------------------*/

#ifndef profiledDieGeometryInterpolator_H
#define profiledDieGeometryInterpolator_H

#include "polyMeshGen.H"
#include "scalarLongList.H"
#include "vectorLongList.H"
#include "triSurf.H"

#include <memory>
#include <tuple>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class cubicBSpline;
class dieGeometryInfo;
class rollingMillPatchNamesHandler;

/*---------------------------------------------------------------------------*\
                Class profiledDieGeometryInterpolator Declaration
\*---------------------------------------------------------------------------*/

class profiledDieGeometryInterpolator
{
    // Private data
        //- reference to dieGeometryInfo
        const dieGeometryInfo& dieProfiles_;

        //- const reference to patch handler
        const rollingMillPatchNamesHandler& patchHandler_;

        //- reference to the volume mesh
        polyMeshGen& mesh_;

        //- list of points at each x coordinate of the inner profile
        VRWGraph sortedProfilePoints_;

        //- points at fixed outer radius
        labelHashSet pointsAtFixedRadius_;

        //- points at y=0 symmetry plane
        labelHashSet pointsInYSymmetry_;

        //- points at z=0 symmetry plane
        labelHashSet pointsInZSymmetry_;

        //- die radius at each axial position. The information is used
        //- to evaluate scaling factors and interpolation
        std::map<scalar, scalar> radiusAtPosition_;

        //- list of all interpolated profiles
        List<std::tuple<scalar, bool, std::shared_ptr<triSurf> > > profiles_;

        //- coordinate of a profile and its id
        std::map<scalar, label> xToId_;

    // Private member functions
        //- calculate average radius at a given position
        scalar radiusAtPosition(const scalar x) const;

        //- calculate outer radius of a die at a given position
        scalar outerRadiusAtPosition(const scalar x) const;

        //- create a cubic spline from a profile
        std::shared_ptr<cubicBSpline> createCubicSpline(const triSurf&) const;

        //- calculate interpolated profile
        std::shared_ptr<triSurf> interpolateProfiles
        (
            const scalar x,
            const triSurf& profile1Surf,
            const std::pair<scalar, std::shared_ptr<cubicBSpline> >& profile1,
            const std::pair<scalar, std::shared_ptr<cubicBSpline> >& profile2
        ) const;

        //- extract points of the inner profile of a die into
        void detectProfilePointsInVolMesh();

        //- analyse axial surface and detect average radius at each position
        void analyseAxialCrossSection();

        //- calculate cross section at every required position
        void calculateProfilesAtRequiredPositions();

        //- create mapping space for inner mesh vertices
        void initializeInterpolationMatrix
        (
            List<List<List<point> > >& interpolationPoints,
            std::map<scalar, label>& indexOfAxialPosition,
            List<std::pair<scalar, scalar> >& minAndMaxRadiusAtAxialPosition
        ) const;

        //- calculate positions of points of the interpolation matrix
        //- using the elliptic meshing scheme
        void calculateInterpolationMatrixElliptic
        (
            List<List<List<point> > >& sectionPoints,
            const List<std::pair<scalar, scalar> >& minAndMaxRadiusAtSection
        ) const;

        //- calculate positions of points of the interpolation matrix
        //- using the surface optimizer
        void calculateInterpolationMatrixSurfOpt
        (
            List<List<List<point> > >& sectionPoints
        ) const;

    // Private member functions

        //- Disallow default bitwise copy construct
        profiledDieGeometryInterpolator(const profiledDieGeometryInterpolator&);

        //- Disallow default bitwise assignment
        void operator=(const profiledDieGeometryInterpolator&);

public:

    // Constructors

        //- Construct from dieSurfaceCreator and polyMeshGen
        profiledDieGeometryInterpolator
        (
            const dieGeometryInfo&,
            const rollingMillPatchNamesHandler&,
            polyMeshGen&
        );

    // Destructor

        ~profiledDieGeometryInterpolator();


    // Member Functions
        //- available axial positions
        void axialPositions(scalarLongList&) const;

        //- returns the surface mesh of the cross-section at the given axial
        //- position
        const triSurf& crossSection(const scalar) const;

        //- calculate displacements of vertices in the patches making
        //- the inner surface of a die
        void calculateDisplacements
        (
            std::map<label, vector>& pointDisplacements
        ) const;

        //- calculate displacement of all vertices
        void calculateDisplacementsAll(vectorLongList&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
