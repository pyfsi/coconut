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

#include "dieSurfaceCreatorProfiledConical.H"
#include "rollingMillMesh.H"
#include "triSurfModifier.H"
#include "triSurfaceExtrude2DEdges.H"
#include "triSurfaceChecks.H"
#include "addToRunTimeSelectionTable.H"
#include "dieCrossSectionSurfaceCreator.H"

//#define DEBUGDie

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(dieSurfaceCreatorProfiledConical, 0);
addToRunTimeSelectionTable
(
    dieSurfaceCreator,
    dieSurfaceCreatorProfiledConical,
    dictionary
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dieSurfaceCreatorProfiledConical::createTriangulatedSurface()
{
    triSurf surf;
    triSurfModifier sMod(surf);

    pointField& pts = sMod.pointsAccess();

    pts.setSize(4);

    pts[0] = point(0.0, outletRadius_, 0.0);
    pts[1] = point(0.0, outerRadius_, 0.0);
    pts[2] = point(-axialLength_, outerRadius_, 0.0);
    pts[3] = point(-axialLength_, inletRadius_, 0.0);

    //- generate feature edge that will extruded into a ribbon
    surf.appendFeatureEdge(edge(0, 1));
    surf.appendFeatureEdge(edge(1, 2));
    surf.appendFeatureEdge(edge(2, 3));
    surf.appendFeatureEdge(edge(3, 0));

    const label downstreamId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEDOWNSTREAM
            )
        );

    const label outerId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEHOUSING
            )
        );

    const label inletId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEUPSTREAM
            )
        );

    const label contactId =
        surf.addEdgeSubset
        (
            patchHandler_.patchNameForDie
            (
                rollingMillPatchNamesHandler::DIEWIRE
            )
        );

    //- add edges to subset
    surf.addEdgeToSubset(downstreamId, 0);
    surf.addEdgeToSubset(outerId, 1);
    surf.addEdgeToSubset(inletId, 2);
    surf.addEdgeToSubset(contactId, 3);

    surfPtr_ = new triSurf();

    triSurfaceExtrude2DEdges extruder(surf);
    extruder.extrudeSurface(*surfPtr_);

    geometricSurfacePatchList& sPatches =
        triSurfModifier(*surfPtr_).patchesAccess();

    forAll(sPatches, patchI)
        sPatches[patchI].geometricType() = "patch";

    # ifdef DEBUGDie
    Info << "Number of surface triangles " << surfPtr_->size() << endl;
    # endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dieSurfaceCreatorProfiledConical::dieSurfaceCreatorProfiledConical
(
    const rollingMillPatchNamesHandler& patchHandler,
    const dictionary& dict,
    const scalar tol
)
:
    dieSurfaceCreator(patchHandler, dict, tol),
    inletRadius_(0.0),
    outletRadius_(0.0),
    axialLength_(0.0),
    outerRadius_(0.0)
{
    //- read axial length
    if( !dict.readIfPresent("axialLength", axialLength_) )
    {
        FatalError << "axialLength not found in dictionary "
                << dict << exit(FatalError);
    }

    //- read outer diameter of a die
    if( !dict.readIfPresent("outerDiameter", outerRadius_) )
    {
        FatalError << "outerDiameter not found in dictionary "
                << dict << exit(FatalError);
    }

    //- get the contact cell size to determine discreatisation of the profile
    scalar contactCellSize = 0.01 * outerRadius_;
    dict.readIfPresent("contactCellSize", contactCellSize);

    if( dict.found("minContactCellSize") )
    {
        scalar mccs;
        dict.readIfPresent("minContactCellSize", mccs);

        contactCellSize = min(contactCellSize, mccs);
    }

    contactCellSize *= 0.25;

    outerRadius_ *= 0.5;

    //- read the profile at the inlet cross-section
    fileName inletProfileName;

    if( dict.found("inletDiameter") )
    {
        scalar d;
        dict.readIfPresent("inletDiameter", d);

        std::shared_ptr<triSurf> surfPtr = createCircularProfile(d);

        crossSections_.append
        (
            std::pair<scalar, std::shared_ptr<triSurf> >
            (
                -axialLength_,
                surfPtr
            )
        );

        inletRadius_ = 0.5 * d;

        if( inletRadius_ >= outerRadius_ )
        {
            FatalErrorIn
            (
                "dieSurfaceCreatorProfiledConical::"
                "dieSurfaceCreatorProfiledConical("
                "const rollingMillPatchNamesHandler&,"
                " const dictionary&, const scalar tol)"
            ) << "inletDiameter of a profile"
              << " is greater than the outer diameter of a die "
              << exit(FatalError);
        }
    }
    else if( dict.readIfPresent("inletProfile", inletProfileName) )
    {
        dieCrossSectionSurfaceCreator inletCreator
        (
            inletProfileName,
            geometryTol_,
            contactCellSize
        );

        std::shared_ptr<triSurf> surfPtr =
            std::make_shared<triSurf>(inletCreator.surface());

        positionCrossSectionAtOrigin(*surfPtr);

        if
        (
            triSurfaceChecks::checkSurfaceManifolds(*surfPtr) != 1 ||
            triSurfaceChecks::checkForNonManifoldEdges(*surfPtr)
        )
        {
            surfPtr->writeSurface(inletProfileName.lessExt()+"_invalid.fms");

            FatalError << "Inlet profile is not valid. "
                << "It is written to "
                << inletProfileName.lessExt()+"_invalid.fms"
                << exit(FatalError);
        }

        crossSections_.append
        (
            std::pair<scalar, std::shared_ptr<triSurf> >
            (
                -axialLength_,
                surfPtr
            )
        );

        //- get the average radius of the inlet profile
        const pointField& pts = surfPtr->points();
        scalar maxRadius(0.0), minRadius(VGREAT);
        forAll(pts, pI)
        {
            point p = pts[pI];
            p.x() = 0.0;

            scalar r = mag(p);

            maxRadius = max(maxRadius, r);
            minRadius = min(minRadius, r);
        }

        if( maxRadius >= outerRadius_ )
        {
            FatalErrorIn
            (
                "dieSurfaceCreatorProfiledConical::"
                "dieSurfaceCreatorProfiledConical("
                "const rollingMillPatchNamesHandler&, const dictionary&,"
                " const scalar tol)"
            ) << "Maximum diameter of a profile given in " << inletProfileName
              << " is greater than the outer radius of a die "
              << exit(FatalError);
        }

        inletRadius_ = 0.5 * (maxRadius + minRadius);
    }
    else
    {
        FatalError << "inletProfile or inletDiameter not found in dictionary "
                << dict << exit(FatalError);
    }

    //- read intermediate profiles
    if
    (
        dict.found("intermediateProfiles") &&
        dict.isDict("intermediateProfiles")
    )
    {
        const dictionary& profilesDict = dict.subDict("intermediateProfiles");

        const wordList names = profilesDict.toc();

        forAll(names, i)
        {
            if( profilesDict.isDict(names[i]) )
            {
                const dictionary& dict = profilesDict.subDict(names[i]);

                scalar axialPosition;
                if( !dict.readIfPresent("axialPosition", axialPosition) )
                {
                    FatalErrorIn
                    (
                        "dieSurfaceCreatorProfiledConical::"
                        "dieSurfaceCreatorProfiledConical("
                        "const rollingMillPatchNamesHandler&,"
                        " const dictionary&, const scalar)"
                    ) << "axialPosition is not available for the intermediate"
                      << " profile " << names[i] << exit(FatalError);
                }

                fileName fName;
                if( !dict.readIfPresent("profile" , fName) )
                {
                    FatalErrorIn
                    (
                        "dieSurfaceCreatorProfiledConical::"
                        "dieSurfaceCreatorProfiledConical("
                        "const rollingMillPatchNamesHandler&,"
                        " const dictionary&, const scalar)"
                    ) << "profile is not available for the intermediate"
                      << " profile " << names[i] << exit(FatalError);
                }

                //- read the profile from the file
                dieCrossSectionSurfaceCreator profileCreator
                (
                    fName,
                    geometryTol_,
                    contactCellSize
                );

                std::shared_ptr<triSurf> surfPtr =
                    std::make_shared<triSurf>(profileCreator.surface());

                positionCrossSectionAtOrigin(*surfPtr);

                if
                (
                    triSurfaceChecks::checkSurfaceManifolds(*surfPtr) != 1 ||
                    triSurfaceChecks::checkForNonManifoldEdges(*surfPtr)
                )
                {
                    surfPtr->writeSurface(names[i]+"_invalid.fms");

                    FatalError << "Intermediate profile " << names[i]
                        << " is not valid. "
                        << "It is written to " << names[i]+"_invalid.fms"
                        << exit(FatalError);
                }

                crossSections_.append
                (
                    std::pair<scalar, std::shared_ptr<triSurf> >
                    (
                        axialPosition-axialLength_,
                        surfPtr
                    )
                );

                //- get the average radius of the inlet profile
                const pointField& pts = surfPtr->points();
                scalar maxRadius(0.0);
                forAll(pts, pI)
                {
                    point p = pts[pI];
                    p.x() = 0.0;

                    scalar r = mag(p);

                    maxRadius = max(maxRadius, r);
                }

                if( maxRadius >= outerRadius_ )
                {
                    FatalErrorIn
                    (
                        "dieSurfaceCreatorProfiledConical::"
                        "dieSurfaceCreatorProfiledConical("
                        "const rollingMillPatchNamesHandler&,"
                        " const dictionary&, const scalar tol)"
                    ) << "Maximum diameter of a profile given in " << fName
                      << " is greater than the outer radius of a die "
                      << exit(FatalError);
                }
            }
            else
            {
                FatalErrorIn
                (
                    "dieSurfaceCreatorProfiledConical::"
                    "dieSurfaceCreatorProfiledConical("
                    "const rollingMillPatchNamesHandler&, const dictionary&,"
                    " const scalar tol)"
                ) << "Intermediate profile " << names[i] << " is not given "
                  << "in a dictionary" << exit(FatalError);
            }
        }
    }

    //- read the profile at the outlet side
    fileName outletProfileName;
    if( !dict.readIfPresent("outletProfile", outletProfileName) )
    {
        FatalError << "outletProfile not found in dictionary "
                << dict << exit(FatalError);
    }
    else
    {
        dieCrossSectionSurfaceCreator outletCreator
        (
            outletProfileName,
            geometryTol_,
            contactCellSize
        );

        std::shared_ptr<triSurf> surfPtr =
            std::make_shared<triSurf>(outletCreator.surface());

        positionCrossSectionAtOrigin(*surfPtr);

        if
        (
            triSurfaceChecks::checkSurfaceManifolds(*surfPtr) != 1 ||
            triSurfaceChecks::checkForNonManifoldEdges(*surfPtr)
        )
        {
            surfPtr->writeSurface(outletProfileName.lessExt()+"_invalid.fms");

            FatalError << "Outlet profile " << outletProfileName
                << " is not valid. "
                << "It is written to "
                << outletProfileName.lessExt()+"_invalid.fms"
                << exit(FatalError);
        }

        crossSections_.append
        (
            std::pair<scalar, std::shared_ptr<triSurf> >
            (
                0.0,
                surfPtr
            )
        );

        //- get the average radius of the inlet profile
        const pointField& pts = surfPtr->points();
        scalar maxRadius(0.0), minRadius(VGREAT);
        forAll(pts, pI)
        {
            point p = pts[pI];
            p.x() = 0.0;

            scalar r = mag(p);

            maxRadius = max(maxRadius, r);
            minRadius = min(minRadius, r);
        }

        if( maxRadius >= outerRadius_ )
        {
            FatalErrorIn
            (
                "dieSurfaceCreatorProfiledConical::"
                "dieSurfaceCreatorProfiledConical("
                "const rollingMillPatchNamesHandler&, const dictionary&,"
                " const scalar tol)"
            ) << "Maximum diameter of a profile given in " << outletProfileName
              << " is greater than the outer radius of a die "
              << exit(FatalError);
        }

        outletRadius_ = 0.5 * (maxRadius + minRadius);
    }

    //- sort axial positions of profiles in the ascending order
    bool changed;
    do
    {
        changed = false;

        for(label i=0;i<crossSections_.size()-1;++i)
        {
            if( crossSections_[i].first > crossSections_[i+1].first )
            {
                const auto copy = crossSections_[i].first;
                crossSections_[i].first = crossSections_[i+1].first;
                crossSections_[i+1].first = copy;
            }
        }
    } while( changed );

    # ifdef DEBUGDie
    forAll(crossSections_, i)
    {
        crossSections_[i].second->writeSurface
        (
            "section_"+std::to_string(i)+".stl"
        );
    }
    # endif

    createTriangulatedSurface();
}

dieSurfaceCreatorProfiledConical::dieSurfaceCreatorProfiledConical
(
    const dieSurfaceCreatorProfiledConical& die
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
    axialLength_(die.axialLength_),
    outerRadius_(die.outerRadius_)
{
    //- copy cross sections if there are any
    forAll(die.crossSections_, i)
    {
        crossSections_.append
        (
            std::pair<scalar, std::shared_ptr<triSurf> >
            (
                die.crossSections_[i].first,
                std::make_shared<triSurf>(*die.crossSections_[i].second)
            )
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar dieSurfaceCreatorProfiledConical::inletDiameter() const
{
    return 2. * inletRadius_;
}

scalar dieSurfaceCreatorProfiledConical::outletDiameter() const
{
    return 2. * outletRadius_;
}

scalar dieSurfaceCreatorProfiledConical::outerDiameter() const
{
    return 2.0 * outerRadius_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
