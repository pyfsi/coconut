/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    cfMeshSmoother

Description
    Mesh smoothing using the cfMesh improveMeshQuality utility.

    The volFields can be advected using an explicit advection transport equation
    where the mesh Courant number can be limited.

    Current issue for parallel in 2-D: a processor may not have cells in the
    subSet, causing it to give an error when calculating the empty directions.

SourceFiles
    cfMeshSmoother.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef cfMeshSmoother_H
#define cfMeshSmoother_H

#include "meshSmoother.H"
#include "newFvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMeshGen;

/*---------------------------------------------------------------------------* \
                         Class cfMeshSmoother Declaration
\*---------------------------------------------------------------------------*/

class cfMeshSmoother
:
    public meshSmoother
{
    // Private Data

        //- Sub-mesh containing the cells to be smoothed
        mutable autoPtr<newFvMeshSubset> subMeshToSmooth_;

        //- cfMesh improveMeshQuality smoothing settings
        class improveMeshQualityCoeffs
        {
            public:

            const int nLoops_;
            const int nIterations_;
            const int nSurfaceIterations_;
            const scalar qualityThreshold_;
            const word constrainedCellsSet_;
            const int nIterationsNearBoundaries_;
            const int nLayers_;
            const Switch smoothSurface_;
            const Switch optimizeMeshFV_;
            const Switch noSurfaceProjectionOnFinalIteration_;

            improveMeshQualityCoeffs
            (
                const int nLoops,
                const int nIterations,
                const int nSurfaceIterations,
                const scalar qualityThreshold,
                const word& constrainedCellsSet,
                const int nIterationsNearBoundaries_,
                const int nLayers_,
                const Switch& smoothSurface,
                const Switch& optimizeMeshFV = false,
                const Switch& noSurfaceProjectionOnFinalIteration = false
            )
            :
                nLoops_(nLoops),
                nIterations_(nIterations),
                nSurfaceIterations_(nSurfaceIterations),
                qualityThreshold_(qualityThreshold),
                constrainedCellsSet_(constrainedCellsSet),
                nIterationsNearBoundaries_(nIterationsNearBoundaries_),
                nLayers_(nLayers_),
                smoothSurface_(smoothSurface),
                optimizeMeshFV_(optimizeMeshFV),
                noSurfaceProjectionOnFinalIteration_
                (
                    noSurfaceProjectionOnFinalIteration
                )
            {}
        } iMQCoeffs_;


    // Private Member Functions

        //- Return a reference to the sub mesh to be smoothed
        fvMesh& subMeshToSmooth();

        //- Return a const reference to the sub mesh to be smoothed
        const fvMesh& subMeshToSmooth() const;

        //- Return a reference to the subset mesh to be smoothed
        newFvMeshSubset& subsetMeshToSmooth();

        //- Return a const reference to the subset mesh to be smoothed
        const newFvMeshSubset& subsetMeshToSmooth() const;

        //- Return a reference to the sub mesh to be smoothed
        void makeSubMeshToSmooth() const;

        //- Calculate mesh point motion Courant number
        void calculatePointMotionCourantNo
        (
            const pointVectorField& pointMotionD,
            const fvMesh& mesh
        );

        //- Smooth surface of mesh
        void smoothSurface(polyMeshGen& pmg, fvMesh& mesh) const;

        //- Calculate areas swept by the edges on the patch
        // void calculateSweptAreas
        // (
        //     PtrList<edgeScalarField>& sweptAreas,
        //     const PtrList<faMesh>& aMeshes,
        //     const fvMesh& mesh,
        //     const pointField& oldPoints
        // ) const;

        //- Disallow default bitwise copy construct
        cfMeshSmoother
        (
            const cfMeshSmoother&
        );

        //- Disallow default bitwise assignment
        void operator=(const cfMeshSmoother&);

public:

    //- Runtime type information
    TypeName("cfMeshImproveMeshQuality");

    // Static data members


    // Constructors

        //- Construct from dictionary
        cfMeshSmoother
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~cfMeshSmoother();


    // Member Functions

        //- Clear out demand driven data
        virtual void clearOut();

        //- Smoothing function
        virtual scalar smooth();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
