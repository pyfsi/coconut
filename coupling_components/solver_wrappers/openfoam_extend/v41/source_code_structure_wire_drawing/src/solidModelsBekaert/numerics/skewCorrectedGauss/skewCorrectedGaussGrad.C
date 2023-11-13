/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "skewCorrectedGaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "zeroGradientFvPatchField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "processorFvPatch.H"
#include "ggiFvPatch.H"
#include "regionCoupleFvPatch.H"
#include "fvc.H"
#include "skewCorrectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
skewCorrectedGaussGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vf.mesh();
    
    GeometricField<Type, fvsPatchField, surfaceMesh> sf =
        linear<Type>(mesh).interpolate(vf);
    
    tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
    (
        gaussGrad<Type>::gradf(sf, name)
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();

    // Skew corretion
    surfaceVectorField scv = skewCorrectionVectors::New(vf.mesh())();
    
    forAll(scv.boundaryField(), patchI)
    {
        if ( isA<ggiFvPatch>(mesh.boundary()[patchI]) )
        {
            scv.boundaryField()[patchI] = vector::zero;
        }
    }

    sf += (scv & linear<GradType>(mesh).interpolate(gGrad));

    gGrad = gaussGrad<Type>::gradf(sf, name);
    
    gGrad.rename("grad(" + vf.name() + ')');
    this->correctBoundaryConditions(vf, gGrad);

    return tgGrad;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
