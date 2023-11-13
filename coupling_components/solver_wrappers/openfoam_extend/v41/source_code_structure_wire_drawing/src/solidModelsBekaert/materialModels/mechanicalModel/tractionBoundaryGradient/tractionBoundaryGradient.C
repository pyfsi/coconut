/*---------------------------------------------------------------------------*\
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

Class
    tractionBoundaryGradient

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "tractionBoundaryGradient.H"
#include "fvPatch.H"
#include "Switch.H"
#include "mechanicalModel.H"
#include "thermalModel.H"


// * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::tractionBoundaryGradient::traction
(
    const tensorField& gradField,
    const word& workingFieldName,
    const word& integralFieldName,
    const fvPatch& patch,
    const bool incremental
)
{
    FatalErrorIn
    (
        "tmp<vectorField> tractionBoundaryGradient::traction\n"
        "(\n"
        "    const tensorField& gradField,\n"
        "    const word& workingFieldName,\n"
        "    const word& integralFieldName,\n"
        "    const fvPatch& patch,\n"
        "    const bool orthotropic,\n"
        "    const nonLinearGeometry::nonLinearType& nonLinear\n"
        ") const"
    )   << "tractionBoundaryGradient::traction is deprecated"
        << exit(FatalError);

    // Create result to keep compiler happy
    tmp<vectorField> ttraction
    (
        new vectorField(gradField.size(), vector::zero)
    );

    return ttraction;
}


// * * * * * * * * * * * * * * * * Operators  * * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::tractionBoundaryGradient::snGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const word& workingFieldName,
    const word& integralFieldName,
    const fvPatch& patch,
    const bool incremental
)
{
    // Lookup nonLinear type from solver
    nonLinearGeometry::nonLinearType nonLinear =
        nonLinearGeometry::nonLinearNames_.read
        (
            patch.boundaryMesh().mesh().solutionDict().subDict
            (
                "solidMechanics"
            ).lookup("nonLinear")
        );

    // Lookup if solver uses an orthotropic approach
    const Switch orthotropic =
        patch.boundaryMesh().mesh().solutionDict().subDict
        (
            "solidMechanics"
        ).lookupOrDefault<Switch>("orthotropic", false);

    // Create result
    tmp<vectorField> tgradient(new vectorField(traction.size(), vector::zero));
    vectorField& gradient = tgradient();

    if (orthotropic)
    {
        FatalErrorIn
        (
            "tmp<vectorField> tractionBoundaryGradient::traction\n"
            "(\n"
            "    const tensorField& gradField,\n"
            "    const word& workingFieldName,\n"
            "    const word& integralFieldName,\n"
            "    const fvPatch& patch,\n"
            "    const bool orthotropic,\n"
            "    const nonLinearGeometry::nonLinearType& nonLinear\n"
            ") const"
        )   << "To be re-implemented for orthotropic elastic"
            << exit(FatalError);
    }
    else
    {
        // Patch unit normals
        const vectorField n = patch.nf();

        // Lookup shear modulus from the solver
        const fvPatchScalarField& mu =
            patch.lookupPatchField<volScalarField, scalar>("mu");

        // Lookup gradient of the working field
        const fvPatchTensorField& gradField =
            patch.lookupPatchField<volTensorField, tensor>
            (
                "grad(" + workingFieldName + ")"
            );

	//Info << "workingFieldName " << workingFieldName << nl << endl;
        if (nonLinear == nonLinearGeometry::OFF)
        {
            // Lookup sigma stress field in the solver; this should be kept
            // up-to-date in the outer iteration loop
            const fvPatchField<symmTensor>& sigma =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "sigma"
                );

            // Calculate traction to be applied
            const vectorField totalTraction = (traction - n*pressure);

            if (patch.boundaryMesh().mesh().foundObject<volScalarField>("p"))
            {
                gradient =
                    (
                        totalTraction
                      - (n & sigma)
                      + (2.0*mu)*(n & gradField)
                    )/(2.0*mu);
            }
            else
            {
                // Lookup first Lame parameter
                const fvPatchScalarField& lambda =
                    patch.lookupPatchField<volScalarField, scalar>("lambda");

                // Set boundary normal gradient such that the user specified
                // traction is enforced
                gradient =
                    (
                        totalTraction
                      - (n & sigma)
                      + (2.0*mu + lambda)*(n & gradField)
                    )/(2.0*mu + lambda);
            }
        }
        else if (nonLinear == nonLinearGeometry::UPDATED_LAGRANGIAN_KIRCHHOFF)
        {
            // Lookup tau (Kirchhoff) stress field in the solver; this should be
            // kept up-to-date in the outer iteration loop
            const fvPatchField<symmTensor>& tau =
                patch.lookupPatchField<volSymmTensorField, symmTensor>
                (
                    "tauKirchhoff"
                );
	//Info << "Tau_beginning_Langrangian " << tau << nl << endl;

            // Relative deformation gradient inverse
            const fvPatchField<tensor>& Finv =
                 patch.lookupPatchField<volTensorField, tensor>("relFinv");

	  //Info << "Finv " << Finv << nl << endl;

            // Total Jacobian relates tau and sigmaCauchy
            const fvPatchField<scalar>& J =
                patch.lookupPatchField<volScalarField, scalar>("J");

            // Calculate current deformed face normals
            // nCurrent should be unit vectors but we will force normalisation
            // to remove any errors
            vectorField nCurrent = Finv.T() & n;
            nCurrent /= mag(nCurrent);

	   //Info << "n" << n << nl << endl;
	 //  Info << "nCurrent " << nCurrent << nl << endl;

            // Calculate Cauchy traction to be applied
            const vectorField tractionCauchy = (traction - nCurrent*pressure);

            // Set boundary normal grdient such that the user specified traction
            // is enforced
            if (patch.boundaryMesh().mesh().foundObject<volScalarField>("p"))
            {
                FatalErrorIn("solidTraction")
	                    << "stop" << abort(FatalError);

                // lambda is not used in pressure approach as it could be
                // infinity
                gradient =
                    (
                        tractionCauchy
                      - (nCurrent & tau/J)
                      + (2.0*mu)*(n & gradField)
                    )/(2.0*mu);
            }
            else
            {
                // Lookup first Lame parameter
                const fvPatchScalarField& lambda =
                    patch.lookupPatchField<volScalarField, scalar>("lambda"); 
		//Info << "gradField " << gradField << endl;
		//Info << "n & gradField " << (n & gradField) << endl;
                gradient =
                    (
                        tractionCauchy
                      - (nCurrent & tau/J)
                      + (2.0*mu + lambda)*(n & gradField)
                    )/(2.0*mu + lambda);
		//Info << "Gradient in tractionBoundaryGradient ########### " << gradient << nl << endl;
            }
        }
        else
        {
            FatalErrorIn
            (
                "tmp<vectorField> tractionBoundaryGradient::snGrad\n"
                "(\n"
                "    const vectorField& traction,\n"
                "    const scalarField& pressure,\n"
                "    const word& workingFieldName,\n"
                "    const word& integralFieldName,\n"
                "    const fvPatch& patch,\n"
                "    const bool incremental\n"
                ") const"
            )   << "unknown nonlinear type: "
                << nonLinearGeometry::nonLinearNames_[nonLinear]
                << abort(FatalError);
        }
    }

    return tgradient;
}


// ************************************************************************* //
