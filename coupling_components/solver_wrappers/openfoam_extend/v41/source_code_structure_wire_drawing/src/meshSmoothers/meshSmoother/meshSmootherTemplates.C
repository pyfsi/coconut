/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "meshSmoother.H"
//#include "conservativeMeshToMesh.H"
#include "IOobjectList.H"
#include "patchToPatchInterpolation.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// template<class Type>
// void meshSmoother::MapConservativeVolFields
// (
//     const conservativeMeshToMesh& meshToMeshInterp,
//     const label method
// )
// {
//     Info<< "    MapConservativeVolFields" << endl;

//     const fvMesh& meshSource = meshToMeshInterp.srcMesh();
//     const fvMesh& meshTarget = meshToMeshInterp.tgtMesh();

//     word fieldClassName
//     (
//         GeometricField<Type, fvPatchField, volMesh>::typeName
//     );

//     HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields =
//     (
//         mesh_.thisDb().objectRegistry::template lookupClass
//         <GeometricField<Type, fvPatchField, volMesh> >()
//     );

//     for
//     (
//      typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
//             iterator fieldIter = fields.begin();
//         fieldIter != fields.end();
//         ++fieldIter
//     )
//     {
//         //if (debug)
//         {
//             Info<< "    Interpolating " << fieldIter()->name();
//         }

//         // Reference to the field
//         GeometricField<Type, fvPatchField, volMesh>& field =
//             const_cast<GeometricField<Type, fvPatchField, volMesh>&>
//             (*fieldIter());

//         // Take a copy of the patch fields as we perform a separate boundary
//         // interpolation below after interpolating the internal field
//         List<Field<Type> > prevPatchField(field.boundaryField().size());
//         forAll(prevPatchField, patchI)
//         {
//             prevPatchField[patchI] = field.boundaryField()[patchI];
//         }

//         // Compute integral of source field
//         Type intSource = gSum(meshSource.V()*field.internalField());
//         //Info<< "Integral source: " << intSource << endl;

//         // Update the field and make sure the name remains unchanged
//         const word fieldName = field.name();
//         field = meshToMeshInterp.interpolate(field, method);
//         field.rename(fieldName);

//         // Compute integral of mapped field
//         Type intTarget = gSum(meshTarget.V()*field.internalField());
//         //Info<< "Integral target: " << intTarget << endl;

//         //if (debug)
//         {
//             Info<< ", error = "
//                 << 100.0*mag(intSource - intTarget)/(mag(intSource) + SMALL)
//                 << "%" << endl;
//         }

//         // Update the boundary fields because conservativeMeshToMesh does not
//         // interpolate on the boundaries; it justs performs a simple map.
//         // We will use the first order inverse distance method in
//         // the patchToPatchInterpolation class

//         if
//         (
//             meshSource.boundaryMesh().size()
//          != meshTarget.boundaryMesh().size()
//         )
//         {
//             FatalErrorIn
//             (
//                 "template<class Type>\n"
//                 "void meshSmoother::MapConservativeVolFields\n"
//                 "(\n"
//                 "    const conservativeMeshToMesh& meshToMeshInterp,\n"
//                 "    const label method\n"
//                 ")"
//             )   << "The source and target meshes have a different number "
//                 << "of patches!" << abort(FatalError);
//         }

//         forAll(meshSource.boundaryMesh(), patchI)
//         {
//             if (meshSource.boundaryMesh()[patchI].type() != "empty")
//             {
//                 // Create the patch interpolator
//                 // Note: it would be better here to use a globalPolyPatch to
//                 // ensure consistent treatment in parallel: actually the
//               // internal field mapping is currently not parallelised so this
//                 // is fine for serial
//                 patchToPatchInterpolation patchInterp
//                 (
//                     meshSource.boundaryMesh()[patchI], // fromPatch
//                     meshTarget.boundaryMesh()[patchI]  // toPatch
//                 );

//                 field.boundaryField()[patchI] =
//                     patchInterp.faceInterpolate(prevPatchField[patchI]);
//             }
//         }
//     }
// }


template<class Type>
void meshSmoother::AdvectVolFields
(
    const surfaceScalarField& sweptVol,
    const volScalarField& volOld,
    const volScalarField& volNew,
    // const PtrList<areaScalarField>& areasOld,
    // const PtrList<areaScalarField>& areasNew,
    // const PtrList<edgeScalarField>& sweptAreas,
    const wordList& ignoreFields
)
{
    word fieldClassName
    (
        GeometricField<Type, fvPatchField, volMesh>::typeName
    );

    Info<< "    AdvectVolFields: " << fieldClassName << endl;

    // Get a list of the Type volFields
    HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields =
    (
        mesh_.thisDb().objectRegistry::template lookupClass
        <GeometricField<Type, fvPatchField, volMesh> >()
    );

    // Lookup options for treating the boundary patch values
    const Switch extrapolateBoundaries =
        Switch(dict().lookup("extrapolateBoundaries"));
    const int maxCorr =
        dict().lookupOrDefault<int>("extrapolateBoundariesMaxCorr", 10);
     const scalar weight =
        dict().lookupOrDefault<scalar>("extrapolateBoundariesWeight", 0.0);
     const Switch zeroGradientBoundaries =
        Switch(dict().lookup("zeroGradientBoundaries"));
    const Switch frozenBoundaries =
        Switch(dict().lookup("frozenBoundaries"));
    const Switch finiteAreaBoundaryAdvection =
        Switch(dict().lookup("finiteAreaBoundaryAdvection"));

    if
    (
        (
            int(bool(extrapolateBoundaries))
          + int(bool(zeroGradientBoundaries))
          + int(bool(frozenBoundaries))
          + int(bool(finiteAreaBoundaryAdvection))
        ) > 1
    )
    {
        FatalErrorIn(type() + "::AdvectVolFields(...)")
            << "At most, only one of extrapolateBoundaries, "
            << "zeroGradientBoundaries, frozenBoundaries and "
            << "finiteAreaBoundaryAdvection should be set to 'on'"
            << abort(FatalError);
    }

    Info<< "extrapolateBoundaries = " << extrapolateBoundaries << nl
        << "zeroGradientBoundaries = " << zeroGradientBoundaries << nl
        << "frozenBoundaries = " << frozenBoundaries << nl
        << "finiteAreaBoundaryAdvection = " << finiteAreaBoundaryAdvection
        << endl;

    for
    (
        typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        // Check if we should skip this field
        // It might be more efficient to use a hashSet but assuming there are
        // only a small number of fields to ignore, then this is fine
        bool skipField = false;
        forAll(ignoreFields, fI)
        {
            if (ignoreFields[fI] == fieldIter()->name())
            {
                skipField = true;
                break;
            }
        }

        if (skipField)
        {
            continue;
        }

        //if (debug)
        {
            Info<< "        " << fieldIter()->name() << endl;
        }

        // Reference to the field: we use const_cast to update it
        GeometricField<Type, fvPatchField, volMesh>& field =
            const_cast<GeometricField<Type, fvPatchField, volMesh>&>
            (*fieldIter());

        // Take a copy of required boundary fields
        PtrList< Field<Type> > bPrev(field.boundaryField().size());
        if (frozenBoundaries)
        {
            // A copy of the boundary values
            forAll(bPrev, patchI)
            {
                bPrev.set
                (
                    patchI, new Field<Type>(field.boundaryField()[patchI])
                );
            }
        }

        // Explicitly advect the field
        // Note: we are assuming that the field is an intrinsic property i.e.
        // per unit volume
        // Note2: fvc::div divides by the current volume so we multiply by it to
        // cancel it out
        field =
            (1.0/volNew)
           *(
               field*volOld
             - volNew*fvc::div(-sweptVol, field, "advectFields")
           );

        // Extrapolate to the boundaries using the previous snGrad
        if (extrapolateBoundaries)
        {
            // Iteratively extrapolate from the internal domain to the boundary
            // We will perform this component by component
            for
            (
                direction cmptI = 0;
                cmptI < pTraits<Type>::nComponents;
                cmptI++
            )
            {
                // Create scalar field from current component
                volScalarField cmptField("cmptField", field.component(cmptI));

                // Calculate the gradient of the field
                volVectorField gradCmptField = fvc::grad(cmptField);

                // Iteratively extrapolate to the boundary patches and
                // re-calculate the gradient
                int iCorr = 0;
                scalar res = 0.0;
                const scalar tol = 1e-10;
                Info<< "        Iteratively extrapolating boundary values for"
                    << " component = " << cmptI << endl;
                do
                {
                    cmptField.storePrevIter();

                    forAll(cmptField.boundaryField(), patchI)
                    {
                        if (!cmptField.boundaryField()[patchI].coupled())
                        {
                            const fvPatch& ppatch =
                                cmptField.mesh().boundary()[patchI];

                            // Vector from cell-centre to face-centre
                            const vectorField delta = ppatch.delta();

                            // Unit normals
                            const vectorField n = ppatch.nf();

                            // Decompose delta in orthogonal and non-orthogonal
                            // components
                            const vectorField nnd = (sqr(n) & delta);
                            const vectorField k = ((I - sqr(n)) & delta);

                            // Field value at closest cell centre
                            const scalarField pif =
                                cmptField.boundaryField()
                                [
                                    patchI
                                ].patchInternalField();

                            // Gradient value at closest cell centre
                            const vectorField gradPif =
                                gradCmptField.boundaryField()
                                [
                                    patchI
                                ].patchInternalField();

                            // Extrapolate to the boundary face centre with
                            // optional scaling of normal component
                            // extrapolation
                            // For weight = 0, it results in zero-gradient with
                            // non-orthogonal correction
                            cmptField.boundaryField()[patchI] =
                                pif + ((weight*nnd + k) & gradPif);
                        }
                    }

                    // Update gradient field, as it depends on the boundary
                    // values
                    gradCmptField = fvc::grad(cmptField);

                    // Calculate a residual
                    const dimensionedScalar maxPrev =
                        max(mag(cmptField.prevIter()));
                    res =
                    (
                        mag(max(mag(cmptField)) - maxPrev)
                       /(maxPrev.value() + SMALL)
                    ).value();

                    Info<< "    i = " << iCorr << ", res = " << res << endl;
                }
                while (res > tol && ++iCorr < maxCorr);

                // Copy extrapolated cmptField back to the field
                forAll(field.boundaryField(), patchI)
                {
                    if (!field.boundaryField()[patchI].coupled())
                    {
                        field.boundaryField()[patchI].replace
                        (
                            cmptI,
                            cmptField.boundaryField()[patchI]
                        );
                    }
                }
            }
        }
        else if (zeroGradientBoundaries)
        {
            forAll(field.boundaryField(), patchI)
            {
                if (!field.boundaryField()[patchI].coupled())
                {
                    // Zero-gradient extrapolation to the boundary
                    field.boundaryField()[patchI] =
                        field.boundaryField()[patchI].patchInternalField();
                }
            }
        }
        else if (frozenBoundaries)
        {
            forAll(bPrev, patchI)
            {
                if (!field.boundaryField()[patchI].coupled())
                {
                    field.boundaryField()[patchI] = bPrev[patchI];
                }
            }
        }
        else if (finiteAreaBoundaryAdvection)
        {
            FatalErrorIn(type() + "::AdvectVolFields(...)")
                << "finiteAreaBoundaryAdvection is not correct so it is "
                << "disabled!" << abort(FatalError);

            // forAll(field.boundaryField(), patchI)
            // {
            //     if (!field.boundaryField()[patchI].coupled())
            //     {
            //         // Create a faMesh for the patch
            //         const faMesh& aMesh = areasOld[patchI].mesh();

            //         // Create a faMesh field for the current field
            //         GeometricField<scalar, faPatchField, areaMesh> aField
            //         (
            //             IOobject
            //             (
            //                 "area_" + field.name(),
            //                 field.mesh().time().timeName(),
            //                 field.mesh(),
            //                 IOobject::NO_READ,
            //                 IOobject::NO_WRITE
            //             ),
            //             aMesh,
            //             dimensioned<scalar>
            //             (
            //                 "zero", field.dimensions(), pTraits<scalar>::zero
            //             )
            //         );

            //         // Take a reference to the area fields
            //         const areaScalarField& areaOld = areasOld[patchI];
            //         const areaScalarField& areaNew = areasNew[patchI];
            //         const edgeScalarField& sweptArea = sweptAreas[patchI];

            //         // Advect each component
            //         // Note: we do this component-wise because the normal
            //         // component of tensors will not be adveccted correctly
            //         // otherwise (fac::div returns zero)
            //         for
            //         (
            //             direction cmptI = 0;
            //             cmptI < pTraits<Type>::nComponents;
            //             cmptI++
            //         )
            //         {
            //             // Transfer field boundary values into aField
            //             aField.internalField() =
            //                 field.boundaryField()[patchI].component(cmptI);

            //             // Advect aField
            //             aField =
            //                 (1.0/areaNew)
            //                *(
            //                     aField*areaOld
            //            - areaNew*fac::div(-sweptArea, aField, "advectFields")
            //                 );

            //             field.boundaryField()[patchI].replace
            //             (
            //                 cmptI,
            //                 aField.internalField()
            //             );
            //         }
            //     }
            //}
        }

        // Force symmetry plane values for tensors and symmTensors as they are
        // not mapped correctly
        if
        (
            fieldClassName
         == GeometricField<tensor, fvPatchField, volMesh>::typeName
         || fieldClassName
         == GeometricField<symmTensor, fvPatchField, volMesh>::typeName
         || fieldClassName
         == GeometricField<sphericalTensor, fvPatchField, volMesh>::typeName
        )
        {
            forAll(field.boundaryField(), patchI)
            {
                if
                (
                    field.mesh().boundaryMesh()[patchI].type()
                 == "symmetryPlane"
                )
                {
                    //Info<< "        Corrrecting symmetryPlane = "
                    //    << field.mesh().boundaryMesh()[patchI].name() << endl;
                    field.boundaryField()[patchI] =
                        field.boundaryField()[patchI].patchInternalField();
                }
            }
        }
    }
}


// template<class Type>
// void meshSmoother::CorrectTensorNormalComponentVolFields
// (
//     const surfaceScalarField& sweptVol,
//     const volScalarField& volOld,
//     const volScalarField& volNew,
//     const PtrList<areaScalarField>& areasOld,
//     const PtrList<areaScalarField>& areasNew,
//     const PtrList<edgeScalarField>& sweptAreas,
//     const wordList& ignoreFields
// )
// {
//     const Switch finiteAreaBoundaryAdvection =
//         Switch(dict().lookup("finiteAreaBoundaryAdvection"));

//     if (!finiteAreaBoundaryAdvection)
//     {
//         return;
//     }

//     word fieldClassName
//     (
//         GeometricField<Type, fvPatchField, volMesh>::typeName
//     );

//     Info<< "    CorrectTensorNormalComponentVolFields: " << fieldClassName << endl;

//     // Get a list of the Type volFields
//     HashTable<const GeometricField<Type, fvPatchField, volMesh>*> fields =
//     (
//         mesh_.thisDb().objectRegistry::template lookupClass
//         <GeometricField<Type, fvPatchField, volMesh> >()
//     );

//     for
//     (
//         typename HashTable<const GeometricField<Type, fvPatchField, volMesh>*>::
//             iterator fieldIter = fields.begin();
//         fieldIter != fields.end();
//         ++fieldIter
//     )
//     {
//         // Check if we should skip this field
//         // It might be more efficient to use a hashSet but assuming there are
//         // only a small number of fields to ignore, then this is fine
//         bool skipField = false;
//         forAll(ignoreFields, fI)
//         {
//             if (ignoreFields[fI] == fieldIter()->name())
//             {
//                 skipField = true;
//                 break;
//             }
//         }

//         if (skipField)
//         {
//             continue;
//         }

//         //if (debug)
//         {
//             Info<< "        " << fieldIter()->name() << endl;
//         }

//         // Reference to the field: we use const_cast to update it
//         GeometricField<Type, fvPatchField, volMesh>& field =
//             const_cast<GeometricField<Type, fvPatchField, volMesh>&>
//             (*fieldIter());

//         // For all patches, correct normal component
//         forAll(field.boundaryField(), patchI)
//         {
//             const tensorField sqrN =
//                 tensorField
//                 (
//                     sqr(field.mesh().boundary()[patchI].nf())
//                   & tensor(1,0,0,0,1,0,0,0,1)
//                 );

//             field.boundaryField()[patchI] -=
//                 transform(sqrN, field.boundaryField()[patchI]);
//             field.boundaryField()[patchI] +=
//                 transform
//                 (
//                     sqrN, field.boundaryField()[patchI].patchInternalField()
//                 );
//         }
//     }
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
