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

\*---------------------------------------------------------------------------*/

#include "readFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Mesh, class GeoField>
void Foam::readFields
(
    const Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeoField>& fields
)
{
    // Search list of objects for volScalarFields
    IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));

    // Remove the cellDist field
    IOobjectList::iterator celDistIter = fieldObjects.find("cellDist");
    if (celDistIter != fieldObjects.end())
    {
        fieldObjects.erase(celDistIter);
    }

    // Remove fields not needed for next pass of rolling or drawing
    wordList ignoreFields(21, "null");
    ignoreFields[0] = "DEpsilonPEq_mat0";
    ignoreFields[1] = "activeYield_mat0";
    ignoreFields[2] = "DEpsilonPEqf";
    ignoreFields[3] = "DU_0";
    ignoreFields[4] = "U_0";
    ignoreFields[5] = "U_0_0";
    ignoreFields[6] = "activeYieldf";
    ignoreFields[7] = "pointDU";
    ignoreFields[8] = "rho_0";
    ignoreFields[9] = "rho_0_0";
    ignoreFields[10] = "DEpsilonP_mat0";
    ignoreFields[11] = "bEbarTrial_mat0";
    ignoreFields[12] = "pointVelocity";
    ignoreFields[13] = "DEpsilonPf";
    ignoreFields[14] = "Velocity";
    ignoreFields[15] = "bEbarTrialf";
    ignoreFields[16] = "epsilonTrue";
    ignoreFields[17] = "sigmaCauchy";
    ignoreFields[18] = "U";
    ignoreFields[19] = "pointU";
    ignoreFields[20] = "pointT";
//    ignoreFields[21] = "T";

    forAll(ignoreFields, fieldI)
    {
        IOobjectList::iterator fieldIter =
            fieldObjects.find(ignoreFields[fieldI]);

        if (fieldIter != fieldObjects.end())
        {
            Info<< "Not mapping field " << ignoreFields[fieldI] << endl;
            fieldObjects.erase(fieldIter);
        }
    }

    // Construct the vol scalar fields
    fields.setSize(fieldObjects.size());

    label fieldi=0;
    for
    (
        IOobjectList::iterator iter = fieldObjects.begin();
        iter != fieldObjects.end();
        ++iter
    )
    {
        fields.set
        (
            fieldi++,
            new GeoField
            (
                *iter(),
                mesh
            )
        );

        // Set all fields to write
        // except for oldTimes of DU and U
        fields[fieldi - 1].writeOpt() = IOobject::AUTO_WRITE;

        if
        (
            fields[fieldi - 1].name() == "U"
         || fields[fieldi - 1].name() == "DU"
         || fields[fieldi - 1].name() == "U_0"
         || fields[fieldi - 1].name() == "U_0_0"
         || fields[fieldi - 1].name() == "DU_0"
         || fields[fieldi - 1].name() == "rho"
         || fields[fieldi - 1].name() == "rho_0"
        )
        {
            fields[fieldi - 1].oldTime().oldTime().writeOpt() =
                IOobject::NO_WRITE;
            fields[fieldi - 1].oldTime().writeOpt() = IOobject::NO_WRITE;
        }
    }
}


// ************************************************************************* //
