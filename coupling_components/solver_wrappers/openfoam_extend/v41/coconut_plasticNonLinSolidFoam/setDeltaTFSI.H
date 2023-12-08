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

Global
    setDeltaT

Description
    Set the time-step based on the plastic strain increment.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep && runTime.timeIndex() > 5)
{
    Info<< "Adjusting the time step" << endl;
    const scalar minDeltaTHere = 1e-8;

    // Calculate the time-step scaling factor and limit it so as it does not
    // grow or shrink too quickly
    // MATHIEU modifications
    const scalar deltaTFact =
        max
        (
            min
            (
                mechanical.newDeltaT()/runTime.deltaTValue(),
                0.5
            ),
            0.5
        );

    // Set the time-step and limit it to the specified max and min values

    scalar newDeltaT = deltaTFact*runTime.deltaTValue();

    if (newDeltaT > maxDeltaT)
    {
        Info<< "**Setting deltaT to the maxDeltaT" << endl;
        newDeltaT = maxDeltaT;
    }
    else if (newDeltaT < minDeltaTHere)
    {
        Info<< "**Setting deltaT to the minDeltaT" << endl;
        newDeltaT = minDeltaTHere;
    }

    runTime.setDeltaT(newDeltaT);

    Info<< "Time step = " <<  runTime.deltaT().value() << nl
        << "Minimum time step = "<< minDeltaTHere << nl
        << "Maximum time step = "<< maxDeltaT << nl
        << "Desired time step scaling factor = "<< deltaTFact << nl
        << nl << nl << endl;
}

// ************************************************************************* //
