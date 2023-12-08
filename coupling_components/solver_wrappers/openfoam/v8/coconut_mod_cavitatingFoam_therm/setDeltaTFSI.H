/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    setDeltaT

Description
    Reset the timestep to maintain a constant maximum courant Number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

if (adjustTimeStep && runTime.timeIndex() > 5)
{
    Info<< "Adjusting the time step" <<endl;

    const scalar minDeltaT = 1e-8;
    const scalar deltaTFact = 0.5;

    scalar newDeltaT = deltaTFact*runTime.deltaTValue();
    
    if (newDeltaT < minDeltaT)
    
    {
	Info<< "**Setting deltaT to the minDeltaT" <<endl;
	newDeltaT = minDeltaT;
    }
   
   runTime.setDeltaT(newDeltaT);

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
    Info<< "Desired time step scaling factor = "<< deltaTFact << endl;
}

// ************************************************************************* //
