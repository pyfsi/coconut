/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

Class
    Foam::pimpleControl

Description
    Coconut class based on the PIMPLE control class to supply convergence information/checks for
    the PIMPLE loop.

    May also be used to for PISO-based algorithms as PISO controls are a
    sub-set of PIMPLE controls.

\*---------------------------------------------------------------------------*/

#ifndef coconutPimpleControl_H
#define coconutPimpleControl_H

#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class coconutPimpleControl Declaration
\*---------------------------------------------------------------------------*/

class coconutPimpleControl: public pimpleControl
{
    // Private member functions
private:
    //- No copy construct
    coconutPimpleControl(const coconutPimpleControl&) = delete;

    //- No copy assignment
    void operator=(const coconutPimpleControl&) = delete;

public:
    using pimpleControl::pimpleControl;

    bool callCriteriaSatisfied() {
        // modified version of pimpleControl::criteriaSatisfied()
        // no checks on first iteration - nothing has been calculated yet
        // if ((corr_ == 1) || residualControl_.empty() || finalIter())
        // {
        //     return false;
        // }

        const bool storeIni = this->storeInitialResiduals();

        bool achieved = true;
        bool checked = false;    // safety that some checks were indeed performed

        const dictionary& solverDict = mesh_.data().solverPerformanceDict();
        for (const entry& solverPerfDictEntry : solverDict)
        {
            const word& fieldName = solverPerfDictEntry.keyword();
            const label fieldi = applyToField(fieldName);

            if (fieldi != -1)
            {
                Pair<scalar> residuals = maxResidual(solverPerfDictEntry);

                checked = true;

                scalar relative = 0.0;
                bool relCheck = false;

                const bool absCheck =
                    (residuals.last() < residualControl_[fieldi].absTol);

                if (storeIni)
                {
                    residualControl_[fieldi].initialResidual = residuals.first();
                }
                else
                {
                    const scalar iniRes =
                        (residualControl_[fieldi].initialResidual + ROOTVSMALL);

                    relative = residuals.last() / iniRes;
                    relCheck = (relative < residualControl_[fieldi].relTol);
                }

                achieved = achieved && (absCheck || relCheck);

                if (debug)
                {
                    Info<< algorithmName_ << " loop:" << endl;

                    Info<< "    " << fieldName
                        << " PIMPLE iter " << corr_
                        << ": ini res = "
                        << residualControl_[fieldi].initialResidual
                        << ", abs tol = " << residuals.last()
                        << " (" << residualControl_[fieldi].absTol << ")"
                        << ", rel tol = " << relative
                        << " (" << residualControl_[fieldi].relTol << ")"
                        << endl;
                }
            }
        }

        return checked && achieved;
    }

    // Constructors
    //- Construct from mesh and the name of control sub-dictionary
    coconutPimpleControl(fvMesh& mesh): pimpleControl(mesh, "PIMPLE", true) {};

    //- Destructor
    virtual ~coconutPimpleControl() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
