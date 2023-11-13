/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

#include "bicgStabUzawa.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bicgStabUzawa, 0);

    //HJ, temporary
    lduSolver::addsymMatrixConstructorToTable<bicgStabUzawa>
        addbicgStabUzawaSymMatrixConstructorToTable_;

    lduSolver::addasymMatrixConstructorToTable<bicgStabUzawa>
        addbicgStabUzawaAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
Foam::bicgStabUzawa::bicgStabUzawa
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    lduSolver
    (
        fieldName,
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        dict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::lduSolverPerformance Foam::bicgStabUzawa::solve
(
    scalarField& x1,
    const scalarField& b1,
    const direction cmpt
) const
{
    // No constraint in the x direction
    if (cmpt == 0)
    {
        return Foam::bicgStabUzawa::regularSolve(x1, b1, cmpt);
    }
    // system is A*x1 = b1
    // have some constraints B*x2 = b2
    // note the input x corresponds to x1
    // and b corresponds to b1

    // example constraint which sets cellx = celly
    // scalarField B(x.size(),0);

    // int cellx = 50;
    // int celly = 3000;
    // B[cellx] = 1.0;    B[celly] = -1.0;
    // scalar b2=0.0;

    // Create constraint coefficients
    scalarField B(x1.size(), 0.0);

    //  Create constraint source
    //scalarField b2(x.size(), 0.0);
    const scalar b2 = 0.0;

    // To-do:
    // what about preconditioning? The constraint coefficients may be a
    // different order of magnitude to the other coefficients.
    // Also, we could create a run-time selectable constraint class to pass
    // the constraint coeffs and source.

    // Testing
    // Let's try implement a basic contact contraint
    // Bottom cells: 70-79
    // Top cells: 80-89
    // Initially we will force the bottom cells to have the same displacement
    // as the top cels

    //for (int i = 70; i < 80; i++)
    for (int i = 73; i < 77; i++)
    {
        // top U == bottom U
        // Hmnn why do we get different behaviour when we scale the
        // coefficients?!
        B[i] = -0.5;
        B[i + 10] = 0.5;
    }

    Foam::lduSolverPerformance solverPerf;

    // UZAWA iteration procedure
    // need to have an initial guess for the lagrange multipliers
    // which are x2 on the wikipedia page
    scalar x2 = 0.0;


    // We start the conjugate gradient iteration by computing the residual

    Foam::bicgStabUzawa::regularSolve(x1, (b1 - B*x2), cmpt);
    //Info<< "x1 is " << x1 << endl;

    scalar r2 = gSumProd(B, x1) - b2;
    // the first search direction is then chosen
    scalar p2 = r2;
    Info<< nl << "cmpt: " << cmpt << nl
        << "p2: " << p2 << endl;

    // scalarField p1(x1.size(), 0.0);

    // It seems to take quite a lot of iterations and also the convergence seems
    // to stall
    const int maxIts = 10;
    for (int iteration = 0; iteration < maxIts; iteration++)
    {
        // Should p1 be reinitialised to zero each iteration?
        scalarField p1(x1.size(), 0.0);
        const scalar a2 = gSumProd(B, p1);

        solverPerf = Foam::bicgStabUzawa::regularSolve(p1, B*p2, cmpt);

        // scaling factor
        const scalar alpha = (p2*r2)/(p2*a2 + SMALL);

        // updates
        x2 += alpha*p2;
        r2 -= alpha*a2;
        x1 -= alpha*p1;

        // new search direction
        const scalar beta = (r2*a2)/(p2*a2 + SMALL);
        p2 = r2 - beta*p2;

        // Info<< nl << "iteration " << iteration << nl
        //     << "p2: " << p2 << nl
        //     << "r2: " << r2 << endl;

        if ((iteration > 1) && mag(r2) < 1e-6)
        {
            Info<< "Uzawa iters: " << iteration << endl;
            break;
        }
    }

    solverPerf.initialResidual() = r2;

    return solverPerf;
}


Foam::lduSolverPerformance Foam::bicgStabUzawa::regularSolve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // --- Setup class containing solver performance data
    lduSolverPerformance solverPerf
    (
        //lduMatrix::preconditioner::getName(dict()) + typeName,
        typeName,
        fieldName()
    );


    register label nCells = x.size();

    scalar* __restrict__ xPtr = x.begin();

    scalarField pA(nCells);
    scalar* __restrict__ pAPtr = pA.begin();

    scalarField pT(nCells, 0.0);
    scalar* __restrict__ pTPtr = pT.begin();

    scalarField wA(nCells);
    scalar* __restrict__ wAPtr = wA.begin();

    scalarField wT(nCells);
    scalar* __restrict__ wTPtr = wT.begin();

    scalar wArT = matrix_.great_;
    scalar wArTold = wArT;

    // Calculate A.x and T.x
    matrix_.Amul(wA, x, coupleBouCoeffs_, interfaces_, cmpt);
    matrix_.Tmul(wT, x, coupleIntCoeffs_, interfaces_, cmpt);

    // Calculate initial residual and transpose residual fields
    scalarField rA(b - wA);
    scalarField rT(b - wT);
    scalar* __restrict__ rAPtr = rA.begin();
    scalar* __restrict__ rTPtr = rT.begin();

    // Calculate normalisation factor
    scalar normFactor = this->normFactor(x, b, wA, pA, cmpt);

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate normalised residual norm
    solverPerf.initialResidual() = gSumMag(rA)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Check convergence, solve if not converged
    if (!stop(solverPerf))
    {
        // Select and construct the preconditioner
        autoPtr<lduPreconditioner> preconPtr;

        preconPtr =
            lduPreconditioner::New
            (
                matrix_,
                coupleBouCoeffs_,
                coupleIntCoeffs_,
                interfaces_,
                dict()
            );

        // Solver iteration
        do
        {
            // Store previous wArT
            wArTold = wArT;

            // Precondition residuals
            preconPtr->precondition(wA, rA, cmpt);
            preconPtr->preconditionT(wT, rT, cmpt);

            // Update search directions:
            wArT = gSumProd(wA, rT);

            if (solverPerf.nIterations() == 0)
            {
                for (register label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell];
                    pTPtr[cell] = wTPtr[cell];
                }
            }
            else
            {
                scalar beta = wArT/wArTold;

                for (register label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
                    pTPtr[cell] = wTPtr[cell] + beta*pTPtr[cell];
                }
            }


            // Update preconditioned residuals
            matrix_.Amul(wA, pA, coupleBouCoeffs_, interfaces_, cmpt);
            matrix_.Tmul(wT, pT, coupleIntCoeffs_, interfaces_, cmpt);

            scalar wApT = gSumProd(wA, pT);


            // Test for singularity
            if (solverPerf.checkSingularity(mag(wApT)/normFactor)) break;


            // Update solution and residual:

            scalar alpha = wArT/wApT;

            for (register label cell=0; cell<nCells; cell++)
            {
                xPtr[cell] += alpha*pAPtr[cell];
                rAPtr[cell] -= alpha*wAPtr[cell];
                rTPtr[cell] -= alpha*wTPtr[cell];
            }

            solverPerf.finalResidual() = gSumMag(rA)/normFactor;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
