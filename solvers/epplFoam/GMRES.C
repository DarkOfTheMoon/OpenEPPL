/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "GMRES.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GMRES, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<GMRES>
        addGMRESSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<GMRES>
        addGMRESAsymMatrixConstructorToTable_;

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::GMRES::givensRotation
(
    const scalar& h,
    const scalar& beta,
    scalar& c,
    scalar& s
) const
{
    if (beta == 0)
    {
        c = 1;
        s = 0;
    }
    else if (mag(beta) > mag(h))
    {
        scalar tau = -h/beta;
        s = 1.0/Foam::sqrt(1.0 + sqr(tau));
        c = s*tau;
    }
    else
    {
        scalar tau = -beta/h;
        c = 1.0/Foam::sqrt(1.0 + sqr(tau));
        s = c*tau;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GMRES::GMRES
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    ),
    nDirs_(readLabel(solverControls.lookup("nDirections")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::GMRES::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_
    );

    const label nCells = psi.size();

    scalar* __restrict__ psiPtr = psi.begin();

    scalarField pA(nCells);
    scalar* __restrict__ pAPtr = pA.begin();

    scalarField wA(nCells);
    scalar* __restrict__ wAPtr = wA.begin();

    // --- Calculate A.psi
    matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // --- Calculate initial residual field
    scalarField rA(source - wA);
    scalar* __restrict__ rAPtr = rA.begin();

    // --- Calculate normalisation factor
    const scalar normFactor = this->normFactor(psi, source, wA, pA);

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() =
        gSumMag(rA, matrix().mesh().comm())/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Note: GMRES cannot be forced to do minIter sweeps
    // if the residual is zero, due to algorithmic reasons
    // HJ, 22/Aug/2012
    if (!solverPerf.checkConvergence(tolerance_, relTol_))
    {
        // Create the Hesenberg matrix
        scalarSquareMatrix H(nDirs_, 0);

        // Create y and b for Hessenberg matrix
        scalarField yh(nDirs_, 0);
        scalarField bh(nDirs_ + 1, 0);

        // Givens rotation vectors
        scalarField c(nDirs_, 0);
        scalarField s(nDirs_, 0);

        // Allocate Krylov space vectors
        FieldField<Field, scalar> V(nDirs_ + 1);

        forAll (V, i)
        {
            V.set(i, new scalarField(psi.size(), 0));
        }

        // --- Select and construct the preconditioner
        autoPtr<lduMatrix::preconditioner> preconPtr =
        lduMatrix::preconditioner::New
        (
            *this,
            controlDict_
        );

        // --- Solver iteration
        do
        {
            // --- Precondition residuals
            preconPtr->precondition(wA, rA, cmpt);

            // Calculate beta and scale first vector
            scalar beta = max(Foam::sqrt(gSumSqr(wA)), SMALL);

            // Set initial rhs and bh[0] = beta
            bh = 0;
            bh[0] = beta;

            for (label i = 0; i < nDirs_; i++)
            {
                // Set search direction
                V[i] = wA;
                V[i] /= beta;

                // Arnoldi's method
                matrix_.Amul(rA, V[i], interfaceBouCoeffs_, interfaces_, cmpt);

                // Execute preconditioning
                preconPtr->precondition(wA, rA, cmpt);

                for (label j = 0; j <= i; j++)
                {
                    beta = gSumProd(wA, V[j]);

                    H[j][i] = beta;

                    forAll (wA, wI)
                    {
                        wA[wI] -= beta*V[j][wI];
                    }
                }

                beta = max(Foam::sqrt(gSumSqr(wA)), SMALL);

                // Apply previous Givens rotations to new column of H.
                for (label j = 0; j < i; j++)
                {
                    const scalar Hji = H[j][i];
                    H[j][i] = c[j]*Hji - s[j]*H[j + 1][i];
                    H[j + 1][i] = s[j]*Hji + c[j]*H[j + 1][i];
                }

                // Apply Givens rotation to current row.
                givensRotation(H[i][i], beta, c[i], s[i]);

                const scalar bhi = bh[i];
                bh[i] = c[i]*bhi - s[i]*bh[i + 1];
                bh[i + 1] = s[i]*bhi + c[i]*bh[i + 1];
                H[i][i] = c[i]*H[i][i] - s[i]*beta;
            }

            // Back substitute to solve Hy = b
            for (label i = nDirs_ - 1; i >= 0; i--)
            {
                scalar sum = bh[i];

                for (label j = i + 1; j < nDirs_; j++)
                {
                    sum -= H[i][j]*yh[j];
                }

                yh[i] = sum/stabilise(H[i][i], VSMALL);
            }

            // --- Update solution and residual:

            for (label i = 0; i < nDirs_; i++)
            {
                const scalarField& Vi = V[i];
                const scalar& yi = yh[i];

                for(label cell=0; cell<nCells; cell++)
                {
                    psiPtr[cell] += yi*Vi[cell];
                }
            }

            // Re-calculate the residual
            matrix_.Amul(wA, psi, interfaceBouCoeffs_, interfaces_, cmpt);

            rA = source - wA;

            solverPerf.finalResidual() =
                gSumMag(rA, matrix().mesh().comm())/normFactor;
        } while
        (
            (
              ++solverPerf.nIterations() < maxIter_
            && !solverPerf.checkConvergence(tolerance_, relTol_)
            )
        );
    }

    return solverPerf;
}


// ************************************************************************* //
