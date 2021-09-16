volScalarField rAU("rAU", 1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

MRF.makeRelative(phiHbyA);

adjustPhi(phiHbyA, U, p);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAU, MRF);

// Non-orthogonal pressure corrector loop
while (eppl.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(eppl.finalInnerIter())));

    if (eppl.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
    else if (pressureRDC < (1. - SMALL))
    {
        p.relax(pressureRDC);
        p.storePrevIter();
    }
}

#include "continuityErrs.H"

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
