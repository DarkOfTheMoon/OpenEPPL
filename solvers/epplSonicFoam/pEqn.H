rho = thermo.rho();

volScalarField rAU("rAU", 1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::interpolate(rho)*fvc::flux(HbyA));

MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAUf, MRF);

surfaceScalarField phid
(
    "phid",
    (fvc::interpolate(psi)/fvc::interpolate(rho))*phiHbyA
);

// Non-orthogonal pressure corrector loop
while (eppl.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p)
      + fvm::div(phid, p)
      - fvm::laplacian(rhorAUf, p)
     ==
        fvOptions(psi, p, rho.name())
    );

    pEqn.solve();

    if (eppl.finalNonOrthogonalIter())
    {
        phi = pEqn.flux();
    }
    else if(pressureRDC < (1. - SMALL))
    {
        p.relax(pressureRDC);
        p.storePrevIter();

        rhorAUf = fvc::interpolate(thermo.rho()*rAU);
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

{
    volScalarField::Boundary& rhoBnd = rho.boundaryFieldRef();

    forAll(rhoBnd, patchi)
    {
        rhoBnd[patchi] = thermo.rho(patchi);
    }
}

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
