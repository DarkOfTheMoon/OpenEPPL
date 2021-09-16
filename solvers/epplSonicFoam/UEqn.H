// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + deferred::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(U)
 ==
    fvOptions(rho, U)
);

//here is no need to relax the UEqn

fvOptions.constrain(UEqn);

U.storePrevIter();

UInitRes = solve(relax(UEqn) == -fvc::grad(p));

if (mag(UInitRes.initialResidual())>momentumResCritia)
{
    U = U.prevIter();
}

fvOptions.correct(U);
K = 0.5*magSqr(U);
