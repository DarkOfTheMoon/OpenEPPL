// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);

//here is no need to relax the UEqn

fvOptions.constrain(UEqn);

if (eppl.momentumPredictor())
{
    U.storePrevIter();

    UInitRes = solve(relax(UEqn) == -fvc::grad(p));

    if (cmptMax(UInitRes.initialResidual())>momentumResCritia)
    {
        U = U.prevIter();
    }
}
