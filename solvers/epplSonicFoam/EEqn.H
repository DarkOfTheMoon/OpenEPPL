{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + deferred::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - (
            fvc::laplacian(turbulence->kappaEff(), T)
          + correction(fvm::laplacian(turbulence->alphaEff(), he))
        )
      + (turbulence->devRhoReff() && fvc::grad(U))
     ==
        fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
