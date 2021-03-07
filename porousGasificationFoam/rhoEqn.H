// eqZx2uHGn042
{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(porosityF,rho)
      + fvc::div(phi)
      ==
        Srho*(1.0 - porosityF)
    );

    rhoEqn.solve();
}

