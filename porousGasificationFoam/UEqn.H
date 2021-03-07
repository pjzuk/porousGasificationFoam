// eqZx2uHGn043

    fvVectorMatrix UEqn
    (
	    fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff()
	==
	    rho*g
    );

    pZone.addResistance(UEqn,Df);  // eqZx2uHGn015

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
	    K = 0.5*magSqr(U);
    }

