// eqZx2uHGn043

// Solve the Momentum equation

fvVectorMatrix UEqn
(   
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff()
 ==
    rho*g
);
    
UEqn.relax();

pZone.addResistance(UEqn,Df);  // eqZx2uHGn015

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
}
else
{   
    volScalarField rUA = 1.0/UEqn.A();
    U = rUA*(UEqn.H() - fvc::grad(p));
    U.correctBoundaryConditions();
}

