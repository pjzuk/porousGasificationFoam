rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();
U = rUA*UEqn.H();

phi =
    fvc::interpolate(rho)*
    (
        (fvc::interpolate(U) & mesh.Sf())
    );

while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p)
      + fvc::div(phi)
      - fvm::laplacian(rho*rUA, p)
    == 
        Srho*(1.0 - porosityF)  
    );

    pEqn.solve
    (
        mesh.solutionDict().solver(p.select(pimple.finalInnerIter()))
    );

    if (pimple.finalNonOrthogonalIter())
    {
        phi += pEqn.flux();
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

//if (oCorr != nOuterCorr-1)
{
    // Explicitly relax pressure for momentum corrector
    p.relax();

    rho = thermo.rho();
    rho.relax();
    Info<< "rho max/min : " << max(rho).value()
        << " " << min(rho).value() << endl;
}

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();


