rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();
U = rUA*UEqn.H();

phi =
    fvc::interpolate(rho)
   *(
        (fvc::interpolate(U) & mesh.Sf())
    );

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(porosityF*psi,p)
      + fvc::div(phi)
      - fvm::laplacian(rho*rUA, p)
    == 
        Srho*(1.0 - porosityF)  
    );

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi += pEqn.flux();
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

if (!pimple.finalIter())
{
    p.relax();
    rho.relax();
}

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();

