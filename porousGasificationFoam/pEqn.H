volScalarField rAU(1.0/UEqn.A());
U = rAU*UEqn.H();

phi =
    fvc::interpolate(rho)
   *(
        (fvc::interpolate(U) & mesh.Sf())
    );

while (pimple.correctNonOrthogonal())
{

    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p)
      + fvc::div(phi)
      - fvm::laplacian(rho*rAU, p)
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

U = U - rAU*fvc::grad(p);

U.correctBoundaryConditions();
K = 0.5*magSqr(U);
