// eqZx2uHGn044

if (solveEnergy)
{
    
    if (pyrolysisZone.equilibrium())
    {}
    else
    {
        volScalarField heatUpGas = pyrolysisZone.heatUpGas()();
        heatTransfer = pyrolysisZone.heatTransfer()();

        Info<< " rho gas min/max = " << gMin(rho) << ", " << gMax(rho) << endl;
        Info<< " Sh gas min/max = " << gMin(chemistrySh) << ", " << gMax(chemistrySh) << endl;
        Info<< " porosity min/max = " << gMin(porosityF) << ", " << gMax(porosityF) << endl;
        Info<< " heatUpGas min/max = " << gMin(heatUpGas) << ", " << gMax(heatUpGas) << endl;
        Info<< " heatTransfer min/max = " << gMin(heatTransfer) << ", " << gMax(heatTransfer) << endl;


        fvScalarMatrix hsEqn
        (
            fvm::ddt(porosityF*rho, hs)
          + fvm::div(phi, hs)
          - fvm::laplacian(turbulence->alphaEff(), hs)
         ==
           chemistrySh*porosityF
         + heatTransfer  
         + heatUpGas     
         + radiation->Shs(thermo)
        );

        hsEqn.relax();
        hsEqn.solve();

        Info<< " hs min/max = " << gMin(hs) << ", " << gMax(hs) << endl;
        Info<< " T gas min/max   = " << gMin(T) << ", "
            << gMax(T) << endl;

        thermo.correct();
    }
}
