
 
if (chemistry.chemistry())
{
    Info<< "Solving gas phase chemistry" << endl;
    const scalar TgasChemistry = chemistry.solve
    (
        runTime.value() - runTime.deltaTValue(),
        runTime.deltaTValue()
    );

    if (adjustTimeStep and gasChemistryTimeStepControl)
    {
        const scalar dt0 = runTime.deltaTValue();
        runTime.setDeltaT
        (
            min
    	    (
     	        TgasChemistry,
    	        dt0 
    	    )
        );
    }

    // turbulent time scale
    if (turbulentReaction)
    {
        tmp<volScalarField> tepsilon(turbulence->epsilon());
        const volScalarField& epsilon = tepsilon();
        tmp<volScalarField> tmuEff(turbulence->muEff());
        const volScalarField& muEff = tmuEff();
        tmp<volScalarField> ttc(chemistry.tc());
        const volScalarField& tc = ttc();

        forAll(epsilon, i)
        {
            if (epsilon[i] > 0)
            {
                // Chalmers PaSR model
                scalar tk = Cmix.value()*Foam::sqrt(muEff[i]/rho[i]/epsilon[i]);
                kappa[i] =
                    (runTime.deltaTValue() + tc[i])
                   /(runTime.deltaTValue() + tc[i] + tk);
            }
            else
            {
                // Return to laminar combustion
                kappa[i] = 1.0;
            }
        }
    }
    else
    {
        kappa = 1.0;
    }

    chemistrySh = kappa*chemistry.Sh()();
}

Srho = pyrolysisZone.Srho()();
