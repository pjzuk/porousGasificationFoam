// eqZx2uHGn041
{
    volScalarField Yt = 0.0*Y[0];
    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];
            fvScalarMatrix YiEqn
            (
                fvm::ddt(porosityF*rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
             ==
                kappa*chemistry.RR(i)
	          + pyrolysisZone.Srho(i)()
            );
            if (mesh.solutionDict().relaxEquation("Yi")) 
            {
                YiEqn.relax(mesh.solutionDict().equationRelaxationFactor("Yi"));
            }
            YiEqn.solve(mesh.solutionDict().solver("Yi"));
            Yi.max(0.0);
	        Yt += Yi;
            if (gMax(Yi) > 2.0) Info << "warning: mass fraction exceeded 2" << endl;
        }
        else
        {
            volScalarField& Yi = Y[i];
        
            fvScalarMatrix YiEqn
            (
                fvm::ddt(porosityF*rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
            );

            if (mesh.solutionDict().relaxEquation("Yi")) 
            {
                YiEqn.relax(mesh.solutionDict().equationRelaxationFactor("Yi"));
            }
            YiEqn.solve(mesh.solutionDict().solver("Yi"));
            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Info << "srho min/max " <<  gMin(Srho) << " " << gMax(Srho) << endl;

    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie)
        {
	        Info << " gas  "<<Y[i].name()<< " Srho min/max   = " << gMin(pyrolysisZone.Srho(i)()) << ", " ;
            Info << gMax(pyrolysisZone.Srho(i)()) << " RR min/max   = " << min(chemistry.RR(i)).value() << ", " ;
            Info << max(chemistry.RR(i)).value() << " min Y = " << gMin(Y[i]) <<" max Y = " << gMax(Y[i]);
            Y[i] = Y[i]/Yt;
            Info<< " min Ynorm = " << gMin(Y[i]) <<" max Ynorm = " << gMax(Y[i]) << endl;		
        }
        else
        {
        	Info << " gas  "<<Y[i].name()<< " is inert specie, min Y = " << gMin(Y[i]) <<" max Y = " << gMax(Y[i]);
            Y[i] = Y[i]/Yt;
            Info<< " min Ynorm = " << gMin(Y[i]) <<" max Ynorm = " << gMax(Y[i]) << endl;			
        }
    } 
}
