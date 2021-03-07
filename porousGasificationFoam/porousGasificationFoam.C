/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    porousGasificationFoam

Description
    Solver for reactive flow through porous medium
    created by Pawel Jan Zuk

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulenceModel.H"
#include "psiChemistryModel.H"
#include "chemistrySolver.H"
#include "porousReactingZone.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "basicHGSSolidThermo.H"
#include "solidChemistryModel.H"
#include "heterogeneousPyrolysisModel.H"
#include "heterogeneousRadiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "readChemistryProperties.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    #include "createPyrolysisModel.H"
    #include "readPyrolysisTimeControls.H"
    #include "createPorosity.H"
    #include "createHeterogeneousRadiationModel.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {

        #include "readTimeControls.H"

        Info<< "Time = " << runTime.timeName() << nl;
        #include "compressibleCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"
        #include "updateChemistryTimeStep.H"
        Info<< "deltaT = " <<  runTime.deltaT().value() << nl << endl;

        #include "radiation.H"

        pyrolysisZone.evolve();

        #include "chemistry.H"

        // --- PIMPLE loop
        while (pimple.loop())
        {

            if (pimple.nCorrPIMPLE() > 0)
            {
                p.storePrevIter();
                rho.storePrevIter();
                U.storePrevIter();
            }

            #include "YEqn.H"
            #include "hsEqn.H"
            #include "rhoEqn.H"
            #include "UEqn.H"

            // --- PISO loop
	        while (pimple.correct())            
            {
                #include "pEqn.H"
            }
            turbulence->correct();

            rho = thermo.rho();
        }

        if (runTime.write())
        {
            chemistry.dQ()().write();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        runTime++;


    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
