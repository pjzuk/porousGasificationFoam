/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "solidOde.H"

#include "ODESolidHeterogeneousChemistryModel.H"

#include "HGSSolidThermo.H"
#include "psiReactionThermo.H"

#include "solidThermoPhysicsTypes.H"

#include "makeChemistrySolver.H"


// TMP, Filip: this should be using mostly OF2406 methods instead of carrying defined from FE. 
// in OF2406 there is src\thermophysicalModels\basic\fluidThermo\makeThermo.H
// which partly include OF8 OpenFOAM-8\src\thermophysicalModels\chemistryModel\chemistrySolver\makeChemistrySolver.H
// and all the below flatten from OF8.
// ------------------------------------------------------------------------------

#include "incompressiblePerfectGas.H"
#include "perfectGas.H"

#include "eConstThermo.H"
#include "hConstThermo.H"
#include "janafThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "thermo.H"


#define typedefThermo(Transport, Energy, Thermo, Equation, Specie)             \
                                                                               \
    typedef                                                                    \
        Transport                                                              \
        <                                                                      \
            Foam::species::thermo                                                    \
            <                                                                  \
                Thermo                                                         \
                <                                                              \
                    Equation                                                   \
                    <                                                          \
                        Specie                                                 \
                    >                                                          \
                >,                                                             \
                Energy                                                         \
            >                                                                  \
        >                                                                      \
        Transport##Energy##Thermo##Equation##Specie

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forThermo(Transport, Energy, Thermo, Equation, Specie, Macro, Args...) \
                                                                               \
    typedefThermo                                                              \
    (                                                                          \
        Transport,                                                             \
        Energy,                                                                \
        Thermo,                                                                \
        Equation,                                                              \
        Specie                                                                 \
    );                                                                         \
                                                                               \
    Macro(Args, Transport##Energy##Thermo##Equation##Specie)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forCommonGasEqns(Mu, He, Cp, Macro, Args...)                           \
    forThermo(Mu, He, Cp, incompressiblePerfectGas, specie, Macro, Args);      \
    forThermo(Mu, He, Cp, perfectGas, specie, Macro, Args)

#define forCommonGasEnergiesAndThermos(Mu, Macro, Args...)                     \
    forCommonGasEqns(Mu, sensibleEnthalpy, hConstThermo, Macro, Args);         \
    forCommonGasEqns(Mu, sensibleEnthalpy, janafThermo, Macro, Args);          \
    forCommonGasEqns(Mu, sensibleInternalEnergy, eConstThermo, Macro, Args);   \
    forCommonGasEqns(Mu, sensibleInternalEnergy, hConstThermo, Macro, Args);   \
    forCommonGasEqns(Mu, sensibleInternalEnergy, janafThermo, Macro, Args)

#define forCommonGasTransports(Macro, Args...)                                 \
    forCommonGasEnergiesAndThermos(constTransport, Macro, Args);               \
    forCommonGasEnergiesAndThermos(sutherlandTransport, Macro, Args)

#define forCommonGases(Macro, Args...)                                         \
    forCommonGasTransports(Macro, Args)


// ------------------------------------------------------------------------------

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistrySolvers(SolidThermo, SolidThermoPhysics, GasThermoPhysics) \
    defineChemistrySolver                                                      \
    (                                                                          \
        ODESolidHeterogeneousChemistryModel,                                   \
        SolidThermo,                                                           \
        SolidThermoPhysics,                                                    \
        GasThermoPhysics                                                       \
    );                                                                         

#define makeChemistrySolvers(Solver, SolidThermo, SolidThermoPhysics, GasThermoPhysics) \
    makeChemistrySolver                                                        \
    (                                                                          \
        Solver,                                                                \
        ODESolidHeterogeneousChemistryModel,                                   \
        SolidThermo,                                                           \
        SolidThermoPhysics,                                                    \
        GasThermoPhysics                                                       \
    );                                                                         
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCommonGases(defineChemistrySolvers, HGSSolidThermo, constSolidThermoPhysics);


    forCommonGases(makeChemistrySolvers, solidOde, HGSSolidThermo, constSolidThermoPhysics);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
