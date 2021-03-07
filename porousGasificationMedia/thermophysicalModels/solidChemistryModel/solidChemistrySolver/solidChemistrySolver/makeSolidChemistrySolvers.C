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

\*---------------------------------------------------------------------------*/

#include "thermoPhysicsTypes.H"
#include "chemistrySolver.H"

#include "psiChemistryModel.H"

#include "EulerImplicit.H"
#include "ode.H"
#include "sequential.H"

#include "solidThermoPhysicsTypes.H"
#include "thermoPhysicsTypes.H"

#include "solidChemistrySolver.H"

#include "ODESolidHeterogeneousChemistryModel.H"
#include "solidChemistryModel.H"

#include "solidOde.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeSolidChemistrySolver(solidChemistryModel,constSolidThermoPhysics,gasThermoPhysics)

    makeSolidChemistrySolverType
    (
        solidOde,
        solidChemistryModel,
        constSolidThermoPhysics,
        gasThermoPhysics
    )

    makeSolidChemistrySolver(solidChemistryModel,expoSolidThermoPhysics,gasThermoPhysics)

    makeSolidChemistrySolverType
    (
        solidOde,
        solidChemistryModel,
        expoSolidThermoPhysics,
        gasThermoPhysics
    )

    makeSolidChemistrySolver(solidChemistryModel,linearSolidThermoPhysics,gasThermoPhysics)

    makeSolidChemistrySolverType
    (
        solidOde,
        solidChemistryModel,
        linearSolidThermoPhysics,
        gasThermoPhysics
    )

}

// ************************************************************************* //
