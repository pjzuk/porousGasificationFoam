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
    totalMassBiomassGasificationFoam

Description
    Integrates solid state mass over all domain in each timestep
    and writes it down to totalMass.txt file
    created by Pawel Jan Zuk

\*---------------------------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <string>

#include "fvc.H"
#include "OSspecific.H"
#include "fixedValueFvPatchFields.H"
#include "calc.H"

#include "timeSelector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

std::ofstream masa("totalMass.txt", ios_base::trunc);

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{

#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

if (rhoSHeader.headerOk() and porosityFHeader.headerOk())
{

volScalarField rhos(rhoSHeader,mesh);
volScalarField porosityF(porosityFHeader,mesh);

dimensionedScalar totalMass = fvc::domainIntegrate(rhos*(1.-porosityF));

masa << runTime.timeName() << " " << totalMass.value() << std::endl;

}

}

// ************************************************************************* //
