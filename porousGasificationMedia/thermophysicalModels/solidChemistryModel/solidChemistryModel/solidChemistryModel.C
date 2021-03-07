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

#include "solidChemistryModel.H"
#include "fvMesh.H"
#include "foamTime.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidChemistryModel, 0);
//    defineRunTimeSelectionTable(solidChemistryModel, fvMesh);
    defineRunTimeSelectionTable(solidChemistryModel, fvMeshGasPhaseGases);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidChemistryModel::solidChemistryModel
(
    const fvMesh& mesh,
    const objectRegistry& obj,
    const word& compTypeName,
    const word& solidThermoTypeName
)
:
    basicChemistryModel(mesh,obj),
    solidThermo_(basicHGSSolidThermo::New(mesh))
{
}



Foam::solidChemistryModel::solidChemistryModel
(
    const fvMesh& mesh,
    const objectRegistry& obj,
    const word& compTypeName,
    const word& solidThermoTypeName,
    PtrList<volScalarField>& gasPhaseGases
)
:
    basicChemistryModel(mesh,obj),
    solidThermo_(basicHGSSolidThermo::New(
	mesh,
	IOdictionary
        (
            IOobject
            (
                "solidThermophysicalProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ),
	gasPhaseGases))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidChemistryModel::~solidChemistryModel()
{}


// ************************************************************************* //
