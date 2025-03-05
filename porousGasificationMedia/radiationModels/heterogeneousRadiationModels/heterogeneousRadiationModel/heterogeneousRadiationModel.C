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

#include "heterogeneousRadiationModel.H"
#include "heterogeneousAbsorptionEmissionModel.H"
#include "scatterModel.H"
#include "sootModel.H"
#include "fvmSup.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heterogeneousRadiationModel, 0);
    defineRunTimeSelectionTable(heterogeneousRadiationModel, T);
    defineRunTimeSelectionTable(heterogeneousRadiationModel, porosity);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::heterogeneousRadiationModel::initialise()
{
    if (radiation_)
    {   
        solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));

        heterogeneousAbsorptionEmission_.reset
        (
            radiationModels::heterogeneousAbsorptionEmissionModel::New(*this, mesh_).ptr()
        );
        
        scatter_.reset(Foam::radiation::scatterModel::New(*this, mesh_).ptr());

        soot_.reset(Foam::radiation::sootModel::New(*this, mesh_).ptr());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heterogeneousRadiationModel::heterogeneousRadiationModel(const volScalarField& T)
:
    radiationModel(type(), T),
    Ts_(T),
    heterogeneousAbsorptionEmission_(nullptr)
{
    initialise();
}


Foam::heterogeneousRadiationModel::heterogeneousRadiationModel
(
    const word& type,
    const volScalarField& T
)
:
    radiationModel(type, T),
    Ts_(T),
    heterogeneousAbsorptionEmission_(nullptr)
{
    initialise();
}


Foam::heterogeneousRadiationModel::heterogeneousRadiationModel
(
    const volScalarField& T,
    const volScalarField& porosityF,
    const List<label>& surfF,
    const volScalarField& Ts
)
:
    radiationModel(type(), T),
    Ts_(Ts),
    heterogeneousAbsorptionEmission_(nullptr)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::heterogeneousRadiationModel::~heterogeneousRadiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::radiationModels::heterogeneousAbsorptionEmissionModel&
Foam::heterogeneousRadiationModel::heterogeneousAbsorptionEmission() const
{
    if (!heterogeneousAbsorptionEmission_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation heterogeneousAbsorptionEmission model, but model is "
            << "not active" << abort(FatalError);
    }

    return heterogeneousAbsorptionEmission_;
}


Foam::label Foam::heterogeneousRadiationModel::nBands() const
{
    if (!heterogeneousAbsorptionEmission_.valid())
    {
        FatalErrorInFunction
            << "Requested nBands, but heterogeneousAbsorptionEmission model is "
            << "not active" << abort(FatalError);
    }
    return heterogeneousAbsorptionEmission_->nBands();
}

// ************************************************************************* //
