/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
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

#include "heterogeneousAbsorptionEmissionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiationModels
    {
        defineTypeNameAndDebug(heterogeneousAbsorptionEmissionModel, 0);
        defineRunTimeSelectionTable(heterogeneousAbsorptionEmissionModel, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousAbsorptionEmissionModel::heterogeneousAbsorptionEmissionModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh)
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiationModels::heterogeneousAbsorptionEmissionModel::~heterogeneousAbsorptionEmissionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::as(const label bandI) const
{
    return asCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::asCont(const label bandI) const
{
    return volScalarField::New
    (
        "asCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::borderAs(const label bandI) const
{
    return borderAsCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::borderAsCont(const label bandI) const
{
    return volScalarField::New
    (
        "borderAsCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::es(const label bandI) const
{
    return esCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::esCont(const label bandI) const
{
    return volScalarField::New
    (
        "esCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::borderEs(const label bandI) const
{
    return borderEsCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::borderEsCont(const label bandI) const
{
    return volScalarField::New
    (
        "borderEsCont",
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    );
}

Foam::dimensionedScalar
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::borderL(const label bandI) const
{
    return  dimensionedScalar("borderL", dimLength, 0.0);
}

Foam::tmp<Foam::volScalarField>
Foam::radiationModels::heterogeneousAbsorptionEmissionModel::addIntensity
(
    const label rayI,
    const volScalarField& ILambda
) const
{
    return ILambda;
}


// ************************************************************************* //
